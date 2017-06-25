//#define TBB_USE_DEBUG 1
#include "pch.h"
#include <thread>
#include <condition_variable>
#include <sys/stat.h>
#include <papi.h>
#include <ittnotify.h>
#include <sys/syscall.h>
#include <boost/graph/topological_sort.hpp>

#include "vpr_types.h"

#include "log.h"
#include "graph.h"
#include "route.h"
#include "misr.h"
#include "utility.h"
#include "init.h"
#include "delta_stepping.h"
#include "route_tree.h"
#include "congestion.h"
#include "router.h"
#include "geometry.h"
#include "clock.h"
#include "new_partitioner.h"
#include "barrier.h"
#include "dijkstra.h"
#include "det_mutex.h"

using namespace std;

//using timer = std::chrono::high_resolution_clock;
using timer = myclock;
using std::chrono::duration_cast;
using std::chrono::nanoseconds;

//#define PRINT_BBS
//#define TEST

#if 0
#define PRINT(msg, ...) printf((msg), __VA_ARGS__)
#else
#define PRINT(msg, ...) 
#endif

//#define PTHREAD_BARRIER

__itt_domain* domain = __itt_domain_create("MyTraces.MyDomain");
__itt_string_handle* shMyTask = __itt_string_handle_create("Route Task");
__itt_string_handle* shAnalysisTask = __itt_string_handle_create("Analysis Task");

typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS> UndirectedGraph;
typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::directedS> DirectedGraph;

struct extra_route_state_t {
	int same;
	int ortho;
	bool existing;
	float upstream_R;
	float delay;
};

struct dijkstra_stats_t {
	unsigned long num_heap_pops;
	unsigned long num_heap_pushes;
	unsigned long num_neighbor_visits;
};

struct global_route_state_t {
	congestion_t *congestion;
	route_tree_t *route_trees;
	//route_tree_t *back_route_trees;
	vector<RRNode> *added_rr_nodes;	
	vector<RRNode> *back_added_rr_nodes;	
	t_net_timing *net_timing;
};

struct runtime_t {
	timer::duration route_time;
	timer::duration wait_time;
	timer::duration get_net_wait_time;
	timer::duration pure_route_time;
	timer::duration last_barrier_wait_time;
	timer::duration last_sync_time;
};

struct new_virtual_net_t {
	int local_index;
	int global_index;
	const net_t *net;
	const source_t *source;
	vector<const sink_t *> sinks;
	//box bounding_box;
	vector<box> overlap_g_boxes;
	vector<vector<box>> bounding_boxes;
};

struct fpga_tree_node_t {
	box bb;

	vector<net_t *> nets;

	vector<new_virtual_net_t> vnets;

	/* coloring */
	vector<vector<new_virtual_net_t *>> scheduled;

	/* undirected to directed */
	vector<new_virtual_net_t *> roots;
	vector<vector<int>> net_g;
	DirectedGraph net_g_2;
	vector<int> num_incoming_edges;
	vector<std::atomic<int> *> current_num_incoming_edges;

	timer::duration route_time;

	fpga_tree_node_t *next;

	fpga_tree_node_t *left;
	fpga_tree_node_t *right;
};

struct /*alignas(64)*/ router_t {
	vector<net_t> all_nets;
	vector<net_t *> current_nets;
	vector<new_virtual_net_t *> all_vnets;
	RRGraph g;
	std::atomic<int> iter;
	global_route_state_t state;
	bool phase_two;

	fpga_tree_node_t net_root;
	int num_levels;

	/* params */
	int num_threads;
	route_parameters_t params;
	float pres_fac;
	int pmc_overflow;

	/* debug */
	//tbb::spin_mutex *locks;

	/* per thread */
	//vector<timer::duration> total_time;

	/* per level, per thread */
	//vector<vector<runtime_t>> time;
	//vector<vector<unsigned long>> last_clock;
	//vector<vector<unsigned long>> num_nets_routed;
	//vector<vector<unsigned long>> num_locks;
	//vector<vector<dijkstra_stats_t>> stats;

	/* per nets */
	vector<bool> routed;
	vector<timer::duration> route_time;
	//vector<timer::duration> acc_net_route_time;
	//vector<unsigned long> acc_num_clock_ticks;
	vector<timer::duration> update_time;
	vector<dijkstra_stats_t> net_stats;
	vector<vector<dijkstra_stats_t>> sink_stats;
};

//static congestion_t **g_congestion;
//static tbb::spin_mutex g_lock;

static char g_dirname[256];

static class Manager *g_manager = nullptr;

double boxes_area(const vector<box> &boxes)
{
	double intersect_area = 0;
	assert(boxes.size() <= 2);
	//printf("start boxes area\n");
	for (int i = 0; i < boxes.size(); ++i) {
		for (int j = i + 1; j < boxes.size(); ++j) {
			if (bg::intersects(boxes[i], boxes[j])) {
				box intersection;
				bg::intersection(boxes[i], boxes[j], intersection);
				//printf("%d %d intersects area %g\n", i, j, bg::area(intersection));
				intersect_area += bg::area(intersection);
			}
		}
	}
	double area = 0;
	//printf("start all boxes %lu\n", boxes.size());
	for (int i = 0; i < boxes.size(); ++i) {
		//printf("box i=%d %d %d, %d %d\n", i, bg::get<bg::min_corner, 0>(boxes[i]), bg::get<bg::max_corner, 0>(boxes[i]), bg::get<bg::min_corner, 1>(boxes[i]), bg::get<bg::max_corner, 1>(boxes[i]));
		area += bg::area(boxes[i]);
	}

	double actual_area = area-intersect_area;
	//printf("actual area %g area %g intersect_area %g\n", actual_area, area, intersect_area);
	assert(actual_area >= 0);

	return actual_area;
}

void flatten(const vector<vector<box>> &boxes, vector<box> &flat)
{
	for (int i = 0; i < boxes.size(); ++i) {
		for (int j = 0; j < boxes[i].size(); ++j) {
			flat.push_back(boxes[i][j]);
		}
	}
}

bool intersects(const vector<vector<box>> &a, const vector<vector<box>> &b)
{
	vector<box> flat_a;
	vector<box> flat_b;

	flatten(a, flat_a);
	flatten(b, flat_b);

	bool inters = false;
	for (int i = 0; i < flat_a.size() && !inters; ++i) {
		const auto &box_a = flat_a[i];
		for (int j = 0; j < flat_b.size() && !inters; ++j) {
			const auto &box_b = flat_b[j];
			inters = bg::intersects(box_a, box_b);
		}
	}
	return inters;
}

double all_boxes_area(const vector<vector<box>> &boxes)
{
	vector<box> flat;

	flatten(boxes, flat);

	return boxes_area(flat);
}

box make_bb(const vector<point> &points, int expand)
{
	box bb;
	bg::assign_inverse(bb);

	for (const auto &p : points) {
		bg::expand(bb, p);
	}

	bg::add_value(bb.max_corner(), expand);
	bg::subtract_value(bb.min_corner(), 1+expand);

	extern int nx, ny;
	clip(bb, make_pair(nx+1, ny+1));

	return bb;
}

static bool inside_bb(const rr_node_property_t &node, const box &bb)
{
	//int xlow, xhigh, ylow, yhigh;

	//switch (node.type) {
		//case CHANX:
		//case CHANY:
			//if (node.inc_direction) {
				//xlow = node.xlow;
				//ylow = node.ylow;
			//} else {
				//xlow = node.xhigh;
				//ylow = node.yhigh;
			//}
			//xhigh = xlow;
			//yhigh = ylow;
			//break;

		//case IPIN:
		//case OPIN:
			//xlow = node.xlow;
			//ylow = node.ylow;
			//xhigh = node.xhigh;
			//yhigh = node.ylow + node.pin_height;
			//assert(xlow == xhigh);
			//assert(ylow <= yhigh);
			//break;

		//case SOURCE:
		//case SINK:
			//xlow = node.xlow;
			//ylow = node.ylow;
			//xhigh = node.xhigh;
			//yhigh = node.yhigh;
			//break;

		//default:
			//assert(false);
			//break;
	//}

	bool inside;

	if (node.real_xhigh < bg::get<bg::min_corner, 0>(bb)
			|| node.real_xlow > bg::get<bg::max_corner, 0>(bb)
			|| node.real_yhigh < bg::get<bg::min_corner, 1>(bb)
			|| node.real_ylow > bg::get<bg::max_corner, 1>(bb)) {
		inside = false;
	} else {
		inside = true;
	}

	return inside;
}

static bool inside_bbs(const rr_node_property_t &node, const vector<box> &boxes)
{
	bool inside = false;
	for (int i = 0; i < boxes.size() && !inside; ++i) {
		inside = inside_bb(node, boxes[i]);
	}
	return inside;
}

static bool inside_region(const rr_node_property_t &node, const box &bb)
{
	int x, y;
	if (node.type == CHANX || node.type == CHANY) {
		if (node.inc_direction) {
			x = node.xlow;
			y = node.ylow;
		} else {
			x = node.xhigh;
			y = node.yhigh;
		}
	} else {
		x = node.xlow;
		y = node.ylow;
	}

	bool inside;
	if (x >= bg::get<bg::min_corner, 0>(bb)
			&& x <= bg::get<bg::max_corner, 0>(bb)
			&& y >= bg::get<bg::min_corner, 1>(bb)
			&& y <= bg::get<bg::max_corner, 1>(bb)) {
		inside = true;
	} else {
		inside = false;
	}

	return inside;
}

struct dijkstra_state_t {
	float *known_distance;
	float *distance;
	RREdge *prev_edge;
	extra_route_state_t *state;
};

class SinkRouter {
	private:
		const RRGraph &_g;

		float _astar_fac;

		float *_known_distance;
		float *_distance;
		RREdge *_prev_edge;
		extra_route_state_t *_state;

		congestion_t *_congestion;
		const float &_pres_fac;

		vector<int> _modified_nodes;
		vector<bool> _modified_node_added;

		RRNode _existing_opin;

		const sink_t *_current_sink;
		const vector<box> *m_bounding_boxes;

		dijkstra_stats_t _stats;

	private:
		void popped_node(const heap_node_t<RREdge, extra_route_state_t> &node)
		{
			char buffer[256];
			sprintf_rr_node(node.node, buffer);
			zlog_level(delta_log, ROUTER_V3, "Current: %s Existing: %d [same %d ortho %d] [kd: %g %a okd: %g %a] [d: %g %a od: %g %a] prev: %d\n",
					buffer, node.extra.existing, node.extra.same, node.extra.ortho,
					node.known_distance, node.known_distance,
					_known_distance[node.node], _known_distance[node.node],
					node.distance, node.distance,
					_distance[node.node], _distance[node.node],
					get_source(_g, node.prev_edge));

			++_stats.num_heap_pops;
		}

		void relax_node(const heap_node_t<RREdge, extra_route_state_t> &node)
		{
			//RouteTreeNode rt_node = route_tree_get_rt_node(*_current_rt, v);

			//if (rt_node != RouteTree::null_vertex()) {
				//const auto &rt_node_p = get_vertex_props(_current_rt->graph, rt_node);
				//if (!valid(_prev_edge[v])) {
					//assert(!valid(rt_node_p.rt_edge_to_parent));
				//} else {
					//int rr = get_vertex_props(_current_rt->graph, get_source(_current_rt->graph, rt_node_p.rt_edge_to_parent)).rr_node;
					//assert(rr == get_source(_g, _prev_edge[v]));
				//}
			//}

			//int v = node.node;
			//const auto &v_p = get_vertex_props(_g, v);

			//if (valid(node.prev_edge)) {
				//assert(v == get_target(_g, node.prev_edge));

				//int u = get_source(_g, node.prev_edge);

				//float u_delay;
				//float u_upstream_R;
				//if (_state[u].upstream_R != std::numeric_limits<float>::max()) {
					//assert(_state[u].delay != std::numeric_limits<float>::max());
					//u_upstream_R = _state[u].upstream_R;
					//u_delay = _state[u].delay;
				//} else {
					//auto rt_node = route_tree_get_rt_node(*_current_rt, u);
					//assert(rt_node != RouteTree::null_vertex());
					//const auto &u_rt_node_p = get_vertex_props(_current_rt->graph, rt_node);
					//u_upstream_R = u_rt_node_p.upstream_R;
					//u_delay = u_rt_node_p.delay;
				//}

				//const auto &e_p = get_edge_props(_g, node.prev_edge);

				//extern struct s_switch_inf *switch_inf;
				//const struct s_switch_inf *sw = &switch_inf[e_p.switch_index];

				//_state[v].upstream_R = sw->R + v_p.R;
				//if (!sw->buffered)  {
					//_state[v].upstream_R += u_upstream_R;
				//} 

				//float delay;
				//if (sw->buffered) {
					//delay = sw->Tdel + v_p.C * (sw->R + 0.5 * v_p.R);
				//} else {
					//delay = sw->Tdel + v_p.C * (u_upstream_R + sw->R + 0.5 * v_p.R);
				//}
				//_state[v].delay = u_delay + delay;
			//} else {
				//_state[v].upstream_R = v_p.R;
				//_state[v].delay = v_p.C * 0.5 * v_p.R;
			//}
			_state[node.node] = node.extra;

			if (!_modified_node_added[node.node]) {
				_modified_nodes.push_back(node.node);
				_modified_node_added[node.node] = true;
			}
		}

		bool expand_node(int node)
		{
			++_stats.num_neighbor_visits;
			
			const auto &prop = get_vertex_props(_g, node);

			char buffer[256];
			sprintf_rr_node(node, buffer);

			zlog_level(delta_log, ROUTER_V3, "\tChecking whether to expand %s ", buffer);

			//const auto &sink_prop = get_vertex_props(_g, _current_sink->rr_node);

			// this has to come before the bounding box checking 
			if (prop.type == IPIN) { 
				assert(num_out_edges(_g, node) == 1);
				int target = get_target(_g, *begin(get_out_edges(_g, node)));
				assert(get_vertex_props(_g, target).type == SINK);

				/* this condition is sufficient because the bounding box would have covered
				 * the SINK */
				if (target != _current_sink->rr_node) {
					zlog_level(delta_log, ROUTER_V3, "not target IPIN\n");
					return false;
				} /*else {
					//zlog_level(delta_log, ROUTER_V3, "found target IPIN\n");
					//return true;
				}*/
			}

			if (prop.type == OPIN && _existing_opin != RRGraph::null_vertex() && node != _existing_opin) {
				zlog_level(delta_log, ROUTER_V3, "not expanding other OPIN\n");
				return false;
			}

			if (!inside_bbs(prop, *m_bounding_boxes)) {
				zlog_level(delta_log, ROUTER_V3, "outside of bounding box\n");
				return false;
			}

			zlog_level(delta_log, ROUTER_V3, " OK\n");

			return true;
		}

		void get_edge_weight(const RREdge &e, float &known_distance, float &distance, extra_route_state_t &extra)
		{
			int u = get_source(_g, e);
			int v = get_target(_g, e);

			const auto &v_p = get_vertex_props(_g, v);
			const auto &e_p = get_edge_props(_g, e);

			extern struct s_switch_inf *switch_inf;
			const struct s_switch_inf *sw = &switch_inf[e_p.switch_index];

			assert(_state[u].upstream_R != std::numeric_limits<float>::max());

			float delay;
			if (sw->buffered) {
				delay = sw->Tdel + v_p.C * (sw->R + 0.5 * v_p.R);
			} else {
				delay = sw->Tdel + v_p.C * (_state[u].upstream_R + sw->R + 0.5 * v_p.R);
			}
			assert(_state[u].delay != std::numeric_limits<float>::max());
			extra.delay = _state[u].delay + delay; 

			extra.upstream_R = sw->R + v_p.R;
			if (!sw->buffered) {
				extra.upstream_R += _state[u].upstream_R;
			} 

			extra.existing = false;

			extern t_rr_indexed_data *rr_indexed_data;

			float congestion_cost = rr_indexed_data[v_p.cost_index].base_cost * get_acc_cost(_congestion, v) * get_pres_cost(_congestion, v);

			known_distance = _current_sink->criticality_fac * delay + (1. - _current_sink->criticality_fac) * congestion_cost;

			float expected_cost = _astar_fac * get_timing_driven_expected_cost(v_p, get_vertex_props(_g, _current_sink->rr_node), _current_sink->criticality_fac, extra.upstream_R, &extra.same, &extra.ortho);

			distance = known_distance + expected_cost;

			unsigned int xe, xk;
			//memcpy(&xe, &expected_cost, sizeof (xe));
			//memcpy(&xk, &known_distance, sizeof (xk));

			zlog_level(delta_log, ROUTER_V3, "\t%d -> %d delay %g upstream_R %g congestion %g crit_fac %g known %g %a %X expected %g %a %X same %d ortho %d predicted %g %a\n", 
					u, v, delay, extra.upstream_R, congestion_cost, _current_sink->criticality_fac, known_distance, known_distance, xk, expected_cost, expected_cost, xe, extra.same, extra.ortho, distance, distance);
			//zlog_level(delta_log, ROUTER_V3, "\t[u: upstream %g] [edge: d %g R %g] [v: R %g C %g]\n",
					//_state[u].upstream_R, sw->Tdel, sw->R, v_p.R, v_p.C);
		}

		void push_node(const heap_node_t<RREdge, extra_route_state_t> &node)
		{
			zlog_level(delta_log, ROUTER_V3, "\tPushing neighbor %d [kd %g %a d %g %a] to heap\n",
					node.node,
					node.known_distance, node.known_distance,
					node.distance, node.distance);

			++_stats.num_heap_pushes;
		}

		RREdge get_previous_edge(int rr_node_id, const route_tree_t &rt)
		{
			RREdge previous_edge;

			if (!valid(_prev_edge[rr_node_id])) {
				previous_edge = RRGraph::null_edge();
			} else {
				RouteTreeNode rt_node = route_tree_get_rt_node(rt, rr_node_id);

				char buffer[256];
				sprintf_rr_node(rr_node_id, buffer);

				if (rt_node != RouteTree::null_vertex() && valid(get_vertex_props(rt.graph, rt_node).rt_edge_to_parent)) {

					/* already reach an existing route tree node but the current state suggests there's a possible path of 
					 * lower cost */
					/*if (rt_node->rr_edge_to_parent != -1) {*/
					/*char parent[256];*/
					/*sprintf_rr_node(get_source(g, rt_node->rr_edge_to_parent), parent);*/
					/*zlog_error(delta_log, "Error: Existing route tree node %s has non-null rr_edge_to_parent that connects to %s\n", buffer, parent);*/
					/*assert(false);*/
					/*}*/

					char s_state[256];
					char s_rt[256];
					sprintf_rr_node(get_source(_g, _prev_edge[rr_node_id]), s_state);
					sprintf_rr_node(get_vertex_props(rt.graph, get_source(rt.graph, get_vertex_props(rt.graph, rt_node).rt_edge_to_parent)).rr_node, s_rt);
					zlog_warn(delta_log, "Warning: Existing route tree node %s does not have a matching route state. (state.prev_edge: %s rt_node.rr_edge_to_parent: %s) because we have found a shorter path to that node\n", buffer, s_state, s_rt);

					previous_edge = RRGraph::null_edge();
				} else {
					previous_edge = _prev_edge[rr_node_id];
				}
			} 

			return previous_edge;
		}

		void backtrack(int sink_node, route_tree_t &rt, vector<RRNode> &added_rr_nodes)
		{
			RRNode rr_node = sink_node;
			RRNode child_rr_node = RRGraph::null_vertex();
			RREdge edge = RRGraph::null_edge();

			bool stop = false;
			while (!stop) {
				const auto &rr_node_p = get_vertex_props(_g, rr_node);

				char buffer[256];
				sprintf_rr_node(rr_node, buffer);
				zlog_level(delta_log, ROUTER_V3, "Adding %s to route tree\n", buffer);

				RouteTreeNode rt_node = route_tree_add_rr_node(rt, rr_node);
				if (rt_node != RouteTree::null_vertex()) {
					assert(inside_bbs(rr_node_p, *m_bounding_boxes));

					assert(_state[rr_node].upstream_R != std::numeric_limits<float>::max() && _state[rr_node].delay != std::numeric_limits<float>::max());
					route_tree_set_node_properties(rt, rt_node, rr_node_p.type != IPIN && rr_node_p.type != SINK, _state[rr_node].upstream_R, _state[rr_node].delay);

					if (rr_node_p.type == SOURCE) {
						route_tree_add_root(rt, rr_node);
					}

					added_rr_nodes.push_back(rr_node);
				} else {
					stop = true;
				}

				if (rr_node_p.type == OPIN && _existing_opin == RRGraph::null_vertex()) {
					_existing_opin = rr_node;
				}

				if (child_rr_node != RRGraph::null_vertex()) {
					char buffer2[256];
					sprintf_rr_node(child_rr_node, buffer2);

					zlog_level(delta_log, ROUTER_V3, "Adding edge %s -> %s\n", buffer, buffer2);

					const RouteTreeEdge &rt_edge = route_tree_add_edge_between_rr_node(rt, rr_node, child_rr_node); 
					auto &rt_edge_props = get_edge_props(rt.graph, rt_edge);
					assert(valid(edge));
					assert(get_source(_g, edge) == rr_node && get_target(_g, edge) == child_rr_node);
					rt_edge_props.rr_edge = edge;
				}

				zlog_level(delta_log, ROUTER_V3, "\n");

				child_rr_node = rr_node;
				edge = _prev_edge[rr_node];
				if (valid(edge)) {
					rr_node = get_source(_g, edge);
				} else {
					stop = true;
				}
			}
		}

		void reset_stats()
		{
			_stats.num_heap_pops = 0;
			_stats.num_heap_pushes = 0;
			_stats.num_neighbor_visits = 0;
		}

		//template<typename Graph, typename Edge, typename EdgeWeightFunc, typename Callbacks>
		//friend void delta_stepping(const Graph &g, const vector<heap_node_t<Edge>> &sources, int sink, float delta, float *known_distance, float *distance, Edge *prev_edge, const EdgeWeightFunc &edge_weight, Callbacks &callbacks);

		template<typename Graph, typename Edge, typename EdgeWeightFunc, typename ExpandCheckFunc, typename Callbacks, typename Extra>
		friend void dijkstra(const Graph &g, priority_queue<heap_node_t<Edge, Extra>> &sources, int sink, float *known_distance, float *distance, Edge *prev_edge, const ExpandCheckFunc &expand_node, const EdgeWeightFunc &edge_weight, Callbacks &callbacks);

		//template<typename Edge, typename Callbacks>
		//friend void relax(Buckets &buckets, float delta, vector<bool> &in_bucket, const vector<bool> &vertex_deleted,
					//float *known_distance, float *distance, Edge *predecessor,
					//int v, float new_known_distance, float new_distance, const Edge &edge,
					//Callbacks &callbacks);

	public:
		SinkRouter(const RRGraph &g, congestion_t *congestion, const float &pres_fac)
			: _g(g), _congestion(congestion), _pres_fac(pres_fac), _modified_node_added(num_vertices(g), false) 
		{
			_known_distance = new float[num_vertices(g)];
			_distance = new float[num_vertices(g)];
			_prev_edge = new RREdge[num_vertices(g)];
			_state = new extra_route_state_t[num_vertices(g)];

			for (int i = 0; i < num_vertices(g); ++i) {
				_known_distance[i] = std::numeric_limits<float>::max();
				_distance[i] = std::numeric_limits<float>::max();
				_prev_edge[i] = RRGraph::null_edge();
				_state[i].delay = std::numeric_limits<float>::max();
				_state[i].upstream_R = std::numeric_limits<float>::max();
			}

			reset_stats();
		}

		void route(const source_t *source, const sink_t *sink, const vector<box> &bounding_boxes, float astar_fac, route_tree_t &rt, vector<RRNode> &added_rr_nodes, float &delay, dijkstra_stats_t *stats)
		{
			_current_sink = sink;
			m_bounding_boxes = &bounding_boxes;
			_astar_fac = astar_fac;

			/* TODO: we should not reset the exisiting OPIN here */
			_existing_opin = RRGraph::null_vertex();

			reset_stats();

			char buffer[256];
			sprintf_rr_node(sink->rr_node, buffer);
			//zlog_level(delta_log, ROUTER_V3, "Current sink: %s BB %d->%d, %d->%d\n", buffer, sink->current_bounding_box.xmin, sink->current_bounding_box.xmax, sink->current_bounding_box.ymin, sink->current_bounding_box.ymax);
			const auto &sink_p = get_vertex_props(_g, sink->rr_node);
			//assert(sink_p.real_xlow == sink_p.real_xhigh && sink_p.real_ylow == sink_p.real_yhigh);
			zlog_level(delta_log, ROUTER_V3, "Current sink: %s [real %d %d] Source (%d, %d) BBs ",
					buffer, sink_p.real_xlow, sink_p.real_ylow, source ? source->x : -1, source ? source->y : -1);
#ifdef PRINT_BBS 
			for (int i = 0; i < bounding_boxes.size(); ++i) {
				zlog_level(delta_log, ROUTER_V3, "%d %d, %d %d -- ",
						bg::get<bg::min_corner, 0>(bounding_boxes[i]), bg::get<bg::max_corner, 0>(bounding_boxes[i]),
						bg::get<bg::min_corner, 1>(bounding_boxes[i]), bg::get<bg::max_corner, 1>(bounding_boxes[i]));
			}
#endif
			zlog_level(delta_log, ROUTER_V3, "\n");

			std::priority_queue<heap_node_t<RREdge, extra_route_state_t>> sources;

			if (source) {
				const auto &source_p = get_vertex_props(_g, source->rr_node);
				float delay = 0.5 * source_p.R * source_p.C;
				float kd = sink->criticality_fac * delay;
				int same, ortho;
				float d = kd + _astar_fac * get_timing_driven_expected_cost(source_p, sink_p, sink->criticality_fac, source_p.R, &same, &ortho); 

				sources.push({ source->rr_node, kd, d, RRGraph::null_edge(), { same, ortho, false, source_p.R, delay } });
				
				sprintf_rr_node(source->rr_node, buffer);
				zlog_level(delta_log, ROUTER_V3, "Using %s as initial heap\n", buffer);
			} else {
				for (const auto &rt_node : route_tree_get_nodes(rt)) {
					const auto &rt_node_p = get_vertex_props(rt.graph, rt_node);
					RRNode node = rt_node_p.rr_node;
					const auto &node_p = get_vertex_props(_g, node);

					sprintf_rr_node(node, buffer);

					if (rt_node_p.reexpand) {
						if (!inside_bbs(node_p, *m_bounding_boxes)) {
							zlog_level(delta_log, ROUTER_V3, "Existing %s out of bounding box\n", buffer);
						} else {
							int same, ortho;

							float kd = sink->criticality_fac * rt_node_p.delay; 
							float d = kd + _astar_fac * get_timing_driven_expected_cost(node_p, sink_p, sink->criticality_fac, rt_node_p.upstream_R, &same, &ortho); 

							RREdge prev = RRGraph::null_edge();
							if (valid(rt_node_p.rt_edge_to_parent)) {
								prev = get_edge_props(rt.graph, rt_node_p.rt_edge_to_parent).rr_edge;
							}

							zlog_level(delta_log, ROUTER_V3, "Adding %s back to heap [same %d ortho %d] [delay %g upstream_R %g] [kd: %g d: %g prev: %d]\n", buffer, same, ortho, rt_node_p.delay, rt_node_p.upstream_R, kd, d, get_source(_g, prev));

							//sources.push_back({ node, kd, d, prev, { same, ortho, true, rt_node_p.upstream_R, rt_node_p.delay } });
							sources.push({ node, kd, d, prev, { same, ortho, true, rt_node_p.upstream_R, rt_node_p.delay } });
						}
					} else {
						zlog_level(delta_log, ROUTER_V3, "Not reexpanding %s\n", buffer);
					}
				}
			}

			//delta_stepping(_g, sources, sink->rr_node, delta, _known_distance, _distance, _prev_edge, [this] (const RREdge &e) -> pair<float, float> { return get_edge_weight(e); }, *this);

			dijkstra(_g, sources, sink->rr_node, _known_distance, _distance, _prev_edge, [this] (int rr_node) -> bool { return expand_node(rr_node); }, [this] (const RREdge &e, float &known_distance, float &distance, extra_route_state_t &extra) -> void { return get_edge_weight(e, known_distance, distance, extra); }, *this);

			backtrack(sink->rr_node, rt, added_rr_nodes);

			assert(get_vertex_props(rt.graph, route_tree_get_rt_node(rt, sink->rr_node)).delay == _state[sink->rr_node].delay);
			//net_timing.delay[sink->id+1] = _state[sink->rr_node].delay;
			delay = _state[sink->rr_node].delay;

			if (stats) {
				*stats = _stats;
			}

			for (const auto &n : _modified_nodes) {
				_known_distance[n] = std::numeric_limits<float>::max();
				_distance[n] = std::numeric_limits<float>::max();
				_prev_edge[n] = RRGraph::null_edge();
				_state[n].delay = std::numeric_limits<float>::max();
				_state[n].upstream_R = std::numeric_limits<float>::max();

				assert(_modified_node_added[n]);
				_modified_node_added[n] = false;
			}

			_modified_nodes.clear();

			//for (const auto &n : get_vertices(_g)) {
			//assert(_known_distance[n] == std::numeric_limits<float>::max() &&
			//_distance[n] == std::numeric_limits<float>::max() &&
			//_prev_edge[n] == RRGraph::null_edge() &&
			//_state[n].delay == std::numeric_limits<float>::max() &&
			//_state[n].upstream_R == std::numeric_limits<float>::max());
			//}
		}
};

class Manager {
	private:
		std::queue<SinkRouter *> _free_routers;
		tbb::spin_mutex _lock;
		int _num_allocs;

		RRGraph &_g;
		congestion_t *_congestion;
		const float &_pres_fac;

	public:
		Manager(RRGraph &g, congestion_t *congestion, const float &pres_fac)
			: _num_allocs(0), _g(g), _congestion(congestion), _pres_fac(pres_fac) 
		{
		}

		SinkRouter *get()
		{
			SinkRouter *router;

			_lock.lock();

			if (_free_routers.empty()) {
				router = new SinkRouter(_g, _congestion, _pres_fac);
				++_num_allocs;
			} else {
				router = _free_routers.front();
				_free_routers.pop();
			}

			_lock.unlock();

			return router;
		}

		void reserve(int size)
		{
			int num_to_alloc = size-_free_routers.size();

			for (int i = 0; i < num_to_alloc; ++i) {
				_free_routers.push(new SinkRouter(_g, _congestion, _pres_fac));
				++_num_allocs;
			}
		}

		void free(SinkRouter *router)
		{
			_lock.lock();

			_free_routers.push(router);

			_lock.unlock();
		}
		
		int get_num_free() const
		{
			return _free_routers.size();
		}
		
		int get_num_allocs() const
		{
			return _num_allocs;
		}
};

void merge(const vector<const sink_t *> &sinks, const vector<route_tree_t> &in, route_tree_t &out, vector<RRNode> *added_rr_nodes)
{
	assert(sinks.size() == in.size());

	for (int i = 0; i < sinks.size(); ++i) {
		const auto &rt = in[i];

		zlog_level(delta_log, ROUTER_V3, "Route tree for sink %d\n", sinks[i]->id);

		RouteTreeNode rt_node = route_tree_get_rt_node(rt, sinks[i]->rr_node);
		assert(rt_node != RouteTree::null_vertex());

		RouteTreeNode child_rt_node = RouteTree::null_vertex();
		bool stop = false;

		while (!stop) {
			const auto &rt_node_p = get_vertex_props(rt.graph, rt_node);

			auto new_rt_node = route_tree_add_rr_node(out, rt_node_p.rr_node);

			if (new_rt_node != RouteTree::null_vertex()) {
				if (added_rr_nodes) {
					added_rr_nodes->push_back(rt_node_p.rr_node);
				}

				//update_first_order_congestion(_congestion, rt_node_p.rr_node, 1, get_vertex_props(_g, rt_node_p.rr_node).capacity, _pres_fac);

				const auto &rr_node_p = get_vertex_props(*rt.rrg, rt_node_p.rr_node);

				route_tree_set_node_properties(out, new_rt_node, rt_node_p.reexpand, rt_node_p.upstream_R, rt_node_p.delay);

				if (rr_node_p.type == SOURCE) {
					route_tree_add_root(out, rt_node_p.rr_node);
				}
			} else {
				//new_rt_node = route_tree_get_rt_node(out, rt_node_p.rr_node);
				stop = true;
			}

			if (child_rt_node != RouteTree::null_vertex()) {
				const auto &child_rt_node_p = get_vertex_props(rt.graph, child_rt_node);

				char s_from[256];
				char s_to[256];
				sprintf_rr_node(rt_node_p.rr_node, s_from);
				sprintf_rr_node(child_rt_node_p.rr_node, s_to);
				zlog_level(delta_log, ROUTER_V3, "\tAdding edge from %s to %s\n", s_from, s_to);

				route_tree_add_edge_between_rr_node(out, rt_node_p.rr_node, child_rt_node_p.rr_node);
			}

			child_rt_node = rt_node;

			if (valid(rt_node_p.rt_edge_to_parent)) {
				rt_node = get_source(rt.graph, rt_node_p.rt_edge_to_parent);
			} else {
				stop = true;
			}
		}
	}
	//for (const auto &rt : sinks) {
	//zlog_level(delta_log, ROUTER_V3, "Route tree\n");

	//for (const auto &rt_node : route_tree_get_nodes(rt)) {
	//const auto &rt_node_p = get_vertex_props(rt.graph, rt_node);

	//auto new_rt_node = route_tree_add_rr_node(out, rt_node_p.rr_node);

	//if (new_rt_node != RouteTree::null_vertex()) {
	//auto &new_rt_node_p = get_vertex_props(out.graph, new_rt_node);

	//new_rt_node_p.delay = rt_node_p.delay;
	//new_rt_node_p.upstream_R = rt_node_p.upstream_R;
	//}
	//}

	//for (const auto &rt_node : route_tree_get_nodes(rt)) {
	//for (const auto &e : route_tree_get_branches(rt, rt_node)) {
	//int from = get_vertex_props(rt.graph, get_source(rt.graph, e)).rr_node;
	//int to = get_vertex_props(rt.graph, get_target(rt.graph, e)).rr_node;

	//if (!route_tree_has_edge(out, from, to)) {
	//char s_from[256];
	//char s_to[256];
	//sprintf_rr_node(from, s_from);
	//sprintf_rr_node(to, s_to);
	//zlog_level(delta_log, ROUTER_V3, "\tAdding edge from %s to %s\n", s_from, s_to);

	//route_tree_add_edge_between_rr_node(out, from, to);
	//}
	//}
	//}
	//}
}

class MultiSinkParallelRouter {
	private:
		const RRGraph &_g;
		congestion_t *_congestion;
		const float &_pres_fac;

	public:
		MultiSinkParallelRouter(const RRGraph &g, congestion_t *congestion, const float &pres_fac)
			: _g(g), _congestion(congestion), _pres_fac(pres_fac)
		{
		}

		void route(const source_t *source, const vector<const sink_t *> &sinks, const vector<box> &bounding_boxes, float astar_fac, route_tree_t &rt, vector<RRNode> &added_rr_nodes, t_net_timing &net_timing, dijkstra_stats_t &stats)
		{
#if 1
			if (sinks.size() > 16) {
				vector<route_tree_t> local_route_trees(sinks.size());
				for (auto &tmp : local_route_trees) {
					route_tree_init(tmp, &_g);
				}

				//tbb::parallel_do(std::begin(sinks), std::end(sinks),
				//[&] (const sink_t *sink) -> void
				//{
				//SinkRouter *router = g_manager->get();
				//vector<int> added_rr_nodes; 

				//router->route(source, sink, astar_fac, local_route_trees[sink->id], added_rr_nodes, net_timing.delay[sink->id+1]);

				//g_manager->free(router);
				//});

				tbb::parallel_for(tbb::blocked_range<int>(0, sinks.size()), 
						[&source, &sinks, &bounding_boxes, &astar_fac, &local_route_trees, &net_timing, &stats]
						(const tbb::blocked_range<int> &r) -> void
						{
						thread_local SinkRouter *router = nullptr;
						if (!router) {
						router = g_manager->get();
						}

						for (int i = r.begin(); i != r.end(); ++i) {
						//SinkRouter *router = g_manager->get();
						vector<int> local_added_rr_nodes; 
						dijkstra_stats_t local_stats;

						router->route(source, sinks[i], bounding_boxes, astar_fac, local_route_trees[i], local_added_rr_nodes, net_timing.delay[sinks[i]->id+1], &local_stats);

						stats.num_heap_pops += local_stats.num_heap_pops;
						stats.num_heap_pushes += local_stats.num_heap_pushes;
						stats.num_neighbor_visits += local_stats.num_neighbor_visits;

						//g_manager->free(router);
						}
						}, tbb::auto_partitioner());

				merge(sinks, local_route_trees, rt, &added_rr_nodes);
			} else {
				thread_local SinkRouter *router = nullptr;
				if (!router) {
					router = g_manager->get();
				}

				for (int i = 0; i < sinks.size(); ++i) {
					dijkstra_stats_t local_stats;

					router->route(i == 0 ? source : nullptr, sinks[i], bounding_boxes, astar_fac, rt, added_rr_nodes, net_timing.delay[sinks[i]->id+1], &local_stats);

					stats.num_heap_pops += local_stats.num_heap_pops;
					stats.num_heap_pushes += local_stats.num_heap_pushes;
					stats.num_neighbor_visits += local_stats.num_neighbor_visits;
				}

				//g_manager->free(router);
#endif
				for (const auto &rr : added_rr_nodes) {
					update_first_order_congestion(_congestion, rr, 1, get_vertex_props(_g, rr).capacity, _pres_fac);
				}

				assert(rt.root_rt_nodes.size() == 1);
			}
		}
};

//struct vec_clock_t {
	//vector<int> counts;
//};
//

/*
static void do_work(router_t<net_t *> *router, MultiSinkParallelRouter &net_router)
{
#ifdef PTHREAD_BARRIER
	pthread_barrier_wait(&router->sync.barrier);
#else
	router->sync.barrier->wait();
#endif

	assert(net_router.get_num_instances() == router->num_threads);

	int tid = net_router.get_instance();

	tl_tid = tid;

	char buffer[256];

	sprintf(buffer, "%d", router->iter.load());
	zlog_put_mdc("iter", buffer);

	vector<timer::duration> net_route_time(router->nets.size(), timer::duration::max());
	vector<dijkstra_stats_t> net_stats(router->nets.size(),
			{ std::numeric_limits<int>::max(), std::numeric_limits<int>::max(), std::numeric_limits<int>::max() });
	vector<int> routed_nets;

	auto pure_route_time = timer::duration::zero();
	auto wait_time = timer::duration::zero();

	auto route_start = timer::now();

//#define ROUND_ROBIN

	//int last_inet = -1;
#ifdef ROUND_ROBIN
	int inet = tid;
#endif

#ifndef ROUND_ROBIN
	while ((*router->e_state[tid].inet = router->sync.net_index++) < router->nets.size()) {
		int inet = *router->e_state[tid].inet;
#else
	while (inet < router->nets.size()) {
		*router->e_state[tid].inet = inet;
#endif
		net_t &net = *router->nets[inet];

		//last_inet = inet;

		//router->sync.lock.lock();

		//assert(!router->sync.global_pending_nets.empty());
		//vnet = router->sync.global_pending_nets.front();
		//router->sync.global_pending_nets.pop();

		//router->sync.lock.unlock();

		auto real_net_route_start = timer::now();

		zlog_level(delta_log, ROUTER_V3, "Routing inet %d net %d num sinks %d\n", inet, net.local_id, net.sinks.size());
		//if (router->iter > 0) {
			//printf("Routing inet %d net %d num sinks %d\n", inet, net.local_id, net.sinks.size());
		//}

		if (router->iter > 0) {
			//for (const auto &rt_node : route_tree_get_nodes(router->state.back_route_trees[net.local_id])) {
			//RRNode rr_node = get_vertex_props(router->state.back_route_trees[net.local_id].graph, rt_node).rr_node;
			//const auto &rr_node_p = get_vertex_props(router->g, rr_node);

			//if (rr_node_p.xhigh >= _current_sink->current_bounding_box.xmin
			//&& rr_node_p.xlow <= _current_sink->current_bounding_box.xmax
			//&& rr_node_p.yhigh >= _current_sink->current_bounding_box.ymin
			//&& rr_node_p.ylow <= _current_sink->current_bounding_box.ymax) {
			//update_one_cost_internal(rr_node, router->g, router->state.congestion, -1, router->pres_fac);
			//}
			//}

			vector<Mods> mods(router->fpga_regions.size(), nullptr);

			for (const auto &rr_node : router->state.back_added_rr_nodes[net.local_id]) {
				const auto &rr_node_p = get_vertex_props(router->g, rr_node);
				//assert(rr_node_p.xlow >= vnet->sinks[0]->current_bounding_box.xmin
						//&& rr_node_p.xhigh <= vnet->sinks[0]->current_bounding_box.xmax
						//&& rr_node_p.ylow >= vnet->sinks[0]->current_bounding_box.ymin
						//&& rr_node_p.yhigh <= vnet->sinks[0]->current_bounding_box.ymax);
				//assert(rr_node_p.xhigh >= vnet->sinks[0]->current_bounding_box.xmin
						//&& rr_node_p.xlow <= vnet->sinks[0]->current_bounding_box.xmax
						//&& rr_node_p.yhigh >= vnet->sinks[0]->current_bounding_box.ymin
						//&& rr_node_p.ylow <= vnet->sinks[0]->current_bounding_box.ymax);
				//assert(!outside_bb(rr_node_p, net->bounding_boxes));
				update_first_order_congestion(router->state.congestion[tid], rr_node, -1, rr_node_p.capacity, router->pres_fac);
				//commit(router->regions[inet], router->g, rr_node, 1, tid);
				mods_add(rr_node, -1, router->rr_regions, mods);
			}

			mods_commit(mods, tid, router->fpga_regions);
		}

		vector<const sink_t *> sinks;
		for (int i = 0; i < net.sinks.size(); ++i) {
			sinks.push_back(&net.sinks[i]);
		}

		source_t *source = &net.source;

		net_router.reset_stats();

		net_router.route(source, sinks, router->params.astar_fac, router->state.route_trees[net.local_id], router->state.added_rr_nodes[net.local_id], router->state.net_timing[net.vpr_id]); 

		net_route_time[net.local_id] = timer::now()-real_net_route_start;
		pure_route_time += net_route_time[net.local_id];

		net_stats[net.local_id] = net_router.get_stats();

		routed_nets.push_back(net.local_id);

#ifdef ROUND_ROBIN
		inet += router->num_threads;
#endif
	}

	auto last_barrier_wait_start = timer::now();
#ifdef PTHREAD_BARRIER
	pthread_barrier_wait(&router->sync.barrier);
#else
	router->sync.barrier->wait(
			[&router, &tid] () -> void {
				++(*router->e_state[tid].logical_clock);
			});
#endif
	router->time[tid].last_barrier_wait_time = timer::now() - last_barrier_wait_start;

	zlog_level(delta_log, ROUTER_V3, "Start last sync\n");

	//assert(last_inet != -1);
	auto last_sync_start = timer::now();
	for (auto &r : router->fpga_regions) {
		//region_get_mods(*r, tid, [&router, &tid] (const mod_t &m) -> void {
				//const auto &props = get_vertex_props(router->g, m.node);
				//update_first_order_congestion(router->state.congestion[tid], m.node, m.delta, props.capacity, router->pres_fac);
				//});
	}
	//sync(router->fpga_regions,
			//[] (int id, const region_t *o) -> bool { return true; },
			//[] (int id, const region_t *o) -> void {  },
			//tid, router->g, router->state.congestion[tid], router->pres_fac);
	router->time[tid].last_sync_time = timer::now()-last_sync_start;

	router->time[tid].route_time = timer::now()-route_start;
	router->time[tid].wait_time = wait_time;
	router->time[tid].pure_route_time = pure_route_time; 

	router->stats[tid].num_heap_pops = 0;
	router->stats[tid].num_heap_pushes = 0;
	router->stats[tid].num_neighbor_visits = 0;

	for (const auto &r : routed_nets) {
		router->net_route_time[r] = net_route_time[r];
		router->net_stats[r] = net_stats[r];
		router->net_router[r] = tid;

		router->stats[tid].num_heap_pops += net_stats[r].num_heap_pops;
		router->stats[tid].num_heap_pushes += net_stats[r].num_heap_pushes;
		router->stats[tid].num_neighbor_visits += net_stats[r].num_neighbor_visits;
	}

#ifdef PTHREAD_BARRIER
	pthread_barrier_wait(&router->sync.barrier);
#else
	router->sync.barrier->wait(
			[&router, &tid] () -> void {
				++(*router->e_state[tid].logical_clock);
			});

#endif
}
*/

class SchedRouterTask : public tbb::task {
	private:
		router_t *m_router;
		fpga_tree_node_t *m_node;
		int m_level;

		void route_vnet(const new_virtual_net_t &vnet)
		{
			const net_t &net = *vnet.net;

			zlog_level(delta_log, ROUTER_V3, "Routing vnet [local %d global %d] net [%d]\n", vnet.local_index, vnet.global_index, net.local_id);

			//last_inet = inet;

			//router->sync.lock.lock();

			//assert(!router->sync.global_pending_nets.empty());
			//vnet = router->sync.global_pending_nets.front();
			//router->sync.global_pending_nets.pop();

			//router->sync.lock.unlock();

			//auto real_net_route_start = timer::now();

			//if (router->iter > 0) {
			//printf("Routing inet %d net %d num sinks %d\n", inet, net.local_id, net.sinks.size());
			//}

			if (m_router->iter > 0) {
				//for (const auto &rt_node : route_tree_get_nodes(router->state.back_route_trees[net.local_id])) {
				//RRNode rr_node = get_vertex_props(router->state.back_route_trees[net.local_id].graph, rt_node).rr_node;
				//const auto &rr_node_p = get_vertex_props(router->g, rr_node);

				//if (rr_node_p.xhigh >= _current_sink->current_bounding_box.xmin
				//&& rr_node_p.xlow <= _current_sink->current_bounding_box.xmax
				//&& rr_node_p.yhigh >= _current_sink->current_bounding_box.ymin
				//&& rr_node_p.ylow <= _current_sink->current_bounding_box.ymax) {
				//update_one_cost_internal(rr_node, router->g, router->state.congestion, -1, router->pres_fac);
				//}
				//}

				assert(!m_router->state.back_added_rr_nodes[vnet.global_index].empty());

				for (const auto &rr_node : m_router->state.back_added_rr_nodes[vnet.global_index]) {
					const auto &rr_node_p = get_vertex_props(m_router->g, rr_node);
					//assert(rr_node_p.xlow >= vnet->sinks[0]->current_bounding_box.xmin
					//&& rr_node_p.xhigh <= vnet->sinks[0]->current_bounding_box.xmax
					//&& rr_node_p.ylow >= vnet->sinks[0]->current_bounding_box.ymin
					//&& rr_node_p.yhigh <= vnet->sinks[0]->current_bounding_box.ymax);
					//assert(rr_node_p.xhigh >= vnet->sinks[0]->current_bounding_box.xmin
					//&& rr_node_p.xlow <= vnet->sinks[0]->current_bounding_box.xmax
					//&& rr_node_p.yhigh >= vnet->sinks[0]->current_bounding_box.ymin
					//&& rr_node_p.ylow <= vnet->sinks[0]->current_bounding_box.ymax);
					//assert(!outside_bb(rr_node_p, net->bounding_boxes));
					update_first_order_congestion(m_router->state.congestion, rr_node, -1, rr_node_p.capacity, m_router->pres_fac);
					//commit(router->regions[inet], router->g, rr_node, 1, tid);
				}
			}

			//_net_router->reset_stats();

			assert(m_router->state.added_rr_nodes[vnet.global_index].empty());

			int num_tasks = (1 << m_level);

			if (vnet.sinks.size() > m_router->num_threads /* && m_router->num_threads / num_tasks > 1*/) {
				zlog_level(delta_log, ROUTER_V3, "Parallel sink\n");

				vector<route_tree_t> route_trees(vnet.sinks.size());
				for (auto &tmp : route_trees) {
					route_tree_init(tmp, &m_router->g);
				}
				vector<vector<int>> added_rr_nodes(vnet.sinks.size());

				//std::atomic<unsigned long> num_heap_pushes(0);
				//std::atomic<unsigned long> num_heap_pops(0);
				//std::atomic<unsigned long> num_neighbor_visits(0);

				/* this assertion is no longer true */
				//assert(route_tree_empty(m_router->state.route_trees[net.local_id]));

				tbb::parallel_for(tbb::blocked_range<int>(0, vnet.sinks.size()),
						[&] (const tbb::blocked_range<int> &r) -> void
						{
						thread_local SinkRouter *sink_router = nullptr;
						if (!sink_router) {
						sink_router = g_manager->get();
						}

						for (int i = r.begin(); i != r.end(); ++i) {
						const sink_t *sink = vnet.sinks[i];
						dijkstra_stats_t dstats;

//#ifdef PRINT_BBS 
						//for (int i = 0; i < vnet.bounding_boxes.size(); ++i) {
						//zlog_level(delta_log, ROUTER_V3, "%d %d, %d %d ",
								//bg::get<bg::min_corner, 0>(vnet.bounding_boxes[i]), bg::get<bg::max_corner, 0>(vnet.bounding_boxes[i]),
								//bg::get<bg::min_corner, 1>(vnet.bounding_boxes[i]), bg::get<bg::max_corner, 1>(vnet.bounding_boxes[i]));
						//}
						//zlog_level(delta_log, ROUTER_V3, "\n");
//#endif

						for (int j = 0; j < vnet.bounding_boxes[i].size(); ++j) {
							assert(bg::covered_by(vnet.bounding_boxes[i][j], m_node->bb));
						}

						sink_router->route(vnet.source, sink, vnet.bounding_boxes[i], m_router->params.astar_fac, route_trees[i], added_rr_nodes[i], m_router->state.net_timing[net.vpr_id].delay[sink->id+1], &dstats); 

						//num_heap_pushes += dstats.num_heap_pushes;
						//num_heap_pops += dstats.num_heap_pops;
						//num_neighbor_visits += dstats.num_neighbor_visits;

						m_router->sink_stats[vnet.global_index][i] = dstats;
						}

						}, tbb::auto_partitioner());

				//m_router->net_stats[vnet.global_index].num_heap_pushes += num_heap_pushes;
				//m_router->net_stats[vnet.global_index].num_heap_pops += num_heap_pops;
				//m_router->net_stats[vnet.global_index].num_neighbor_visits += num_neighbor_visits;

				//assert(m_router->locks[net.local_id].try_lock());
				merge(vnet.sinks, route_trees, m_router->state.route_trees[net.local_id], &m_router->state.added_rr_nodes[vnet.global_index]);
				//m_router->locks[net.local_id].unlock();
			} else {
				zlog_level(delta_log, ROUTER_V3, "Serial sink\n");

				thread_local SinkRouter *sink_router = nullptr;
				if (!sink_router) {
					sink_router = g_manager->get();
				}

				const source_t *source;
				//assert(m_router->locks[net.local_id].try_lock());
				if (route_tree_empty(m_router->state.route_trees[net.local_id])) {
					source = vnet.source;
				} else {
					source = nullptr;
				}
				//m_router->locks[net.local_id].unlock();

				for (int i = 0; i < vnet.sinks.size(); ++i) {
					const sink_t *sink = vnet.sinks[i];
					dijkstra_stats_t dstats;

					for (int j = 0; j < vnet.bounding_boxes[i].size(); ++j) {
						assert(bg::covered_by(vnet.bounding_boxes[i][j], m_node->bb));
					}

					//assert(m_router->locks[net.local_id].try_lock());
					sink_router->route(i == 0 ? source : nullptr, sink, vnet.bounding_boxes[i], m_router->params.astar_fac, m_router->state.route_trees[net.local_id], m_router->state.added_rr_nodes[vnet.global_index], m_router->state.net_timing[net.vpr_id].delay[sink->id+1], &dstats); 
					//m_router->locks[net.local_id].unlock();

					m_router->sink_stats[vnet.global_index][i] = dstats;

					//m_router->net_stats[vnet.global_index].num_heap_pushes += dstats.num_heap_pushes;
					//m_router->net_stats[vnet.global_index].num_heap_pops += dstats.num_heap_pops;
					//m_router->net_stats[vnet.global_index].num_neighbor_visits += dstats.num_neighbor_visits;
				}
			}

			for (const auto &node : m_router->state.added_rr_nodes[vnet.global_index]) {
				const auto &rr_node_p = get_vertex_props(m_router->g, node);
				update_first_order_congestion(m_router->state.congestion, node, 1, rr_node_p.capacity, m_router->pres_fac);
			}

			//auto real_net_route_end = timer::now();

			//net_route_time[net.local_id] = real_net_route_end-real_net_route_start;
			//acc_net_route_time[net.local_id] = real_net_route_end-route_start;
			//acc_num_locks[net.local_id] = get_wait_stats(level, tid).size();
			//pure_route_time += net_route_time[net.local_id];

			//net_stats[net.local_id] = _net_router->get_stats();

			//net_level[net.local_id] = level;

			//routed_nets[level].push_back(&net);
		}

	public:
		SchedRouterTask(router_t *router, fpga_tree_node_t *node, int level)
			: m_router(router), m_node(node), m_level(level)
		{
		}

		void route_nets()
		{
			if (false) {
				for (const auto &vnet : m_node->vnets) {
					route_vnet(vnet);
				}
			} else {
			//int num_tasks = (1 << m_level);
			//if (m_router->num_threads / num_tasks > 1) {
				tbb::parallel_do(begin(m_node->roots), end(m_node->roots),
						[&] (new_virtual_net_t *vnet, tbb::parallel_do_feeder<new_virtual_net_t *> &feeder) -> void
						{
							assert(!m_router->routed[vnet->global_index]);
							m_router->routed[vnet->global_index] = true;

							auto vnet_route_start = timer::now();	

							route_vnet(*vnet);

							auto update_start = timer::now();
							m_router->route_time[vnet->global_index] = update_start - vnet_route_start;

							for (const auto &v : m_node->net_g[vnet->local_index]) {
								//zlog_level(delta_log, ROUTER_V3, "Net %d neighbor %d incoming %d\n", vnet->local_index, v, m_node->current_num_incoming_edges[v]->load());

								assert(m_node->vnets[v].local_index == v);

								int after = --(*m_node->current_num_incoming_edges[v]);
								assert(after >= 0);
								if (after == 0) {
									feeder.add(&m_node->vnets[v]);
								}
							}

							m_router->update_time[vnet->global_index] = timer::now() - update_start;
						});
			//} else {
				//for (const auto &vnet : m_node->vnets) {
					//work(vnet);
				//}
			//}
			}
		}

		tbb::task *execute()
		{
			zlog_level(delta_log, ROUTER_V3, "Node BB %d %d, %d %d\n", 
					bg::get<bg::min_corner, 0>(m_node->bb), bg::get<bg::max_corner, 0>(m_node->bb),
					bg::get<bg::min_corner, 1>(m_node->bb), bg::get<bg::max_corner, 1>(m_node->bb));


			__itt_task_begin(domain, __itt_null, __itt_null, shMyTask);

			auto route_start = timer::now();

			//auto route_start = timer::now();

			route_nets();
			
			//time[level].route_time = timer::now()-route_start;
			//time[level].wait_time = wait_time;
			//time[level].get_net_wait_time = get_net_wait_time;
			//time[level].pure_route_time = pure_route_time; 
			m_node->route_time = timer::now() - route_start;
			
			__itt_task_end(domain);

			bool recycle_left = false;
			bool recycle_right = false;
			if (m_node->left || m_node->right) {
				tbb::task &c = *new(allocate_continuation()) tbb::empty_task;

				int num_refs = 0;
				if (m_node->left) {
					/* cannot assign here because we are still accessing m_node later */
					//m_node = m_node->left;
					recycle_as_child_of(c);
					recycle_left = true;
					++num_refs;
					zlog_level(delta_log, ROUTER_V3, "Recycled left BB %d %d, %d %d\n", 
							bg::get<bg::min_corner, 0>(m_node->left->bb), bg::get<bg::max_corner, 0>(m_node->left->bb),
							bg::get<bg::min_corner, 1>(m_node->left->bb), bg::get<bg::max_corner, 1>(m_node->left->bb));
				}

				tbb::task *t = nullptr;
				if (m_node->right) {
					if (!recycle_left) {
						//m_node = m_node->right;
						recycle_as_child_of(c);
						recycle_right = true;
						zlog_level(delta_log, ROUTER_V3, "Recycled right BB %d %d, %d %d\n", 
								bg::get<bg::min_corner, 0>(m_node->bb), bg::get<bg::max_corner, 0>(m_node->bb),
								bg::get<bg::min_corner, 1>(m_node->bb), bg::get<bg::max_corner, 1>(m_node->bb));
					} else {
						t = new(c.allocate_child()) SchedRouterTask(m_router, m_node->right, m_level+1);
						zlog_level(delta_log, ROUTER_V3, "Spawned right BB %d %d, %d %d\n", 
								bg::get<bg::min_corner, 0>(m_node->right->bb), bg::get<bg::max_corner, 0>(m_node->right->bb),
								bg::get<bg::min_corner, 1>(m_node->right->bb), bg::get<bg::max_corner, 1>(m_node->right->bb));
					}
					++num_refs;
				}

				c.set_ref_count(num_refs);

				if (t) {
					spawn(*t);
				}

				zlog_level(delta_log, ROUTER_V3, "Finished task\n");

				assert(!(recycle_left && recycle_right));

				if (recycle_left) {
					m_node = m_node->left; 
					++m_level;
				} else if (recycle_right) {
					m_node = m_node->right;
					++m_level;
				}
			}

			return (recycle_left || recycle_right) ? this : nullptr;
		}
};

class MultiLayerSchedRouterTask : public tbb::task {
	private:
		router_t *m_router;
		fpga_tree_node_t *m_node;
		int m_level;

	public:
		MultiLayerSchedRouterTask(router_t *router, fpga_tree_node_t *node, int level)
			: m_router(router), m_node(node), m_level(level)
		{
		}


		tbb::task *execute()
		{
			zlog_level(delta_log, ROUTER_V3, "Node BB %d %d, %d %d\n", 
					bg::get<bg::min_corner, 0>(m_node->bb), bg::get<bg::max_corner, 0>(m_node->bb),
					bg::get<bg::min_corner, 1>(m_node->bb), bg::get<bg::max_corner, 1>(m_node->bb));


			__itt_task_begin(domain, __itt_null, __itt_null, shMyTask);

			auto route_start = timer::now();

			//auto route_start = timer::now();

			int num_tasks = (1 << m_level);

			if (m_node->next) {
				tbb::task &root = *new(tbb::task::allocate_root()) SchedRouterTask(m_router, m_node->next, m_level);

				tbb::task::spawn_root_and_wait(root);
			} else {
				SchedRouterTask task(m_router, m_node, m_level);
				task.route_nets();
			}
			
			//time[level].route_time = timer::now()-route_start;
			//time[level].wait_time = wait_time;
			//time[level].get_net_wait_time = get_net_wait_time;
			//time[level].pure_route_time = pure_route_time; 
			m_node->route_time = timer::now() - route_start;
			
			__itt_task_end(domain);

			bool recycle_left = false;
			bool recycle_right = false;
			if (m_node->left || m_node->right) {
				tbb::task &c = *new(allocate_continuation()) tbb::empty_task;

				int num_refs = 0;
				if (m_node->left) {
					/* cannot assign here because we are still accessing m_node later */
					//m_node = m_node->left;
					recycle_as_child_of(c);
					recycle_left = true;
					++num_refs;
					zlog_level(delta_log, ROUTER_V3, "Recycled left BB %d %d, %d %d\n", 
							bg::get<bg::min_corner, 0>(m_node->left->bb), bg::get<bg::max_corner, 0>(m_node->left->bb),
							bg::get<bg::min_corner, 1>(m_node->left->bb), bg::get<bg::max_corner, 1>(m_node->left->bb));
				}

				tbb::task *t = nullptr;
				if (m_node->right) {
					if (!recycle_left) {
						//m_node = m_node->right;
						recycle_as_child_of(c);
						recycle_right = true;
						zlog_level(delta_log, ROUTER_V3, "Recycled right BB %d %d, %d %d\n", 
								bg::get<bg::min_corner, 0>(m_node->bb), bg::get<bg::max_corner, 0>(m_node->bb),
								bg::get<bg::min_corner, 1>(m_node->bb), bg::get<bg::max_corner, 1>(m_node->bb));
					} else {
						t = new(c.allocate_child()) MultiLayerSchedRouterTask(m_router, m_node->right, m_level+1);
						zlog_level(delta_log, ROUTER_V3, "Spawned right BB %d %d, %d %d\n", 
								bg::get<bg::min_corner, 0>(m_node->right->bb), bg::get<bg::max_corner, 0>(m_node->right->bb),
								bg::get<bg::min_corner, 1>(m_node->right->bb), bg::get<bg::max_corner, 1>(m_node->right->bb));
					}
					++num_refs;
				}

				c.set_ref_count(num_refs);

				if (t) {
					spawn(*t);
				}

				zlog_level(delta_log, ROUTER_V3, "Finished task\n");

				assert(!(recycle_left && recycle_right));

				if (recycle_left) {
					m_node = m_node->left; 
					++m_level;
				} else if (recycle_right) {
					m_node = m_node->right;
					++m_level;
				}
			}

			return (recycle_left || recycle_right) ? this : nullptr;
		}
};

class VirtualRouterTask : public tbb::task {
	private:
		router_t *_router;
		fpga_tree_node_t *_node;
		int _level;

		void work(int color, int low, int high)
		{
			thread_local SinkRouter *sink_router = nullptr;
			if (!sink_router) {
				sink_router = g_manager->get();
			}

			for (int inet = low; inet != high; ++inet) {
				new_virtual_net_t &vnet = *_node->scheduled[color][inet];
				const net_t &net = *vnet.net;

				zlog_level(delta_log, ROUTER_V3, "Routing vnet %d net %d BB ", vnet.global_index, net.local_id);
#if 0 
				for (int i = 0; i < vnet.bounding_boxes.size(); ++i) {
					zlog_level(delta_log, ROUTER_V3, "%d %d, %d %d ",
							bg::get<bg::min_corner, 0>(vnet.bounding_boxes[i]), bg::get<bg::max_corner, 0>(vnet.bounding_boxes[i]),
							bg::get<bg::min_corner, 1>(vnet.bounding_boxes[i]), bg::get<bg::max_corner, 1>(vnet.bounding_boxes[i]));
				}
#endif
				zlog_level(delta_log, ROUTER_V3, "\n");

				//last_inet = inet;

				//router->sync.lock.lock();

				//assert(!router->sync.global_pending_nets.empty());
				//vnet = router->sync.global_pending_nets.front();
				//router->sync.global_pending_nets.pop();

				//router->sync.lock.unlock();

				auto real_net_route_start = timer::now();

				//if (router->iter > 0) {
				//printf("Routing inet %d net %d num sinks %d\n", inet, net.local_id, net.sinks.size());
				//}

				if (_router->iter > 0) {
					//for (const auto &rt_node : route_tree_get_nodes(router->state.back_route_trees[net.local_id])) {
					//RRNode rr_node = get_vertex_props(router->state.back_route_trees[net.local_id].graph, rt_node).rr_node;
					//const auto &rr_node_p = get_vertex_props(router->g, rr_node);

					//if (rr_node_p.xhigh >= _current_sink->current_bounding_box.xmin
					//&& rr_node_p.xlow <= _current_sink->current_bounding_box.xmax
					//&& rr_node_p.yhigh >= _current_sink->current_bounding_box.ymin
					//&& rr_node_p.ylow <= _current_sink->current_bounding_box.ymax) {
					//update_one_cost_internal(rr_node, router->g, router->state.congestion, -1, router->pres_fac);
					//}
					//}

					assert(!_router->state.back_added_rr_nodes[vnet.global_index].empty());

					for (const auto &rr_node : _router->state.back_added_rr_nodes[vnet.global_index]) {
						const auto &rr_node_p = get_vertex_props(_router->g, rr_node);
						//assert(rr_node_p.xlow >= vnet->sinks[0]->current_bounding_box.xmin
						//&& rr_node_p.xhigh <= vnet->sinks[0]->current_bounding_box.xmax
						//&& rr_node_p.ylow >= vnet->sinks[0]->current_bounding_box.ymin
						//&& rr_node_p.yhigh <= vnet->sinks[0]->current_bounding_box.ymax);
						//assert(rr_node_p.xhigh >= vnet->sinks[0]->current_bounding_box.xmin
						//&& rr_node_p.xlow <= vnet->sinks[0]->current_bounding_box.xmax
						//&& rr_node_p.yhigh >= vnet->sinks[0]->current_bounding_box.ymin
						//&& rr_node_p.ylow <= vnet->sinks[0]->current_bounding_box.ymax);
						//assert(!outside_bb(rr_node_p, net->bounding_boxes));
						update_first_order_congestion(_router->state.congestion, rr_node, -1, rr_node_p.capacity, _router->pres_fac);
						//commit(router->regions[inet], router->g, rr_node, 1, tid);
					}
				}

				const source_t *source;
				if (route_tree_empty(_router->state.route_trees[net.local_id])) {
					source = vnet.source;
				} else {
					source = nullptr;
				}

				//_net_router->reset_stats();

				dijkstra_stats_t lmao;

				for (int i = 0; i < vnet.sinks.size(); ++i) {
					const sink_t *sink = vnet.sinks[i];

					for (int j = 0; j < vnet.bounding_boxes.size(); ++j) {
						assert(bg::covered_by(vnet.bounding_boxes[i][j], _node->bb));
					}

					sink_router->route(source, sink, vnet.bounding_boxes[i], _router->params.astar_fac, _router->state.route_trees[net.local_id], _router->state.added_rr_nodes[vnet.global_index], _router->state.net_timing[net.vpr_id].delay[sink->id+1], &lmao); 
				}

				for (const auto &node : _router->state.added_rr_nodes[vnet.global_index]) {
					const auto &rr_node_p = get_vertex_props(_router->g, node);
					update_first_order_congestion(_router->state.congestion, node, 1, rr_node_p.capacity, _router->pres_fac);
				}

				auto real_net_route_end = timer::now();

				//net_route_time[net.local_id] = real_net_route_end-real_net_route_start;
				//acc_net_route_time[net.local_id] = real_net_route_end-route_start;
				//acc_num_locks[net.local_id] = get_wait_stats(level, tid).size();
				//pure_route_time += net_route_time[net.local_id];

				//net_stats[net.local_id] = _net_router->get_stats();

				//net_level[net.local_id] = level;

				//routed_nets[level].push_back(&net);
			}
		}

	public:
		VirtualRouterTask(router_t *router, fpga_tree_node_t *node, int level)
			: _router(router), _node(node), _level(level)
		{
		}


		tbb::task *execute()
		{
			int num_tasks = (2 << _level);

			zlog_level(delta_log, ROUTER_V3, "Node BB %d %d, %d %d Num tasks %d\n", 
					bg::get<bg::min_corner, 0>(_node->bb), bg::get<bg::max_corner, 0>(_node->bb),
					bg::get<bg::min_corner, 1>(_node->bb), bg::get<bg::max_corner, 1>(_node->bb),
					num_tasks);

			auto pure_route_time = timer::duration::zero();
			auto wait_time = timer::duration::zero();
			auto get_net_wait_time = timer::duration::zero();

			__itt_task_begin(domain, __itt_null, __itt_null, shMyTask);

			auto route_start = timer::now();

			for (int color = 0; color < _node->scheduled.size(); ++color) {
				zlog_level(delta_log, ROUTER_V3, "Color %d\n", color);

				if (_router->num_threads/num_tasks > 1 && _node->scheduled[color].size() > 16) {
					tbb::parallel_for(tbb::blocked_range<int>(0, _node->scheduled[color].size()),
							[&] (const tbb::blocked_range<int> &r) -> void
							{
								work(color, r.begin(), r.end());
							});
				} else {
					work(color, 0, _node->scheduled[color].size());
				}
			}	
			
			__itt_task_end(domain);

			//time[level].route_time = timer::now()-route_start;
			//time[level].wait_time = wait_time;
			//time[level].get_net_wait_time = get_net_wait_time;
			//time[level].pure_route_time = pure_route_time; 

			bool recycle_left = false;
			bool recycle_right = false;
			if (_node->left || _node->right) {
				tbb::task &c = *new(allocate_continuation()) tbb::empty_task;

				int num_refs = 0;
				if (_node->left) {
					/* cannot assign here because we are still accessing _node later */
					//_node = _node->left;
					recycle_as_child_of(c);
					recycle_left = true;
					++num_refs;
					zlog_level(delta_log, ROUTER_V3, "Recycled left BB %d %d, %d %d\n", 
							bg::get<bg::min_corner, 0>(_node->left->bb), bg::get<bg::max_corner, 0>(_node->left->bb),
							bg::get<bg::min_corner, 1>(_node->left->bb), bg::get<bg::max_corner, 1>(_node->left->bb));
				}

				if (_node->right) {
					if (!recycle_left) {
						//_node = _node->right;
						recycle_as_child_of(c);
						recycle_right = true;
						zlog_level(delta_log, ROUTER_V3, "Recycled right BB %d %d, %d %d\n", 
								bg::get<bg::min_corner, 0>(_node->bb), bg::get<bg::max_corner, 0>(_node->bb),
								bg::get<bg::min_corner, 1>(_node->bb), bg::get<bg::max_corner, 1>(_node->bb));
					} else {
						tbb::task &t = *new(c.allocate_child()) VirtualRouterTask(_router, _node->right, _level+1);
						spawn(t);
						zlog_level(delta_log, ROUTER_V3, "Spawned right BB %d %d, %d %d\n", 
								bg::get<bg::min_corner, 0>(_node->right->bb), bg::get<bg::max_corner, 0>(_node->right->bb),
								bg::get<bg::min_corner, 1>(_node->right->bb), bg::get<bg::max_corner, 1>(_node->right->bb));
					}
					++num_refs;
				}

				c.set_ref_count(num_refs);

				zlog_level(delta_log, ROUTER_V3, "Finished task\n");

				assert(!(recycle_left && recycle_right));

				if (recycle_left) {
					_node = _node->left; 
					++_level;
				} else if (recycle_right) {
					_node = _node->right;
					++_level;
				}
			}

			return (recycle_left || recycle_right) ? this : nullptr;
		}
};

class RouterTask : public tbb::task {
	private:
		router_t *_router;
		fpga_tree_node_t *_node;
		MultiSinkParallelRouter *_net_router;

	public:
		RouterTask(router_t *router, fpga_tree_node_t *node)
			: _router(router), _node(node), _net_router(new MultiSinkParallelRouter(router->g, router->state.congestion, router->pres_fac))
		{
		}

		tbb::task *execute()
		{
			zlog_level(delta_log, ROUTER_V3, "\tNode BB %d %d, %d %d\n", 
					bg::get<bg::min_corner, 0>(_node->bb), bg::get<bg::max_corner, 0>(_node->bb),
					bg::get<bg::min_corner, 1>(_node->bb), bg::get<bg::max_corner, 1>(_node->bb));

			auto pure_route_time = timer::duration::zero();
			auto wait_time = timer::duration::zero();
			auto get_net_wait_time = timer::duration::zero();

			__itt_task_begin(domain, __itt_null, __itt_null, shMyTask);

			auto route_start = timer::now();

			for (int inet = 0; inet < _node->nets.size(); ++inet) {
				net_t &net = *_node->nets[inet];

				assert(bg::covered_by(net.bounding_box, _node->bb));

				//last_inet = inet;

				//router->sync.lock.lock();

				//assert(!router->sync.global_pending_nets.empty());
				//vnet = router->sync.global_pending_nets.front();
				//router->sync.global_pending_nets.pop();

				//router->sync.lock.unlock();

				auto real_net_route_start = timer::now();

				zlog_level(delta_log, ROUTER_V3, "Routing inet %d net %d num sinks %lu\n", inet, net.local_id, net.sinks.size());
				//if (router->iter > 0) {
				//printf("Routing inet %d net %d num sinks %d\n", inet, net.local_id, net.sinks.size());
				//}

				if (_router->iter > 0) {
					//for (const auto &rt_node : route_tree_get_nodes(router->state.back_route_trees[net.local_id])) {
					//RRNode rr_node = get_vertex_props(router->state.back_route_trees[net.local_id].graph, rt_node).rr_node;
					//const auto &rr_node_p = get_vertex_props(router->g, rr_node);

					//if (rr_node_p.xhigh >= _current_sink->current_bounding_box.xmin
					//&& rr_node_p.xlow <= _current_sink->current_bounding_box.xmax
					//&& rr_node_p.yhigh >= _current_sink->current_bounding_box.ymin
					//&& rr_node_p.ylow <= _current_sink->current_bounding_box.ymax) {
					//update_one_cost_internal(rr_node, router->g, router->state.congestion, -1, router->pres_fac);
					//}
					//}

					assert(!_router->state.back_added_rr_nodes[net.local_id].empty());

					for (const auto &rr_node : _router->state.back_added_rr_nodes[net.local_id]) {
						const auto &rr_node_p = get_vertex_props(_router->g, rr_node);
						//assert(rr_node_p.xlow >= vnet->sinks[0]->current_bounding_box.xmin
						//&& rr_node_p.xhigh <= vnet->sinks[0]->current_bounding_box.xmax
						//&& rr_node_p.ylow >= vnet->sinks[0]->current_bounding_box.ymin
						//&& rr_node_p.yhigh <= vnet->sinks[0]->current_bounding_box.ymax);
						//assert(rr_node_p.xhigh >= vnet->sinks[0]->current_bounding_box.xmin
						//&& rr_node_p.xlow <= vnet->sinks[0]->current_bounding_box.xmax
						//&& rr_node_p.yhigh >= vnet->sinks[0]->current_bounding_box.ymin
						//&& rr_node_p.ylow <= vnet->sinks[0]->current_bounding_box.ymax);
						//assert(!outside_bb(rr_node_p, net->bounding_boxes));
						update_first_order_congestion(_router->state.congestion, rr_node, -1, rr_node_p.capacity, _router->pres_fac);
						//commit(router->regions[inet], router->g, rr_node, 1, tid);
					}
				}

				assert(route_tree_empty(_router->state.route_trees[net.local_id]));

				vector<const sink_t *> sinks;
				for (int i = 0; i < net.sinks.size(); ++i) {
					sinks.push_back(&net.sinks[i]);
				}

				source_t *source = &net.source;

				//_net_router->reset_stats();

				vector<box> boxes { net.bounding_box };
				dijkstra_stats_t lmao;
				_net_router->route(source, sinks, boxes, _router->params.astar_fac, _router->state.route_trees[net.local_id], _router->state.added_rr_nodes[net.local_id], _router->state.net_timing[net.vpr_id], lmao); 

				auto real_net_route_end = timer::now();

				//net_route_time[net.local_id] = real_net_route_end-real_net_route_start;
				//acc_net_route_time[net.local_id] = real_net_route_end-route_start;
				//acc_num_locks[net.local_id] = get_wait_stats(level, tid).size();
				//pure_route_time += net_route_time[net.local_id];

				//net_stats[net.local_id] = _net_router->get_stats();

				//net_level[net.local_id] = level;

				assert(&net == _node->nets[inet]);
				//routed_nets[level].push_back(&net);
			}

			__itt_task_end(domain);

			//time[level].route_time = timer::now()-route_start;
			//time[level].wait_time = wait_time;
			//time[level].get_net_wait_time = get_net_wait_time;
			//time[level].pure_route_time = pure_route_time; 

			bool recycle_left = false;
			bool recycle_right = false;
			if (_node->left || _node->right) {
				tbb::task &c = *new(allocate_continuation()) tbb::empty_task;

				int num_refs = 0;
				if (_node->left) {
					/* cannot assign here because we are still accessing _node later */
					//_node = _node->left;
					recycle_as_child_of(c);
					recycle_left = true;
					++num_refs;
					zlog_level(delta_log, ROUTER_V3, "Recycled left BB %d %d, %d %d\n", 
							bg::get<bg::min_corner, 0>(_node->left->bb), bg::get<bg::max_corner, 0>(_node->left->bb),
							bg::get<bg::min_corner, 1>(_node->left->bb), bg::get<bg::max_corner, 1>(_node->left->bb));
				}

				if (_node->right) {
					if (!recycle_left) {
						//_node = _node->right;
						recycle_as_child_of(c);
						recycle_right = true;
						zlog_level(delta_log, ROUTER_V3, "Recycled right BB %d %d, %d %d\n", 
								bg::get<bg::min_corner, 0>(_node->bb), bg::get<bg::max_corner, 0>(_node->bb),
								bg::get<bg::min_corner, 1>(_node->bb), bg::get<bg::max_corner, 1>(_node->bb));
					} else {
						tbb::task &t = *new(c.allocate_child()) RouterTask(_router, _node->right);
						spawn(t);
						zlog_level(delta_log, ROUTER_V3, "Spawned right BB %d %d, %d %d\n", 
								bg::get<bg::min_corner, 0>(_node->right->bb), bg::get<bg::max_corner, 0>(_node->right->bb),
								bg::get<bg::min_corner, 1>(_node->right->bb), bg::get<bg::max_corner, 1>(_node->right->bb));
					}
					++num_refs;
				}

				c.set_ref_count(num_refs);

				zlog_level(delta_log, ROUTER_V3, "Finished task\n");

				assert(!(recycle_left && recycle_right));

				if (recycle_left) {
					_node = _node->left; 
				} else if (recycle_right) {
					_node = _node->right;
				}
			}

			return (recycle_left || recycle_right) ? this : nullptr;
		}
};

void write_routes(const char *fname, const route_tree_t *route_trees, int num_nets)
{
	FILE *routes = fopen(fname, "w");

	for (int i = 0; i < num_nets; ++i) {
		fprintf(routes, "Net %d\n", i);
		for (const auto &node : route_tree_get_nodes(route_trees[i])) {
			const auto &props = get_vertex_props(route_trees[i].graph, node);
			fprintf(routes, "%d\n", props.rr_node);
		}
	}

	fclose(routes);
}

//void write_net_stats(const char *fname, const router_t *router)
//{
	//FILE *file = fopen(fname, "w");
	//assert(router->net_stats.size() == router->update_time.size());
	//assert(router->net_stats.size() == router->route_time.size());

	//for (int i = 0; i < router->net_stats.size(); ++i) {
		//const dijkstra_stats_t &stats = router->net_stats[i];
		//const new_virtual_net_t *vnet = router->all_vnets[i];
		//assert(router->all_vnets[i]->global_index == i);

		//fprintf(file, "%d %d %d %d %d\n", i, vnet->sinks.size(), stats.num_heap_pops, stats.num_neighbor_visits, stats.num_heap_pushes);
	//}
	
	//fclose(file);
//}

void write_overlap(const char *fname, const vector<vector<int>> &overlap)
{
	FILE *state = fopen(fname, "w");

	for (int i = 0; i < overlap.size(); ++i) {
		fprintf(state, "%d", i);
		for (int j = 0; j < overlap[i].size(); ++j) {
			fprintf(state, " %d", overlap[i][j]);
		}
		fprintf(state, "\n");
	}
	
	fclose(state);
}

template<typename Graph>
int longest_path_length(const Graph &g)
{
	vector<int> sorted;
	topological_sort(g, back_inserter(sorted));
	assert(num_vertices(g) == sorted.size());
	vector<int> length(sorted.size(), 0);
	int max_length = std::numeric_limits<int>::min();
	for (auto iter = sorted.rbegin(); iter != sorted.rend(); ++iter) {
		int u = *iter;
		auto edges = out_edges(u, g);
		for (auto ei = edges.first; ei != edges.second; ++ei) {
			int v = target(*ei, g);
			if (length[u] + 1 > length[v]) {
				length[v] = length[u]+1;
				max_length = std::max(max_length, length[v]);
			}
		}
	}
	return max_length;
}

void write_nets(const char *dirname, fpga_tree_node_t *node, const router_t *router, int level, int index, int layer)
{
	if (node) {
		if (node->next) {
			write_nets(dirname, node->next, router, level, index, layer+1);
		} else {
			char buffer[256];
			sprintf(buffer, "%s/net_stats_iter_%d_level_%d_node_%d.txt", dirname, router->iter.load(), level, index);
			FILE *file = fopen(buffer, "w");

			for (const auto &vnet : node->vnets) {
				const dijkstra_stats_t &stats = router->net_stats[vnet.global_index];

				//double bb_area = boxes_area(vnet.overlap_g_boxes);
				double bb_area = 0;

				fprintf(file, "%d %lu %lu %g %g %g %lu %lu %lu\n",
						vnet.global_index, vnet.net->sinks.size(), vnet.sinks.size(), bb_area, bb_area/bg::area(node->bb),
						duration_cast<nanoseconds>(router->route_time[vnet.global_index]).count()/1e9,
						stats.num_heap_pops, stats.num_neighbor_visits, stats.num_heap_pushes);
			}

			fclose(file);
		}

		write_nets(dirname, node->left, router, layer > 0 ? level : level+1, layer > 0 ? (index*2) - (1 << (level+1)) : index*2, layer);
		write_nets(dirname, node->right, router, layer > 0 ? level : level+1, layer > 0 ? (index*2+1) - (1 << (level+1)) : index*2+1, layer);
	}
}

void write_sink_stats(const char *dirname, fpga_tree_node_t *node, const router_t *router, int level, int index, int layer)
{
	if (node) {
		if (node->next) {
			write_sink_stats(dirname, node->next, router, level, index, layer+1);
		} else {
			char buffer[256];
			sprintf(buffer, "%s/sink_stats_iter_%d_level_%d_node_%d.txt", dirname, router->iter.load(), level, index);
			FILE *file = fopen(buffer, "w");

			for (const auto &vnet : node->vnets) {
				assert(vnet.sinks.size() == router->sink_stats[vnet.global_index].size());

				for (int i = 0; i < vnet.sinks.size(); ++i) {
					const auto &stats = router->sink_stats[vnet.global_index][i];
					fprintf(file, "%d %d %d %d %d %d %g %lu %lu %lu\n",
							vnet.global_index, i,
							vnet.source->x, vnet.source->y, vnet.sinks[i]->x, vnet.sinks[i]->y, boxes_area(vnet.bounding_boxes[i]),
							stats.num_heap_pops, stats.num_neighbor_visits, stats.num_heap_pushes);
				}
			}

			fclose(file);
		}

		write_sink_stats(dirname, node->left, router, layer > 0 ? level : level+1, layer > 0 ? (index*2) - (1 << (level+1)) : index*2, layer);
		write_sink_stats(dirname, node->right, router, layer > 0 ? level : level+1, layer > 0 ? (index*2+1) - (1 << (level+1)) : index*2+1, layer);
	}
}

void write_net_g(const char *dirname, fpga_tree_node_t *node, int level, int index)
{
	if (node) {
		char buffer[256];

		sprintf(buffer, "%s/net_g_level_%d_node_%d.txt", dirname, level, index);
		write_overlap(buffer, node->net_g);

		write_net_g(dirname, node->left, level+1, index*2);
		write_net_g(dirname, node->right, level+1, index*2+1);
	}
}

void write_congestion_state(const char *fname, const congestion_t *congestion, int num_rr_nodes)
{
	FILE *state = fopen(fname, "w");

	for (int i = 0; i < num_rr_nodes; ++i) {
		fprintf(state, "%d %d\n", i, congestion[i].occ);
	}
	
	fclose(state);
}

void print_tabs(int num)
{
	for (int i = 0; i < num; ++i) {
		printf("\t");
	}
}

void split_edges_lb(const std::vector<std::pair<int, net_t *>> &edges,
		int opt_cut, bool horizontal,
		std::vector<std::pair<int, net_t *>> &left, std::vector<std::pair<int, net_t *>> &right)
{
	for (const auto &edge : edges) {
		const auto &box = edge.second->bounding_box;

		int low;
		int high;

		if (horizontal) {
			low = bg::get<bg::min_corner, 1>(box);
			high = bg::get<bg::max_corner, 1>(box);
		} else {
			low = bg::get<bg::min_corner, 0>(box);
			high = bg::get<bg::max_corner, 0>(box);
		}

		assert(low < high);

		//printf("horizontal %d edge %d low %d high %d\n", horizontal, edge.first, low, high);
		//assert((edge.first == low && edge.first != high) || (edge.first != low && edge.first == high));

		if (high <= opt_cut) {
			left.push_back(edge);
		} else if (low > opt_cut) {
			right.push_back(edge);
		} else {
			assert(low <= opt_cut && high > opt_cut);
			//assert((low < opt_cut && high >= opt_cut) || );
		}
	}
}

template<typename Load, typename Edge>
void bar2(const std::vector<net_t *> &nets, const Edge &edge, const Load &load, vector<vector<net_t *>> &sorted_nets, vector<int> &num_nets, vector<float> &loads)
{
	//int cut = -1;
	//int total_num_nets = 0;
	//for (int i = 0; i < sorted.size(); ++i) {
		//if (cut == -1 || cut != sorted[i].first) {
			//cut = sorted[i].first;
			//assert(cut >= 0 && cut < num_nets.size());
			//num_nets[cut] = total_num_nets;
		//} else {
			//assert(cut == sorted[i].first);
		//} 

		//++num_nets[cut];
		//++total_num_nets;
	//}
	//assert(total_num_nets == sorted.size());
	assert(sorted_nets.size() == num_nets.size());
	assert(num_nets.size() == loads.size());

	for (int i = 0; i < nets.size(); ++i) {
		int cut = edge(nets[i]);

		assert(cut >= 0 && cut < sorted_nets.size());

		sorted_nets[cut].push_back(nets[i]);

		++num_nets[cut];
		loads[cut] += load(nets[i]);
	}
}

template<typename T>
void acc_left(vector<T> &items)
{
	T sum = items.front();
	for (int i = 1; i < items.size(); ++i) {
		items[i] += sum;
		sum = items[i];
	}
}

template<typename T>
void acc_right(vector<T> &items)
{
	T sum = items.back();
	for (int i = items.size()-2; i >= 0; --i) {
		items[i] += sum;
		sum = items[i];
	}
}

template<int Axis, typename Load>
void split_nets(const std::vector<net_t *> &nets, const Load &load,
		int num_all_nets, const pair<int, int> &bounds,
		int &opt_cut, float &min_diff,
		vector<net_t *> &left_nets, vector<net_t *> &right_nets, vector<net_t *> &middle_nets)
{
	assert(bounds.second > bounds.first);

	vector<vector<net_t *>> sorted_end(bounds.second+1);
	vector<int> num_left_nets(bounds.second+1, 0);
	vector<float> left_loads(bounds.second+1, 0);
	bar2(nets, [] (const net_t *net) -> int {
			return bg::get<bg::max_corner, Axis>(net->bounding_box);
			}, load,
			sorted_end, num_left_nets, left_loads);

	acc_left(num_left_nets);
	acc_left(left_loads);

	assert(num_left_nets.back() == nets.size());

	//printf("bounds: %d %d\n", bounds.first, bounds.second);

	//printf("all left: ");
	//for (int i = bounds.first; i <= bounds.second; ++i) {
		//if (i % 8 == 0) {
			//printf("\n");
		//}
		//printf("%d ", num_left_nets[i]);
	//}
	//printf("\n");

	vector<vector<net_t *>> sorted_start(bounds.second+1);
	vector<int> num_right_nets(bounds.second+1, 0);
	vector<float> right_loads(bounds.second+1, 0);
	bar2(nets, [] (const net_t *net) -> int {
			return bg::get<bg::min_corner, Axis>(net->bounding_box);
			}, load,
			sorted_start, num_right_nets, right_loads);

	acc_right(num_right_nets);
	acc_right(right_loads);

	assert(num_right_nets.front() == nets.size());

	//for (int i = num_grid_points-2; i >= 0; --i) {
		//if (num_right_nets[i] == 0) {
			//num_right_nets[i] = num_right_nets[i+1];
		//}
	//}

	//printf("all right: ");
	//for (int i = bounds.first; i <= bounds.second; ++i) {
		//if (i % 8 == 0) {
			//printf("\n");
		//}
		//printf("%d ", num_right_nets[i]);
	//}
	//printf("\n");

	/* checking */
	for (int cut = bounds.first; cut <= bounds.second; ++cut) {
		//printf("cut %d\n", cut);

		vector<bool> added(num_all_nets, false);

		int num_left = 0;
		//printf("num_left %d\n", num_left);
		for (int i = 0; i <= cut; ++i) {
			for (int j = 0; j < sorted_end[i].size(); ++j) {
				assert((bg::get<bg::min_corner, Axis>(sorted_end[i][j]->bounding_box) <= cut));
				assert((bg::get<bg::max_corner, Axis>(sorted_end[i][j]->bounding_box) <= cut));

				int id = sorted_end[i][j]->local_id;

				assert(!added[id]);
				added[id] = true;

				++num_left;

				//printf("left net %d, %d %d\n", id, bg::get<bg::min_corner, 0>(sorted_end[i].second->bounding_box),
				//bg::get<bg::max_corner, 0>(sorted_end[i].second->bounding_box));
			}
		}
		assert(num_left == num_left_nets[cut]);

		int num_right = 0;
		if (cut + 1 <= bounds.second) {
			for (int i = cut+1; i <= bounds.second; ++i) {
				for (int j = 0; j < sorted_start[i].size(); ++j) {
					assert((bg::get<bg::min_corner, Axis>(sorted_start[i][j]->bounding_box) > cut));
					assert((bg::get<bg::max_corner, Axis>(sorted_start[i][j]->bounding_box) > cut));

					int id = sorted_start[i][j]->local_id;

					assert(!added[id]);
					added[id] = true;

					++num_right;
					//printf("right net %d, %d %d\n", id, bg::get<bg::min_corner, 0>(sorted_start[i].second->bounding_box),
					//bg::get<bg::max_corner, 0>(sorted_start[i].second->bounding_box));
				}
			}
			assert(num_right == num_right_nets[cut+1]);
		}
		//printf("num_right %d\n", num_right);

		int num_cross_nets = 0;
		for (int i = 0; i < nets.size(); ++i) {
			int id = nets[i]->local_id;

			if (!added[id]) {
				//printf("cross net %d, %d %d\n", id, bg::get<bg::min_corner, 0>(sorted_start[i].second->bounding_box),
						//bg::get<bg::max_corner, 0>(sorted_start[i].second->bounding_box));

				assert((bg::get<bg::min_corner, Axis>(nets[i]->bounding_box) <= cut));
				assert((bg::get<bg::max_corner, Axis>(nets[i]->bounding_box) > cut));

				++num_cross_nets;
			}
		}

		//printf("cut %d num left/right/cross/total nets: %d/%d/%d/%lu\n",
				//cut,
				//num_left, num_right, num_cross_nets, nets.size());
		//printf("\tleft/right loads %g/%g\n",
				//left_loads[cut], cut+1 <= bounds.second ? right_loads[cut+1] : 0);
		assert(num_left + num_right + num_cross_nets == nets.size());
	}


	min_diff = std::numeric_limits<float>::max();
	opt_cut = -1;
	for (int cut = bounds.first; cut <= bounds.second; ++cut) {
		float diff;
		if (cut+1 <= bounds.second) {
			diff = abs(left_loads[cut] - right_loads[cut+1]);
		} else {
			diff = left_loads[cut];
		}
		assert(diff >= 0);
		if (diff < min_diff) {
			opt_cut = cut;
			min_diff = diff;
		}
	}
	assert(opt_cut != -1);
	assert(opt_cut + 1 <= bounds.second);

	vector<bool> added(num_all_nets, false);
	for (int i = 0; i <= opt_cut; ++i) {
		for (int j = 0; j < sorted_end[i].size(); ++j) {
			int id = sorted_end[i][j]->local_id;

			assert(!added[id]);
			added[id] = true;

			left_nets.push_back(sorted_end[i][j]);

			//printf("left net %d, %d %d\n", id, bg::get<bg::min_corner, 0>(sorted_end[i].second->bounding_box),
			//bg::get<bg::max_corner, 0>(sorted_end[i].second->bounding_box));
		}
	}
	assert(num_left_nets[opt_cut] == left_nets.size());

	for (int i = opt_cut+1; i <= bounds.second; ++i) {
		for (int j = 0; j < sorted_start[i].size(); ++j) {

			int id = sorted_start[i][j]->local_id;

			assert(!added[id]);
			added[id] = true;

			right_nets.push_back(sorted_start[i][j]);

			//printf("right net %d, %d %d\n", id, bg::get<bg::min_corner, 0>(sorted_start[i].second->bounding_box),
			//bg::get<bg::max_corner, 0>(sorted_start[i].second->bounding_box));
		}
	}
	assert(num_right_nets[opt_cut+1] == right_nets.size());

	for (int i = 0; i < nets.size(); ++i) {
		int id = nets[i]->local_id;

		if (!added[id]) {
			middle_nets.push_back(nets[i]);
		}
	}

	assert(left_nets.size() + right_nets.size() + middle_nets.size() == nets.size());
}

void bar(const std::vector<std::pair<int, net_t *>> &sorted, int cut, int type,
		int &cur, int &num_added, vector<net_t *> &nets, vector<bool> &added)
{
	num_added = 0;
	for (cur = 0; cur < sorted.size() && sorted[cur].first <= cut; ++cur) {
		int id = sorted[cur].second->local_id;

		if (type == 1) {
			printf("left net %d, %d %d\n", id, bg::get<bg::min_corner, 0>(sorted[cur].second->bounding_box),
					bg::get<bg::max_corner, 0>(sorted[cur].second->bounding_box));
		} else if (type == 2) {
			printf("cross net %d, %d %d\n", id, bg::get<bg::min_corner, 0>(sorted[cur].second->bounding_box),
					bg::get<bg::max_corner, 0>(sorted[cur].second->bounding_box));
		}

		assert(!added[id]);

		nets.push_back(sorted[cur].second);
		added[id] = true;

		++num_added;
	}
}

void bar_right(const std::vector<std::pair<int, net_t *>> &sorted, int cut, int type,
		int &cur, int &num_added, vector<net_t *> &nets, vector<bool> &added)
{
	num_added = 0;
	for (cur = 0; cur < sorted.size() && sorted[cur].first <= cut; ++cur) {
		int id = sorted[cur].second->local_id;

		if (type == 1) {
			printf("left net %d, %d %d\n", id, bg::get<bg::min_corner, 0>(sorted[cur].second->bounding_box),
					bg::get<bg::max_corner, 0>(sorted[cur].second->bounding_box));
		} else if (type == 2) {
			printf("cross net %d, %d %d\n", id, bg::get<bg::min_corner, 0>(sorted[cur].second->bounding_box),
					bg::get<bg::max_corner, 0>(sorted[cur].second->bounding_box));
		}

		assert(!added[id]);

		nets.push_back(sorted[cur].second);
		added[id] = true;

		++num_added;
	}
}

void foo(const std::vector<std::pair<int, net_t *>> &sorted_start, const std::vector<std::pair<int, net_t *>> &sorted_end,
		int num_all_nets,
		int num_grid_points
	   	)
{
	vector<int> num_left_nets(num_grid_points, 0);
	vector<int> num_right_nets(num_grid_points, 0);

	assert(sorted_start.size() == sorted_end.size());

	int cur_sorted_end = 0;
	int cur_sorted_start = 0;

	for (int cut = 0; cut < num_grid_points; ++cut) {
		vector<net_t *> left_nets;
		/* cross net can change back to left net so we need per 'cut' vector */
		vector<net_t *> right_nets;
		vector<bool> added(num_all_nets, false);

		printf("point %d\n", cut);

		bar(sorted_end, cut, 1,
				cur_sorted_end, num_left_nets[cut], left_nets, added);
		printf("num left nets %lu\n", left_nets.size());
		bar(sorted_start, cut, 2,
				cur_sorted_start, num_right_nets[cut], right_nets, added);
		printf("num right nets %lu\n", right_nets.size());

		for (const auto &net : left_nets) {
			assert((cut >= bg::get<bg::min_corner, 0>(net->bounding_box)));
			assert((cut >= bg::get<bg::max_corner, 0>(net->bounding_box)));
		}
		for (const auto &net : right_nets) {
			assert((cut < bg::get<bg::min_corner, 0>(net->bounding_box)));
			assert((cut < bg::get<bg::max_corner, 0>(net->bounding_box)));
		}
		int num_cross_nets = 0;
		for (int i = 0; i < sorted_start.size(); ++i) {
			int id = sorted_start[i].second->local_id;
			if (!added[id]) {
				printf("cross net %d, %d %d\n", id, bg::get<bg::min_corner, 0>(sorted_start[i].second->bounding_box),
						bg::get<bg::max_corner, 0>(sorted_start[i].second->bounding_box));

				assert((cut >= bg::get<bg::min_corner, 0>(sorted_start[i].second->bounding_box)));
				assert((cut <= bg::get<bg::max_corner, 0>(sorted_start[i].second->bounding_box)));

				++num_cross_nets;
			}
		}
		printf("num cross nets: %d\n", num_cross_nets);
		printf("num nets: %lu\n", sorted_start.size());
		assert(left_nets.size() + right_nets.size() + num_cross_nets == sorted_start.size());
	}
}

template<typename Load>
void split_nets_recursive(
		const std::vector<net_t *> &nets, const Load &load,
		int num_all_nets, const box &bb, int level, bool horizontal,
		fpga_tree_node_t *node)
{
	printf("level %d\n", level);

	print_tabs(level); printf("box %d %d %d %d hor %d\n",
			bg::get<bg::min_corner, 0>(bb), bg::get<bg::max_corner, 0>(bb),
			bg::get<bg::min_corner, 1>(bb), bg::get<bg::max_corner, 1>(bb),
			horizontal);

	assert(!nets.empty());

	if (level > 0) {
		//int prev = -1;
		//for (const auto &e : ver_start) {
			//assert(e.first >= prev);
			//prev = e.first;
		//}
		//prev = -1;
		//for (const auto &e : ver) {
			//assert(e.first >= prev);
			//prev = e.first;
		//}

		//for (const auto &ver_e : ver_edges) {
			//PRINT("ver edge %d box %X\n", ver_e.first, ver_e.second);
		//}
		
		int opt_cut;
		float min_diff;

		vector<net_t *> left_nets;
		vector<net_t *> right_nets;

		if (horizontal) {
			split_nets<1>(nets, load,
					num_all_nets, make_pair(bg::get<bg::min_corner, 1>(bb), bg::get<bg::max_corner, 1>(bb)),
					opt_cut, min_diff,
					left_nets, right_nets, node->nets);
		} else {
			split_nets<0>(nets, load,
					num_all_nets, make_pair(bg::get<bg::min_corner, 0>(bb), bg::get<bg::max_corner, 0>(bb)),
					opt_cut, min_diff,
					left_nets, right_nets, node->nets);
		}

		printf("opt_cut %d\n", opt_cut);
		
		box left_box, right_box;
		if (horizontal) {
			assert((opt_cut >= bg::get<bg::min_corner, 1>(bb) && opt_cut <= bg::get<bg::max_corner, 1>(bb)));

			bg::set<bg::min_corner, 0>(left_box, bg::get<bg::min_corner, 0>(bb));
			bg::set<bg::max_corner, 0>(left_box, bg::get<bg::max_corner, 0>(bb));
			bg::set<bg::min_corner, 1>(left_box, bg::get<bg::min_corner, 1>(bb));
			bg::set<bg::max_corner, 1>(left_box, opt_cut);

			bg::set<bg::min_corner, 0>(right_box, bg::get<bg::min_corner, 0>(bb));
			bg::set<bg::max_corner, 0>(right_box, bg::get<bg::max_corner, 0>(bb));
			bg::set<bg::min_corner, 1>(right_box, bg::get<bg::max_corner, 1>(left_box)+1);
			bg::set<bg::max_corner, 1>(right_box, bg::get<bg::max_corner, 1>(bb));

			//printf("hsplit 0 %d->%d,%d->%d 1 %d->%d,%d->%d\n",
			//bg::get<bg::min_corner, 0>(split_0), bg::get<bg::max_corner, 0>(split_0),
			//bg::get<bg::min_corner, 1>(split_0), bg::get<bg::max_corner, 1>(split_0),
			//bg::get<bg::min_corner, 0>(split_1), bg::get<bg::max_corner, 0>(split_1),
			//bg::get<bg::min_corner, 1>(split_1), bg::get<bg::max_corner, 1>(split_1));

			assert((bg::get<bg::min_corner, 1>(left_box) <= bg::get<bg::max_corner, 1>(left_box)));
			assert((bg::get<bg::min_corner, 1>(right_box) <= bg::get<bg::max_corner, 1>(right_box)));
		} else {
			assert((opt_cut >= bg::get<bg::min_corner, 0>(bb) && opt_cut <= bg::get<bg::max_corner, 0>(bb)));

			bg::set<bg::min_corner, 1>(left_box, bg::get<bg::min_corner, 1>(bb));
			bg::set<bg::max_corner, 1>(left_box, bg::get<bg::max_corner, 1>(bb));
			bg::set<bg::min_corner, 0>(left_box, bg::get<bg::min_corner, 0>(bb));
			bg::set<bg::max_corner, 0>(left_box, opt_cut);

			bg::set<bg::min_corner, 1>(right_box, bg::get<bg::min_corner, 1>(bb));
			bg::set<bg::max_corner, 1>(right_box, bg::get<bg::max_corner, 1>(bb));
			bg::set<bg::min_corner, 0>(right_box, bg::get<bg::max_corner, 0>(left_box)+1);
			bg::set<bg::max_corner, 0>(right_box, bg::get<bg::max_corner, 0>(bb));

			//printf("vsplit 0 %d->%d,%d->%d 1 %d->%d,%d->%d\n",
			//bg::get<bg::min_corner, 0>(split_0), bg::get<bg::max_corner, 0>(split_0),
			//bg::get<bg::min_corner, 1>(split_0), bg::get<bg::max_corner, 1>(split_0),
			//bg::get<bg::min_corner, 0>(split_1), bg::get<bg::max_corner, 0>(split_1),
			//bg::get<bg::min_corner, 1>(split_1), bg::get<bg::max_corner, 1>(split_1));

			assert((bg::get<bg::min_corner, 0>(left_box) <= bg::get<bg::max_corner, 0>(left_box)));
			assert((bg::get<bg::min_corner, 0>(right_box) <= bg::get<bg::max_corner, 0>(right_box)));
		}

		print_tabs(level); printf("\tleft %d %d %d %d\n",
				bg::get<bg::min_corner, 0>(left_box), bg::get<bg::max_corner, 0>(left_box),
				bg::get<bg::min_corner, 1>(left_box), bg::get<bg::max_corner, 1>(left_box));

		print_tabs(level); printf("\tright %d %d %d %d\n",
				bg::get<bg::min_corner, 0>(right_box), bg::get<bg::max_corner, 0>(right_box),
				bg::get<bg::min_corner, 1>(right_box), bg::get<bg::max_corner, 1>(right_box));

		//#ifdef DEBUG_MISR
		//#endif
		
		node->bb = bb;

		//std::vector<std::pair<int, net_t *>> left_ver_start;
		//std::vector<std::pair<int, net_t *>> left_ver_end;
		//std::vector<std::pair<int, net_t *>> right_ver_start;
		//std::vector<std::pair<int, net_t *>> right_ver_end;
		//split_edges_lb(ver_start, 
				//opt_cut, horizontal,
				//left_ver_start, right_ver_start);
		//split_edges_lb(ver_end, 
				//opt_cut, horizontal,
				//left_ver_end, right_ver_end);

		//std::vector<std::pair<int, net_t *>> left_hor_start;
		//std::vector<std::pair<int, net_t *>> left_hor_end;
		//std::vector<std::pair<int, net_t *>> right_hor_start;
		//std::vector<std::pair<int, net_t *>> right_hor_end;
		//split_edges_lb(hor_start, 
				//opt_cut, horizontal,
				//left_hor_start, right_hor_start);
		//split_edges_lb(hor_end, 
				//opt_cut, horizontal,
				//left_hor_end, right_hor_end);

		if (!left_nets.empty()) {
			node->left = new fpga_tree_node_t;

			split_nets_recursive(left_nets, load, 
					num_all_nets, left_box, level-1, !horizontal,
					node->left);
		} else {
			node->left = nullptr;
		}

		if (!right_nets.empty()) {
			node->right = new fpga_tree_node_t;

			split_nets_recursive(right_nets, load,
					num_all_nets, right_box, level-1, !horizontal,
					node->right);
		} else { 
			node->right = nullptr;
		}
	} else {
		printf("last level\n");

		node->bb = bb;
		for (const auto &item : nets) {
			assert(bg::covered_by(item->bounding_box, bb));
			node->nets.push_back(item);
		}
		node->left = nullptr;
		node->right = nullptr;
	}
}

template<typename Load>
void split_nets_recursive_multi_layer(
		const std::vector<net_t *> &nets, const Load &load, int num_all_nets, int num_levels, int num_extra_layers,
		const box &bb, int level, int index, int layer, bool horizontal,
		fpga_tree_node_t *node)
{
	printf("level %d index %d layer %d\n", level, index, layer);

 	printf("\tnum nets %lu box %d %d, %d %d hor %d\n",
			nets.size(),
			bg::get<bg::min_corner, 0>(bb), bg::get<bg::max_corner, 0>(bb),
			bg::get<bg::min_corner, 1>(bb), bg::get<bg::max_corner, 1>(bb),
			horizontal);

	assert(!nets.empty());

	if (level < num_levels-1) {
		//int prev = -1;
		//for (const auto &e : ver_start) {
			//assert(e.first >= prev);
			//prev = e.first;
		//}
		//prev = -1;
		//for (const auto &e : ver) {
			//assert(e.first >= prev);
			//prev = e.first;
		//}

		//for (const auto &ver_e : ver_edges) {
			//PRINT("ver edge %d box %X\n", ver_e.first, ver_e.second);
		//}
		
		int opt_cut;
		float min_diff;

		vector<net_t *> left_nets;
		vector<net_t *> right_nets;

		if (horizontal) {
			split_nets<1>(nets, load,
					num_all_nets, make_pair(bg::get<bg::min_corner, 1>(bb), bg::get<bg::max_corner, 1>(bb)),
					opt_cut, min_diff,
					left_nets, right_nets, node->nets);
		} else {
			split_nets<0>(nets, load,
					num_all_nets, make_pair(bg::get<bg::min_corner, 0>(bb), bg::get<bg::max_corner, 0>(bb)),
					opt_cut, min_diff,
					left_nets, right_nets, node->nets);
		}

		printf("\tleft %lu right %lu middle %lu\n", left_nets.size(), right_nets.size(), node->nets.size());
		printf("\topt_cut %d\n", opt_cut);
		
		box left_box, right_box;
		if (horizontal) {
			assert((opt_cut >= bg::get<bg::min_corner, 1>(bb) && opt_cut <= bg::get<bg::max_corner, 1>(bb)));

			bg::set<bg::min_corner, 0>(left_box, bg::get<bg::min_corner, 0>(bb));
			bg::set<bg::max_corner, 0>(left_box, bg::get<bg::max_corner, 0>(bb));
			bg::set<bg::min_corner, 1>(left_box, bg::get<bg::min_corner, 1>(bb));
			bg::set<bg::max_corner, 1>(left_box, opt_cut);

			bg::set<bg::min_corner, 0>(right_box, bg::get<bg::min_corner, 0>(bb));
			bg::set<bg::max_corner, 0>(right_box, bg::get<bg::max_corner, 0>(bb));
			bg::set<bg::min_corner, 1>(right_box, bg::get<bg::max_corner, 1>(left_box)+1);
			bg::set<bg::max_corner, 1>(right_box, bg::get<bg::max_corner, 1>(bb));

			//printf("hsplit 0 %d->%d,%d->%d 1 %d->%d,%d->%d\n",
			//bg::get<bg::min_corner, 0>(split_0), bg::get<bg::max_corner, 0>(split_0),
			//bg::get<bg::min_corner, 1>(split_0), bg::get<bg::max_corner, 1>(split_0),
			//bg::get<bg::min_corner, 0>(split_1), bg::get<bg::max_corner, 0>(split_1),
			//bg::get<bg::min_corner, 1>(split_1), bg::get<bg::max_corner, 1>(split_1));

			assert((bg::get<bg::min_corner, 1>(left_box) <= bg::get<bg::max_corner, 1>(left_box)));
			assert((bg::get<bg::min_corner, 1>(right_box) <= bg::get<bg::max_corner, 1>(right_box)));
		} else {
			assert((opt_cut >= bg::get<bg::min_corner, 0>(bb) && opt_cut <= bg::get<bg::max_corner, 0>(bb)));

			bg::set<bg::min_corner, 1>(left_box, bg::get<bg::min_corner, 1>(bb));
			bg::set<bg::max_corner, 1>(left_box, bg::get<bg::max_corner, 1>(bb));
			bg::set<bg::min_corner, 0>(left_box, bg::get<bg::min_corner, 0>(bb));
			bg::set<bg::max_corner, 0>(left_box, opt_cut);

			bg::set<bg::min_corner, 1>(right_box, bg::get<bg::min_corner, 1>(bb));
			bg::set<bg::max_corner, 1>(right_box, bg::get<bg::max_corner, 1>(bb));
			bg::set<bg::min_corner, 0>(right_box, bg::get<bg::max_corner, 0>(left_box)+1);
			bg::set<bg::max_corner, 0>(right_box, bg::get<bg::max_corner, 0>(bb));

			//printf("vsplit 0 %d->%d,%d->%d 1 %d->%d,%d->%d\n",
			//bg::get<bg::min_corner, 0>(split_0), bg::get<bg::max_corner, 0>(split_0),
			//bg::get<bg::min_corner, 1>(split_0), bg::get<bg::max_corner, 1>(split_0),
			//bg::get<bg::min_corner, 0>(split_1), bg::get<bg::max_corner, 0>(split_1),
			//bg::get<bg::min_corner, 1>(split_1), bg::get<bg::max_corner, 1>(split_1));

			assert((bg::get<bg::min_corner, 0>(left_box) <= bg::get<bg::max_corner, 0>(left_box)));
			assert((bg::get<bg::min_corner, 0>(right_box) <= bg::get<bg::max_corner, 0>(right_box)));
		}

		printf("\tleft %d %d, %d %d\n",
				bg::get<bg::min_corner, 0>(left_box), bg::get<bg::max_corner, 0>(left_box),
				bg::get<bg::min_corner, 1>(left_box), bg::get<bg::max_corner, 1>(left_box));

		printf("\tright %d %d, %d %d\n",
				bg::get<bg::min_corner, 0>(right_box), bg::get<bg::max_corner, 0>(right_box),
				bg::get<bg::min_corner, 1>(right_box), bg::get<bg::max_corner, 1>(right_box));

		//#ifdef DEBUG_MISR
		//#endif
		
		node->bb = bb;

		//std::vector<std::pair<int, net_t *>> left_ver_start;
		//std::vector<std::pair<int, net_t *>> left_ver_end;
		//std::vector<std::pair<int, net_t *>> right_ver_start;
		//std::vector<std::pair<int, net_t *>> right_ver_end;
		//split_edges_lb(ver_start, 
				//opt_cut, horizontal,
				//left_ver_start, right_ver_start);
		//split_edges_lb(ver_end, 
				//opt_cut, horizontal,
				//left_ver_end, right_ver_end);

		//std::vector<std::pair<int, net_t *>> left_hor_start;
		//std::vector<std::pair<int, net_t *>> left_hor_end;
		//std::vector<std::pair<int, net_t *>> right_hor_start;
		//std::vector<std::pair<int, net_t *>> right_hor_end;
		//split_edges_lb(hor_start, 
				//opt_cut, horizontal,
				//left_hor_start, right_hor_start);
		//split_edges_lb(hor_end, 
				//opt_cut, horizontal,
				//left_hor_end, right_hor_end);

		if (layer < num_extra_layers) {
			node->next = new fpga_tree_node_t;
			split_nets_recursive_multi_layer(node->nets, load, num_all_nets, num_levels, num_extra_layers, 
					bb, level, index, layer+1, !horizontal,
					node->next);
		} else {
			node->next = nullptr;
		}

		if (!left_nets.empty()) {
			node->left = new fpga_tree_node_t;

			split_nets_recursive_multi_layer(left_nets, load, num_all_nets, num_levels, num_extra_layers,
					left_box, level+1, index*2, layer, !horizontal,
					node->left);
		} else {
			node->left = nullptr;
		}

		if (!right_nets.empty()) {
			node->right = new fpga_tree_node_t;

			split_nets_recursive_multi_layer(right_nets, load, num_all_nets, num_levels, num_extra_layers,
					right_box, level+1, index*2+1, layer, !horizontal,
					node->right);
		} else { 
			node->right = nullptr;
		}
	} else {
		printf("\tlast level\n");

		node->bb = bb;
		for (const auto &item : nets) {
			assert(bg::covered_by(item->bounding_box, bb));
			node->nets.push_back(item);
		}
		node->next = nullptr;
		node->left = nullptr;
		node->right = nullptr;
	}
}

void clear_tree(fpga_tree_node_t *node)
{
	if (node) {
		clear_tree(node->left);
		if (node->left) {
			delete node->left;
			node->left = nullptr;
		}

		clear_tree(node->right);
		if (node->right) {
			delete node->right;
			node->right = nullptr;
		}

		node->bb = bg::make_inverse<box>();
		node->nets.clear();
	}
}

void split_box(const box &current, bool horizontal, box &split_0, box &split_1)
{
	if (horizontal) {
		bg::set<bg::min_corner, 0>(split_0, bg::get<bg::min_corner, 0>(current));
		bg::set<bg::max_corner, 0>(split_0, bg::get<bg::max_corner, 0>(current));
		bg::set<bg::min_corner, 1>(split_0, bg::get<bg::min_corner, 1>(current));
		bg::set<bg::max_corner, 1>(split_0, (bg::get<bg::min_corner, 1>(current)+bg::get<bg::max_corner, 1>(current))/2);

		bg::set<bg::min_corner, 0>(split_1, bg::get<bg::min_corner, 0>(current));
		bg::set<bg::max_corner, 0>(split_1, bg::get<bg::max_corner, 0>(current));
		bg::set<bg::min_corner, 1>(split_1, bg::get<bg::max_corner, 1>(split_0)+1);
		bg::set<bg::max_corner, 1>(split_1, bg::get<bg::max_corner, 1>(current));

		//printf("hsplit 0 %d->%d,%d->%d 1 %d->%d,%d->%d\n",
		//bg::get<bg::min_corner, 0>(split_0), bg::get<bg::max_corner, 0>(split_0),
		//bg::get<bg::min_corner, 1>(split_0), bg::get<bg::max_corner, 1>(split_0),
		//bg::get<bg::min_corner, 0>(split_1), bg::get<bg::max_corner, 0>(split_1),
		//bg::get<bg::min_corner, 1>(split_1), bg::get<bg::max_corner, 1>(split_1));

		assert((bg::get<bg::min_corner, 1>(split_0) < bg::get<bg::max_corner, 1>(split_0)));
		assert((bg::get<bg::min_corner, 1>(split_1) < bg::get<bg::max_corner, 1>(split_1)));
	} else {
		bg::set<bg::min_corner, 1>(split_0, bg::get<bg::min_corner, 1>(current));
		bg::set<bg::max_corner, 1>(split_0, bg::get<bg::max_corner, 1>(current));
		bg::set<bg::min_corner, 0>(split_0, bg::get<bg::min_corner, 0>(current));
		bg::set<bg::max_corner, 0>(split_0, (bg::get<bg::min_corner, 0>(current)+bg::get<bg::max_corner, 0>(current))/2);

		bg::set<bg::min_corner, 1>(split_1, bg::get<bg::min_corner, 1>(current));
		bg::set<bg::max_corner, 1>(split_1, bg::get<bg::max_corner, 1>(current));
		bg::set<bg::min_corner, 0>(split_1, bg::get<bg::max_corner, 0>(split_0)+1);
		bg::set<bg::max_corner, 0>(split_1, bg::get<bg::max_corner, 0>(current));

		//printf("vsplit 0 %d->%d,%d->%d 1 %d->%d,%d->%d\n",
		//bg::get<bg::min_corner, 0>(split_0), bg::get<bg::max_corner, 0>(split_0),
		//bg::get<bg::min_corner, 1>(split_0), bg::get<bg::max_corner, 1>(split_0),
		//bg::get<bg::min_corner, 0>(split_1), bg::get<bg::max_corner, 0>(split_1),
		//bg::get<bg::min_corner, 1>(split_1), bg::get<bg::max_corner, 1>(split_1));

		assert((bg::get<bg::min_corner, 0>(split_0) < bg::get<bg::max_corner, 0>(split_0)));
		assert((bg::get<bg::min_corner, 0>(split_1) < bg::get<bg::max_corner, 0>(split_1)));
	}
}

void fpga_bipartition(fpga_tree_node_t *node, int level, bool horizontal)
{
	if (level > 0) {
		node->left = new fpga_tree_node_t;
		node->right = new fpga_tree_node_t;

		split_box(node->bb, horizontal, node->left->bb, node->right->bb);

		fpga_bipartition(node->left, level-1, !horizontal);
		fpga_bipartition(node->right, level-1, !horizontal);
	} else {
		node->left = nullptr;
		node->right = nullptr;
	}
}

void net_tree_to_region_tree(fpga_tree_node_t *net_node, int extra_level, fpga_tree_node_t *region_node)
{
	region_node->bb = net_node->bb;

	if (net_node->left) {
		assert(net_node->right);

		region_node->left = new fpga_tree_node_t;
		region_node->right = new fpga_tree_node_t;
		//region_node->left->bb = net_node->left->bb;
		//region_node->right->bb = net_node->right->bb;

		net_tree_to_region_tree(net_node->left, extra_level, region_node->left);
		net_tree_to_region_tree(net_node->right, extra_level, region_node->right);
	} else {
		int width = bg::get<bg::max_corner, 0>(region_node->bb) - bg::get<bg::min_corner, 0>(region_node->bb);
		int height = bg::get<bg::max_corner, 1>(region_node->bb) - bg::get<bg::min_corner, 1>(region_node->bb);
		assert(width > 0 && height > 0);
		fpga_bipartition(region_node, extra_level, height > width);
	}
}

bool insert_net_into_fpga_tree(fpga_tree_node_t *node, net_t *net, int *total_inserted = nullptr)
{
	if (!node) {
		return false;
	}

	bool covered = bg::covered_by(net->bounding_box, node->bb);
	bool left_inserted = false;
	bool right_inserted = false;

	if (covered) {
		left_inserted = insert_net_into_fpga_tree(node->left, net, total_inserted);
		right_inserted = insert_net_into_fpga_tree(node->right, net, total_inserted);
	}

	//printf("cur %d->%d,%d->%d left %d right %d\n",
				//bg::get<bg::min_corner, 0>(node->bb), bg::get<bg::max_corner, 0>(node->bb),
				//bg::get<bg::min_corner, 1>(node->bb), bg::get<bg::max_corner, 1>(node->bb), left_inserted, right_inserted);

	bool inserted = false;
	if (covered && !left_inserted && !right_inserted) {
		//printf("inserted %d->%d,%d->%d to %d->%d,%d->%d\n",
				//bg::get<bg::min_corner, 0>(net->bounding_box), bg::get<bg::max_corner, 0>(net->bounding_box),
				//bg::get<bg::min_corner, 1>(net->bounding_box), bg::get<bg::max_corner, 1>(net->bounding_box),
				//bg::get<bg::min_corner, 0>(node->bb), bg::get<bg::max_corner, 0>(node->bb),
				//bg::get<bg::min_corner, 1>(node->bb), bg::get<bg::max_corner, 1>(node->bb));
		node->nets.push_back(net);
		inserted = true;
		if (total_inserted) {
			++(*total_inserted);
		}
	}

	assert((inserted && !left_inserted && !right_inserted) ||
			(!inserted && left_inserted && !right_inserted) ||
			(!inserted && !left_inserted && right_inserted) ||
			(!inserted && !left_inserted && !right_inserted));


	return inserted || left_inserted || right_inserted;
}

void fill_fpga_tree(fpga_tree_node_t *root, vector<net_t *> &nets)
{
	for (const auto &net : nets) {
		int total_inserted = 0;
		insert_net_into_fpga_tree(root, net, &total_inserted);
		assert(total_inserted == 1);
	}
}

void print_tree(fpga_tree_node_t *node, int level, int index, int layer)
{
	if (node) {
		if (node->next) {
			print_tree(node->next, level, index, layer+1);
		} else {
			print_tabs(level); printf("level %d index %d layer %d BB %d %d %d %d num nets %lu num vnets %lu max length %d\n",
					level, index, layer,
					bg::get<bg::min_corner, 0>(node->bb), bg::get<bg::max_corner, 0>(node->bb),
					bg::get<bg::min_corner, 1>(node->bb), bg::get<bg::max_corner, 1>(node->bb),
					node->nets.size(), node->vnets.size(), longest_path_length(node->net_g_2));
		}

		print_tree(node->left, layer > 0 ? level : level+1, layer > 0 ? (index*2) - (1 << (level+1)) : index*2, layer);
		print_tree(node->right, layer > 0 ? level : level+1, layer > 0 ? (index*2+1) - (1 << (level+1)) : index*2+1, layer);
	}
}

void print_tree_stats(fpga_tree_node_t *node, router_t *router, int level, int layer)
{
	if (node) {
		if (node->next) {
			print_tree_stats(node->next, router, level, layer+1);
		} else {
			dijkstra_stats_t node_net_stats;
			timer::duration total_update_time = timer::duration::zero();
			timer::duration total_real_route_time = timer::duration::zero();

			node_net_stats.num_heap_pushes = 0;
			node_net_stats.num_heap_pops = 0;
			node_net_stats.num_neighbor_visits = 0;

			for (const auto vnet : node->vnets) {
				node_net_stats.num_heap_pushes += router->net_stats[vnet.global_index].num_heap_pushes;
				node_net_stats.num_heap_pops += router->net_stats[vnet.global_index].num_heap_pops;
				node_net_stats.num_neighbor_visits += router->net_stats[vnet.global_index].num_neighbor_visits;

				total_update_time += router->update_time[vnet.global_index];
				total_real_route_time += router->route_time[vnet.global_index];
			}

			print_tabs(level); printf("num vnets %lu route time %g real route time %g (%g) update time %g (%g) pops/visits/pushes %lu %lu %lu\n",
					node->vnets.size(),
					duration_cast<nanoseconds>(node->route_time).count()/1e9,
					duration_cast<nanoseconds>(total_real_route_time).count()/1e9, 100.0*duration_cast<nanoseconds>(total_real_route_time).count()/duration_cast<nanoseconds>(node->route_time).count(),
					duration_cast<nanoseconds>(total_update_time).count()/1e9, 100.0*duration_cast<nanoseconds>(total_update_time).count()/duration_cast<nanoseconds>(node->route_time).count(),
					node_net_stats.num_heap_pops, node_net_stats.num_neighbor_visits, node_net_stats.num_heap_pushes);
		}

		print_tree_stats(node->left, router, layer > 0 ? level : level+1, layer);
		print_tree_stats(node->right, router, layer > 0 ? level : level+1, layer);
	}
}

void sum_dijkstra_stats(fpga_tree_node_t *node, router_t *router, dijkstra_stats_t &stats)
{
	//if (node) {
		//if (node->next) {
			//sum_dijkstra_stats(node->next, router, stats);
		//} else {
			//for (const auto &vnet : node->vnets) {
				//stats.num_heap_pushes += router->net_stats[vnet.global_index].num_heap_pushes;
				//stats.num_heap_pops += router->net_stats[vnet.global_index].num_heap_pops;
				//stats.num_neighbor_visits += router->net_stats[vnet.global_index].num_neighbor_visits;
			//}
		//}

		//sum_dijkstra_stats(node->left, router, stats);
		//sum_dijkstra_stats(node->right, router, stats);
	//}
	
	for (const auto &net_stats : router->net_stats) {
		stats.num_heap_pushes += net_stats.num_heap_pushes;
		stats.num_heap_pops += net_stats.num_heap_pops;
		stats.num_neighbor_visits += net_stats.num_neighbor_visits;
	}
}

void test_fpga_bipartition()
{
	fpga_tree_node_t root;
	bg::assign_values(root.bb, 1, 2, 10, 13);
	fpga_bipartition(&root, 2, false);
	print_tree(&root, 0, 0, 0);

	net_t net;
	bg::assign_values(net.bounding_box, 1, 2, 10, 13);
	int total_inserted = 0;
	insert_net_into_fpga_tree(&root, &net, &total_inserted);
	assert(total_inserted == 1);
}

void verify_tree_internal(fpga_tree_node_t *node, vector<net_t *> &nets)
{
	for (const auto &net : node->nets) {
		assert(bg::covered_by(net->bounding_box, node->bb));
		nets.push_back(net);
	}

	if (node->left) {
		assert(bg::covered_by(node->left->bb, node->bb));

		verify_tree_internal(node->left, nets);
	}

	if (node->right) {
		assert(bg::covered_by(node->right->bb, node->bb));

		verify_tree_internal(node->right, nets);
	}

	if (node->left && node->right) {
		assert(!bg::intersects(node->left->bb, node->right->bb));
	}
}

void verify_tree(fpga_tree_node_t *root, const vector<net_t *> &nets)
{
	vector<net_t *> net_ptrs;
	verify_tree_internal(root, net_ptrs);

	std::sort(begin(net_ptrs), end(net_ptrs));

	auto nets_copy = nets;
	std::sort(begin(nets_copy), end(nets_copy));

	assert(net_ptrs == nets_copy);
}

//init_region(region_t &r, int gid, )

			//r.gid = gid;
			//bg::assign_values(r.bb, x, y, xmax, ymax);
			//det_mutex_init(r.mutex, e_state, num_threads);
			//r.mods.resize(num_threads);
			//r.num_committed_mods.resize(num_threads);
			//for (auto &n : r.num_committed_mods) {
				//n.resize(num_threads, 0);
			//}

template<typename Load>
void build_net_tree(router_t &router, const t_router_opts *opts, const Load &load)
{
	router.num_levels = opts->num_net_cuts+1;

	extern int nx, ny;
	bool horizontal_cut = ny > nx ? true : false;

	if (opts->net_partitioner == NetPartitioner::Uniform) {
		bg::assign_values(router.net_root.bb, 0, 0, nx+1, ny+1);
		fpga_bipartition(&router.net_root, opts->num_net_cuts, horizontal_cut);

		fill_fpga_tree(&router.net_root, router.current_nets);
	} else {
		assert(opts->net_partitioner == NetPartitioner::Median);

		extern int nx, ny;
		bool horizontal = ny > nx;

		split_nets_recursive_multi_layer(router.current_nets, load, router.all_nets.size(), router.num_levels, 1,
				bg::make<box>(0, 0, nx+1, ny+1), 0, 0, 0, horizontal,
				&router.net_root);
	}

	verify_tree(&router.net_root, router.current_nets);
}

template <class OrderPA, class ColorMap>
	typename boost::property_traits<ColorMap>::value_type
custom_vertex_coloring(const vector<vector<int>> &G, OrderPA &order, 
		ColorMap &color)
{
	typedef typename boost::property_traits<ColorMap>::value_type size_type;

	size_type max_color = 0;
	const size_type V = G.size();

	// We need to keep track of which colors are used by
	// adjacent vertices. We do this by marking the colors
	// that are used. The mark array contains the mark
	// for each color. The length of mark is the
	// number of vertices since the maximum possible number of colors
	// is the number of vertices.

	//Initialize colors 
	for (int v = 0; v < V; ++v)
		put(color, v, V-1);

	std::vector<size_type> mark(V, 
			std::numeric_limits<size_type>::max BOOST_PREVENT_MACRO_SUBSTITUTION());
	//Determine the color for every vertex one by one
	for ( size_type i = 0; i < V; i++) {
		int current = get(order,i);

		/* we are done coloring this node */
		if (get(color, current) != V-1) {
			continue;
		}

		//Mark the colors of vertices adjacent to current.
		//i can be the value for marking since i increases successively
		for (int j = 0; j < G[current].size(); ++j) {
			mark[get(color, G[current][j])] = i; 
		}

		//Next step is to assign the smallest un-marked color
		//to the current vertex.
		size_type j = 0;

		//Scan through all useable colors, find the smallest possible
		//color that is not used by neighbors.  Note that if mark[j]
		//is equal to i, color j is used by one of the current vertex's
		//neighbors.
		while ( j < max_color && mark[j] == i ) 
			++j;

		assert(j <= max_color);
		if (j == max_color) {
			++max_color;
		}
		//while ( max_color <= j)  //All colors are used up. Add one more color
		//++max_color;

		//At this point, j is the smallest possible color
		put(color, current, j);  //Save the color of vertex current
	}

	return max_color;
}

template <class Order, class Degree, 
		 class Marker, class BucketSorter>
void 
custom_smallest_last_vertex_ordering(const vector<vector<int>> & G, Order order, 
		Degree degree, Marker marker,
		BucketSorter& degree_buckets) {
	typedef std::size_t size_type;

	const size_type num = G.size();

	for (int v = 0; v < G.size(); ++v) {
		put(marker, v, num);
		put(degree, v, G[v].size());
		degree_buckets.push(v);
	}

	size_type minimum_degree = 0;
	size_type current_order = num - 1;

	while ( 1 ) {
		typedef typename BucketSorter::stack MDStack;
		MDStack minimum_degree_stack = degree_buckets[minimum_degree];
		while (minimum_degree_stack.empty())
			minimum_degree_stack = degree_buckets[++minimum_degree];

		int node = minimum_degree_stack.top();
		put(order, current_order, node);

		if ( current_order == 0 ) //find all vertices
			break;

		minimum_degree_stack.pop();
		put(marker, node, 0); //node has been ordered.

		for (int v = 0; v < G[node].size(); ++v) {
			int vn = G[node][v];
			if ( get(marker, vn) > current_order ) { //vn is unordered vertex
				put(marker, vn, current_order);  //mark the columns adjacent to node

				//delete vn from the bucket sorter         
				degree_buckets.remove(vn);

				//It is possible minimum degree goes down
				//Here we keep tracking it.
				put(degree, vn, get(degree, vn) - 1); 
				minimum_degree = std::min(minimum_degree, get(degree, vn)); 

				//reinsert vn in the bucket sorter with the new degree
				degree_buckets.push(vn);
			}
		}

		current_order--;
	}

	//at this point, order[i] = v_i;
}

template <class Order, class Degree, class Marker>
void 
custom_smallest_last_vertex_ordering(const vector<vector<int>> & G, Order order, 
		Degree degree, Marker marker) {
	typedef std::size_t size_type;

	const size_type num = G.size();

	typedef boost::bucket_sorter<size_type, int, Degree, boost::identity_property_map> BucketSorter;

	BucketSorter degree_bucket_sorter(num, num, degree);

	custom_smallest_last_vertex_ordering(G, order, degree, marker, degree_bucket_sorter);
}

template <class Order>
void 
custom_smallest_last_vertex_ordering(const vector<vector<int>> & G, Order order) {
	custom_smallest_last_vertex_ordering(G, order,
			make_shared_array_property_map(G.size(), std::size_t(0), boost::identity_property_map()),
			make_shared_array_property_map(G.size(), 0, boost::identity_property_map()));
}

void create_virtual_nets(const vector<net_t *> &nets, int bb_factor, bool all, vector<new_virtual_net_t> &vnets, int &total_num_vnets)
{
	for (const auto &net : nets) {
		double total_area = 0;

		box net_bb;
		bg::assign_inverse(net_bb);
		box net_ortho_bb;
		bg::assign_inverse(net_ortho_bb);

		for (const auto &sink : net->sinks) {
			vector<box> bbs;

#if 0
			box bb;
			bg::assign_inverse(bb);
			bg::expand(bb, point(net->source.x, net->source.y));
			bg::expand(bb, point(sink.x, sink.y));
			bbs.push_back(bb);
#else
			extern int nx, ny;
			make_l_segment(point(net->source.x, net->source.y), point(sink.x, sink.y), true, make_pair(1+bb_factor, bb_factor), make_pair(1+bb_factor, bb_factor), make_pair(nx+1, ny+1), false, bbs);

			/* checking */
			vector<box> ortho_bbs;
			make_l_segment(point(net->source.x, net->source.y), point(sink.x, sink.y), false, make_pair(1+bb_factor, bb_factor), make_pair(1+bb_factor, bb_factor), make_pair(nx+1, ny+1), false, ortho_bbs);

			assert(boxes_area(bbs) == boxes_area(ortho_bbs));

			for (int i = 0; i < ortho_bbs.size(); ++i) {
				bg::expand(net_ortho_bb, ortho_bbs[i]);
			}
#endif
			
			for (int i = 0; i < bbs.size(); ++i) {
				bg::expand(net_bb, bbs[i]);
			}

			total_area += boxes_area(bbs);
		}

		assert(bg::area(net_ortho_bb) == bg::area(net_bb));

		//if (false) {
		if (!all && total_area > bg::area(net_bb)) {
			new_virtual_net_t vnet;

			vnet.local_index = vnets.size();
			vnet.global_index = total_num_vnets++;
			vnet.net = net;
			vnet.source = &net->source;

			//vector<point> points;
			//points.emplace_back(net->source.x, net->source.y);
			//for (const auto &sink : net->sinks) {
				//points.emplace_back(sink.x, sink.y);
			//}
			//box net_bb = make_bb(points, bb_factor);
			//assert(bg::equals(net_bb, net->bounding_box));
			for (int i = 0; i < net->sinks.size(); ++i) {
				vnet.sinks.push_back(&net->sinks[i]);
			}
			//vnet.bounding_boxes.resize(net->sinks.size());

			//assert(vnet.sinks.size() == vnet.bounding_boxes.size());

			vnets.push_back(vnet);
		} else {
			for (const auto &sink : net->sinks) {
				new_virtual_net_t vnet;

				vnet.local_index = vnets.size();
				vnet.global_index = total_num_vnets++;
				vnet.net = net;
				vnet.source = &net->source;

				vnet.sinks.push_back(&sink);
				//vnet.bounding_boxes.resize(1);

				//assert(vnet.sinks.size() == vnet.bounding_boxes.size());
#if 0
				extern int nx, ny;
				make_l_segment(point(net->source.x, net->source.y), point(sink.x, sink.y), true, make_pair(1+bb_factor, bb_factor), make_pair(1+bb_factor, bb_factor), make_pair(nx+1, ny+1), true, vnet.bounding_boxes[0]);

				box normal_bb = make_bb({ point(vnet.source->x, vnet.source->y), point(sink.x, sink.y) }, bb_factor);
				assert(boxes_area(vnet.bounding_boxes[0]) <= bg::area(normal_bb));

				for (int i = 0; i < vnet.bounding_boxes[0].size(); ++i) {
					vnet.overlap_g_boxes.push_back(vnet.bounding_boxes[0][i]);
				}
#endif

				vnets.push_back(vnet);
			}
		}
	}
}

int build_overlap_graph(const vector<new_virtual_net_t> &vnets, bool add_same_net, UndirectedGraph &overlap)
{
	for (int i = 0; i < vnets.size(); ++i) {
		add_vertex(overlap);
	}

	int num_edges = 0;

	for (int i = 0; i < vnets.size(); ++i) {
		for (int j = i+1; j < vnets.size(); ++j) {
			//if (vnets[i].net != vnets[j].net
					//&& bg::intersects(vnets[i].bounding_box, vnets[j].bounding_box)) {
			bool inters = false;
			for (int k = 0; k < vnets[i].overlap_g_boxes.size() && !inters; ++k) {
				const auto &box_a = vnets[i].overlap_g_boxes[k];
				for (int l = 0; l < vnets[j].overlap_g_boxes.size() && !inters; ++l) {
					const auto &box_b = vnets[j].overlap_g_boxes[l];
					inters = bg::intersects(box_a, box_b);
				}
			}
			if (inters) {
			//if (bg::intersects(vnets[i].whole_bounding_box, vnets[j].whole_bounding_box)) {
			//if (intersects(vnets[i].bounding_boxes, vnets[j].bounding_boxes)) {
				if (vnets[i].net != vnets[j].net || add_same_net) {
					assert(i == vnets[i].local_index);
					assert(j == vnets[j].local_index);
					add_edge(vnets[i].local_index, vnets[j].local_index, overlap);
				}
			}
		}
	}

	return num_edges;
}

int build_overlap_graph(const vector<new_virtual_net_t> &vnets, bool add_same_net, vector<vector<int>> &overlap)
{
	overlap.resize(vnets.size());
	for (int i = 0; i < vnets.size(); ++i) {
		assert(overlap[i].empty());
	}

	//printf("building overlap graph\n");

	int num_edges = 0;

	for (int i = 0; i < vnets.size(); ++i) {
		for (int j = i+1; j < vnets.size(); ++j) {
			//if (vnets[i].net != vnets[j].net
					//&& bg::intersects(vnets[i].bounding_box, vnets[j].bounding_box)) {
			bool inters = false;
			for (int k = 0; k < vnets[i].overlap_g_boxes.size() && !inters; ++k) {
				const auto &box_a = vnets[i].overlap_g_boxes[k];
				for (int l = 0; l < vnets[j].overlap_g_boxes.size() && !inters; ++l) {
					const auto &box_b = vnets[j].overlap_g_boxes[l];
					inters = bg::intersects(box_a, box_b);
				}
			}
			if (inters) {
				if (vnets[i].net != vnets[j].net || add_same_net) {
			//if (bg::intersects(vnets[i].whole_bounding_box, vnets[j].whole_bounding_box)) {
			//if (intersects(vnets[i].bounding_boxes, vnets[j].bounding_boxes)) {
					assert(i == vnets[i].local_index);
					assert(j == vnets[j].local_index);
					overlap[vnets[i].local_index].push_back(vnets[j].local_index);
					overlap[vnets[j].local_index].push_back(vnets[i].local_index);

				//printf("adding edge between %d %d\n", i, j);

					num_edges += 2;
				}
			}
		}
	}

	/* TODO: this is very pessimistic. we can remove some edges if we are sure that
	 * the bounding box overlaps with the existing route tree sufficiently */
	//for (auto &virtual_nets : all_virtual_nets) {
		//for (int i = 0; i < virtual_nets.size(); ++i) {
			//for (int j = 0; j < virtual_nets.size(); ++j) {
				//if (i != j) {
					//auto e = edge(virtual_nets[i].v, virtual_nets[j].v, g);
					//if (!e.second) {
						//add_edge(virtual_nets[i].v, virtual_nets[j].v, g);
					//}
				//}
			//}
		//}
	//}

	/* Don't need this here because we are checking the parent based on the virtual net index */
	//for (auto &virtual_nets : all_virtual_nets) {
		//for (int i = virtual_nets.size()-1; i > 0; --i) {
			//auto e = edge(virtual_nets[i].v, virtual_nets[i-1].v, g);
			//assert(!e.second);
			//add_edge(virtual_nets[i].v, virtual_nets[i-1].v, g);
		//}
	//}
	
	//for (auto &virtual_nets : all_virtual_nets) {
		//printf("virtual nets for net %d num sinks %d\n", virtual_nets[0].net->local_id, virtual_nets[0].net->sinks.size());

		//for (int i = 0; i < virtual_nets.size(); ++i) {
			//printf("vnet %d\n", i);

			//if (i == 0) {
				//virtual_nets[i].parent = -1;
			//} else {
				//int num_overlaps = 0;
				//double max_overlap_area = std::numeric_limits<double>::lowest();
				//int best = -1;
				//for (int j = 0; j < i; ++j) {
					//box intersection;
					//bool intersects = bg::intersection(virtual_nets[i].bounding_box, virtual_nets[j].bounding_box, intersection);
					//printf("bb0 %d,%d -> %d,%d bb1 %d,%d -> %d,%d\n", bg::get<0>(virtual_nets[i].bounding_box.min_corner()), bg::get<1>(virtual_nets[i].bounding_box.min_corner()),
							//bg::get<0>(virtual_nets[i].bounding_box.max_corner()), bg::get<1>(virtual_nets[i].bounding_box.max_corner()),
							//bg::get<0>(virtual_nets[j].bounding_box.min_corner()), bg::get<1>(virtual_nets[j].bounding_box.min_corner()),
							//bg::get<0>(virtual_nets[j].bounding_box.max_corner()), bg::get<1>(virtual_nets[j].bounding_box.max_corner()));
					//if (intersects) {
						//assert(bg::intersects(virtual_nets[i].bounding_box, virtual_nets[j].bounding_box));
						//++num_overlaps;
						//double area = bg::area(intersection);
						//printf("\toverlap between vnet %d and %d area %g max area %g\n", i, j, area, max_overlap_area);
						//if (area > max_overlap_area) {
							//max_overlap_area = area;
							//best = virtual_nets[j].global_index;
							//printf("\tsetting max area to %g and best vnet to %d\n", area, best);
						//}
					//} else {
						//assert(!bg::intersects(virtual_nets[i].bounding_box, virtual_nets[j].bounding_box));
						//printf("\tno overlap between vnet %d and %d\n", i, j);
					//}
				//}
				//assert(num_overlaps > 0);
				//assert(best != -1);
				//virtual_nets[i].parent = best;
				
				//overlap[best].push_back(virtual_nets[i].global_index);
			//}
		//}
	//}
	//for (auto &virtual_nets : all_virtual_nets) {
		//for (int i = 0; i < virtual_nets.size(); ++i) {
			//if (i == 0) {
				//virtual_nets[i].parent = -1;
			//} else {
				//virtual_nets[i].parent = virtual_nets[i-1].global_index;

				//overlap[virtual_nets[i-1].global_index].push_back(virtual_nets[i].global_index);
				//++num_edges;
			//}
		//}
	//}
	//int num_orig_directed_edges = 0;
	//for (int i = 0; i < all_virtual_nets_ptr.size(); ++i) {
		//assert(all_virtual_nets_ptr[i]->v == i);
		//assert(all_virtual_nets_ptr[i]->global_index == i);

		//if (all_virtual_nets_ptr[i]->index > 0) {
			//overlap[all_virtual_nets_ptr[i]->parent].push_back(all_virtual_nets_ptr[i]->v);
			//++num_orig_directed_edges;
		//} 
	//}

	return num_edges;
}

void schedule_nets_ind(const vector<net_t *> &nets, vector<vector<const net_t *>> &net_scheduled_at_time)
{
	int max_local_id = 0;
	for (const auto &net : nets) {
		max_local_id = std::max(net->local_id, max_local_id);
	}

	vector<bool> scheduled(max_local_id+1, true);
	for (const auto &net : nets) {
		scheduled[net->local_id] = false;
	}
	net_scheduled_at_time.reserve(nets.size()); //worst case
	int num_scheduled_nets = 0;
	int time = 0;

	using namespace boost::accumulators;
	accumulator_set<decltype(vector<net_t **>().size()), stats<tag::mean, tag::max, tag::min, tag::variance>> acc;

	while (num_scheduled_nets < nets.size()) {
		vector<net_t *> unscheduled_nets;
		for (int i = 0; i < nets.size(); ++i) {
			if (!scheduled[nets[i]->local_id]) {
				unscheduled_nets.push_back(nets[i]);
			}
		}

		vector<net_t **> chosen;
		max_independent_rectangles(unscheduled_nets,
				[] (const net_t *net) -> box { return net->bounding_box; },
				[] (const net_t *net) -> pair<int, int> { return make_pair(bg::get<bg::min_corner, 1>(net->bounding_box), bg::get<bg::max_corner, 1>(net->bounding_box)); },
				chosen);

		net_scheduled_at_time.push_back(vector<const net_t *>());
		for (const auto &c : chosen) {
			scheduled[(*c)->local_id] = true;
			net_scheduled_at_time[time].push_back(*c);
		}
		num_scheduled_nets += chosen.size();
		acc(chosen.size());
		++time;
	}

	printf("NET: nets %lu num colors %lu min %lu max %lu mean %g std dev %g\n", nets.size(), net_scheduled_at_time.size(), boost::accumulators::min(acc), boost::accumulators::max(acc), mean(acc), std::sqrt(variance(acc)));

	assert(all_of(scheduled.begin(), scheduled.end(), [] (const bool &item) -> bool { return item; }));

	for (const auto &at_time : net_scheduled_at_time) {
		verify_ind(at_time, [] (const net_t *net) -> box { return net->bounding_box; });
	}
}

template<typename LoadFunc>
void undirected_to_directed_by_heap(const vector<vector<int>> &g, const vector<new_virtual_net_t> &vnets, const LoadFunc &load, int num_partitions, vector<vector<int>> &directed_g, vector<int> &num_incoming_edges)
{
	vector<int> sorted_nodes(g.size());
	for (int i = 0; i < g.size(); ++i) {
		sorted_nodes[i] = i;
	}
	std::sort(begin(sorted_nodes), end(sorted_nodes), [&g] (int a, int b) -> bool { return g[a].size() > g[b].size(); });

	struct node_item_t {
		int id;
		float start_time;
		unsigned long num_out_edges;

		//bool operator<(const node_item_t &other) const {
			//return start_time > other.start_time || (start_time == other.start_time && num_out_edges < other.num_out_edges);
		//}
		bool operator<(const node_item_t &other) const {
			return start_time > other.start_time;
		}
	};

	boost::heap::binomial_heap<node_item_t, boost::heap::stable<true>> min_node_heap;
	using node_handle_t = typename boost::heap::binomial_heap<node_item_t, boost::heap::stable<true>>::handle_type;

	vector<node_handle_t> node_handles(g.size());

	vector<float> start_time(g.size());

	for (const auto &n : sorted_nodes) {
		start_time[n] = 0;
		node_handles[n] = min_node_heap.push({ n, 0, g[n].size() });
	}

	struct partition_item_t {
		int id;
		float start_time;

		bool operator<(const partition_item_t &other) const {
			return start_time > other.start_time;
		}
	};
	boost::heap::binomial_heap<partition_item_t, boost::heap::stable<true>> min_partition_heap;

	for (int i = 0; i < num_partitions; ++i) {
		min_partition_heap.push({ i, 0 });
	}

	vector<bool> scheduled(g.size(), false);

	directed_g.resize(g.size());
	for (auto &u : directed_g) {
		u.clear();
	}

	num_incoming_edges.resize(g.size());
	std::fill(begin(num_incoming_edges), end(num_incoming_edges), 0);

	while (!min_node_heap.empty()) {
		node_item_t min_node = min_node_heap.top();
		min_node_heap.pop();

		int u = min_node.id;

		assert(!scheduled[u]);
		scheduled[u] = true;

		start_time[u] = min_node.start_time;

		partition_item_t min_partition = min_partition_heap.top();
		min_partition_heap.pop();

		float uload = load(vnets[u]);
		assert(uload > 0);
		float new_start = std::max(min_node.start_time, min_partition.start_time) + uload;
		//float new_start = min_node.start_time + uload;

		min_partition_heap.push({ min_partition.id, new_start });

		PRINT("Current %d [start time %g load %g] sent to partition %d with start time %g\n", u, min_node.start_time, uload, min_partition.id, min_partition.start_time);

		PRINT("New start %g = max(item %g, part %g) + load %g\n", new_start, min_node.start_time, min_partition.start_time, uload);

		/* it is possible that we have less amount of threads than the amount of parallelism available */
		//assert(min_node.start_time >= min_partition.start_time);
		
		for (const auto &v : g[u]) {
			if (scheduled[v]) {
				assert(start_time[v] < start_time[u]);
				continue;
			}

			assert((*node_handles[v]).id == v);
			assert((*node_handles[v]).start_time == start_time[v]);

			PRINT("\tNeighbor %d Current start %g\n", v, (*node_handles[v]).start_time);

			directed_g[u].push_back(v);
			++num_incoming_edges[v];

			if ((*node_handles[v]).start_time == std::numeric_limits<float>::max() || new_start > (*node_handles[v]).start_time) {
				PRINT("\t\tUpdating neighbor %d start time from %g to %g\n", v, (*node_handles[v]).start_time, new_start);

				min_node_heap.update(node_handles[v], { v, new_start, g[v].size() });

				start_time[v] = new_start;
			} 
		}
	}

	assert(std::all_of(begin(scheduled), end(scheduled), [] (const bool &s) -> bool { return s; }));

	int num_edges = 0;
	for (const auto &u : g) {
		num_edges += u.size();
	}

	int num_directed_edges = 0;
	for (const auto &u : directed_g) {
		num_directed_edges += u.size();
	}

	assert((num_edges % 2) == 0);
	assert(num_edges/2 == num_directed_edges);
}

void undirected_to_directed_by_ordering(const vector<vector<int>> &g, vector<vector<int>> &directed_g, vector<int> &num_incoming_edges)
{
	boost::vector_property_map<int> order_map; 

	/* this gives better result than largest first order */
	custom_smallest_last_vertex_ordering(g, order_map); 

	vector<int> order(g.size(), -1);
	for (int i = 0; i < g.size(); ++i) {
		int u = get(order_map, i);
		assert(u >= 0 && u < order.size());
		assert(order[u] == -1);
		order[u] = i;
	}

	//for (int i = 0; i < order.size(); ++i) {
		//printf("order of %d = %d\n", i, order[i]);
	//}

	assert(all_of(begin(order), end(order), [] (const int &val) -> bool { return val != -1; }));

	directed_g.resize(g.size());
	for (auto &u : directed_g) {
		u.clear();
	}
	num_incoming_edges.resize(g.size());
	std::fill(begin(num_incoming_edges), end(num_incoming_edges), 0);

	for (int i = 0; i < g.size(); ++i) {
		int u = get(order_map, i);
		for (int j = 0; j < g[u].size(); ++j) {
			int v = g[u][j];

			//printf("u %d uo %d v %d vo %d\n", u, order[u], v, order[v]);

			assert(order[v] != order[u]);
			if (order[v] > order[u]) {
				directed_g[u].push_back(v);
				++num_incoming_edges[v];
			}
		}
	}

	int num_edges = 0;
	for (const auto &u : g) {
		num_edges += u.size();
	}

	int num_directed_edges = 0;
	for (const auto &u : directed_g) {
		num_directed_edges += u.size();
	}

	assert((num_edges % 2) == 0);
	assert(num_edges/2 == num_directed_edges);
}

void undirected_to_directed_by_coloring(const vector<vector<int>> &g, vector<vector<int>> &directed_g, vector<int> &num_incoming_edges)
{
	boost::vector_property_map<int> order_map;

	/* this gives better result than largest first order */
	custom_smallest_last_vertex_ordering(g, order_map); 

	boost::vector_property_map<int> color_map;

	int num_colors = custom_vertex_coloring(g, order_map, color_map);

	directed_g.resize(g.size());
	for (auto &u : directed_g) {
		u.clear();
	}
	num_incoming_edges.resize(g.size());
	std::fill(begin(num_incoming_edges), end(num_incoming_edges), 0);

	for (int u = 0; u < g.size(); ++u) {
		int uc = get(color_map, u);

		for (int j = 0; j < g[u].size(); ++j) {
			int v = g[u][j];
			int vc = get(color_map, v);

			if (vc > uc) {
				directed_g[u].push_back(v);
				++num_incoming_edges[v];
			}
		}
	}

	int num_edges = 0;
	for (const auto &u : g) {
		num_edges += u.size();
	}

	int num_directed_edges = 0;
	for (const auto &u : directed_g) {
		num_directed_edges += u.size();
	}

	assert((num_edges % 2) == 0);
	assert(num_edges/2 == num_directed_edges);
}

template<typename ColorMap>
void sort_color_map(const ColorMap &color_map, int num_colors, const vector<new_virtual_net_t> &vnets, const vector<net_t *> &nets, vector<int> &color_map_map)
{
	set<int> added;
	vector<int> color_num_sinks(num_colors, 0);

	for (int i = 0; i < vnets.size(); ++i) {
		if (added.find(vnets[i].net->local_id) == added.end()) {
			int c = get(color_map, i);
			assert(c >= 0 && c < color_num_sinks.size());
			color_num_sinks[c] += vnets[i].net->sinks.size();

			added.insert(vnets[i].net->local_id);
		}
	}

	/* checking */
	int total1 = 0;
	for (int val : color_num_sinks) {
		total1 += val;
	}
	int total2 = 0;
	for (const auto &net : nets) {
		total2 += net->sinks.size();
	}
	assert(total1 == total2);

	color_map_map.resize(num_colors);
	std::iota(begin(color_map_map), end(color_map_map), 0);
	std::sort(begin(color_map_map), end(color_map_map), [&color_num_sinks] (int a, int b) -> bool {
			return color_num_sinks[a] > color_num_sinks[b];
			});
}

template<typename ColorMap>
void sort_color_map_2(const ColorMap &color_map, int num_colors, const vector<new_virtual_net_t> &vnets, const vector<net_t *> &nets, vector<int> &color_num_sinks)
{
	set<int> added;

	color_num_sinks.resize(num_colors);
	std::fill(begin(color_num_sinks), end(color_num_sinks), 0);

	for (int i = 0; i < vnets.size(); ++i) {
		if (added.find(vnets[i].net->local_id) == added.end()) {
			int c = get(color_map, i);
			assert(c >= 0 && c < color_num_sinks.size());
			color_num_sinks[c] += vnets[i].net->sinks.size();

			added.insert(vnets[i].net->local_id);
		}
	}

	/* checking */
	int total1 = 0;
	for (int val : color_num_sinks) {
		total1 += val;
	}
	int total2 = 0;
	for (const auto &net : nets) {
		total2 += net->sinks.size();
	}
	assert(total1 == total2);
}

template<typename ColorMap>
void sort_color_map_3(const ColorMap &color_map, int num_colors, const vector<new_virtual_net_t> &vnets, const vector<net_t *> &nets, vector<int> &color_map_map)
{
	vector<set<int>> added(num_colors);
	vector<unsigned long> color_num_sinks(num_colors, 0);

	for (int i = 0; i < vnets.size(); ++i) {
		int c = get(color_map, i);

		assert(c >= 0 && c < color_num_sinks.size());
		assert(added[c].find(vnets[i].net->local_id) == added[c].end());

		//color_num_sinks[c] += vnets[i].net->sinks.size();
		color_num_sinks[c] = std::max(color_num_sinks[c], vnets[i].net->sinks.size());

		added[c].insert(vnets[i].net->local_id);
	}

	/* checking */
	//int total1 = 0;
	//for (int val : color_num_sinks) {
		//total1 += val;
	//}
	//int total2 = 0;
	//for (const auto &net : nets) {
		//total2 += net->sinks.size();
	//}
	//printf("t1t2 %d %d\n", total1, total2);
	//assert(total1 == total2);

	vector<int> sorted_colors(num_colors);
	std::iota(begin(sorted_colors), end(sorted_colors), 0);

	//for (int i = 0; i < sorted_colors.size(); ++i) {
		//printf("iota color %d\n", sorted_colors[i]);
	//}

	std::sort(begin(sorted_colors), end(sorted_colors), [&color_num_sinks] (int a, int b) -> bool {
			return color_num_sinks[a] < color_num_sinks[b];
			});

	/* using value as color and index as weight */
	color_map_map.resize(num_colors);
	for (int i = 0; i < color_map_map.size(); ++i) {
		color_map_map[sorted_colors[i]] = i;
	}
}

template<typename HasEdge>
void undirected_to_directed_internal(const vector<vector<int>> &g, const HasEdge &has_edge, vector<vector<int>> &directed_g, vector<int> &num_incoming_edges)
{
	directed_g.resize(g.size());
	for (auto &u : directed_g) {
		u.clear();
	}
	num_incoming_edges.resize(g.size());
	std::fill(begin(num_incoming_edges), end(num_incoming_edges), 0);

	for (int u = 0; u < g.size(); ++u) {
		for (int j = 0; j < g[u].size(); ++j) {
			int v = g[u][j];

			if (has_edge(u, v)) {
				directed_g[u].push_back(v);
				++num_incoming_edges[v];
			}
		}
	}

	int num_edges = 0;
	for (const auto &u : g) {
		num_edges += u.size();
	}

	int num_directed_edges = 0;
	for (const auto &u : directed_g) {
		num_directed_edges += u.size();
	}

	assert((num_edges % 2) == 0);
	assert(num_edges/2 == num_directed_edges);
}

void undirected_to_directed_by_coloring_sorted(const vector<vector<int>> &g, const vector<new_virtual_net_t> &vnets, const vector<net_t *> &nets, vector<vector<int>> &directed_g, vector<int> &num_incoming_edges)
{
	boost::vector_property_map<int> order_map;

	/* this gives better result than largest first order */
	custom_smallest_last_vertex_ordering(g, order_map); 

	boost::vector_property_map<int> color_map;

	int num_colors = custom_vertex_coloring(g, order_map, color_map);

	vector<int> color_map_map;
	sort_color_map(color_map, num_colors, vnets, nets, color_map_map);

	directed_g.resize(g.size());
	for (auto &u : directed_g) {
		u.clear();
	}
	num_incoming_edges.resize(g.size());
	std::fill(begin(num_incoming_edges), end(num_incoming_edges), 0);

	for (int u = 0; u < g.size(); ++u) {
		int uc = get(color_map, u);
		int sorted_uc = color_map_map[uc];

		for (int j = 0; j < g[u].size(); ++j) {
			int v = g[u][j];
			int vc = get(color_map, v);
			int sorted_vc = color_map_map[vc];

			assert(uc != vc);
			assert(sorted_uc != sorted_vc);

			if (sorted_vc > sorted_uc) {
				directed_g[u].push_back(v);
				++num_incoming_edges[v];
			}
		}
	}

	int num_edges = 0;
	for (const auto &u : g) {
		num_edges += u.size();
	}

	int num_directed_edges = 0;
	for (const auto &u : directed_g) {
		num_directed_edges += u.size();
	}

	assert((num_edges % 2) == 0);
	assert(num_edges/2 == num_directed_edges);
}

void undirected_to_directed_by_coloring_sorted_2(const vector<vector<int>> &g, const vector<new_virtual_net_t> &vnets, const vector<net_t *> &nets, vector<vector<int>> &directed_g, vector<int> &num_incoming_edges)
{
	boost::vector_property_map<int> order_map;

	/* this gives better result than largest first order */
	custom_smallest_last_vertex_ordering(g, order_map); 

	boost::vector_property_map<int> color_map;

	int num_colors = custom_vertex_coloring(g, order_map, color_map);

	vector<int> color_num_sinks;
	sort_color_map_2(color_map, num_colors, vnets, nets, color_num_sinks);

	directed_g.resize(g.size());
	for (auto &u : directed_g) {
		u.clear();
	}
	num_incoming_edges.resize(g.size());
	std::fill(begin(num_incoming_edges), end(num_incoming_edges), 0);

	vector<vector<bool>> edge_added(g.size(), vector<bool>(g.size(), false));

	for (int u = 0; u < g.size(); ++u) {
		int uc = get(color_map, u);

		for (int j = 0; j < g[u].size(); ++j) {
			int v = g[u][j];
			int vc = get(color_map, v);

			assert(uc != vc);

			//if (!edge_added[u][v] && !edge_added[v][u] && color_num_sinks[uc] >= color_num_sinks[vc]) {
			if (uc < vc) {
				directed_g[u].push_back(v);
				//edge_added[u][v] = true;
				++num_incoming_edges[v];
			}
		}
	}

	int num_edges = 0;
	for (const auto &u : g) {
		num_edges += u.size();
	}

	int num_directed_edges = 0;
	for (const auto &u : directed_g) {
		num_directed_edges += u.size();
	}

	assert((num_edges % 2) == 0);
	assert(num_edges/2 == num_directed_edges);
}

void undirected_to_directed_by_coloring_sorted_3(const vector<vector<int>> &g, const vector<new_virtual_net_t> &vnets, const vector<net_t *> &nets, vector<vector<int>> &directed_g, vector<int> &num_incoming_edges)
{
	boost::vector_property_map<int> order_map;

	/* this gives better result than largest first order */
	custom_smallest_last_vertex_ordering(g, order_map); 

	boost::vector_property_map<int> color_map;

	int num_colors = custom_vertex_coloring(g, order_map, color_map);

	vector<int> color_map_map;
	sort_color_map_3(color_map, num_colors, vnets, nets, color_map_map);

	undirected_to_directed_internal(g, [&color_map, &color_map_map] (int u, int v) -> bool {
			int uc = get(color_map, u);
			int vc = get(color_map, v);
			assert(uc != vc);
			int sorted_uc = color_map_map[uc];
			int sorted_vc = color_map_map[vc];
			assert(sorted_uc != sorted_vc);
			return sorted_uc > sorted_vc;
			},
			directed_g, num_incoming_edges);
}

template<typename LoadFunc>
void undirected_to_directed_by_heap_2(const vector<vector<int>> &g, const vector<new_virtual_net_t> &vnets, const LoadFunc &load, int num_partitions, vector<vector<int>> &directed_g, vector<int> &num_incoming_edges)
{
	PRINT("Converting g of size %d\n", g.size());

	boost::vector_property_map<int> order_map;

	custom_smallest_last_vertex_ordering(g, order_map); 

	struct node_item_t {
		int id;
		float start_time;
		unsigned long num_out_edges;

		//bool operator<(const node_item_t &other) const {
			//return start_time > other.start_time || (start_time == other.start_time && num_out_edges < other.num_out_edges);
		//}
		bool operator<(const node_item_t &other) const {
			return start_time > other.start_time;
		}
	};

	boost::heap::binomial_heap<node_item_t, boost::heap::stable<true>> min_node_heap;
	using node_handle_t = typename boost::heap::binomial_heap<node_item_t, boost::heap::stable<true>>::handle_type;

	vector<node_handle_t> node_handles(g.size());

	vector<float> start_time(g.size());

	for (int i = 0; i < g.size(); ++i) {
		int n = get(order_map, i);

		PRINT("Pushing %d with degree %lu\n", n, g[n].size());

		start_time[n] = 0;
		node_handles[n] = min_node_heap.push({ n, 0, g[n].size() });
	}

	struct partition_item_t {
		int id;
		float start_time;

		bool operator<(const partition_item_t &other) const {
			return start_time > other.start_time;
		}
	};
	boost::heap::binomial_heap<partition_item_t, boost::heap::stable<true>> min_partition_heap;

	for (int i = 0; i < num_partitions; ++i) {
		min_partition_heap.push({ i, 0 });
	}

	vector<bool> scheduled(g.size(), false);

	directed_g.resize(g.size());
	for (auto &u : directed_g) {
		u.clear();
	}

	num_incoming_edges.resize(g.size());
	std::fill(begin(num_incoming_edges), end(num_incoming_edges), 0);

	while (!min_node_heap.empty()) {
		node_item_t min_node = min_node_heap.top();
		min_node_heap.pop();

		int u = min_node.id;

		assert(!scheduled[u]);
		scheduled[u] = true;

		start_time[u] = min_node.start_time;

		partition_item_t min_partition = min_partition_heap.top();
		min_partition_heap.pop();

		float uload = load(vnets[u]);
		assert(uload > 0);
		//float new_start = std::max(min_node.start_time, min_partition.start_time) + uload;
		float new_start = min_node.start_time + uload;

		min_partition_heap.push({ min_partition.id, new_start });

		PRINT("Current %d [start time %g load %g] sent to partition %d with start time %g\n", u, min_node.start_time, uload, min_partition.id, min_partition.start_time);

		PRINT("New start %g = max(item %g, part %g) + load %g\n", new_start, min_node.start_time, min_partition.start_time, uload);

		/* it is possible that we have less amount of threads than the amount of parallelism available */
		//assert(min_node.start_time >= min_partition.start_time);
		
		for (const auto &v : g[u]) {
			if (scheduled[v]) {
				assert(start_time[v] < start_time[u]);
				continue;
			}

			assert((*node_handles[v]).id == v);
			assert((*node_handles[v]).start_time == start_time[v]);

			PRINT("\tNeighbor %d Current start %g\n", v, (*node_handles[v]).start_time);

			directed_g[u].push_back(v);
			++num_incoming_edges[v];

			if ((*node_handles[v]).start_time == std::numeric_limits<float>::max() || new_start > (*node_handles[v]).start_time) {
				PRINT("\t\tUpdating neighbor %d start time from %g to %g\n", v, (*node_handles[v]).start_time, new_start);

				min_node_heap.update(node_handles[v], { v, new_start, g[v].size() });

				start_time[v] = new_start;
			} 
		}
	}

	assert(std::all_of(begin(scheduled), end(scheduled), [] (const bool &s) -> bool { return s; }));

	int num_edges = 0;
	for (const auto &u : g) {
		num_edges += u.size();
	}

	int num_directed_edges = 0;
	for (const auto &u : directed_g) {
		num_directed_edges += u.size();
	}

	assert((num_edges % 2) == 0);
	assert(num_edges/2 == num_directed_edges);
}

template<typename LoadFunc>
void schedule_nets_by_convert(fpga_tree_node_t *node, const LoadFunc &load, int num_partitions)
{
	vector<vector<int>> overlap;
	build_overlap_graph(node->vnets, true, overlap);

	//undirected_to_directed(overlap, node->vnets, load, num_partitions, node->net_g, node->num_incoming_edges);
	//undirected_to_directed_4(overlap, node->vnets, load, num_partitions, node->net_g, node->num_incoming_edges);
	//undirected_to_directed_by_ordering(overlap, node->net_g, node->num_incoming_edges);
	//undirected_to_directed_by_coloring(overlap, node->net_g, node->num_incoming_edges);
	undirected_to_directed_by_coloring_sorted(overlap, node->vnets, node->nets, node->net_g, node->num_incoming_edges);

	for (int i = 0; i < node->net_g.size(); ++i) {
		add_vertex(node->net_g_2);
	}

	for (int i = 0; i < node->net_g.size(); ++i) {
		for (int j = 0; j < node->net_g[i].size(); ++j) {
			add_edge(i, node->net_g[i][j], node->net_g_2);
		}
	}

	for (int i = 0; i < node->num_incoming_edges.size(); ++i) {
		if (node->num_incoming_edges[i] == 0) {
			node->roots.push_back(&node->vnets[i]);
		}
	}
	assert(!node->roots.empty());
}

void add_new_merge_vnet(const new_virtual_net_t &vnet, vector<new_virtual_net_t> &merged_vnets)
{
	merged_vnets.push_back(vnet);

	auto &mvnet = merged_vnets.back();

	//assert(mvnet.bounding_boxes.size() == 1);
	//assert(mvnet.overlap_g_boxes.size() == 2);
	//assert(mvnet.bounding_boxes[0].size() == mvnet.overlap_g_boxes.size());
	//assert(bg::equals(mvnet.bounding_boxes[0][0], mvnet.overlap_g_boxes[0]));
	//assert(bg::equals(mvnet.bounding_boxes[0][1], mvnet.overlap_g_boxes[1]));

	//mvnet.overlap_g_boxes.resize(1);
	//bg::assign_inverse(mvnet.overlap_g_boxes[0]);
	//for (const auto &b : mvnet.bounding_boxes[0]) {
		//bg::expand(mvnet.overlap_g_boxes[0], b);
	//}
}

void schedule_nets_by_merge_and_convert(fpga_tree_node_t *node, int &total_num_vnets)
{
	vector<vector<int>> overlap;
	build_overlap_graph(node->vnets, false, overlap);

	//undirected_to_directed(overlap, node->vnets, load, num_partitions, node->net_g, node->num_incoming_edges);
	//undirected_to_directed_4(overlap, node->vnets, load, num_partitions, node->net_g, node->num_incoming_edges);
	//undirected_to_directed_by_ordering(overlap, node->net_g, node->num_incoming_edges);
	//undirected_to_directed_by_coloring(overlap, node->net_g, node->num_incoming_edges);
	for (int i = 0; i < node->vnets.size(); ++i) {
		assert(node->vnets[i].local_index == i);
	}
	
	boost::vector_property_map<int> order_map;

	/* this gives better result than largest first order */
	custom_smallest_last_vertex_ordering(overlap, order_map); 

	boost::vector_property_map<int> color_map;

	int num_colors = custom_vertex_coloring(overlap, order_map, color_map);

	vector<new_virtual_net_t> sorted_vnets = node->vnets;
	sort(begin(sorted_vnets), end(sorted_vnets), [&] (const new_virtual_net_t &a, const new_virtual_net_t &b) -> bool {
			return make_pair(get(color_map, a.local_index), a.net->local_id) < make_pair(get(color_map, b.local_index), b.net->local_id);
			});

	vector<new_virtual_net_t> merged_vnets;
	vector<int> vnet_to_merged_vnet(sorted_vnets.size(), -1);
	vector<int> color_start;

	add_new_merge_vnet(sorted_vnets[0], merged_vnets);
	vnet_to_merged_vnet[sorted_vnets[0].local_index] = 0;
	color_start.push_back(0);

	int current_net_local_id = sorted_vnets[0].net->local_id;
	int current_color = get(color_map, sorted_vnets[0].local_index);

	for (int i = 1; i < sorted_vnets.size(); ++i) {
		bool is_new_net = sorted_vnets[i].net->local_id != current_net_local_id;
		bool is_new_color = get(color_map, sorted_vnets[i].local_index) != current_color;

		//printf("current net local id %d current color %d vnet local index %d color %d net local id %d\n", current_net_local_id, current_color, sorted_vnets[i].local_index, get(color_map, sorted_vnets[i].local_index), sorted_vnets[i].net->local_id);

		if (is_new_net || is_new_color) {
			if (is_new_net) {
				current_net_local_id = sorted_vnets[i].net->local_id;
				//printf("\tnew net local id %d\n", current_net_local_id);
			}
			if (is_new_color) {
				current_color = get(color_map, sorted_vnets[i].local_index);
				color_start.push_back(merged_vnets.size());
				//printf("\tnew color %d color start %d\n", current_color, color_start.back());
			}

			vnet_to_merged_vnet[sorted_vnets[i].local_index] = merged_vnets.size();

			add_new_merge_vnet(sorted_vnets[i], merged_vnets);

			//printf("\tvnet to merged vnet %d -> %d\n", sorted_vnets[i].local_index, vnet_to_merged_vnet[sorted_vnets[i].local_index]);
		} else {
			assert(sorted_vnets[i].sinks.size() == 1);
			assert(sorted_vnets[i].sinks[0]->net == merged_vnets.back().net);

			auto &mvnet = merged_vnets.back();

			for (int j = 0; j < sorted_vnets[i].sinks.size(); ++j) {
				mvnet.sinks.push_back(sorted_vnets[i].sinks[j]);
				mvnet.bounding_boxes.push_back(sorted_vnets[i].bounding_boxes[j]);

				//for (const auto &b : sorted_vnets[i].bounding_boxes[0]) {
				//bg::expand(mvnet.overlap_g_boxes[0], b);
				//}
				for (const auto &b : sorted_vnets[i].bounding_boxes[j]) {
					mvnet.overlap_g_boxes.push_back(b);
				}
			}

			vnet_to_merged_vnet[sorted_vnets[i].local_index] = merged_vnets.size()-1;

			//printf("\tadding to existing mvnet\n");
			//printf("\tvnet to merged vnet %d -> %d\n", sorted_vnets[i].local_index, vnet_to_merged_vnet[sorted_vnets[i].local_index]);
		}
	}
	color_start.push_back(merged_vnets.size());

	for (int i = 0; i < merged_vnets.size(); ++i) {
		merged_vnets[i].local_index = i;
		merged_vnets[i].global_index = total_num_vnets++;
	}

	/* checking */
	int total_sinks_0 = 0;
	for (const auto &vnet : sorted_vnets) {
		total_sinks_0 += vnet.sinks.size();
	}
	int total_sinks_1 = 0;
	for (const auto &vnet : merged_vnets) {
		total_sinks_1 += vnet.sinks.size();
	}
	assert(total_sinks_0 == total_sinks_1);
	for (const auto &vnet : merged_vnets) {
		assert(vnet.sinks.size() == vnet.bounding_boxes.size());
		for (const auto &sink : vnet.sinks) {
			assert(sink->net == vnet.net);
		}
	}
	int debug_num_colors = 0;
	int debug_num_mvnets = 0;
	for (int i = 0; i < color_start.size()-1; ++i) {
		++debug_num_colors;	
		for (int k = color_start[i]; k < color_start[i+1]; ++k) {
			++debug_num_mvnets;
			for (int l = k+1; l < color_start[i+1]; ++l) {
				const auto &vnet_a = merged_vnets[k];
				const auto &vnet_b = merged_vnets[l];

				assert(vnet_a.net != vnet_b.net);
				//if (vnet_a.net != vnet_b.net) {
					bool dis = true;
					for (int m = 0; m < vnet_a.overlap_g_boxes.size() && dis; ++m) {
						const auto &box_a = vnet_a.overlap_g_boxes[m];
						for (int n = 0; n < vnet_b.overlap_g_boxes.size() && dis; ++n) {
							const auto &box_b = vnet_b.overlap_g_boxes[n];
							dis = bg::disjoint(box_a, box_b);
						}
					}
					assert(dis);
				//}
			}
		}
	}
	assert(debug_num_colors == num_colors);
	assert(debug_num_mvnets == merged_vnets.size());
	for(int i = 0; i < vnet_to_merged_vnet.size(); ++i) {
		const auto &mvnet = merged_vnets[vnet_to_merged_vnet[i]];
		assert(mvnet.net == node->vnets[i].net);
	}
	vector<set<int>> added_local_id(num_colors);
	for (int i = 0; i < color_start.size()-1; ++i) {
		for (int k = color_start[i]; k < color_start[i+1]; ++k) {
			assert(added_local_id[i].find(merged_vnets[k].net->local_id) == added_local_id[i].end());
			added_local_id[i].insert(merged_vnets[k].net->local_id);
		}
	}

	/* build graph */
	//assert(node->net_g.empty());
	//node->net_g.resize(merged_vnets.size());
	//node->num_incoming_edges.resize(merged_vnets.size());
	//std::fill(begin(node->num_incoming_edges), end(node->num_incoming_edges), 0);

	vector<vector<int>> new_overlap(merged_vnets.size());
	int num_edges = 0;
	int num_actual_edges = 0;
	vector<vector<bool>> edge_added(merged_vnets.size(), vector<bool>(merged_vnets.size(), false));
	vector<vector<bool>> should_have(merged_vnets.size(), vector<bool>(merged_vnets.size(), false));
	for (int i = 0; i < overlap.size(); ++i) {
		int mvi_a = vnet_to_merged_vnet[i];
		const auto &mvnet_a = merged_vnets[mvi_a];
		assert(mvnet_a.local_index == mvi_a);
		int mvnet_color_a = get(color_map, i);

		for (int j = 0; j < overlap[i].size(); ++j) {
			int mvi_b = vnet_to_merged_vnet[overlap[i][j]];
			const auto &mvnet_b = merged_vnets[mvi_b];
			assert(mvnet_b.local_index == mvi_b);
			int mvnet_color_b = get(color_map, overlap[i][j]);

			assert(mvnet_color_a != mvnet_color_b);

			if (!edge_added[mvi_a][mvi_b]) {
				//printf("should have index %d %d color %d %d net local id %d %d\n", mvi_a, mvi_b, color_num_sinks[mvnet_color_a], color_num_sinks[mvnet_color_b], merged_vnets[mvi_a].net->local_id, merged_vnets[mvi_b].net->local_id);
				if (!should_have[mvi_a][mvi_b]) {
				   	should_have[mvi_a][mvi_b] = true;
					++num_actual_edges;
				}
			}

			if (!edge_added[mvi_a][mvi_b]) {
			//if (!edge_added[mvi_a][mvi_b] && !edge_added[mvi_b][mvi_a] && mvnet_color_a < mvnet_color_b) {
			//if (!edge_added[mvi_a][mvi_b] && !edge_added[mvi_b][mvi_a] && color_num_sinks[mvnet_color_a] >= color_num_sinks[mvnet_color_b]) {
				new_overlap[mvi_a].push_back(mvi_b);
				//++node->num_incoming_edges[mvi_b];

				//printf("\tshould have, added index %d %d color %d %d net local id %d %d\n", mvi_a, mvi_b, color_num_sinks[mvnet_color_a], color_num_sinks[mvnet_color_b], merged_vnets[mvi_a].net->local_id, merged_vnets[mvi_b].net->local_id);

				edge_added[mvi_a][mvi_b] = true;
				++num_edges;
			}
		}
	}
	//assert(num_edges * 2 == num_actual_edges);
	assert(num_edges == num_actual_edges);

	vector<vector<int>> debug_overlap;
	int debug_num_edges = build_overlap_graph(merged_vnets, false, debug_overlap);
	assert(debug_num_edges == num_actual_edges);

	int max_net_local_id = 0;
	for (const auto &net : node->nets) {
		max_net_local_id = std::max(net->local_id, max_net_local_id);
	}

	vector<vector<int>> net_mvnets(max_net_local_id+1);
	for (const auto &mvnet : merged_vnets) {
		net_mvnets[mvnet.net->local_id].push_back(mvnet.local_index);
	}

	int num_extra_edges = 0;
	for (int i = 0; i < net_mvnets.size(); ++i) {
		for (int k = 0; k < net_mvnets[i].size(); ++k) {
			int u = net_mvnets[i][k];
			assert(merged_vnets[u].local_index == u);
			for (int l = k+1; l < net_mvnets[i].size(); ++l) {
				int v = net_mvnets[i][l];
				assert(merged_vnets[v].local_index == v);

				assert(!edge_added[u][v] && !edge_added[v][u]);

				new_overlap[u].push_back(v);
				new_overlap[v].push_back(u);

				num_extra_edges += 2;
			}
		}
	}

	boost::vector_property_map<int> new_order_map;

	/* this gives better result than largest first order */
	custom_smallest_last_vertex_ordering(new_overlap, new_order_map); 

	boost::vector_property_map<int> new_color_map;

	int new_num_colors = custom_vertex_coloring(new_overlap, new_order_map, new_color_map);

	vector<int> color_map_map;
	sort_color_map_3(new_color_map, new_num_colors, merged_vnets, node->nets, color_map_map);

	node->net_g.clear();
	undirected_to_directed_internal(new_overlap, [&new_color_map, &color_map_map] (int u, int v) -> bool {
			int uc = get(new_color_map, u);
			int vc = get(new_color_map, v);
			assert(uc != vc);
			int sorted_uc = color_map_map[uc];
			int sorted_vc = color_map_map[vc];
			assert(sorted_uc != sorted_vc);
			return sorted_uc > sorted_vc;
			//return sorted_vc > sorted_uc;
			},
			node->net_g, node->num_incoming_edges);

	printf("vnets size %lu merged vnets size %lu num colors %d num edges %d num actual edges %d num extra edges %d new num colors %d\n", node->vnets.size(), merged_vnets.size(), num_colors, num_edges, num_actual_edges, num_extra_edges, new_num_colors);

#if 0
	vector<int> mvnet_color(merged_vnets.size(), -1);
	for (int i = 0; i < node->vnets.size(); ++i) {
		int mvi = vnet_to_merged_vnet[i];
		if (mvnet_color[mvi] == -1) {
			mvnet_color[mvi] = get(color_map, i);
		} else {
			assert(mvnet_color[mvi] == get(color_map, i));
		}
	}
	printf("starting check\n");
	bool correct = true;
	for (int i = 0; i < merged_vnets.size(); ++i) {
		for (int j = i+1; j < merged_vnets.size(); ++j) {
			if (should_have[i][j]) {
				assert(should_have[j][i]);
				//assert((edge_added[i][j] && !edge_added[j][i]) || (!edge_added[i][j] && edge_added[j][i]));
				assert(!(edge_added[i][j] && edge_added[j][i]));
				if (!edge_added[i][j] && !edge_added[j][i]) {
					printf("failed check! index %d %d color %d %d net local id %d %d\n", i, j, color_num_sinks[mvnet_color[i]], color_num_sinks[mvnet_color[j]], merged_vnets[i].net->local_id, merged_vnets[j].net->local_id);
					correct = false;
				}
			}
		}
	}
	assert(correct);
#endif

	/* assign back */
	node->vnets = std::move(merged_vnets);

	for (int i = 0; i < node->net_g.size(); ++i) {
		add_vertex(node->net_g_2);
	}

	for (int i = 0; i < node->net_g.size(); ++i) {
		for (int j = 0; j < node->net_g[i].size(); ++j) {
			add_edge(i, node->net_g[i][j], node->net_g_2);
		}
	}

	assert(num_vertices(node->net_g_2) == node->net_g.size());

	sorted_vnets = node->vnets;
	sort(begin(sorted_vnets), end(sorted_vnets), [&] (const new_virtual_net_t &a, const new_virtual_net_t &b) -> bool {
			return make_pair(color_map_map[get(new_color_map, a.local_index)], a.net->local_id) < make_pair(color_map_map[get(new_color_map, b.local_index)], b.net->local_id);
			});

	printf("mvnet sorted start\n");
	for (int i = 0; i < sorted_vnets.size(); ++i) {
		printf("mvnet color %d num sinks %d\n", color_map_map[get(new_color_map, sorted_vnets[i].local_index)], sorted_vnets[i].net->sinks.size());
	}

	vector<int> sorted;
	topological_sort(node->net_g_2, back_inserter(sorted));
	assert(node->vnets.size() == sorted.size());
	//printf("mvnet sorted start\n");
	//for (auto iter = sorted.rbegin(); iter != sorted.rend(); ++iter) {
		//printf("mvnet color %d num sinks %d\n", color_map_map[get(new_color_map, *iter)], node->vnets[*iter].net->sinks.size());
	//}

	for (int i = 0; i < node->num_incoming_edges.size(); ++i) {
		if (node->num_incoming_edges[i] == 0) {
			node->roots.push_back(&node->vnets[i]);
			printf("roots num sinks %lu\n", node->vnets[i].net->sinks.size());
		}
	}
	assert(!node->roots.empty());
}

void schedule_nets_by_coloring(fpga_tree_node_t *node)
{
	//int num_non_decomposed_nets = 0;

	//vector<vector<const net_t*>> blah;
	//schedule_nets_ind(nets, blah);

	//printf("num non %d/%lu (%g)\n", num_non_decomposed_nets, nets.size(), 100.0*num_non_decomposed_nets/nets.size());

	vector<vector<int>> overlap;
	int num_edges = build_overlap_graph(node->vnets, true, overlap);

	vector<int> order(overlap.size());
	auto order_map = make_iterator_property_map(begin(order), boost::identity_property_map());

	/* this gives better result than largest first order */
	custom_smallest_last_vertex_ordering(overlap, order_map); 

	boost::vector_property_map<int> color_map;

	int num_colors = custom_vertex_coloring(overlap, order_map, color_map);

	//OverlapGraph g;
	//boost::vector_property_map<int> order_map_2, color_map_2;

	//build_overlap_graph(vnets, g);
	//smallest_last_vertex_ordering(g, order_map_2);
	//sequential_vertex_coloring(g, order_map_2, color_map_2);

	//int num_colors_2 = 0;
	//for (int i = 0; i < vnets.size(); ++i) {
		//num_colors_2 = std::max(num_colors_2, get(color_map_2, i));
	//}
	//++num_colors_2;

	//for (int i = 0; i< overlap.size(); ++i) {
		//assert(get(order_map, i) == get(order_map_2, i));
		//assert(get(color_map_2, i) == get(color_map_2, i));
	//}

	node->scheduled.resize(num_colors);
	for (int i = 0; i < node->vnets.size(); ++i) {
		int c = get(color_map, i);
		assert(c >= 0 && c < node->scheduled.size());
		node->scheduled[c].push_back(&node->vnets[i]);
	}

	using namespace boost::accumulators;
	accumulator_set<decltype(vector<int>().size()), stats<tag::mean, tag::max, tag::min, tag::variance > > acc;
	for (const auto &color : node->scheduled) {
		acc(color.size());
	}
	printf("num nets %lu num vnets %lu num overlap edges %d num colors %d min %lu max %lu mean %g std dev %g\n", node->nets.size(), node->vnets.size(), num_edges, num_colors, boost::accumulators::min(acc), boost::accumulators::max(acc), mean(acc), std::sqrt(variance(acc)));

	for (const auto &color : node->scheduled) {
		for (int i = 0; i < color.size(); ++i) {
			for (int j = i + 1; j < color.size(); ++j) {
				//if (color[i]->net != color[j]->net) {
				bool dis = true;
				for (int k = 0; k < color[i]->overlap_g_boxes.size() && dis; ++k) {
					const auto &box_a = color[i]->overlap_g_boxes[k];
					for (int l = 0; l < color[j]->overlap_g_boxes.size() && dis; ++l) {
						const auto &box_b = color[j]->overlap_g_boxes[l];
						dis = bg::disjoint(box_a, box_b);
					}
				}
				assert(dis);
					//assert(bg::disjoint(color[i]->whole_bounding_box, color[j]->whole_bounding_box));
					//assert(!intersects(color[i]->bounding_boxes, color[j]->bounding_boxes));
				//}
			}	
		}
	}

	for (const auto &color : node->scheduled) {
		double area = 0;
		for (int i = 0; i < color.size(); ++i) {

		}
	}
}

template<typename LoadFunc>
void schedule_net_tree_by_convert(fpga_tree_node_t *node, int level, int num_logical_threads, const LoadFunc &load)
{
	if (node) {
		int num_tasks = (1 << level);
		int num_partitions = num_logical_threads/num_tasks;

		if (node->next) {
			schedule_net_tree_by_convert(node->next, level, num_logical_threads, load);
		} else {
			schedule_nets_by_convert(node, load, num_partitions);
		}

		assert(node->net_g.size() == node->num_incoming_edges.size());
		assert(node->net_g.size() == node->vnets.size());

		node->current_num_incoming_edges.resize(node->num_incoming_edges.size());
		for (int i = 0; i < node->current_num_incoming_edges.size(); ++i) {
			node->current_num_incoming_edges[i] = new std::atomic<int>;
		}

		schedule_net_tree_by_convert(node->left, level+1, num_logical_threads, load);
		schedule_net_tree_by_convert(node->right, level+1, num_logical_threads, load);
	}
}

template<typename LoadFunc>
void schedule_net_tree_by_merge_and_convert(fpga_tree_node_t *node, int level, int num_logical_threads, const LoadFunc &load, int &total_num_vnets)
{
	if (node) {
		int num_tasks = (1 << level);
		int num_partitions = num_logical_threads/num_tasks;

		if (node->next) {
			schedule_net_tree_by_merge_and_convert(node->next, level, num_logical_threads, load, total_num_vnets);
		} else {
			schedule_nets_by_merge_and_convert(node, total_num_vnets);
		}

		assert(node->net_g.size() == node->num_incoming_edges.size());
		assert(node->net_g.size() == node->vnets.size());

		node->current_num_incoming_edges.resize(node->num_incoming_edges.size());
		for (int i = 0; i < node->current_num_incoming_edges.size(); ++i) {
			node->current_num_incoming_edges[i] = new std::atomic<int>;
		}

		schedule_net_tree_by_merge_and_convert(node->left, level+1, num_logical_threads, load, total_num_vnets);
		schedule_net_tree_by_merge_and_convert(node->right, level+1, num_logical_threads, load, total_num_vnets);
	}
}

void schedule_net_tree_by_coloring(fpga_tree_node_t *node)
{
	if (node) {
		if (node->next) {
			schedule_net_tree_by_coloring(node->next);
		} else {
			schedule_nets_by_coloring(node);
		}

		schedule_net_tree_by_coloring(node->left);
		schedule_net_tree_by_coloring(node->right);
	}
}

void create_virtual_net_tree(fpga_tree_node_t *node, int bb_factor, bool all, int level, int index, int layer, int &total_num_vnets)
{
	if (node) {
		if (node->next) {
			create_virtual_net_tree(node->next, bb_factor, all, level, index, layer+1, total_num_vnets);
		} else {
			printf("Creating virtual nets for level %d index %d layer %d\n", level, index, layer);

			create_virtual_nets(node->nets, bb_factor, all, node->vnets, total_num_vnets);
		}

		create_virtual_net_tree(node->left, bb_factor, all, level+1, index*2, layer, total_num_vnets);
		create_virtual_net_tree(node->right, bb_factor, all, level+1, index*2+1, layer, total_num_vnets);
	}
}

void reset_incoming_edges(fpga_tree_node_t *node)
{
	if (node) {
		if (node->next) {
			reset_incoming_edges(node->next);
		} else {
			for (int i = 0; i < node->num_incoming_edges.size(); ++i) {
				*node->current_num_incoming_edges[i] = node->num_incoming_edges[i];
			}
		}
		//std::copy(begin(node->num_incoming_edges), end(node->num_incoming_edges), begin(node->current_num_incoming_edges));

		reset_incoming_edges(node->left);
		reset_incoming_edges(node->right);
	}
}

void add_sorted(fpga_tree_node_t *node, vector<new_virtual_net_t *> &sorted_vnets)
{
	vector<int> sorted;
	topological_sort(node->net_g_2, back_inserter(sorted));
	assert(node->vnets.size() == sorted.size());
	for (auto iter = sorted.rbegin(); iter != sorted.rend(); ++iter) {
		int u = *iter;
		sorted_vnets.push_back(&node->vnets[u]);
	}
	printf("adding vnets %lu\n", node->vnets.size());
}

void add_pre_schedule_sorted(fpga_tree_node_t *node, vector<new_virtual_net_t *> &sorted_vnets)
{
	int old_size = sorted_vnets.size();

	for (int i = 0; i < node->vnets.size(); ++i) {
		sorted_vnets.push_back(&node->vnets[i]);
	}

	assert(std::distance(begin(sorted_vnets)+old_size, end(sorted_vnets)) == node->vnets.size());

	sort(begin(sorted_vnets)+old_size, end(sorted_vnets), [] (const new_virtual_net_t *a, const new_virtual_net_t *b) -> bool { return a->net->sinks.size() > b->net->sinks.size(); });

	printf("adding vnets %lu\n", node->vnets.size());
}

void get_sorted_vnets_to_route(fpga_tree_node_t *root, vector<new_virtual_net_t *> &sorted_vnets)
{
	struct node_bfs_t {
		fpga_tree_node_t *node;
		bool is_next;
	};
	std::queue<node_bfs_t> nodes;
	nodes.push({ root, false });

	while (!nodes.empty()) {
		node_bfs_t current = nodes.front();
		nodes.pop();

		fpga_tree_node_t *node = current.node;

		//printf("popped vnets %lu next vnets %lu\n", node->vnets.size(), node->next ? node->next->vnets.size() : 0);

		if (node->next) {
			add_pre_schedule_sorted(node->next, sorted_vnets);

			if (node->next->left) {
				add_pre_schedule_sorted(node->next->left, sorted_vnets);
			}

			if (node->next->right) {
				add_pre_schedule_sorted(node->next->right, sorted_vnets);
			}
		} else {
			add_pre_schedule_sorted(node, sorted_vnets);
			//printf("num nets %d num vnets %d\n", node->nets.size(), node->vnets.size());
		}

		if (node->left) {
			nodes.push({ node->left, false });
			//printf("pushing normal left vnets %lu\n", node->left->vnets.size());
		}

		if (node->right) {
			nodes.push({ node->right, false });
			//printf("pushing normal right vnets %lu\n", node->right->vnets.size());
		}
	}
}

#if 0
void l_segment_to_polygon(const vector<box> &segments, polygon &p)
{
	polygon sp1, sp2;
	vector<polygon> spu;

	assert(segments.size() == 2);

	box_to_polygon(segments[0], sp1);
	box_to_polygon(segments[1], sp2);
	bg::union_(sp1, sp2, spu);

	assert(spu.size() == 1);
	assert(boxes_area(segments) == bg::area(spu[0]));
	p = spu[0];
}

//bool operator<(const box &b1, const box &b2)
//{
	//return make_tuple(bg::get<0>(b1.min_corner()), bg::get<1>(b1.min_corner()), bg::get<0>(b1.max_corner()), bg::get<1>(b1.max_corner())) < make_tuple(bg::get<0>(b2.min_corner()), bg::get<1>(b2.min_corner()), bg::get<0>(b2.max_corner()), bg::get<1>(b2.max_corner()));
//}

struct rtree_value_equal {
	template<typename GeometryPtrPair>
		bool operator()(const GeometryPtrPair &a, const GeometryPtrPair &b) const
		{
			return bg::equals(a.first, b.first) && a.second == b.second;
		}
};

void set_vnet_bounding_boxes(vector<new_virtual_net_t *> &sorted_vnets, int bb_factor)
{
	using virtual_rtree_value = pair<box, new_virtual_net_t *>;

	//struct virtual_net_to_rtree_item {
		//virtual_rtree_value operator()(new_virtual_net_t *vnet) const
		//{
			//return make_pair(vnet->bounding_box, vnet);
		//}
	//};

	bgi::rtree<virtual_rtree_value, bgi::rstar<64>, bgi::indexable<virtual_rtree_value>, rtree_value_equal> tree;

	for (int i = 0; i < sorted_vnets.size(); ++i) {
		new_virtual_net_t *vnet = sorted_vnets[i];

		//printf("checking vnet %d\n", i);

		if (vnet->sinks.size() > 1) {
			box vnet_bb;
			bg::assign_inverse(vnet_bb);

			for (int i = 0; i < vnet->sinks.size(); ++i) {
#if 1
				vector<box> bb { make_bb({ point(vnet->source->x, vnet->source->y), point(vnet->sinks[i]->x, vnet->sinks[i]->y) }, bb_factor) };
#else
				vector<box> bb { vnet_bb };
#endif
				vnet->bounding_boxes.emplace_back(bb);
				for (int i = 0; i < bb.size(); ++i) {
					bg::expand(vnet_bb, bb[i]);
				}
			}

			assert(bg::equals(vnet_bb, vnet->net->bounding_box));

			vnet->overlap_g_boxes.push_back(vnet_bb);

			vector<box> flat;
			flatten(vnet->bounding_boxes, flat);
			for (int i = 0; i < flat.size(); ++i) {
				bool covered = false;
				for (int j = 0; j < vnet->overlap_g_boxes.size() && !covered; ++j) {
					covered = bg::covered_by(flat[i], vnet->overlap_g_boxes[j]);
				}
				assert(covered);
			}

			assert(vnet->bounding_boxes.size() == vnet->sinks.size());
			assert(vnet->overlap_g_boxes.size() == 1);

			tree.insert(make_pair(vnet->overlap_g_boxes[0], vnet));
		} else {
			assert(vnet->sinks.size() == 1);

			vector<box> segments, ortho_segments;
			//polygon p, ortho_p;
			//vector<box> overlaps, ortho_overlaps;
			set<int> overlaps, ortho_overlaps;

			extern int nx, ny;

			make_l_segment(point(vnet->source->x, vnet->source->y), point(vnet->sinks[0]->x, vnet->sinks[0]->y), true, make_pair(1+bb_factor, bb_factor), make_pair(1+bb_factor, bb_factor), make_pair(nx+1, ny+1), true, segments);
			//l_segment_to_polygon(segments, p);
			//tree.query(bgi::intersects(p), back_inserter(overlaps)); 
			for (const auto &s : segments) {
				for (auto iter = tree.qbegin(bgi::intersects(s)); iter != tree.qend(); ++iter) {
					assert(iter->second->global_index != vnet->global_index);
					if (iter->second->net != vnet->net) {
						overlaps.insert(iter->second->global_index);
					}
				}
			}

			make_l_segment(point(vnet->source->x, vnet->source->y), point(vnet->sinks[0]->x, vnet->sinks[0]->y), false, make_pair(1+bb_factor, bb_factor), make_pair(1+bb_factor, bb_factor), make_pair(nx+1, ny+1), true, ortho_segments);
			//l_segment_to_polygon(ortho_segments, ortho_p);
			//tree.query(bgi::intersects(ortho_p), back_inserter(ortho_overlaps)); 
			for (const auto &s : ortho_segments) {
				for (auto iter = tree.qbegin(bgi::intersects(s)); iter != tree.qend(); ++iter) {
					assert(iter->second->global_index != vnet->global_index);
					if (iter->second->net != vnet->net) {
						ortho_overlaps.insert(iter->second->global_index);
					}
				}
			}

			const vector<box> &actual_segments = overlaps.size() < ortho_overlaps.size() ? segments : ortho_segments;

			vnet->bounding_boxes.push_back(actual_segments);
			vnet->overlap_g_boxes = actual_segments;
			for (const auto &b : actual_segments) {
				tree.insert(make_pair(b, vnet));
			}

			assert(vnet->bounding_boxes.size() == 1);
			assert(vnet->overlap_g_boxes.size() == 2);
			
			box normal_bb = make_bb({ point(vnet->source->x, vnet->source->y), point(vnet->sinks[0]->x, vnet->sinks[0]->y) }, bb_factor);
			assert(boxes_area(vnet->bounding_boxes[0]) <= bg::area(normal_bb));
		}
	}
}
#endif

void blit(const box &b, int val, vector<vector<int>> &all_net_map)
{
	for (int x = bg::get<bg::min_corner, 0>(b); x <= bg::get<bg::max_corner, 0>(b); ++x) {
		for (int y = bg::get<bg::min_corner, 1>(b); y <= bg::get<bg::max_corner, 1>(b); ++y) {
			assert(x >= 0 && x < all_net_map.size());
			assert(y >= 0 && y < all_net_map[x].size());

			all_net_map[x][y] += val;
		}
	}
}

void blit(const vector<box> &boxes, vector<vector<int>> &all_net_map)
{
	assert(boxes.size() == 2);

	blit(boxes[0], 1, all_net_map);
	blit(boxes[1], 1, all_net_map);
	box inters;
	bg::intersection(boxes[0], boxes[1], inters);
	blit(inters, -1, all_net_map);
}

int get_blit(const box &b, vector<vector<int>> &all_net_map)
{
	int total = 0;
	for (int x = bg::get<bg::min_corner, 0>(b); x <= bg::get<bg::max_corner, 0>(b); ++x) {
		for (int y = bg::get<bg::min_corner, 1>(b); y <= bg::get<bg::max_corner, 1>(b); ++y) {
			assert(x >= 0 && x < all_net_map.size());
			assert(y >= 0 && y < all_net_map[x].size());

			total += all_net_map[x][y];
		}
	}
	return total;
}

int get_blit(const vector<box> &boxes, vector<vector<int>> &all_net_map)
{
	assert(boxes.size() == 2);

	int total = 0;
	total += get_blit(boxes[0], all_net_map);
	total += get_blit(boxes[1], all_net_map);
	box inters;
	bg::intersection(boxes[0], boxes[1], inters);
	total -= get_blit(inters, all_net_map);
	return total;
}

void set_vnet_bounding_boxes_2(vector<new_virtual_net_t *> &sorted_vnets, int bb_factor)
{
	extern int nx, ny;

	vector<vector<int>> all_net_map(nx+2, vector<int>(ny+2, 0));
	vector<vector<vector<int>>> net_map;
	net_map.resize(sorted_vnets.size(), vector<vector<int>>(nx+2, vector<int>(ny+2, 0)));

	for (int i = 0; i < sorted_vnets.size(); ++i) {
		new_virtual_net_t *vnet = sorted_vnets[i];

		//printf("checking vnet %d\n", i);

		if (vnet->sinks.size() > 1) {
			box vnet_bb;
			bg::assign_inverse(vnet_bb);

			for (int i = 0; i < vnet->sinks.size(); ++i) {
#if 1
				vector<box> bb { make_bb({ point(vnet->source->x, vnet->source->y), point(vnet->sinks[i]->x, vnet->sinks[i]->y) }, bb_factor) };
#else
				vector<box> bb { vnet_bb };
#endif
				vnet->bounding_boxes.emplace_back(bb);
				for (int i = 0; i < bb.size(); ++i) {
					bg::expand(vnet_bb, bb[i]);
				}
			}

			assert(bg::equals(vnet_bb, vnet->net->bounding_box));

			vnet->overlap_g_boxes.push_back(vnet_bb);

			vector<box> flat;
			flatten(vnet->bounding_boxes, flat);
			for (int i = 0; i < flat.size(); ++i) {
				bool covered = false;
				for (int j = 0; j < vnet->overlap_g_boxes.size() && !covered; ++j) {
					covered = bg::covered_by(flat[i], vnet->overlap_g_boxes[j]);
				}
				assert(covered);
			}

			assert(vnet->bounding_boxes.size() == vnet->sinks.size());
			assert(vnet->overlap_g_boxes.size() == 1);

			blit(vnet->overlap_g_boxes[0], 1, all_net_map);
			blit(vnet->overlap_g_boxes[0], 1, net_map[vnet->global_index]);
		} else {
			assert(vnet->sinks.size() == 1);

			vector<box> segments, ortho_segments;
			//polygon p, ortho_p;
			//vector<box> overlaps, ortho_overlaps;
			set<int> overlaps, ortho_overlaps;

			make_l_segment(point(vnet->source->x, vnet->source->y), point(vnet->sinks[0]->x, vnet->sinks[0]->y), true, make_pair(1+bb_factor, bb_factor), make_pair(1+bb_factor, bb_factor), make_pair(nx+1, ny+1), true, segments);
			//l_segment_to_polygon(segments, p);
			//tree.query(bgi::intersects(p), back_inserter(overlaps)); 
			int cong = get_blit(segments, all_net_map) - get_blit(segments, net_map[vnet->global_index]);

			make_l_segment(point(vnet->source->x, vnet->source->y), point(vnet->sinks[0]->x, vnet->sinks[0]->y), false, make_pair(1+bb_factor, bb_factor), make_pair(1+bb_factor, bb_factor), make_pair(nx+1, ny+1), true, ortho_segments);
			//l_segment_to_polygon(ortho_segments, ortho_p);
			//tree.query(bgi::intersects(ortho_p), back_inserter(ortho_overlaps)); 
			int ortho_cong = get_blit(ortho_segments, all_net_map) - get_blit(ortho_segments, net_map[vnet->global_index]);

			const vector<box> &actual_segments = cong < ortho_cong ? segments : ortho_segments;

			vnet->bounding_boxes.push_back(actual_segments);
			vnet->overlap_g_boxes = actual_segments;
			blit(actual_segments, all_net_map);
			blit(actual_segments, net_map[vnet->global_index]);

			assert(vnet->bounding_boxes.size() == 1);
			assert(vnet->overlap_g_boxes.size() == 2);
			
			box normal_bb = make_bb({ point(vnet->source->x, vnet->source->y), point(vnet->sinks[0]->x, vnet->sinks[0]->y) }, bb_factor);
			assert(boxes_area(vnet->bounding_boxes[0]) <= bg::area(normal_bb));
		}
	}
}

//void reset_tree_net_stats(fpga_tree_node_t *node)
//{
	//if (node) {
		//node->net_stats.num_heap_pops = 0;
		//node->net_stats.num_heap_pushes = 0;
		//node->net_stats.num_neighbor_visits = 0;

		//reset_tree_net_stats(node->left);
		//reset_tree_net_stats(node->right);
	//}
//}

class my_observer : public tbb::task_scheduler_observer {
	private:
		std::atomic<int> current_core;

	public:
		my_observer()
			: current_core(0)
		{
		}

		virtual void on_scheduler_entry( bool is_worker )
		{
			cpu_set_t cpuset;
			CPU_ZERO(&cpuset);
			for (int i = 0; i < 4; ++i) {
				CPU_SET(i, &cpuset);
			}
			assert(CPU_COUNT(&cpuset) == 4);
			assert(sched_setaffinity(0, sizeof(cpu_set_t), &cpuset) == 0);
		}

		//virtual void on_scheduler_entry( bool is_worker )
		//{
			////int arena_tid = tbb::this_task_arena::current_thread_index();
			////int tid = syscall(SYS_gettid);

			//cpu_set_t cpuset;
			//assert(sched_getaffinity(0, sizeof(cpu_set_t), &cpuset) == 0);

			////printf("[%d,%d] entry, worker %d count %d\n", arena_tid, tid, is_worker, CPU_COUNT(&cpuset));

			//if (CPU_COUNT(&cpuset) > 1) {
				//CPU_ZERO(&cpuset);
				////int setto = current_core++;
				//int setto = current_core.fetch_add(1);
				//CPU_SET(setto, &cpuset);
				//assert(CPU_COUNT(&cpuset) == 1);
				////printf("[%d,%d] affinity set to %d\n", arena_tid, tid, setto);
				//assert(sched_setaffinity(0, sizeof(cpu_set_t), &cpuset) == 0);
			//} else {
				//assert(CPU_COUNT(&cpuset) == 1);
				////printf("[%d,%d] affinity already set\n", arena_tid, tid);
			//}
		//}

		virtual void on_scheduler_exit( bool is_worker )
		{
			//printf("[%d,%d] exit, worker %d\n", tbb::this_task_arena::current_thread_index(), syscall(SYS_gettid), is_worker);
		}
};

void fill_all_vnets(fpga_tree_node_t *node, vector<new_virtual_net_t *> &all_vnets)
{
	if (node) {
		if (node->next) {
			fill_all_vnets(node->next, all_vnets);
		} else {
			for (auto &vnet : node->vnets) {
				all_vnets.push_back(&vnet);
			}
		}

		fill_all_vnets(node->left, all_vnets);
		fill_all_vnets(node->right, all_vnets);
	}
}

void congestion_map(const vector<new_virtual_net_t *> &vnets)
{
	extern int nx, ny;

	vector<vector<int>> cmap(nx+2, vector<int>(ny+2, 0));

	for (const auto &vnet : vnets) {
		for (int i = 0; i < vnet->sinks.size(); ++i) {
			for (int j = 0; j < vnet->bounding_boxes[i].size(); ++j) {
				const box &b = vnet->bounding_boxes[i][j];

				for (int x = bg::get<bg::min_corner, 0>(b); x <= bg::get<bg::max_corner, 0>(b); ++x) {
					for (int y = bg::get<bg::min_corner, 1>(b); y <= bg::get<bg::max_corner, 1>(b); ++y) {
						assert(x >= 0 && x < cmap.size());
						assert(y >= 0 && y < cmap[x].size());

						++cmap[x][y];
					}
				}
			}
		}
	}
	using namespace boost::accumulators;
	accumulator_set<int, stats<tag::mean, tag::max, tag::min, tag::variance>> acc;
	for (int x = 0; x < nx+2; ++x) {
		for (int y = 0; y < ny+2; ++y) {
			acc(cmap[x][y]);
		}
	}
	printf("congestion map size %dx%d min %d max %d mean %g std dev %g\n", nx+2, ny+2, boost::accumulators::min(acc), boost::accumulators::max(acc), mean(acc), std::sqrt(variance(acc)));
}

bool partitioning_multi_sink_delta_stepping_route(const t_router_opts *opts)
{
#if TBB_USE_DEBUG == 1
#error Using debug version of TBB
#endif
	my_observer obs;
	obs.observe();

	tbb::task_scheduler_init init(opts->num_threads);

	char buffer[256];

	int dir_index;
	if (!opts->log_dir) {
		bool dir_exists;
		dir_index = 0;
		do {
			extern char *s_circuit_name;
			sprintf(g_dirname, "%s_stats_%d", s_circuit_name, dir_index);
			struct stat s;
			dir_exists = stat(g_dirname, &s) == 0 ? true : false;
			if (dir_exists) {
				++dir_index;
			}
		} while (dir_exists);
	} else {
		sprintf(g_dirname, opts->log_dir);

		struct stat s;
		if (stat(g_dirname, &s) == 0) {
			printf("log dir %s already exists\n", g_dirname);
			dir_index = -1;
		} else {
			dir_index = 0;
		}
	}

	if (dir_index >= 0) {
		if (mkdir(g_dirname, 0700) == -1) {
			printf("failed to create dir %s\n", g_dirname);
			exit(-1);
		}
	}

	printf("sizeof exec_state_t %lu\n", sizeof(exec_state_t));

	//test_fpga_bipartition();
	//exit(0);

//#if defined (PMC)
	//int retval;
	/* papi library initialization */
	//if ((retval=PAPI_library_init(PAPI_VER_CURRENT)) != PAPI_VER_CURRENT) {
		//printf("Library initialization error! \n");
		//exit(1);
	//}

	/* thread initialization */
    //retval=PAPI_thread_init((unsigned long(*)(void))(pthread_self));
    //if (retval != PAPI_OK)
        //ERROR_RETURN(retval);
//#endif

	init_logging(opts->num_threads);
    zlog_set_record("custom_output", concurrent_log_impl_2);
	sprintf(buffer, "%d", 0);
	zlog_put_mdc("iter", buffer);
	zlog_put_mdc("tid", buffer);
	zlog_put_mdc("log_dir", g_dirname);

	router_t router;

	//assert(((unsigned long)&router.sync.stop_routing & 63) == 0);
	//assert(((unsigned long)&router.sync.net_index & 63) == 0);
#ifdef __linux__
	//assert(((unsigned long)&router.sync.barrier & 63) == 0);
#endif

	init_graph(router.g);

	init_sprintf_rr_node(&router.g);

	vector<net_t> global_nets;
	init_nets(router.all_nets, global_nets, opts->bb_factor, opts->large_bb);

	//for (const auto &net : nets) {
		//const auto &src_prop = get_vertex_props(router.g, net.source.rr_node);
		//assert(net.source.y == src_prop.ylow + src_prop.pin_height);
		//for (const auto &sink : net.sinks) {
			//const auto &sink_prop = get_vertex_props(router.g, sink.rr_node);
			//assert(sink.y == sink_prop.ylow + sink_prop.pin_height);
		//}
	//}

#if defined(TEST)
	//net_t *largest_net = &nets[0];
	//for (auto &net : nets) {
		//if (net.sinks.size() > largest_net->sinks.size()) {
			//largest_net = &net;
		//}
	//}

	//vector<net_t> test_nets;
	//for (int i = 0; i < opts->num_threads; ++i) {
		//test_nets.push_back(*largest_net);
		//test_nets.back().local_id = i;
	//}

	//for (auto &net : test_nets) {
		//router.nets.push_back(&net);
	//}
	sprintf(buffer, "/tmp/nets_inject.txt", g_dirname);
	FILE *nets_in = fopen(buffer, "r");
	while (!feof(nets_in)) {
		if (buffer != fgets(buffer, sizeof(buffer), nets_in)) {
			break;
		}
		int id = atoi(buffer);
		router.current_nets.push_back(&router.all_nets[id]);
		printf("Inject net %d\n", id);
	}
	fclose(nets_in);
#else
	for (int i = 0; i < router.all_nets.size(); ++i) {
		router.current_nets.push_back(&router.all_nets[i]);
	}
	std::sort(begin(router.current_nets), end(router.current_nets), [] (const net_t *a, const net_t *b) -> bool {
			return a->sinks.size() > b->sinks.size();
			});
#endif

	//router.e_state.resize(opts->num_threads);
	//for (auto &e : router.e_state) {
		//e.logical_clock = new std::atomic<unsigned long>(0);
		//e.inet = new std::atomic<int>(0);
	//}
	//g_e_state = router.e_state.data();
	//vector<exec_state_t, tbb::cache_aligned_allocator<exec_state_t>> e_states(4);
	//router.e_state = e_states.data();
	//router.e_state = new exec_state_t[opts->num_threads];

	router.num_threads = opts->num_threads;

#if 0
	fpga_tree_node_t test_root;
	max_independent_rectangles(nets, 3, &test_root);

	verify_tree(&test_root, nets);
	print_tree(&test_root, 0);

	fpga_tree_node_t test_region_root;

	net_tree_to_region_tree(&test_root, &test_region_root, 2);

	print_tree(&test_region_root, 0);
#endif

	//extern int nx, ny;
	//split_nets_recursive(router.current_nets, [] (const net_t *net) -> float { return net->sinks.size(); }, router.all_nets.size(), bg::make<box>(0, 0, nx+1, ny+1), opts->num_net_cuts, 1, true, &router.net_root);

	//exit(0);

	build_net_tree(router, opts, [] (const net_t *net) -> float { return net->sinks.size(); });

	int total_num_vnets = 0;
	create_virtual_net_tree(&router.net_root, opts->bb_factor, false, 0, 0, 0, total_num_vnets);

	printf("total num nets %lu total num vnets %d\n", router.all_nets.size(), total_num_vnets);

	vector<new_virtual_net_t *> sorted_vnets;
	get_sorted_vnets_to_route(&router.net_root, sorted_vnets);

	assert(sorted_vnets.size() == total_num_vnets);

	vector<bool> vnets_added(total_num_vnets, false);
	for (const auto &vnet : sorted_vnets) {
		vnets_added[vnet->global_index] = true;
	}
	assert(all_of(begin(vnets_added), end(vnets_added), [] (const bool &val) -> bool { return val; }));

	/* the check below is incorrect because vnets in different levels are 
	 * sorted independently of one another */
	//int prev_num_sinks = std::numeric_limits<int>::max();
	//for (const auto &vnet : sorted_vnets) {
		//printf("%lu %d\n", vnet->net->sinks.size(), prev_num_sinks);
		//assert(vnet->net->sinks.size() <= prev_num_sinks);
		//prev_num_sinks = vnet->net->sinks.size();
	//}

	//printf("vnets start\n");
	//for (const auto &vnet : sorted_vnets) {
		//printf("vnet %d net %d vnet sinks %d net sinks %d\n", vnet->global_index, vnet->net->local_id, vnet->sinks.size(), vnet->net->sinks.size());
	//}

	set_vnet_bounding_boxes_2(sorted_vnets, opts->bb_factor);

	fill_all_vnets(&router.net_root, router.all_vnets);

	congestion_map(router.all_vnets);

	printf("Scheduling by coloring\n");

	schedule_net_tree_by_coloring(&router.net_root);

	printf("Scheduling by convert\n");

#if 0
	schedule_net_tree_by_convert(&router.net_root, 0, pow(2, opts->num_net_cuts), [] (const new_virtual_net_t &net) ->float { return net.sinks.size(); });
#else
	total_num_vnets = 0;
	schedule_net_tree_by_merge_and_convert(&router.net_root, 0, pow(2, opts->num_net_cuts), [] (const new_virtual_net_t &net) ->float { return net.sinks.size(); }, total_num_vnets);
#endif

	printf("total num nets %lu total num vnets %d\n", router.all_nets.size(), total_num_vnets);

	//write_net_g(g_dirname, &router.net_root, 0, 0);

	printf("Net tree\n");
	print_tree(&router.net_root, 0, 0, 0);

	router.all_vnets.clear();

	fill_all_vnets(&router.net_root, router.all_vnets);

	congestion_map(router.all_vnets);

	//init_fpga_regions(4, router.e_state, opts->num_threads, fpga_regions);

	//vector<net_t> nets_copy = nets;
	//std::sort(begin(nets_copy), end(nets_copy), [] (const net_t &a, const net_t &b) -> bool {
			//return a.sinks.size() > b.sinks.size();
			//});
	//vector<net_t> only;
	//boost::copy(nets_copy | boost::adaptors::sliced(0, 5000), std::back_inserter(only));
	//boost::copy(nets_copy, std::back_inserter(only));

	//best_case(nets, router.nets);
	
	//vector<vector<int>> overlap2;
	//vector<new_virtual_net_t *> all_virtual_nets_ptr2;

	//auto build_start = timer::now();

	//build_overlap_graph_3(all_virtual_nets, all_virtual_nets_ptr2, overlap2);

	//auto build_time = timer::now()-build_start;

	//printf("overlap graph build 3 time %g\n", duration_cast<nanoseconds>(build_time).count() / 1e9);

	//vector<vector<int>> overlap1;
	//vector<new_virtual_net_t *> all_virtual_nets_ptr1;
	//build_start = timer::now();
	//build_overlap_graph_2(all_virtual_nets, overlap1, all_virtual_nets_ptr1);
	//build_time = timer::now()-build_start;

	//printf("overlap graph build 2 time %g\n", duration_cast<nanoseconds>(build_time).count() / 1e9);

	//for (int i = 0; i < overlap1.size(); ++i) {
		//vector<int> o1_temp = overlap1[i];
		//vector<int> o2_temp = overlap2[i];

		//std::sort(begin(o1_temp), end(o1_temp));
		//std::sort(begin(o2_temp), end(o2_temp));

		//assert(o1_temp == o2_temp);
	//}

	//int num_orig_directed_edges = 0;
	//for (int i = 0; i < all_virtual_nets_ptr.size(); ++i) {
		//assert(all_virtual_nets_ptr[i]->v == i);
		//assert(all_virtual_nets_ptr[i]->global_index == i);

		//if (all_virtual_nets_ptr[i]->index > 0) {
			//overlap[all_virtual_nets_ptr[i]->parent].push_back(all_virtual_nets_ptr[i]->v);
			//++num_orig_directed_edges;
		//} 
	//}

	//schedule_virtual_nets_5(overlap, all_virtual_nets_ptr, [&all_virtual_nets_ptr] (int v) -> float { return all_virtual_nets_ptr[v]->sinks.size(); }, opts->num_threads, directed, num_incoming_edges);

	//assert(directed.size() == router.directed.size());
	//for (int i = 0; i < directed.size(); ++i) {
		//auto sorted_directed = directed[i];
		//auto sorted_router_directed = router.directed[i];
		//std::sort(begin(sorted_directed), end(sorted_directed));
		//std::sort(begin(sorted_router_directed), end(sorted_router_directed));
		//assert(sorted_directed == sorted_router_directed);
	//}

	//assert(num_incoming_edges.size() == router.num_incoming_edges.size());
	//for (int i = 0; i < num_incoming_edges.size(); ++i) {
		//assert(num_incoming_edges[i] == router.num_incoming_edges[i]);
	//}

	//int num_edges = 0;
	//for (int i = 0; i < overlap.size(); ++i) {
		//num_edges += overlap[i].size();
	//}
	
	//int num_directed_edges = 0;
	//for (int i = 0; i < router.directed.size(); ++i) {
		//num_directed_edges += router.directed[i].size();
	//}

	//assert(((num_edges-num_orig_directed_edges) % 2) == 0);
	//assert(num_directed_edges-num_orig_directed_edges == (num_edges-num_orig_directed_edges) / 2);


	//schedule_virtual_nets_3(overlap, all_virtual_nets_ptr, order, opts->num_threads, nullptr, router.nets);
	//schedule_virtual_nets_4(overlap, all_virtual_nets_ptr, order, opts->num_threads, nullptr, router.directed, router.topo_nets);
	//test_misr(nets);
	//write_nets(nets);

	//sprintf(buffer, "%s/virtual_nets_bb_%g.txt", dirname, opts->bb_area_threshold_scale);
	//write_net_stats(all_virtual_nets_ptr, buffer, [] (const new_virtual_net_t *vnet) -> box { return vnet->bounding_box; }, [] (const new_virtual_net_t *vnet) -> int { return vnet->net->sinks.size(); });

	//congestion_t *congestion_aligned;
//#ifdef __linux__
	//assert(posix_memalign((void **)&congestion_aligned, 64, sizeof(congestion_t)*num_vertices(router.g)) == 0);
	//assert(((unsigned long)congestion_aligned & 63) == 0);
//#else
	//congestion_aligned = (congestion_t *)malloc(sizeof(congestion_t)*num_vertices(router.g));
//#endif 
	//router.state.congestion = new(congestion_aligned) congestion_t[num_vertices(router.g)];
	router.phase_two = false;

	router.state.congestion = new congestion_t[num_vertices(router.g)];
	for (int i = 0; i < num_vertices(router.g); ++i) {
		router.state.congestion[i].acc_cost = 1;
		router.state.congestion[i].pres_cost = 1;
		router.state.congestion[i].occ = 0;
	}

	//g_congestion = router.state.congestion;

	router.state.route_trees = new route_tree_t[router.all_nets.size()];
	//router.state.back_route_trees = new route_tree_t[nets.size()];
	for (int i = 0; i < router.all_nets.size(); ++i) {
		route_tree_init(router.state.route_trees[i], &router.g);
		//route_tree_init(router.state.back_route_trees[i]);
	}

	router.state.added_rr_nodes = new vector<RRNode>[total_num_vnets];
	router.state.back_added_rr_nodes = new vector<RRNode>[total_num_vnets];

    router.state.net_timing = new t_net_timing[router.all_nets.size()+global_nets.size()];
    init_net_timing(router.all_nets, global_nets, router.state.net_timing);

    router.params.criticality_exp = opts->criticality_exp;
    router.params.astar_fac = opts->astar_fac;
    router.params.max_criticality = opts->max_criticality;
    router.params.bend_cost = opts->bend_cost;

	router.pres_fac = opts->first_iter_pres_fac;

	router.pmc_overflow = opts->pmc_overflow;

	//void *aligned_current_num_incoming_edges;
	//assert(posix_memalign(&aligned_current_num_incoming_edges, 64, sizeof(aligned_atomic<int, 64>)*nets.size()) == 0);
	//assert(((unsigned long)aligned_current_num_incoming_edges & 63) == 0);
	//router.sync.current_num_incoming_edges = new(aligned_current_num_incoming_edges) aligned_atomic<int, 64>[num_virtual_nets];
	//router.sync.current_num_incoming_edges = new std::atomic<int>[num_virtual_nets];

	router.routed.resize(total_num_vnets);
	router.route_time.resize(total_num_vnets);
	router.update_time.resize(total_num_vnets);
	router.net_stats.resize(total_num_vnets);
	router.sink_stats.resize(total_num_vnets);
	for (int i = 0; i < router.sink_stats.size(); ++i) {
		router.sink_stats[i].resize(router.all_vnets[i]->sinks.size());
	}

	//router.locks = new tbb::spin_mutex[router.all_nets.size()];

	router.iter = 0;

	g_manager = new Manager(router.g, router.state.congestion, router.pres_fac);
	g_manager->reserve(opts->num_threads * 4);
	printf("num allocs %d\n", g_manager->get_num_allocs());
	assert(g_manager->get_num_allocs() == opts->num_threads*4);

	dijkstra_stats_t total_stats { 0, 0, 0 };
	timer::duration total_time = timer::duration::zero();

	init_wait_stats(router.num_levels, opts->num_threads);
	init_inc_time(opts->num_threads);

	int num_repartitions = 0;

	bool routed = false;
	unsigned long prev_num_overused_nodes = std::numeric_limits<unsigned long>::max();
	for (; router.iter < opts->max_router_iterations && !routed; ++router.iter) {
		/* hack */
		//tbb::parallel_for(tbb::blocked_range<int>(0, opts->num_threads, 1),
				//[&] (const tbb::blocked_range<int> &r) -> void {
				//sprintf(buffer, "%d", router.iter.load());
				//zlog_put_mdc("iter", buffer);
				//sprintf(buffer, "%d", 0);
				//zlog_put_mdc("tid", buffer);
				//zlog_put_mdc("log_dir", g_dirname);
				//},
				//tbb::simple_partitioner());

		printf("Routing iteration: %d Num repartitions: %d\n", router.iter.load(), num_repartitions);

		//assert(router.sync.num_incoming_edges.empty());

		//router.sync.handles.resize(router.topo_nets.size());

		//for (const auto &vnet : router.topo_nets) {
			//router.sync.handles[vnet->global_index] = router.sync.num_incoming_edges.emplace(make_pair(router.num_incoming_edges[vnet->global_index], vnet));
		//}
		
		//assert(router.sync.fast_global_pending_nets.size_approx() > 0);
		sprintf(buffer, "%s/nets_iter_%d.txt", g_dirname, router.iter.load());

		FILE *nets_out = fopen(buffer, "w");
		for (const auto &net : router.current_nets) {
			fprintf(nets_out, "%d\n", net->local_id);
		}

		for (int i = 0; i < total_num_vnets; ++i) {
			router.routed[i] = false;
			router.route_time[i] = timer::duration::zero();
			router.update_time[i] = timer::duration::zero();
		}

		for (auto &net : router.current_nets) {
			update_sink_criticalities(*net, router.state.net_timing[net->vpr_id], router.params);
		}

		clear_wait_stats();
		clear_inc_time();

		reset_incoming_edges(&router.net_root);
		//reset_tree_net_stats(&router.net_root);

		//do_virtual_work(&router, *net_router);
		tbb::task &root = *new(tbb::task::allocate_root()) MultiLayerSchedRouterTask(&router, &router.net_root, 0);

		auto route_start = timer::now();

		tbb::task::spawn_root_and_wait(root);
		
		auto route_time = timer::now()-route_start;

		//printf("Resched time [has_resched %d]: %g\n", has_resched ? 1 : 0, duration_cast<nanoseconds>(resched_time).count()/1e9);
		
		__itt_task_begin(domain, __itt_null, __itt_null, shAnalysisTask);

		for (const auto &vnet : router.all_vnets) {
			router.net_stats[vnet->global_index].num_heap_pops = 0;
			router.net_stats[vnet->global_index].num_heap_pushes = 0;
			router.net_stats[vnet->global_index].num_neighbor_visits = 0;

			for (int i = 0; i < vnet->sinks.size(); ++i) {
				router.net_stats[vnet->global_index].num_heap_pops += router.sink_stats[vnet->global_index][i].num_heap_pops;
				router.net_stats[vnet->global_index].num_heap_pushes += router.sink_stats[vnet->global_index][i].num_heap_pushes;
				router.net_stats[vnet->global_index].num_neighbor_visits += router.sink_stats[vnet->global_index][i].num_neighbor_visits;
			}
		}

		printf("Free sink routers: %d Num allocs %d\n", g_manager->get_num_free(), g_manager->get_num_allocs());
		assert(g_manager->get_num_free() > 0);
		
		printf("Iteration route time: %g\n", duration_cast<nanoseconds>(route_time).count()/1e9);

		dijkstra_stats_t iteration_stats { 0, 0, 0 };
		sum_dijkstra_stats(&router.net_root, &router, iteration_stats);

		printf("Iteration num heap pops: %lu\n", iteration_stats.num_heap_pops);
		printf("Iteration num neighbor visits: %lu\n", iteration_stats.num_neighbor_visits);
		printf("Iteration num heap pushes: %lu\n", iteration_stats.num_heap_pushes);

		total_time += route_time;

		total_stats.num_heap_pops += iteration_stats.num_heap_pops;
		total_stats.num_heap_pushes += iteration_stats.num_heap_pushes;
		total_stats.num_neighbor_visits += iteration_stats.num_neighbor_visits;

		print_tree_stats(&router.net_root, &router, 0, 0);

		//printf("Route time:\n");
		//for (int level = 0; level < router.num_levels; ++level) {
			//for (int i = 0; i < opts->num_threads; ++i) {
				//printf("%g (%g) ", duration_cast<nanoseconds>(router.time[level][i].route_time).count()/1e9, 100.0*router.time[level][i].route_time.count()/router.total_time[i].count());
			//}
			//printf("\n");
		//}

		//printf("Route time debug: ");
		//for (auto iter = begin(router.route_time); iter != end(router.route_time); ++iter) {
			//printf("%g ", duration_cast<nanoseconds>(*iter).count()/1e9);
		//}
		//printf("\n");

		//printf("Stats wait time:\n");
		//for (int level = 0; level < router.num_levels; ++level) {
			//for (int i = 0; i < opts->num_threads; ++i) {
				//const auto &stats = get_wait_stats(level, i);
				//timer::duration total = timer::duration::zero();
				//for (const auto &s : stats) {
					//total += s.wait_time;
				//}
				//printf("%g (%g) ", duration_cast<nanoseconds>(total).count()/1e9, 100.0*total.count()/router.total_time[i].count());
			//}
			//printf("\n");
		//}

		//printf("Wait time:\n");
		//for (int level = 0; level < router.num_levels; ++level) {
			//for (int i = 0; i < opts->num_threads; ++i) {
				//printf("%g (%g) ", duration_cast<nanoseconds>(router.time[level][i].wait_time).count()/1e9, 100.0*router.time[level][i].wait_time.count()/router.total_time[i].count());
			//}
			//printf("\n");
		//}

		//printf("Get net wait time:\n");
		//for (int level = 0; level < router.num_levels; ++level) {
			//for (int i = 0; i < opts->num_threads; ++i) {
				//printf("%g (%g) ", duration_cast<nanoseconds>(router.time[level][i].get_net_wait_time).count()/1e9, 100.0*router.time[level][i].get_net_wait_time.count()/router.total_time[i].count());
			//}
			//printf("\n");
		//}

		//printf("Pure route time:\n");
		//for (int level = 0; level < router.num_levels; ++level) {
			//for (int i = 0; i < opts->num_threads; ++i) {
				//printf("%g (%g) ", duration_cast<nanoseconds>(router.time[level][i].pure_route_time).count()/1e9, 100.0*router.time[level][i].pure_route_time.count()/router.total_time[i].count());
			//}
			//printf("\n");
		//}

		//printf("Last barrier wait time:\n");
		//for (int level = 0; level < router.num_levels; ++level) {
			//for (int i = 0; i < opts->num_threads; ++i) {
				//printf("%g (%g) ", duration_cast<nanoseconds>(router.time[level][i].last_barrier_wait_time).count()/1e9, 100.0*router.time[level][i].last_barrier_wait_time.count()/router.total_time[i].count());
			//}
			//printf("\n");
		//}

		//printf("Last sync time:\n");
		//for (int level = 0; level < router.num_levels; ++level) {
			//for (int i = 0; i < opts->num_threads; ++i) {
				//printf("%g (%g) ", duration_cast<nanoseconds>(router.time[level][i].last_sync_time).count()/1e9, 100.0*router.time[level][i].last_sync_time.count()/router.total_time[i].count());
			//}
			//printf("\n");
		//}

		//sprintf(buffer, "%s/net_route_time_%d.txt", dirname, router.iter);
		//FILE *nrt = fopen(buffer, "w");

		//for (int i = 0; i < router.nets.size(); ++i) {
			//for (int j = 0; j < router.nets[i].size(); ++j) {
				//const auto *vnet = router.nets[i][j];
				//fprintf(nrt, "%d %d %d %g %d %d %d %g\n", i, router.nets[i].size(), router.net_router[vnet->global_index], duration_cast<nanoseconds>(router.net_route_time[vnet->global_index]).count()/1e9, vnet->index, vnet->net->sinks.size(), vnet->sinks.size(), bg::area(vnet->bounding_box));
			//}
		//}

		//fclose(nrt);

		if (routed) {
			sprintf(buffer, "%s/routes_final.txt", g_dirname);
		} else {
			sprintf(buffer, "%s/routes_iter_%d.txt", g_dirname, router.iter.load());
		}
		write_routes(buffer, router.state.route_trees, router.all_nets.size());

		//if (routed) {
			//sprintf(buffer, "%s/net_stats_final.txt", g_dirname);
		//} else {
			//sprintf(buffer, "%s/net_stats_%d.txt", g_dirname, router.iter.load());
		//}
		//write_net_stats(buffer, &router);

		if (routed) {
			sprintf(buffer, "%s/congestion_state_final.txt", g_dirname);
		} else {
			sprintf(buffer, "%s/congestion_state_%d.txt", g_dirname, router.iter.load());
		}
		write_congestion_state(buffer, router.state.congestion, num_vertices(router.g));

		write_nets(g_dirname, &router.net_root, &router, 0, 0, 0);
		write_sink_stats(g_dirname, &router.net_root, &router, 0, 0, 0);

		//write_mutex_stats(g_dirname, router.iter, router.num_levels, opts->num_threads);

		//for (int i = 0; i < opts->num_threads; ++i) {
			//sprintf(buffer, "%s/inc_iter_%d_tid_%d.txt", g_dirname, router.iter.load(), i);
			//write_inc_time(buffer, get_inc_time(i));
		//}

		/* checking */
		all_of(begin(router.routed), end(router.routed), [] (const bool &v) -> bool { return v; });
		for (int i = 0; i < num_vertices(router.g); ++i) {
			router.state.congestion[i].recalc_occ = 0; 
		}

		for (const auto &net : router.all_nets) {
			check_route_tree(router.state.route_trees[net.local_id], net, router.g);
		}

		for (const auto &net : router.all_nets) {
			recalculate_occ(router.state.route_trees[net.local_id], router.g, router.state.congestion);
		}

		bool valid = true;
		for (int i = 0; i < num_vertices(router.g); ++i) {
			const auto &props = get_vertex_props(router.g, i);
			if (router.state.congestion[i].recalc_occ != router.state.congestion[i].occ) {
				sprintf_rr_node_impl(i, buffer);
				printf("Node %s occ mismatch, recalc: %d original: %d\n", buffer, router.state.congestion[i].recalc_occ, router.state.congestion[i].occ);
				valid = false;
			}
		}
		assert(valid);

		if (router.iter == 0) {
			//lmao2(router, opts, nets.size(), [&router] (const net_t *net) -> float { return router.net_stats[net->local_id].num_heap_pops; });
		}

		unsigned long num_overused_nodes = 0;
		vector<int> overused_nodes_by_type(NUM_RR_TYPES, 0);

		if (feasible_routing(router.g, router.state.congestion)) {
			//dump_route(*current_traces_ptr, "route.txt");
			routed = true;
		} else {
			for (int i = 0; i < num_vertices(router.g); ++i) {
				if (router.state.congestion[i].occ > get_vertex_props(router.g, i).capacity) {
					++num_overused_nodes;
					const auto &v_p = get_vertex_props(router.g, i);
					++overused_nodes_by_type[v_p.type];
				}
			}

			//bool switch_to_phase_two = (router.iter > 10 && num_overused_nodes > prev_num_overused_nodes);
			bool switch_to_phase_two = false;
			if (switch_to_phase_two) {
				vector<net_t *> congested_nets;
				for (const auto &net : router.current_nets) {
					if (route_tree_is_congested(router.state.route_trees[net->local_id], router.g, router.state.congestion)) {
						congested_nets.push_back(net);
					}
				}
				assert(congested_nets.size() <= router.current_nets.size());

				printf("Num congested nets: %lu\n", congested_nets.size());

				router.current_nets = std::move(congested_nets);

				clear_tree(&router.net_root);
				build_net_tree(router, opts, [&router] (const net_t *net) -> float { return router.net_stats[net->local_id].num_heap_pops; });

				router.phase_two = true;

				++num_repartitions;

				prev_num_overused_nodes = std::numeric_limits<unsigned long>::max();
			} else {
				prev_num_overused_nodes = num_overused_nodes;
			}

			//auto update_cost_start = timer::now();

			if (router.iter == 0) {
				router.pres_fac = opts->initial_pres_fac;

				update_costs(router.g, router.state.congestion, router.pres_fac, 0);
			} else {
				router.pres_fac *= opts->pres_fac_mult;

				/* Avoid overflow for high iteration counts, even if acc_cost is big */
				router.pres_fac = std::min(router.pres_fac, static_cast<float>(HUGE_POSITIVE_FLOAT / 1e5));

				update_costs(router.g, router.state.congestion, router.pres_fac, opts->acc_fac);
			}

			//update_cost_time = timer::now()-update_cost_start;
		}

		//auto analyze_timing_start = timer::now();

		float crit_path_delay = analyze_timing(router.state.net_timing);

		//analyze_timing_time = timer::now()-analyze_timing_start;
		
		printf("Overused: %lu/%d (%g) Crit path delay: %g\n", num_overused_nodes, num_vertices(router.g), 100.0*num_overused_nodes/num_vertices(router.g), crit_path_delay);
		printf("\n");

		//for (int i = 0; i < nets.size(); ++i) {
			//route_tree_clear(router.state.back_route_trees[i]);
		//}
		//std::swap(router.state.route_trees, router.state.back_route_trees);
		for (int i = 0; i < total_num_vnets; ++i) {
			router.state.back_added_rr_nodes[i].clear();
		}
		std::swap(router.state.added_rr_nodes, router.state.back_added_rr_nodes);
		for (int i = 0; i < router.current_nets.size(); ++i) {
			int id = router.current_nets[i]->local_id;
			route_tree_clear(router.state.route_trees[id]);
		}

		__itt_task_end(domain);
	}

	if (routed) {
		printf("Routed in %d iterations\n", router.iter.load());
	}

	printf("Total num heap pops: %lu\n", total_stats.num_heap_pops);
	printf("Total num neighbor visits: %lu\n", total_stats.num_neighbor_visits);
	printf("Total num heap pushes: %lu\n", total_stats.num_heap_pushes);
	printf("Total route time: %g\n", duration_cast<nanoseconds>(total_time).count() / 1e9);

	return false;
}
