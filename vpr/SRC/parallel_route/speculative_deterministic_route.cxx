#include "pch.h"
#include <thread>
#include <condition_variable>
#include <sys/stat.h>

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

using namespace std;

using timer = std::chrono::high_resolution_clock;
using std::chrono::duration_cast;
using std::chrono::nanoseconds;

typedef struct extra_route_state_t {
	float upstream_R;
	float delay;
} extra_route_state_t;

typedef struct dijkstra_stats_t {
	int num_heap_pops;
	int num_heap_pushes;
	int num_neighbor_visits;
} dijkstra_stats_t;

typedef struct global_route_state_t {
	congestion_local_t *congestion;
	route_tree_t *route_trees;
	//route_tree_t *back_route_trees;
	vector<RRNode> *added_rr_nodes;	
	vector<RRNode> *back_added_rr_nodes;	
	t_net_timing *net_timing;
} global_route_state_t;

typedef struct new_virtual_net_t {
	int global_index;
	int index;
	net_t *net;
	vector<sink_t *> sinks;
	box bounding_box;
	bool l_segment;
	vector<box> bounding_boxes;
	point last_point;
	/* set when build the overlap graph */
	int parent;
	/* set by the scheduler */
	//int num_incoming_edges;
} new_virtual_net_t;

static char dirname[256];

static bool outside_bb(const rr_node_property_t &node, const vector<box> &boxes)
{
	bool outside = true;

	for (int i = 0; i < boxes.size() && outside; ++i) {
		if (node.type == CHANX || node.type == CHANY) {
			int x, y;
			if (node.inc_direction) {
				x = node.xlow;
				y = node.ylow;
			} else {
				x = node.xhigh;
				y = node.yhigh;
			}

			if (x >= bg::get<bg::min_corner, 0>(boxes[i])
					&& x <= bg::get<bg::max_corner, 0>(boxes[i])
					&& y >= bg::get<bg::min_corner, 1>(boxes[i])
					&& y <= bg::get<bg::max_corner, 1>(boxes[i])) {
				outside = false;
			}
		} else {
			if (node.xhigh >= bg::get<bg::min_corner, 0>(boxes[i])
					&& node.xlow <= bg::get<bg::max_corner, 0>(boxes[i])
					&& node.yhigh >= bg::get<bg::min_corner, 1>(boxes[i])
					&& node.ylow <= bg::get<bg::max_corner, 1>(boxes[i])) {
				outside = false;
			}
		}
	}

	return outside;
}

class SpeculativeDeterministicRouter {
	private:
		RRGraph &_g;

		float _astar_fac;

		float *_known_distance;
		float *_distance;
		RREdge *_prev_edge;
		extra_route_state_t *_state;

		congestion_local_t *_congestion;
		const float &_pres_fac;

		vector<int> _modified_nodes;
		vector<bool> _modified_node_added;

		const sink_t *_current_sink;
		RRNode _existing_opin;

		route_tree_t *_current_rt;

		dijkstra_stats_t _stats;

		static tbb::atomic<int> _num_instances;
		int _instance;

	private:
		void popped_node(const heap_node_t<RREdge> &node)
		{
			char buffer[256];
			sprintf_rr_node(node.node, buffer);
			zlog_level(delta_log, ROUTER_V3, "%s\n", buffer);
			++_stats.num_heap_pops;

			//rand_delay();
		}

		void relax_node(const heap_node_t<RREdge> &node)
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

			int v = node.node;
			const auto &v_p = get_vertex_props(_g, v);

			if (valid(node.prev_edge)) {
				assert(v == get_target(_g, node.prev_edge));

				int u = get_source(_g, node.prev_edge);

				float u_delay;
				float u_upstream_R;
				if (_state[u].upstream_R != std::numeric_limits<float>::max()) {
					assert(_state[u].delay != std::numeric_limits<float>::max());
					u_upstream_R = _state[u].upstream_R;
					u_delay = _state[u].delay;
				} else {
					auto rt_node = route_tree_get_rt_node(*_current_rt, u);
					assert(rt_node != RouteTree::null_vertex());
					const auto &u_rt_node_p = get_vertex_props(_current_rt->graph, rt_node);
					u_upstream_R = u_rt_node_p.upstream_R;
					u_delay = u_rt_node_p.delay;
				}

				const auto &e_p = get_edge_props(_g, node.prev_edge);

				extern struct s_switch_inf *switch_inf;
				const struct s_switch_inf *sw = &switch_inf[e_p.switch_index];

				_state[v].upstream_R = sw->R + v_p.R;
				if (!sw->buffered)  {
					_state[v].upstream_R += u_upstream_R;
				} 

				float delay;
				if (sw->buffered) {
					delay = sw->Tdel + v_p.C * (sw->R + 0.5f * v_p.R);
				} else {
					delay = sw->Tdel + v_p.C * (u_upstream_R + sw->R + 0.5f * v_p.R);
				}
				_state[v].delay = u_delay + delay;
			} else {
				_state[v].upstream_R = v_p.R;
				_state[v].delay = v_p.C * 0.5f * v_p.R;
			}

			if (!_modified_node_added[v]) {
				_modified_nodes.push_back(v);
				_modified_node_added[v] = true;
			}
		}

		bool expand_node(int node)
		{
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
				}
			}

			if (_existing_opin != RRGraph::null_vertex() && prop.type == OPIN && node != _existing_opin) {
				zlog_level(delta_log, ROUTER_V3, "not expanding other OPIN\n");
				return false;
			}

			//if (prop.xlow >= _current_sink->current_bounding_box.xmin
					//&& prop.xhigh <= _current_sink->current_bounding_box.xmax
					//&& prop.ylow >= _current_sink->current_bounding_box.ymin
					//&& prop.yhigh <= _current_sink->current_bounding_box.ymax) {

			//if (prop.xlow < _current_sink->current_bounding_box.xmin
					//|| prop.xhigh > _current_sink->current_bounding_box.xmax
					//|| prop.ylow < _current_sink->current_bounding_box.ymin
					//|| prop.yhigh > _current_sink->current_bounding_box.ymax) {
				//zlog_level(delta_log, ROUTER_V3, "outside of bounding box\n");
				//return false;
			//}
			
			//if (prop.xhigh < _current_sink->current_bounding_box.xmin
					//|| prop.xlow > _current_sink->current_bounding_box.xmax
					//|| prop.yhigh < _current_sink->current_bounding_box.ymin
					//|| prop.ylow > _current_sink->current_bounding_box.ymax) {
				//zlog_level(delta_log, ROUTER_V3, "outside of bounding box\n");
				//return false;
			//}


			if (outside_bb(prop, _current_sink->vnet->bounding_boxes)) {
				zlog_level(delta_log, ROUTER_V3, "outside of bounding box\n");
				return false;
			}

			zlog_level(delta_log, ROUTER_V3, " OK\n");

			++_stats.num_neighbor_visits;

			return true;
		}

		void push_node(const heap_node_t<RREdge> &node)
		{
			++_stats.num_heap_pushes;
		}

		pair<float, float> get_edge_weight(const RREdge &e)
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
				delay = sw->Tdel + v_p.C * (sw->R + 0.5f * v_p.R);
			} else {
				delay = sw->Tdel + v_p.C * (_state[u].upstream_R + sw->R + 0.5f * v_p.R);
			}
			extern t_rr_indexed_data *rr_indexed_data;
			float congestion_cost = rr_indexed_data[v_p.cost_index].base_cost * get_acc_cost(_congestion, v) * get_pres_cost(_congestion, v);
			float known_cost = _current_sink->criticality_fac * delay + (1 - _current_sink->criticality_fac) * congestion_cost;

			float upstream_R = sw->R + v_p.R;
			if (!sw->buffered) {
				upstream_R += _state[u].upstream_R;
			}
			float expected_cost = get_timing_driven_expected_cost(v_p, get_vertex_props(_g, _current_sink->rr_node), _current_sink->criticality_fac, upstream_R);

			zlog_level(delta_log, ROUTER_V3, "\t%d -> %d delay %g congestion %g crit_fac %g expected %g expected_hex %X known %g predicted %g\n", 
					u, v, delay, congestion_cost, _current_sink->criticality_fac, expected_cost, *(unsigned int *)&expected_cost, known_cost, known_cost + _astar_fac * expected_cost);
			//zlog_level(delta_log, ROUTER_V3, "\t[u: upstream %g] [edge: d %g R %g] [v: R %g C %g]\n",
					//_state[u].upstream_R, sw->Tdel, sw->R, v_p.R, v_p.C);

			return make_pair(known_cost, known_cost + _astar_fac * expected_cost);
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
			char buffer[256];
			sprintf_rr_node(sink_node, buffer);
			zlog_level(delta_log, ROUTER_V3, "Adding %s to route tree\n", buffer);

			RouteTreeNode rt_node = route_tree_add_rr_node(rt, sink_node, _g);
			assert(rt_node != RouteTree::null_vertex());
			assert(_state[sink_node].upstream_R != std::numeric_limits<float>::max() && _state[sink_node].delay != std::numeric_limits<float>::max());
			route_tree_set_node_properties(rt, rt_node, false, _state[sink_node].upstream_R, _state[sink_node].delay);

			const auto &sink_node_p = get_vertex_props(_g, sink_node);
			//assert(sink_node_p.xlow >= _current_sink->current_bounding_box.xmin
					//&& sink_node_p.xhigh <= _current_sink->current_bounding_box.xmax
					//&& sink_node_p.ylow >= _current_sink->current_bounding_box.ymin
					//&& sink_node_p.yhigh <= _current_sink->current_bounding_box.ymax);
			//assert(sink_node_p.xhigh >= _current_sink->current_bounding_box.xmin
					//&& sink_node_p.xlow <= _current_sink->current_bounding_box.xmax
					//&& sink_node_p.yhigh >= _current_sink->current_bounding_box.ymin
					//&& sink_node_p.ylow <= _current_sink->current_bounding_box.ymax);
			assert(!outside_bb(sink_node_p, _current_sink->vnet->bounding_boxes));

			update_one_cost_internal(sink_node, _g, _congestion, 1, _pres_fac);
			added_rr_nodes.push_back(sink_node);

			zlog_level(delta_log, ROUTER_V3, "\n");

			RREdge edge = get_previous_edge(sink_node, rt);

			RRNode child_rr_node = sink_node;

			while (valid(edge)) {
				RRNode parent_rr_node = get_source(_g, edge);

				const auto &parent_rr_node_p = get_vertex_props(_g, parent_rr_node);

				sprintf_rr_node(parent_rr_node, buffer);
				zlog_level(delta_log, ROUTER_V3, "Adding %s to route tree\n", buffer);

				RouteTreeNode parent_rt_node = route_tree_add_rr_node(rt, parent_rr_node, _g);
				if (parent_rt_node != RouteTree::null_vertex()) {
					assert(_state[parent_rr_node].upstream_R != std::numeric_limits<float>::max() && _state[parent_rr_node].delay != std::numeric_limits<float>::max());
					route_tree_set_node_properties(rt, parent_rt_node, parent_rr_node_p.type != IPIN && parent_rr_node_p.type != SINK, _state[parent_rr_node].upstream_R, _state[parent_rr_node].delay);
				} 

				assert(child_rr_node == get_target(_g, edge));

				char buffer2[256];
				sprintf_rr_node(child_rr_node, buffer2);

				zlog_level(delta_log, ROUTER_V3, "Adding edge %s -> %s\n", buffer, buffer2);

				const RouteTreeEdge &rt_edge = route_tree_add_edge_between_rr_node(rt, parent_rr_node, child_rr_node); 
				auto &rt_edge_props = get_edge_props(rt.graph, rt_edge);
				rt_edge_props.rr_edge = edge;

				if (parent_rr_node_p.type == OPIN && _existing_opin == RRGraph::null_vertex()) {
					_existing_opin = parent_rr_node;
				}

				if ((parent_rt_node != RouteTree::null_vertex()) || (get_vertex_props(_g, parent_rr_node).type == SOURCE)) {
					const auto &rr_node_p = get_vertex_props(_g, parent_rr_node);
					//assert(rr_node_p.xlow >= _current_sink->current_bounding_box.xmin
							//&& rr_node_p.xhigh <= _current_sink->current_bounding_box.xmax
							//&& rr_node_p.ylow >= _current_sink->current_bounding_box.ymin
							//&& rr_node_p.yhigh <= _current_sink->current_bounding_box.ymax);
					//assert(rr_node_p.xhigh >= _current_sink->current_bounding_box.xmin
							//&& rr_node_p.xlow <= _current_sink->current_bounding_box.xmax
							//&& rr_node_p.yhigh >= _current_sink->current_bounding_box.ymin
							//&& rr_node_p.ylow <= _current_sink->current_bounding_box.ymax);
					assert(!outside_bb(rr_node_p, _current_sink->vnet->bounding_boxes));

					update_one_cost_internal(parent_rr_node, _g, _congestion, 1, _pres_fac);
					added_rr_nodes.push_back(parent_rr_node);
				}

				zlog_level(delta_log, ROUTER_V3, "\n");

				child_rr_node = parent_rr_node;
				edge = get_previous_edge(parent_rr_node, rt);
			}
		}

		//template<typename Graph, typename Edge, typename EdgeWeightFunc, typename Callbacks>
		//friend void delta_stepping(const Graph &g, const vector<heap_node_t<Edge>> &sources, int sink, float delta, float *known_distance, float *distance, Edge *prev_edge, const EdgeWeightFunc &edge_weight, Callbacks &callbacks);

		template<typename Graph, typename Edge, typename EdgeWeightFunc, typename ExpandCheckFunc, typename Callbacks>
		friend void dijkstra(const Graph &g, const vector<heap_node_t<Edge>> &sources, int sink, float *known_distance, float *distance, Edge *prev_edge, const ExpandCheckFunc &expand_node, const EdgeWeightFunc &edge_weight, Callbacks &callbacks);

		//template<typename Edge, typename Callbacks>
		//friend void relax(Buckets &buckets, float delta, vector<bool> &in_bucket, const vector<bool> &vertex_deleted,
					//float *known_distance, float *distance, Edge *predecessor,
					//int v, float new_known_distance, float new_distance, const Edge &edge,
					//Callbacks &callbacks);

	public:
		SpeculativeDeterministicRouter(RRGraph &g, congestion_local_t *congestion, const float &pres_fac)
			: _g(g), _pres_fac(pres_fac), _modified_node_added(num_vertices(g), false)
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

			_instance = _num_instances++;
			_congestion = &congestion[_instance];
		}

		int get_instance() const
		{
			return _instance;
		}

		int get_num_instances() const
		{
			return _num_instances;
		}

		void reset_stats()
		{
			_stats.num_heap_pops = 0;
			_stats.num_heap_pushes = 0;
			_stats.num_neighbor_visits = 0;
		}

		const dijkstra_stats_t &get_stats() const
		{
			return _stats;
		}

		void route(const source_t *source, const vector<const sink_t *> &sinks, float astar_fac, route_tree_t &rt, vector<RRNode> &added_rr_nodes, t_net_timing &net_timing)
		{
			char buffer[256];

			_astar_fac = astar_fac;

			/* TODO: we should not reset the exisiting OPIN here */
			_existing_opin = RRGraph::null_vertex();

			if (source) {
				const auto &source_p = get_vertex_props(_g, source->rr_node);
				RouteTreeNode root_rt_node = route_tree_add_rr_node(rt, source->rr_node, _g);
				route_tree_set_node_properties(rt, root_rt_node, true, source_p.R, 0.5f * source_p.R * source_p.C);
				route_tree_add_root(rt, source->rr_node);
			} else {
				assert(rt.root_rt_nodes.size() == 1);
			}

			vector<const sink_t *> sorted_sinks = sinks;

			std::sort(begin(sorted_sinks), end(sorted_sinks), [] (const sink_t *a, const sink_t *b) -> bool {
					return a->criticality_fac > b->criticality_fac;
					});

			_current_rt = &rt;

			int isink = 0;
			for (const auto &sink : sorted_sinks) {
				_current_sink = sink;

				sprintf_rr_node(sink->rr_node, buffer);
				//zlog_level(delta_log, ROUTER_V3, "Current sink: %s BB %d->%d, %d->%d\n", buffer, sink->current_bounding_box.xmin, sink->current_bounding_box.xmax, sink->current_bounding_box.ymin, sink->current_bounding_box.ymax);
				zlog_level(delta_log, ROUTER_V3, "Current sink: %s BB ", buffer);
				for (int i = 0; i < sink->vnet->bounding_boxes.size(); ++i) {
				   	zlog_level(delta_log, ROUTER_V3, "%d->%d,%d->%d ", bg::get<bg::min_corner, 0>(sink->vnet->bounding_boxes[i]), bg::get<bg::max_corner, 0>(sink->vnet->bounding_boxes[i]), bg::get<bg::min_corner, 1>(sink->vnet->bounding_boxes[i]), bg::get<bg::max_corner, 1>(sink->vnet->bounding_boxes[i]));
				}
				zlog_level(delta_log, ROUTER_V3, "\n");

				vector<heap_node_t<RREdge>> sources;

				for (const auto &rt_node : route_tree_get_nodes(rt)) {
					const auto &rt_node_p = get_vertex_props(rt.graph, rt_node);
					RRNode node = rt_node_p.rr_node;
					const auto &node_p = get_vertex_props(_g, node);

					sprintf_rr_node(node, buffer);

					if (rt_node_p.reexpand) {
						//if (node_p.xhigh < sink->current_bounding_box.xmin
								//|| node_p.xlow > sink->current_bounding_box.xmax
								//|| node_p.yhigh < sink->current_bounding_box.ymin
								//|| node_p.ylow > sink->current_bounding_box.ymax) {
						if (outside_bb(node_p, _current_sink->vnet->bounding_boxes)) {
							zlog_level(delta_log, ROUTER_V3, "Existing %s out of bounding box\n", buffer);
						} else {
							float kd = sink->criticality_fac * rt_node_p.delay; 
							float d = kd + _astar_fac * get_timing_driven_expected_cost(node_p, get_vertex_props(_g, sink->rr_node), sink->criticality_fac, rt_node_p.upstream_R); 

							RREdge prev = RRGraph::null_edge();
							if (valid(rt_node_p.rt_edge_to_parent)) {
								prev = get_edge_props(rt.graph, rt_node_p.rt_edge_to_parent).rr_edge;
							}

							zlog_level(delta_log, ROUTER_V3, "Adding %s back to heap [delay %g upstream_R %g] [kd=%g d=%g prev=%d]\n", buffer, rt_node_p.delay, rt_node_p.upstream_R, kd, d, get_source(_g, prev));

							sources.push_back({ node, kd, d, prev });
						}
					} else {
						zlog_level(delta_log, ROUTER_V3, "Not reexpanding %s\n", buffer);
					}
				}

				//delta_stepping(_g, sources, sink->rr_node, delta, _known_distance, _distance, _prev_edge, [this] (const RREdge &e) -> pair<float, float> { return get_edge_weight(e); }, *this);
				dijkstra(_g, sources, sink->rr_node, _known_distance, _distance, _prev_edge, [this] (int rr_node) -> bool { return expand_node(rr_node); }, [this] (const RREdge &e) -> pair<float, float> { return get_edge_weight(e); }, *this);

				backtrack(sink->rr_node, rt, added_rr_nodes);

				assert(get_vertex_props(rt.graph, route_tree_get_rt_node(rt, sink->rr_node)).delay == _state[sink->rr_node].delay);
				net_timing.delay[sink->id+1] = _state[sink->rr_node].delay;

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

				++isink;
			}
		}
};

tbb::atomic<int> SpeculativeDeterministicRouter::_num_instances = 0;

typedef struct worker_sync_t {
	/*alignas(64)*/ tbb::atomic<int> stop_routing;
	/*alignas(64)*/ tbb::atomic<int> net_index;
#ifdef PTHREAD_BARRIER
	/*alignas(64)*/ pthread_barrier_t barrier;
#else
	//alignas(64) boost::barrier *barrier;
	//alignas(64) SpinningBarrier *barrier;
	SpinningBarrier *barrier;
#endif
	//std::atomic<int> nums_net_to_route;
	//std::atomic<bool> has_global;
} worker_sync_t;

template<typename Net>
struct /*alignas(64)*/ router_t {
	vector<Net> nets;
	tbb::atomic<int> iter;
	worker_sync_t sync;
	RRGraph g;
	global_route_state_t state;
	route_parameters_t params;
	float pres_fac;
	int num_threads;
	vector<dijkstra_stats_t> stats;
	vector<timer::duration> route_time;
	vector<timer::duration> wait_time;
	vector<timer::duration> pure_route_time;
	vector<timer::duration> update_time;
	vector<timer::duration> last_barrier_wait_time;
	vector<vector<timer::duration>> phase_time;
	vector<vector<timer::duration>> pure_phase_time;
	vector<timer::duration> net_route_time;
	vector<int> net_router;
	vector<dijkstra_stats_t> net_stats;
};

void do_work(router_t<net_t *> *router, SpeculativeDeterministicRouter &net_router)
{
#ifdef PTHREAD_BARRIER
	pthread_barrier_wait(&router->sync.barrier);
#else
	router->sync.barrier->wait();
#endif

	assert(net_router.get_num_instances() == router->num_threads);

	char buffer[256];

	sprintf(buffer, "%d", router->iter);
	zlog_put_mdc("iter", buffer);

	vector<timer::duration> net_route_time(router->nets.size(), timer::duration::max());
	vector<dijkstra_stats_t> net_stats(router->nets.size(),
			{ std::numeric_limits<int>::max(), std::numeric_limits<int>::max(), std::numeric_limits<int>::max() });
	vector<int> routed_vnets;

	auto pure_route_time = timer::duration::zero();
	auto wait_time = timer::duration::zero();
	auto update_time = timer::duration::zero();

	auto route_start = timer::now();

	int inet;

	while ((inet = router->sync.net_index++) < router->nets.size()) {
		net_t &net = *router->nets[inet];

		//router->sync.lock.lock();

		//assert(!router->sync.global_pending_nets.empty());
		//vnet = router->sync.global_pending_nets.front();
		//router->sync.global_pending_nets.pop();

		//router->sync.lock.unlock();

		auto real_net_route_start = timer::now();

		zlog_level(delta_log, ROUTER_V3, "Routing net %d num sinks %d\n", net.local_id, net.sinks.size());

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
				update_one_cost_internal(rr_node, router->g, router->state.congestion, -1, router->pres_fac);
			}

			router->sync.barrier->wait();
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

		routed_vnets.push_back(net.local_id);
	}

	auto last_barrier_wait_start = timer::now();
#ifdef PTHREAD_BARRIER 
	pthread_barrier_wait(&router->sync.barrier);
#else
	router->sync.barrier->wait();
#endif
	router->last_barrier_wait_time[net_router.get_instance()] = timer::now() - last_barrier_wait_start;

	router->route_time[net_router.get_instance()] = timer::now()-route_start;
	router->wait_time[net_router.get_instance()] = wait_time;
	router->pure_route_time[net_router.get_instance()] = pure_route_time; 
	router->update_time[net_router.get_instance()] = update_time; 

	router->stats[net_router.get_instance()].num_heap_pops = 0;
	router->stats[net_router.get_instance()].num_heap_pushes = 0;
	router->stats[net_router.get_instance()].num_neighbor_visits = 0;

	for (const auto &r : routed_vnets) {
		router->net_route_time[r] = net_route_time[r];
		router->net_router[r] = net_router.get_instance();
		router->net_stats[r] = net_stats[r];

		router->stats[net_router.get_instance()].num_heap_pops += net_stats[r].num_heap_pops;
		router->stats[net_router.get_instance()].num_heap_pushes += net_stats[r].num_heap_pushes;
		router->stats[net_router.get_instance()].num_neighbor_visits += net_stats[r].num_neighbor_visits;
	}

#ifdef PTHREAD_BARRIER 
	pthread_barrier_wait(&router->sync.barrier);
#else
	router->sync.barrier->wait();
#endif
}

static void *worker_thread(void *args)
{
	router_t<net_t *> *router = static_cast<router_t<net_t *> *>(args);

#ifdef PTHREAD_BARRIER
	pthread_barrier_wait(&router->sync.barrier);
#else
	router->sync.barrier->wait();
#endif

	SpeculativeDeterministicRouter net_router(router->g, router->state.congestion, router->pres_fac);

	char buffer[256];
	sprintf(buffer, "%d", net_router.get_instance());
	zlog_put_mdc("tid", buffer);
	zlog_put_mdc("log_dir", dirname);

	while (!router->sync.stop_routing) {
		do_work(router, net_router);
	}

	return nullptr;
}

void write_routes(const char *fname, const route_tree_t *route_trees, int num_nets);
void write_net_dijkstra_stats(const char *fname, const vector<dijkstra_stats_t> &net_stats);
void write_congestion_state(const char *fname, const congestion_t *congestion, int num_rr_nodes);

bool speculative_deterministic_route(t_router_opts *opts)
{
	char buffer[256];

	int dir_index;
	if (!opts->log_dir) {
		bool dir_exists;
		dir_index = 0;
		do {
			extern char *s_circuit_name;
			sprintf(dirname, "%s_stats_%d", s_circuit_name, dir_index);
			struct stat s;
			dir_exists = stat(dirname, &s) == 0 ? true : false;
			if (dir_exists) {
				++dir_index;
			}
		} while (dir_exists);
	} else {
		sprintf(dirname, opts->log_dir);

		struct stat s;
		if (stat(dirname, &s) == 0) {
			printf("log dir %s already exists\n", dirname);
			dir_index = -1;
		} else {
			dir_index = 0;
		}
	}

	if (dir_index >= 0) {
		if (mkdir(dirname, 0700) == -1) {
			printf("failed to create dir %s\n", dirname);
			exit(-1);
		}
	}

	init_logging(opts->num_threads);
    zlog_set_record("custom_output", concurrent_log_impl_2);
	sprintf(buffer, "%d", 0);
	zlog_put_mdc("iter", buffer);
	zlog_put_mdc("tid", buffer);
	zlog_put_mdc("log_dir", dirname);

	router_t<net_t *> router;

	//assert(((unsigned long)&router.sync.stop_routing & 63) == 0);
	//assert(((unsigned long)&router.sync.net_index & 63) == 0);
#ifdef __linux__
	//assert(((unsigned long)&router.sync.barrier & 63) == 0);
#endif

	router.num_threads = opts->num_threads;

	init_graph(router.g);

	init_sprintf_rr_node(&router.g);

	vector<net_t> nets;
	vector<net_t> global_nets;
	init_nets(nets, global_nets, opts->bb_factor, opts->large_bb);

	for (auto &net : nets) {
		router.nets.push_back(&net);
	}

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
	congestion_t *global = new congestion_t[num_vertices(router.g)];
    for (int i = 0; i < num_vertices(router.g); ++i) {
        global[i].acc_cost = 1;
        global[i].pres_cost = 1;
        global[i].occ = 0;
    }
	tbb::spin_mutex *locks = new tbb::spin_mutex[num_vertices(router.g)];

	router.state.congestion = new congestion_local_t[opts->num_threads];
	for (int i = 0; i < opts->num_threads; ++i) {
		//router.state.congestion[i].tid = i;
		router.state.congestion[i].global = global;
		router.state.congestion[i].locks = locks;
		router.state.congestion[i].local.resize(num_vertices(router.g), { 0, std::numeric_limits<float>::max() });
	}

	router.state.route_trees = new route_tree_t[nets.size()];
	//router.state.back_route_trees = new route_tree_t[nets.size()];
	for (int i = 0; i < nets.size(); ++i) {
		route_tree_init(router.state.route_trees[i]);
		//route_tree_init(router.state.back_route_trees[i]);
	}

	router.state.added_rr_nodes = new vector<RRNode>[nets.size()];
	router.state.back_added_rr_nodes = new vector<RRNode>[nets.size()];

    router.state.net_timing = new t_net_timing[nets.size()+global_nets.size()];
    init_net_timing(nets, global_nets, router.state.net_timing);

    router.params.criticality_exp = opts->criticality_exp;
    router.params.astar_fac = opts->astar_fac;
    router.params.max_criticality = opts->max_criticality;
    router.params.bend_cost = opts->bend_cost;

	//std::sort(begin(nets), end(nets), [] (const net_t &a, const net_t &b) -> bool {
			//return a.sinks.size() > b.sinks.size();
			//});

	router.pres_fac = opts->first_iter_pres_fac;

#ifdef PTHREAD_BARRIER
	assert(pthread_barrier_init(&router.sync.barrier, nullptr, opts->num_threads) == 0);
#else
	//router.sync.barrier = new boost::barrier(opts->num_threads);
	router.sync.barrier = new SpinningBarrier(opts->num_threads);
#endif
	router.sync.stop_routing = 0;

	//void *aligned_current_num_incoming_edges;
	//assert(posix_memalign(&aligned_current_num_incoming_edges, 64, sizeof(aligned_atomic<int, 64>)*nets.size()) == 0);
	//assert(((unsigned long)aligned_current_num_incoming_edges & 63) == 0);
	//router.sync.current_num_incoming_edges = new(aligned_current_num_incoming_edges) aligned_atomic<int, 64>[num_virtual_nets];
	//router.sync.current_num_incoming_edges = new std::atomic<int>[num_virtual_nets];

	router.stats.resize(opts->num_threads);
	router.route_time.resize(opts->num_threads);
	router.wait_time.resize(opts->num_threads);
	router.pure_route_time.resize(opts->num_threads);
	router.update_time.resize(opts->num_threads);
	router.last_barrier_wait_time.resize(opts->num_threads);
	router.pure_phase_time.resize(opts->num_threads);
	router.phase_time.resize(opts->num_threads);

	router.net_route_time.resize(nets.size());
	router.net_router.resize(nets.size());
	router.net_stats.resize(nets.size());

	router.iter = 0;

	vector<pthread_t> tids(opts->num_threads);
	for (int i = 1; i < opts->num_threads; ++i) {
		pthread_attr_t attr;
		assert(pthread_attr_init(&attr) == 0);

#ifdef __linux__
		cpu_set_t cpuset;
		CPU_ZERO(&cpuset);
		CPU_SET(4+i, &cpuset);
		assert(pthread_attr_setaffinity_np(&attr, sizeof(cpu_set_t), &cpuset) == 0);
#endif

		//assert(pthread_create(&tids[i], &attr, virtual_worker_thread, (void *)&router) == 0);
		assert(pthread_create(&tids[i], &attr, worker_thread, &router) == 0);
	}

#ifdef PTHREAD_BARRIER
	pthread_barrier_wait(&router.sync.barrier);
#else
	router.sync.barrier->wait();
#endif

	timer::duration total_route_time = timer::duration::zero();

	/* just to make sure that we are allocated on the correct NUMA node */
	SpeculativeDeterministicRouter *net_router = new SpeculativeDeterministicRouter(router.g, router.state.congestion, router.pres_fac);

	sprintf(buffer, "%d", net_router->get_instance());
	zlog_put_mdc("tid", buffer);

	bool routed = false;
	bool has_resched = false;
	for (; router.iter < opts->max_router_iterations && !routed; ++router.iter) {
		sprintf(buffer, "%d", router.iter);
		zlog_put_mdc("iter", buffer);

		printf("Routing iteration: %d\n", router.iter);

		//assert(router.sync.num_incoming_edges.empty());

		//router.sync.handles.resize(router.topo_nets.size());

		//for (const auto &vnet : router.topo_nets) {
			//router.sync.handles[vnet->global_index] = router.sync.num_incoming_edges.emplace(make_pair(router.num_incoming_edges[vnet->global_index], vnet));
		//}
		
		//assert(router.sync.fast_global_pending_nets.size_approx() > 0);

		router.sync.net_index = 0;

		for (int i = 0; i < nets.size(); ++i) {
			router.net_stats[i].num_heap_pops = std::numeric_limits<int>::max();
			router.net_stats[i].num_heap_pushes = std::numeric_limits<int>::max();
			router.net_stats[i].num_neighbor_visits = std::numeric_limits<int>::max();
		}

		for (int i = 0; i < opts->num_threads; ++i) {
			router.stats[i].num_heap_pops = 0;
			router.stats[i].num_heap_pushes = 0;
			router.stats[i].num_neighbor_visits = 0;
			router.route_time[i] = timer::duration::zero();
			router.wait_time[i] = timer::duration::zero();
			router.pure_route_time[i] = timer::duration::zero();
			router.update_time[i] = timer::duration::zero();
			router.last_barrier_wait_time[i] = timer::duration::zero();
		}

		for (auto &net : nets) {
			update_sink_criticalities(net, router.state.net_timing[net.vpr_id], router.params);
		}

		//do_virtual_work(&router, *net_router);
		do_work(&router, *net_router);

		/* checking */
		for (int i = 0; i < num_vertices(router.g); ++i) {
			global[i].recalc_occ = 0; 
		}

		for (const auto &net : nets) {
			check_route_tree(router.state.route_trees[net.local_id], net, router.g);
			recalculate_occ(router.state.route_trees[net.local_id], router.g, global);
		}

		//printf("Resched time [has_resched %d]: %g\n", has_resched ? 1 : 0, duration_cast<nanoseconds>(resched_time).count()/1e9);

		printf("Route time: ");
		for (int i = 0; i < opts->num_threads; ++i) {
			printf("%g ", duration_cast<nanoseconds>(router.route_time[i]).count()/1e9);
		}
		printf("\n");

		//printf("Route time debug: ");
		//for (auto iter = begin(router.route_time); iter != end(router.route_time); ++iter) {
			//printf("%g ", duration_cast<nanoseconds>(*iter).count()/1e9);
		//}
		//printf("\n");

		auto max_route_time_i = std::max_element(begin(router.route_time), end(router.route_time));
		assert(max_route_time_i != end(router.route_time));
		timer::duration max_route_time = *max_route_time_i;
		total_route_time += max_route_time;

		printf("Max route time: %lu %g\n", max_route_time.count(), max_route_time.count() / 1e9);

		printf("Wait time: ");
		for (int i = 0; i < opts->num_threads; ++i) {
			printf("%g (%g) ", duration_cast<nanoseconds>(router.wait_time[i]).count()/1e9, 100.0*router.wait_time[i].count()/router.route_time[i].count());
		}
		printf("\n");

		printf("Pure route time: ");
		for (int i = 0; i < opts->num_threads; ++i) {
			printf("%g (%g) ", duration_cast<nanoseconds>(router.pure_route_time[i]).count()/1e9, 100.0*router.pure_route_time[i].count()/router.route_time[i].count());
		}
		printf("\n");

		printf("Update time: ");
		for (int i = 0; i < opts->num_threads; ++i) {
			printf("%g (%g) ", duration_cast<nanoseconds>(router.update_time[i]).count()/1e9, 100.0*router.update_time[i].count()/router.route_time[i].count());
		}
		printf("\n");

		printf("Last barrier wait time: ");
		for (int i = 0; i < opts->num_threads; ++i) {
			printf("%g (%g) ", duration_cast<nanoseconds>(router.last_barrier_wait_time[i]).count()/1e9, 100.0*router.last_barrier_wait_time[i].count()/router.route_time[i].count());
		}
		printf("\n");

		unsigned long total_num_heap_pops = 0;
		printf("Num heap pops: ");
		for (int i = 0; i < opts->num_threads; ++i) {
			printf("%d ", router.stats[i].num_heap_pops);
			total_num_heap_pops += router.stats[i].num_heap_pops;
		}
		printf("\n");
		printf("Total num heap pops: %lu\n", total_num_heap_pops);

		unsigned long total_num_neighbor_visits = 0;
		printf("Num neighbor visits: ");
		for (int i = 0; i < opts->num_threads; ++i) {
			printf("%d ", router.stats[i].num_neighbor_visits);
			total_num_neighbor_visits += router.stats[i].num_neighbor_visits;
		}
		printf("\n");
		printf("Total num neighbor visits: %lu\n", total_num_neighbor_visits);

		unsigned long total_num_heap_pushes = 0;
		printf("Num heap pushes: ");
		for (int i = 0; i < opts->num_threads; ++i) {
			printf("%d ", router.stats[i].num_heap_pushes);
			total_num_heap_pushes += router.stats[i].num_heap_pushes;
		}
		printf("\n");
		printf("Total num heap pushes: %lu\n", total_num_heap_pushes);

		sprintf(buffer, "%s/stats_%d.txt", dirname, router.iter);
		FILE *stats = fopen(buffer, "w");

		for (int i = 0; i < router.nets.size(); ++i) {
			for (int j = 0; j < opts->num_threads; ++j) {
				fprintf(stats, "%g %g %g ", duration_cast<nanoseconds>(router.pure_phase_time[j][i]).count()/1e9, duration_cast<nanoseconds>(router.phase_time[j][i]).count()/1e9,
						100.0*router.pure_phase_time[j][i].count()/router.phase_time[j][i].count());
			}
			fprintf(stats, "\n");
		}

		fclose(stats);

		//sprintf(buffer, "%s/net_route_time_%d.txt", dirname, router.iter);
		//FILE *nrt = fopen(buffer, "w");

		//for (int i = 0; i < router.nets.size(); ++i) {
			//for (int j = 0; j < router.nets[i].size(); ++j) {
				//const auto *vnet = router.nets[i][j];
				//fprintf(nrt, "%d %d %d %g %d %d %d %g\n", i, router.nets[i].size(), router.net_router[vnet->global_index], duration_cast<nanoseconds>(router.net_route_time[vnet->global_index]).count()/1e9, vnet->index, vnet->net->sinks.size(), vnet->sinks.size(), bg::area(vnet->bounding_box));
			//}
		//}

		//fclose(nrt);

		unsigned long num_overused_nodes = 0;
		vector<int> overused_nodes_by_type(NUM_RR_TYPES, 0);

		if (feasible_routing(router.g, global)) {
			//dump_route(*current_traces_ptr, "route.txt");
			routed = true;
		} else {
			for (int i = 0; i < num_vertices(router.g); ++i) {
				if (global[i].occ > get_vertex_props(router.g, i).capacity) {
					++num_overused_nodes;
					const auto &v_p = get_vertex_props(router.g, i);
					++overused_nodes_by_type[v_p.type];
				}
			}

			//auto update_cost_start = timer::now();

			if (router.iter == 0) {
				router.pres_fac = opts->initial_pres_fac;
				update_costs(router.g, global, router.pres_fac, 0);
			} else {
				router.pres_fac *= opts->pres_fac_mult;

				/* Avoid overflow for high iteration counts, even if acc_cost is big */
				router.pres_fac = std::min(router.pres_fac, static_cast<float>(HUGE_POSITIVE_FLOAT / 1e5));

				update_costs(router.g, global, router.pres_fac, opts->acc_fac);
			}

			//update_cost_time = timer::now()-update_cost_start;
		}

		if (router.iter == 0) {
			sprintf(buffer, "%s/routes_iter_%d.txt", dirname, router.iter);
			write_routes(buffer, router.state.route_trees, nets.size());
		} else if (routed) {
			sprintf(buffer, "%s/routes_final.txt", dirname);
			write_routes(buffer, router.state.route_trees, nets.size());
		}

		if (routed) {
			sprintf(buffer, "%s/net_dijkstra_stats_final.txt", dirname);
			write_net_dijkstra_stats(buffer, router.net_stats);
			sprintf(buffer, "%s/congestion_state_final.txt", dirname);
			write_congestion_state(buffer, global, num_vertices(router.g));
		} else {
			sprintf(buffer, "%s/net_dijkstra_stats_%d.txt", dirname, router.iter);
			write_net_dijkstra_stats(buffer, router.net_stats);
			sprintf(buffer, "%s/congestion_state_%d.txt", dirname, router.iter);
			write_congestion_state(buffer, global, num_vertices(router.g));
		}

		//auto analyze_timing_start = timer::now();

		float crit_path_delay = analyze_timing(router.state.net_timing);

		//analyze_timing_time = timer::now()-analyze_timing_start;
		
		printf("Overused: %lu Crit path delay: %g\n", num_overused_nodes, crit_path_delay);
		printf("\n");

		//for (int i = 0; i < nets.size(); ++i) {
			//route_tree_clear(router.state.back_route_trees[i]);
		//}
		//std::swap(router.state.route_trees, router.state.back_route_trees);
		for (int i = 0; i < nets.size(); ++i) {
			router.state.back_added_rr_nodes[i].clear();
		}
		std::swap(router.state.added_rr_nodes, router.state.back_added_rr_nodes);
		for (int i = 0; i < nets.size(); ++i) {
			route_tree_clear(router.state.route_trees[i]);
		}
	}

	if (routed) {
		printf("Routed in %d iterations\n", router.iter);
		printf("Total route time: %g\n", duration_cast<nanoseconds>(total_route_time).count() / 1e9);
	}

	return false;
}
