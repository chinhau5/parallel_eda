#include "pch.h"

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
#include <boost/range/adaptor/transformed.hpp>

float get_timing_driven_expected_cost(const rr_node_property_t &current, const rr_node_property_t &target, float criticality_fac, float R_upstream);
void update_sink_criticalities(net_t &net, const t_net_timing &net_timing, const route_parameters_t &params);

using rtree_value = pair<box, net_t *>;

struct net_to_rtree_item {
	rtree_value operator()(net_t *net) const
	{
		return make_pair(box(point(net->bounding_box.xmin, net->bounding_box.ymin), point(net->bounding_box.xmax, net->bounding_box.ymax)), net);
	}
};

struct rtree_value_equal {
	bool operator()(const rtree_value &a, const rtree_value &b) const
	{
		return bg::equals(a.first, b.first) && a.second == b.second;
	}
};

void test_rtree(vector<net_t> &nets)
{
	using clock = std::chrono::high_resolution_clock;
	using std::chrono::duration_cast;
	using std::chrono::nanoseconds;

	vector<net_t *> sorted_nets;
	for (auto &net : nets) {
		sorted_nets.push_back(&net);
	}
	std::sort(begin(sorted_nets), end(sorted_nets), [] (const net_t *a, const net_t *b) -> bool {
			return a->sinks.size() > b->sinks.size();
			});

	auto rtree_start = clock::now();
	//bgi::rtree<rtree_value, bgi::rstar<16>> tree(nets | boost::adaptors::transformed([] (const net_t &net) -> pair<box, int> {
				//return make_pair(box(point(net.bounding_box.xmin, net.bounding_box.ymin), point(net.bounding_box.xmax, net.bounding_box.ymax)), net.local_id);
				//}
				//));

	net_to_rtree_item to_rtree_item;
	bgi::rtree<rtree_value, bgi::rstar<16>, bgi::indexable<rtree_value>, rtree_value_equal> tree(sorted_nets | boost::adaptors::transformed(to_rtree_item));

	auto rtree_time = clock::now()-rtree_start;

	printf("rtree build time %g\n", duration_cast<nanoseconds>(rtree_time).count() / 1e9);

	auto to_box = [] (const net_t &net) -> box { return box(point(net.bounding_box.xmin, net.bounding_box.ymin), point(net.bounding_box.xmax, net.bounding_box.ymax)); };

	//auto lol = tree.qbegin(bgi::satisfies([] (const rtree_value &val) -> bool { return true; }));
	vector<bool> net_scheduled(sorted_nets.size(), false);
	int inet = 0;
	int num_rounds = 0;
	int max_con = 0;
	auto sched_start = clock::now();
	while (!tree.empty()) {
		net_t *net = sorted_nets[inet];

		if (!net_scheduled[net->local_id]) {
			vector<rtree_value> disjoint_nets;
			vector<rtree_value> dispatched;

			dispatched.push_back(make_pair(to_box(*net), net));

			if (tree.query(bgi::disjoint(to_box(*net)), std::back_inserter(disjoint_nets)) > 0) {
				/* we could do a brute force search here instead of using MISR if sizeof disjoin_nets is small enough to consider all
				 * possible combinations */
				/* run graph coloring on disjoint_nets */
				/* for each item in disjoint_nets, find out their disjoint nets and do a set union or set covering */
				std::sort(begin(disjoint_nets), end(disjoint_nets), [] (const rtree_value &a, const rtree_value &b) -> bool {
						return bg::area(a.first) > bg::area(b.first);
						});

				for (int i = 0; i < disjoint_nets.size(); ++i) {
					const auto &cur_dis = disjoint_nets[i];
					bool intersect = std::any_of(begin(dispatched), end(dispatched), [&cur_dis] (const rtree_value &n) -> bool {
							return bg::intersects(cur_dis.first, n.first);
							});
					if (!intersect) {
						dispatched.push_back(cur_dis);
					}
				}


				//int max_area = -1;
				//int max = -1;
				//for (int i = 0; i < disjoint_nets.size(); ++i) {
					//int area = bg::area(to_box(disjoint_nets[i]));
					//if (area > max_area) {
						//max_area = area;
						//max = i;
					//}
				//}

				//assert(max >= 0);
				//tree.remove(make_pair(to_box(disjoint_nets[max]), disjoint_nets[max].local_id));
			}

			assert(tree.remove(begin(dispatched), end(dispatched)) == dispatched.size());

			//printf("Num dispatched: %d\n", dispatched.size());
			verify_ind(dispatched, [] (const rtree_value &a) -> box { return a.first; });

			for (const auto &d : dispatched) {
				assert(!net_scheduled[d.second->local_id]);
				net_scheduled[d.second->local_id] = true;
			}

			++num_rounds;
			max_con = std::max(max_con, (int)dispatched.size());
		}

		++inet;
	}
	auto sched_time = clock::now() - sched_start;

	printf("sched_time: %g max_con: %d ave_con: %g\n", duration_cast<nanoseconds>(sched_time).count()/1e9, max_con, sorted_nets.size()*1.0/num_rounds);

	assert(std::all_of(begin(net_scheduled), end(net_scheduled), [] (bool a) -> bool { return a; }));
}

void test_misr(const vector<net_t> &nets)
{
	set<int> all_nets;

	auto to_box = [] (const net_t &net) -> box { return box(point(net.bounding_box.xmin, net.bounding_box.ymin), point(net.bounding_box.xmax, net.bounding_box.ymax)); };

	int subiter = 0;
	int max_con = 0;

	vector<net_t> nets_copy = nets;

	while (!nets_copy.empty()) {
		vector<net_t *> chosen;

		max_independent_rectangles(nets_copy, to_box, [] (const net_t &net) -> pair<int, int> { return make_pair(net.bounding_box.ymin, net.bounding_box.ymax); }, chosen);

		verify_ind(chosen, [&to_box] (const net_t *net) -> box { return to_box(*net); });

		extern int nx, ny;

		FILE *file;
		if (subiter < 10) {
			char buffer[256];
			sprintf(buffer, "%d_boxes.txt", subiter);
			file = fopen(buffer, "w");
			fprintf(file, "0 0 %d %d 0\n", nx+2, ny+2);
		}

		int total_area = 0;
		for (const auto &c : chosen) {
			assert(all_nets.find(c->vpr_id) == all_nets.end());
			all_nets.insert(c->vpr_id);

			total_area += bg::area(to_box(*c));

			if (subiter < 10) {
				fprintf(file, "%d %d %d %d 0\n", c->bounding_box.xmin, c->bounding_box.ymin, c->bounding_box.xmax-c->bounding_box.xmin, c->bounding_box.ymax-c->bounding_box.ymin);
			}
		}

		vector<net_t *> new_chosen;
		for (auto &net : nets_copy) {
			const auto &box = to_box(net);
			bool dis = std::all_of(begin(chosen), end(chosen), [&to_box, &box] (const net_t *c) -> bool {
						return bg::disjoint(to_box(*c), box);
					}) && std::all_of(begin(new_chosen), end(new_chosen), [&to_box, &box] (const net_t *c) -> bool {
						return bg::disjoint(to_box(*c), box);
					});
			if (dis) {
				if (subiter < 10) {
					fprintf(file, "%d %d %d %d 1\n", net.bounding_box.xmin, net.bounding_box.ymin, net.bounding_box.xmax-net.bounding_box.xmin, net.bounding_box.ymax-net.bounding_box.ymin);
				}
				new_chosen.push_back(&net);
			}
		}

		if (subiter < 10) {
			fclose(file);
		}

		++subiter;

		printf("chosen size: %d %d %g\n", chosen.size(), total_area, 100.0*total_area/((nx+2)*(ny+2)));
		max_con = std::max((int)chosen.size(), max_con);

		nets_copy.erase(std::remove_if(begin(nets_copy), end(nets_copy), [&chosen] (const net_t &other) -> bool {
					return std::any_of(begin(chosen), end(chosen), [&other] (const net_t *c) -> bool {
							return c->vpr_id == other.vpr_id;
							});
					}), end(nets_copy));
	}

	printf("max_con: %d ave_con: %g\n", max_con, nets.size()*1.0/subiter);
}

template<typename EdgeProperties>
bool operator<(const existing_source_t<EdgeProperties> &a, const existing_source_t<EdgeProperties> &b)
{
	return a.distance > b.distance;
}

template<typename VertexProperties, typename EdgeProperties, typename EdgeWeightFunc, typename Callbacks>
void dijkstra(fast_graph_t<VertexProperties, EdgeProperties> &g, const vector<existing_source_t<EdgeProperties>> &sources, int sink, float *known_distance, float *distance, fast_edge_t<EdgeProperties> *prev_edge, const EdgeWeightFunc &edge_weight, Callbacks &callbacks)
{
	using Item = existing_source_t<EdgeProperties>;
	std::priority_queue<Item> heap;

	for (const auto &s : sources) {
		heap.push(s);
	}

	bool found = false;
	while (!heap.empty() && !found) {
		Item item = heap.top();
		heap.pop();

		zlog_level(delta_log, ROUTER_V3, "Current: %d [kd=%g okd=%g] [d=%g od=%g] prev=%d\n",
					item.node, item.known_distance, known_distance[item.node], item.distance, distance[item.node], get_source(g, item.prev_edge));
		callbacks.popped_node(item.node);

		//if (!callbacks.expand_node(item.node)) {
			//continue;
		//}

		if (item.distance < distance[item.node]) {
			assert(item.known_distance <= known_distance[item.node]);

			known_distance[item.node] = item.known_distance;
			distance[item.node] = item.distance;
			prev_edge[item.node] = item.prev_edge;

			zlog_level(delta_log, ROUTER_V3, "Relaxing %d\n", item.node);
			callbacks.relax_node(item.node, item.prev_edge);

			if (item.node != sink) {
				for (const auto &e : get_out_edges(g, item.node)) {
					int v = get_target(g, e);

					zlog_level(delta_log, ROUTER_V3, "\tNeighbor: %d\n", v);

					if (!callbacks.expand_node(v)) {
						continue;
					}

					const auto &weight = edge_weight(e);
					float kd = known_distance[item.node] + weight.first;
					float d = known_distance[item.node] + weight.second;

					zlog_level(delta_log, ROUTER_V3, "\t[w1 %X w2 %X] [kd=%X okd=%X] [d=%X od=%X] [kd=%g okd=%g] [d=%g od=%g]\n",
							*(unsigned int *)&weight.first, *(unsigned int *)&weight.second, *(unsigned int *)&kd, *(unsigned int *)&known_distance[v], *(unsigned int *)&d, *(unsigned int *)&distance[v], kd, known_distance[v], d, distance[v]);

					if (d < distance[v]) {
						assert(kd <= known_distance[v]);

						zlog_level(delta_log, ROUTER_V3, "\tPushing neighbor %d to heap\n", v);

						heap.push({ v, kd, d, e });
					}
				}
			} else {
				found = true;
			}
		}
	}

	assert(found);
}

typedef struct extra_route_state_t {
	float upstream_R;
	float delay;
} extra_route_state_t;

class DeltaSteppingRouter {
	private:
		RRGraph &_g;

		float _astar_fac;

		float *_known_distance;
		float *_distance;
		RREdge *_prev_edge;
		extra_route_state_t *_state;

		congestion_t *_congestion;
		float *_pres_fac;

		vector<int> _modified_nodes;
		vector<bool> _modified_node_added;

		const sink_t *_current_sink;
		RRNode _existing_opin;

		route_tree_t *_current_rt;

		int _num_heap_pops;

	private:
		void popped_node(int v)
		{
			char buffer[256];
			sprintf_rr_node(v, buffer);
			zlog_level(delta_log, ROUTER_V3, "%s\n", buffer);
			++_num_heap_pops;
		}

		void relax_node(int v, const RREdge &e)
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

			const auto &v_p = get_vertex_props(_g, v);

			if (valid(e)) {
				assert(v == get_target(_g, e));

				int u = get_source(_g, e);

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

				const auto &e_p = get_edge_props(_g, e);

				_state[v].upstream_R = e_p.R + v_p.R;
				if (!e_p.buffered)  {
					_state[v].upstream_R += u_upstream_R;
				} 

				float delay;
				if (e_p.buffered) {
					delay = e_p.switch_delay + v_p.C * (e_p.R + 0.5f * v_p.R);
				} else {
					delay = e_p.switch_delay + v_p.C * (u_upstream_R + e_p.R + 0.5f * v_p.R);
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

			if (_existing_opin != RRGraph::null_vertex() && prop.type == OPIN && node != _existing_opin) {
				zlog_level(delta_log, ROUTER_V3, "not expanding other OPIN\n");
				return false;
			}

			if (prop.xhigh < _current_sink->current_bounding_box.xmin
					|| prop.xlow > _current_sink->current_bounding_box.xmax
					|| prop.yhigh < _current_sink->current_bounding_box.ymin
					|| prop.ylow > _current_sink->current_bounding_box.ymax) {
				zlog_level(delta_log, ROUTER_V3, "outside of bounding box\n");
				return false;
			}

			const auto &sink_prop = get_vertex_props(_g, _current_sink->rr_node);

			if (prop.type == IPIN
					&& (prop.xhigh != sink_prop.xhigh 
						|| prop.yhigh != sink_prop.yhigh)) {
				zlog_level(delta_log, ROUTER_V3, "not target IPIN\n");
				return false;
			}

			zlog_level(delta_log, ROUTER_V3, "\n");

			return true;
		}

		pair<float, float> get_edge_weight(const RREdge &e)
		{
			int u = get_source(_g, e);
			int v = get_target(_g, e);

			const auto &v_p = get_vertex_props(_g, v);
			const auto &e_p = get_edge_props(_g, e);

			assert(_state[u].upstream_R != std::numeric_limits<float>::max());
			float delay;
			if (e_p.buffered) {
				delay = e_p.switch_delay + v_p.C * (e_p.R + 0.5f * v_p.R);
			} else {
				delay = e_p.switch_delay + v_p.C * (_state[u].upstream_R + e_p.R + 0.5f * v_p.R);
			}
			extern t_rr_indexed_data *rr_indexed_data;
			float congestion_cost = rr_indexed_data[v_p.cost_index].base_cost * _congestion[v].acc_cost * _congestion[v].pres_cost;
			float known_cost = _current_sink->criticality_fac * delay + (1 - _current_sink->criticality_fac) * congestion_cost;

			float upstream_R = e_p.R + v_p.R;
			if (!e_p.buffered) {
				upstream_R += _state[u].upstream_R;
			}
			float expected_cost = get_timing_driven_expected_cost(v_p, get_vertex_props(_g, _current_sink->rr_node), _current_sink->criticality_fac, upstream_R);

			zlog_level(delta_log, ROUTER_V3, "\t%d -> %d delay %g congestion %g crit_fac %g expected %g expected_hex %X known %g predicted %g\n", 
					u, v, delay, congestion_cost, _current_sink->criticality_fac, expected_cost, *(unsigned int *)&expected_cost, known_cost, known_cost + _astar_fac * expected_cost);
			zlog_level(delta_log, ROUTER_V3, "\t[u: upstream %g] [edge: d %g R %g] [v: R %g C %g]\n",
					_state[u].upstream_R, e_p.switch_delay, e_p.R, v_p.R, v_p.C);

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

		void backtrack(int sink_node, route_tree_t &rt)
		{
			char buffer[256];
			sprintf_rr_node(sink_node, buffer);
			zlog_level(delta_log, ROUTER_V3, "Adding %s to route tree\n", buffer);

			RouteTreeNode rt_node = route_tree_add_rr_node(rt, sink_node, _g);
			assert(rt_node != RouteTree::null_vertex());
			assert(_state[sink_node].upstream_R != std::numeric_limits<float>::max() && _state[sink_node].delay != std::numeric_limits<float>::max());
			route_tree_set_node_properties(rt, rt_node, false, _state[sink_node].upstream_R, _state[sink_node].delay);
			update_one_cost_internal(sink_node, _g, _congestion, 1, *_pres_fac);

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
					update_one_cost_internal(parent_rr_node, _g, _congestion, 1, *_pres_fac);
				}

				zlog_level(delta_log, ROUTER_V3, "\n");

				child_rr_node = parent_rr_node;
				edge = get_previous_edge(parent_rr_node, rt);
			}
		}

		template<typename VertexProperties, typename EdgeProperties, typename EdgeWeightFunc, typename Callbacks>
		friend void delta_stepping(fast_graph_t<VertexProperties, EdgeProperties> &g, const vector<existing_source_t<EdgeProperties>> &sources, int sink, float delta, float *known_distance, float *distance, fast_edge_t<EdgeProperties> *prev_edge, const EdgeWeightFunc &edge_weight, Callbacks &callbacks);

		template<typename VertexProperties, typename EdgeProperties, typename EdgeWeightFunc, typename Callbacks>
		friend void dijkstra(fast_graph_t<VertexProperties, EdgeProperties> &g, const vector<existing_source_t<EdgeProperties>> &sources, int sink, float *known_distance, float *distance, fast_edge_t<EdgeProperties> *prev_edge, const EdgeWeightFunc &edge_weight, Callbacks &callbacks);

		template<typename EdgeProperties, typename Callbacks>
		friend void relax(Buckets &buckets, float delta, vector<bool> &in_bucket, const vector<bool> &vertex_deleted,
					float *known_distance, float *distance, fast_edge_t<EdgeProperties> *predecessor,
					int v, float new_known_distance, float new_distance, const fast_edge_t<EdgeProperties> &edge,
					Callbacks &callbacks);

	public:
		DeltaSteppingRouter(RRGraph &g, congestion_t *congestion, float *pres_fac)
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
		}

		void route(const source_t *source, const vector<const sink_t *> &sinks, int net_local_id, float delta, float astar_fac, route_tree_t &rt, t_net_timing &net_timing)
		{
			char buffer[256];

			_astar_fac = astar_fac;

			_existing_opin = RRGraph::null_vertex();

			const auto &source_p = get_vertex_props(_g, source->rr_node);
			RouteTreeNode root_rt_node = route_tree_add_rr_node(rt, source->rr_node, _g);
			route_tree_set_node_properties(rt, root_rt_node, true, source_p.R, 0.5f * source_p.R * source_p.C);
			route_tree_add_root(rt, source->rr_node);

			vector<const sink_t *> sorted_sinks = sinks;

			//std::sort(begin(sorted_sinks), end(sorted_sinks), [] (const sink_t *a, const sink_t *b) -> bool {
					//return a->criticality_fac > b->criticality_fac;
					//});

			_current_rt = &rt;

			int isink = 0;
			for (const auto &sink : sorted_sinks) {
				_current_sink = sink;

				sprintf_rr_node(sink->rr_node, buffer);
				zlog_level(delta_log, ROUTER_V3, "Current sink: %s\n", buffer);

				vector<existing_source_t<rr_edge_property_t>> sources;

				for (const auto &rt_node : route_tree_get_nodes(rt)) {
					const auto &rt_node_p = get_vertex_props(rt.graph, rt_node);
					RRNode node = rt_node_p.rr_node;
					const auto &node_p = get_vertex_props(_g, node);

					sprintf_rr_node(node, buffer);

					if (rt_node_p.reexpand) {
						if (node_p.xhigh < sink->current_bounding_box.xmin
								|| node_p.xlow > sink->current_bounding_box.xmax
								|| node_p.yhigh < sink->current_bounding_box.ymin
								|| node_p.ylow > sink->current_bounding_box.ymax) {
							zlog_level(delta_log, ROUTER_V3, "Existing %s out of bounding box\n", buffer);

						} else {
							float kd = sink->criticality_fac * rt_node_p.delay; 
							float d = kd + _astar_fac * get_timing_driven_expected_cost(node_p, get_vertex_props(_g, sink->rr_node), sink->criticality_fac, rt_node_p.upstream_R); 

							zlog_level(delta_log, ROUTER_V3, "Adding %s back to heap [delay %g upstream_R %g] [kd=%g d=%g]\n", buffer, rt_node_p.delay, rt_node_p.upstream_R, kd, d);

							RREdge prev = RRGraph::null_edge();
							if (valid(rt_node_p.rt_edge_to_parent)) {
								prev = get_edge_props(rt.graph, rt_node_p.rt_edge_to_parent).rr_edge;
							}
							sources.push_back({ node, kd, d, prev });
						}
					} else {
						zlog_level(delta_log, ROUTER_V3, "Not reexpanding %s\n", buffer);
					}
				}

				_num_heap_pops = 0;
				//delta_stepping(_g, sources, sink->rr_node, delta, _known_distance, _distance, _prev_edge, [this] (const RREdge &e) -> pair<float, float> { return get_edge_weight(e); }, *this);
				dijkstra(_g, sources, sink->rr_node, _known_distance, _distance, _prev_edge, [this] (const RREdge &e) -> pair<float, float> { return get_edge_weight(e); }, *this);

				zlog_level(delta_log, ROUTER_V3, "num_heap_pop %d\n", _num_heap_pops);

				backtrack(sink->rr_node, rt);

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

				for (const auto &n : get_vertices(_g)) {
					assert(_known_distance[n] == std::numeric_limits<float>::max() &&
							_distance[n] == std::numeric_limits<float>::max() &&
							_prev_edge[n] == RRGraph::null_edge() &&
							_state[n].delay == std::numeric_limits<float>::max() &&
							_state[n].upstream_R == std::numeric_limits<float>::max());
				}

				++isink;
			}
		}
};

bool partitioning_delta_stepping_deterministic_route(t_router_opts *opts)
{
	init_logging();
    zlog_set_record("custom_output", concurrent_log_impl);

	zlog_put_mdc("iter", "0");
	zlog_put_mdc("tid", "0");

	RRGraph g;
	init_graph(g);

	init_sprintf_rr_node(&g);

	vector<net_t> nets;
	vector<net_t> global_nets;
	init_nets(nets, global_nets, opts->bb_factor, opts->large_bb);

	test_rtree(nets);
	test_misr(nets);
	return false;

	congestion_t *congestion = new congestion_t[num_vertices(g)];
    for (int i = 0; i < num_vertices(g); ++i) {
        congestion[i].acc_cost = 1;
        congestion[i].pres_cost = 1;
        congestion[i].occ = 0;
    }

	route_tree_t *route_trees = new route_tree_t[nets.size()];
	for (int i = 0; i < nets.size(); ++i) {
		route_tree_init(route_trees[i]);
	}

    t_net_timing *net_timing = new t_net_timing[nets.size()+global_nets.size()];
    init_net_timing(nets, global_nets, net_timing);

    route_parameters_t params;
    params.criticality_exp = opts->criticality_exp;
    params.astar_fac = opts->astar_fac;
    params.max_criticality = opts->max_criticality;
    params.bend_cost = opts->bend_cost;

	//std::sort(begin(nets), end(nets), [] (const net_t &a, const net_t &b) -> bool {
			//return a.sinks.size() > b.sinks.size();
			//});

	float pres_fac = opts->first_iter_pres_fac;

	float delta;
	for (const auto &n : get_vertices(g)) {
		const auto &prop = get_vertex_props(g, n);
		if (prop.type == CHANX) {
			extern t_rr_indexed_data *rr_indexed_data;
			delta = rr_indexed_data[prop.cost_index].base_cost;
			break;
		}
	}

	DeltaSteppingRouter router(g, congestion, &pres_fac);

	bool routed = false;
	for (int iter = 0; iter < opts->max_router_iterations && !routed; ++iter) {
		int inet = 0;
		for (auto &net : nets) {
			update_sink_criticalities(net, net_timing[net.vpr_id], params);

			route_tree_mark_all_nodes_to_be_ripped(route_trees[net.local_id], g);

			route_tree_rip_up_marked(route_trees[net.local_id], g, congestion, pres_fac);

			vector<const sink_t *> sinks;
			for (int i = 0; i < net.sinks.size(); ++i) {
				sinks.push_back(&net.sinks[i]);
			}

			zlog_level(delta_log, ROUTER_V3, "Routing net %d\n", net.local_id);
			router.route(&net.source, sinks, net.local_id, delta, opts->astar_fac, route_trees[net.local_id], net_timing[net.vpr_id]); 

			++inet;
		}

		/* checking */
		for (int i = 0; i < num_vertices(g); ++i) {
			congestion[i].recalc_occ = 0; 
		}

		for (const auto &net : nets) {
			check_route_tree(route_trees[net.local_id], net, g);
			recalculate_occ(route_trees[net.local_id], g, congestion);
		}

		unsigned long num_overused_nodes = 0;
		vector<int> overused_nodes_by_type(NUM_RR_TYPES, 0);

		if (feasible_routing(g, congestion)) {
			//dump_route(*current_traces_ptr, "route.txt");
			routed = true;
		} else {
			for (int i = 0; i < num_vertices(g); ++i) {
				if (congestion[i].occ > get_vertex_props(g, i).capacity) {
					++num_overused_nodes;
					const auto &v_p = get_vertex_props(g, i);
					++overused_nodes_by_type[v_p.type];
				}
			}

			//auto update_cost_start = clock::now();

			if (iter == 0) {
				pres_fac = opts->initial_pres_fac;
				update_costs(g, congestion, pres_fac, 0);
			} else {
				pres_fac *= opts->pres_fac_mult;

				/* Avoid overflow for high iteration counts, even if acc_cost is big */
				pres_fac = std::min(pres_fac, static_cast<float>(HUGE_POSITIVE_FLOAT / 1e5));

				update_costs(g, congestion, pres_fac, opts->acc_fac);
			}

			//update_cost_time = clock::now()-update_cost_start;
		}

		//auto analyze_timing_start = clock::now();

		float crit_path_delay = analyze_timing(net_timing);

		//analyze_timing_time = clock::now()-analyze_timing_start;
		
		printf("Overused: %lu Crit path delay: %g\n", num_overused_nodes, crit_path_delay);
	}

	return false;
}
