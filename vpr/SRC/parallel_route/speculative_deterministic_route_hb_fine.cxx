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
#include "det_mutex.h"

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
	congestion_t **congestion;
	route_tree_t *route_trees;
	//route_tree_t *back_route_trees;
	vector<RRNode> *added_rr_nodes;	
	vector<RRNode> *back_added_rr_nodes;	
	t_net_timing *net_timing;
} global_route_state_t;

typedef struct worker_sync_t {
	/*alignas(64)*/ std::atomic<int> stop_routing;
	/*alignas(64)*/ std::atomic<int> net_index;
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

struct runtime_t {
	timer::duration route_time;
	timer::duration wait_time;
	timer::duration pure_route_time;
	timer::duration last_barrier_wait_time;
	timer::duration last_sync_time;
};

struct mod_t {
	int node;
	int delta;
};

using Mods = vector<mod_t> *;

struct region_t {
	int gid;
	box bb;
	det_mutex_t mutex;
	//vector<vector<mod_t>>> mods; [> [tid][mod_id] <]
	vector<vector<Mods>, tbb::cache_aligned_allocator<vector<Mods>>> mods; /* [tid][mod_id] */
	vector<vector<int>, tbb::cache_aligned_allocator<vector<int>>> num_committed_mods; /* [tid][other_tid] */
};

template<typename Net>
struct /*alignas(64)*/ router_t {
	vector<Net> nets;
	vector<region_t *> fpga_regions;
	vector<region_t *> rr_regions;
	std::atomic<int> iter;
	worker_sync_t sync;
	RRGraph g;
	global_route_state_t state;
	route_parameters_t params;
	float pres_fac;

	vector<exec_state_t, tbb::cache_aligned_allocator<exec_state_t>> e_state;
	int num_threads;
	int window_size;

	/* per thread */
	vector<runtime_t> time;
	vector<dijkstra_stats_t> stats;

	/* per nets */
	vector<timer::duration> net_route_time;
	vector<dijkstra_stats_t> net_stats;
	vector<int> net_router;
};

void rand_delay();

thread_local int tl_tid;

static char dirname[256];

static bool inside_bb(const rr_node_property_t &node, const box &bb)
{
	bool inside = false;

	if (node.type == CHANX || node.type == CHANY) {
		int x, y;
		if (node.inc_direction) {
			x = node.xlow;
			y = node.ylow;
		} else {
			x = node.xhigh;
			y = node.yhigh;
		}

		if (x >= bg::get<bg::min_corner, 0>(bb)
				&& x <= bg::get<bg::max_corner, 0>(bb)
				&& y >= bg::get<bg::min_corner, 1>(bb)
				&& y <= bg::get<bg::max_corner, 1>(bb)) {
			inside = true;
		}
	} else {
		/* partially inside the bounding box */
		if (node.xhigh >= bg::get<bg::min_corner, 0>(bb)
				&& node.xlow <= bg::get<bg::max_corner, 0>(bb)
				&& node.yhigh >= bg::get<bg::min_corner, 1>(bb)
				&& node.ylow <= bg::get<bg::max_corner, 1>(bb)) {
			inside = true;
		}
		/* fully inside */
		//if (node.xlow >= bg::get<bg::min_corner, 0>(bb)
				//&& node.xhigh <= bg::get<bg::max_corner, 0>(bb)
				//&& node.ylow >= bg::get<bg::min_corner, 1>(bb)
				//&& node.yhigh <= bg::get<bg::max_corner, 1>(bb)) {
			//inside = true;
		//}
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

static router_t<net_t *> *g_router;

void region_acquire(region_t &o, int tid)
{
	det_mutex_lock(o.mutex, tid);
}

void region_release(region_t &o, int tid)
{
	det_mutex_unlock(o.mutex, tid);
}

//void region_add_mod(region_t &o, int rr_node, int delta, int tid)
//{
	//region_acquire(o, tid);

	//mod_t m;
	//m.node = rr_node;
	//m.delta = rr_node;
	//o.mods[tid].emplace_back(m);

	//region_release(o, tid);
//}

void region_add_mod(region_t &o, int tid, vector<mod_t> *mods)
{
	region_acquire(o, tid);

	zlog_level(delta_log, ROUTER_V3, "[Net %d] Adding mods (size %d) to region %d (old size %d)\n", g_router->e_state[tid].inet->load(), mods->size(), o.gid, o.mods[tid].size());

	o.mods[tid].emplace_back(mods);

	region_release(o, tid);
}

//template<typename Callback>
//void region_get_mods(region_t &o, int tid, const Callback &f)
//{
	//region_acquire(o, tid);

	//for (int i = 0; i < o.mods.size(); ++i) {
		//if (i == tid) {
			//continue;
		//}

		//assert(o.num_committed_mods[tid][i] <= o.mods[i].size());

		//for (int j = o.num_committed_mods[tid][i]; j < o.mods[i].size(); ++j) {
			//f(o.mods[i][j]);
		//}

		//o.num_committed_mods[tid][i] = o.mods[i].size();
	//}

	//region_release(o, tid);
//}

template<typename Callback>
void region_get_mods(region_t &o, int tid, const Callback &f)
{
	region_acquire(o, tid);

	zlog_level(delta_log, ROUTER_V3, "[Net %d] Thread %d Getting mods in region %d\n", g_router->e_state[tid].inet->load(), tl_tid, o.gid);

	for (int i = 0; i < o.mods.size(); ++i) {
		if (i == tid) {
			zlog_level(delta_log, ROUTER_V3, "Skipping mods from thread %d\n", i);
			continue;
		}

		assert(o.num_committed_mods[tid][i] <= o.mods[i].size());

		zlog_level(delta_log, ROUTER_V3, "[Net %d region %d] Mods from thread %d. Commited: %d Size: %d\n", g_router->e_state[tid].inet->load(), o.gid, i, o.num_committed_mods[tid][i], o.mods[i].size());

		for (int j = o.num_committed_mods[tid][i]; j < o.mods[i].size(); ++j) {
			zlog_level(delta_log, ROUTER_V3, "Mods %d size %d\n", j, o.mods[i][j]->size());
			for (int k = 0; k < o.mods[i][j]->size(); ++k) {
				f((*o.mods[i][j])[k]);
			}
		}

		o.num_committed_mods[tid][i] = o.mods[i].size();
	}

	region_release(o, tid);
}

void region_clear_mods(region_t &o)
{
	for (int i = 0; i < o.mods.size(); ++i) {
		for (int j = 0; j < o.mods[i].size(); ++j) {
			delete o.mods[i][j];
		}
		o.mods[i].clear();
	}

	for (int i = 0; i < o.num_committed_mods.size(); ++i) {
		for (int j = 0; j < o.num_committed_mods[i].size(); ++j) {
			o.num_committed_mods[i][j] = 0;
		}
	}
}

void mods_add(int rr_node, int delta, vector<region_t *> &rr_regions, vector<Mods> &mods)
{
	mod_t m;
	m.node = rr_node;
	m.delta = delta;

	region_t *r = rr_regions[rr_node];

	assert(r->gid >= 0 && r->gid < mods.size());

	if (!mods[r->gid]) {
		mods[r->gid] = new vector<mod_t>;
	}

	mods[r->gid]->emplace_back(m);
}

void mods_commit(vector<Mods> &mods, int tid, vector<region_t *> &fpga_regions)
{
	for (int i = 0; i < mods.size(); ++i) {
		if (mods[i]) {
			region_add_mod(*fpga_regions[i], tid, mods[i]);
		}
	}
}

//template<typename ToSync, typename DoneSync>
//void sync(vector<region_t *> &regions, const ToSync &to_sync, const DoneSync &done_sync, int tid, const RRGraph &g, congestion_t *congestion, float pres_fac)
//{
	//for (int i = 0; i < regions.size(); ++i) {
		//region_t *o = regions[i];

		//assert(o);

		//if (to_sync(i, o)) {
			//int num_updates = 0;
			//region_get_mods(*o, tid, [&g, &congestion, &pres_fac, &num_updates] (const mod_t &m) -> void {
					//const auto &props = get_vertex_props(g, m.node);
					//update_first_order_congestion(congestion, m.node, m.delta, props.capacity, pres_fac);
					//++num_updates;
					//});
			////region_set_commited(*o, _instance);
			////_acquired_regions[o->id] = true;
			//done_sync(i, o);
			
			//zlog_level(delta_log, ROUTER_V3, "Num updates in region %d = %d\n", o->gid, num_updates);
		//} else {
			//zlog_level(delta_log, ROUTER_V3, "Not getting mods in region %d\n", o->gid);
		//}
	//}
//}

//void region_set_commited(region_t &o, int tid)
//{
	//region_acquire(o, tid);

	//for (int i = 0; i < o.mods.size(); ++i) {
		//if (i == tid) {
			//continue;
		//}

		//o.num_committed_mods[tid][i] = o.mods[i].size();
	//}

	//region_release(o, tid);
//}

//void commit(vector<region_t *> &regions, const RRGraph &g, int node, int delta, int tid)
//{
	//[>TODO: cache the in_region checking<]
	//for (auto &o : regions) {
		//if (o == nullptr) {
			//continue;
		//}

		//assert(o->id >= 0);

		//if (inside_bb(get_vertex_props(g, node), o->bb)) {
			//region_add_mod(*o, node, delta, tid);
		//}
	//}
//}

class SpeculativeDeterministicRouter {
	private:
		RRGraph &_g;
		vector<region_t *> &_fpga_regions;
		vector<region_t *> &_rr_regions;

		float _astar_fac;

		float *_known_distance;
		float *_distance;
		RREdge *_prev_edge;
		extra_route_state_t *_state;

		congestion_t *_congestion;
		const float &_pres_fac;

		exec_state_t *_e_state;

		vector<int> _modified_nodes;
		vector<bool> _modified_node_added;

		RRNode _existing_opin;

		route_tree_t *_current_rt;

		vector<bool> _acquired_regions;

		const sink_t *_current_sink;

		dijkstra_stats_t _stats;

		static std::atomic<int> _num_instances;
		int _instance;

	private:
		void popped_node(int v)
		{
			char buffer[256];
			sprintf_rr_node(v, buffer);
			zlog_level(delta_log, ROUTER_V3, "%s\n", buffer);
			++_stats.num_heap_pops;

			//rand_delay();
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
			++(*_e_state->logical_clock);
			
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

			if (prop.type == OPIN && _existing_opin != RRGraph::null_vertex() && node != _existing_opin) {
				zlog_level(delta_log, ROUTER_V3, "not expanding other OPIN\n");
				return false;
			}

			if (!inside_bb(prop, _current_sink->bounding_box)) {
				zlog_level(delta_log, ROUTER_V3, "outside of bounding box\n");
				return false;
			}

			zlog_level(delta_log, ROUTER_V3, " OK\n");

			++_stats.num_neighbor_visits;

			return true;
		}

		void push_node(int node)
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

			region_t *r = _rr_regions[v];

			if (!_acquired_regions[r->gid]) {
				region_get_mods(*r, _instance, [this] (const mod_t &m) -> void {
						const auto &props = get_vertex_props(_g, m.node);
						update_first_order_congestion(_congestion, m.node, m.delta, props.capacity, _pres_fac);
						});
				_acquired_regions[r->gid] = true;
			}

			//sync(_rr_regions,
					//[this, &v_p] (int id, const region_t *o) -> bool { return !_acquired_regions[id];  },
					////[this, &v_p] (int id, const region_t *o) -> bool { return true;  },
					//[this] (int id, const region_t *o) -> void { _acquired_regions[id] = true; },
					//_instance, _g, _congestion, _pres_fac);

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
			assert(inside_bb(sink_node_p, _current_sink->bounding_box));

			vector<Mods> mods(_fpga_regions.size(), nullptr);
			update_first_order_congestion(_congestion, sink_node, 1, sink_node_p.capacity, _pres_fac);
			//commit(*_current_regions, _g, sink_node, 1, _instance);
			mods_add(sink_node, 1, _rr_regions, mods);

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
					assert(inside_bb(rr_node_p, _current_sink->bounding_box));

					update_first_order_congestion(_congestion, parent_rr_node, 1, rr_node_p.capacity, _pres_fac);
					//commit(*_current_regions, _g, parent_rr_node, 1, _instance);
					mods_add(parent_rr_node, 1, _rr_regions, mods);

					added_rr_nodes.push_back(parent_rr_node);
				}

				zlog_level(delta_log, ROUTER_V3, "\n");

				child_rr_node = parent_rr_node;
				edge = get_previous_edge(parent_rr_node, rt);
			}

			mods_commit(mods, _instance, _fpga_regions);
		}

		//template<typename Graph, typename Edge, typename EdgeWeightFunc, typename Callbacks>
		//friend void delta_stepping(const Graph &g, const vector<existing_source_t<Edge>> &sources, int sink, float delta, float *known_distance, float *distance, Edge *prev_edge, const EdgeWeightFunc &edge_weight, Callbacks &callbacks);

		template<typename Graph, typename Edge, typename EdgeWeightFunc, typename ExpandCheckFunc, typename Callbacks>
		friend void dijkstra(const Graph &g, const vector<existing_source_t<Edge>> &sources, int sink, float *known_distance, float *distance, Edge *prev_edge, const ExpandCheckFunc &expand_node, const EdgeWeightFunc &edge_weight, Callbacks &callbacks);

		//template<typename Edge, typename Callbacks>
		//friend void relax(Buckets &buckets, float delta, vector<bool> &in_bucket, const vector<bool> &vertex_deleted,
					//float *known_distance, float *distance, Edge *predecessor,
					//int v, float new_known_distance, float new_distance, const Edge &edge,
					//Callbacks &callbacks);

	public:
		SpeculativeDeterministicRouter(RRGraph &g, congestion_t **congestion, vector<region_t *> &fpga_regions, vector<region_t *> &rr_regions, exec_state_t *e_state, const float &pres_fac)
			: _g(g), _fpga_regions(fpga_regions), _rr_regions(rr_regions), _pres_fac(pres_fac), _modified_node_added(num_vertices(g), false) 
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
			
			_instance = _num_instances++;

			_congestion = congestion[_instance];
			_e_state = &e_state[_instance];

			_acquired_regions.resize(fpga_regions.size());

			reset_stats();
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

				std::fill(begin(_acquired_regions), end(_acquired_regions), false);

				sprintf_rr_node(sink->rr_node, buffer);
				//zlog_level(delta_log, ROUTER_V3, "Current sink: %s BB %d->%d, %d->%d\n", buffer, sink->current_bounding_box.xmin, sink->current_bounding_box.xmax, sink->current_bounding_box.ymin, sink->current_bounding_box.ymax);
				zlog_level(delta_log, ROUTER_V3, "Current sink: %s BB %d->%d,%d->%d\n",
						buffer, 
						sink->current_bounding_box.xmin, sink->current_bounding_box.xmax, sink->current_bounding_box.ymin, sink->current_bounding_box.ymax);

				vector<existing_source_t<RREdge>> sources;

				for (const auto &rt_node : route_tree_get_nodes(rt)) {
					const auto &rt_node_p = get_vertex_props(rt.graph, rt_node);
					RRNode node = rt_node_p.rr_node;
					const auto &node_p = get_vertex_props(_g, node);

					sprintf_rr_node(node, buffer);

					if (rt_node_p.reexpand) {
						if (!inside_bb(node_p, _current_sink->bounding_box)) {
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

std::atomic<int> SpeculativeDeterministicRouter::_num_instances{ 0 };

//struct vec_clock_t {
	//vector<int> counts;
//};
//

int get_next_net(det_mutex_t &lock, int tid, std::atomic<int> &net_index)
{
	int inet;
	det_mutex_lock(lock, tid);	
	inet = net_index++;
	det_mutex_unlock(lock, tid);
	return inet;
}

static void do_work_2(router_t<net_t *> *router, SpeculativeDeterministicRouter &net_router)
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

	det_mutex_t lock;

	det_mutex_init(lock, router->e_state.data(), router->num_threads);

	auto pure_route_time = timer::duration::zero();
	auto wait_time = timer::duration::zero();

	auto route_start = timer::now();

	int inet;

	while ((inet = get_next_net(lock, tid, router->sync.net_index)) < router->nets.size()) {
		*router->e_state[tid].inet = inet;

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
		region_get_mods(*r, tid, [&router, &tid] (const mod_t &m) -> void {
				const auto &props = get_vertex_props(router->g, m.node);
				update_first_order_congestion(router->state.congestion[tid], m.node, m.delta, props.capacity, router->pres_fac);
				});
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

static void do_work(router_t<net_t *> *router, SpeculativeDeterministicRouter &net_router)
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
		region_get_mods(*r, tid, [&router, &tid] (const mod_t &m) -> void {
				const auto &props = get_vertex_props(router->g, m.node);
				update_first_order_congestion(router->state.congestion[tid], m.node, m.delta, props.capacity, router->pres_fac);
				});
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

static void *worker_thread(void *args)
{
	router_t<net_t *> *router = static_cast<router_t<net_t *> *>(args);

#ifdef PTHREAD_BARRIER
	pthread_barrier_wait(&router->sync.barrier);
#else
	router->sync.barrier->wait();
#endif

	SpeculativeDeterministicRouter net_router(router->g, router->state.congestion, router->fpga_regions, router->rr_regions, router->e_state.data(), router->pres_fac);

	char buffer[256];
	sprintf(buffer, "%d", net_router.get_instance());
	zlog_put_mdc("tid", buffer);
	zlog_put_mdc("log_dir", dirname);

	while (!router->sync.stop_routing) {
		do_work_2(router, net_router);
	}

	return nullptr;
}

void write_routes(const char *fname, const route_tree_t *route_trees, int num_nets);
void write_net_dijkstra_stats(const char *fname, const vector<dijkstra_stats_t> &net_stats);
void write_congestion_state(const char *fname, const congestion_t *congestion, int num_rr_nodes);

static void init_fpga_regions(int split, exec_state_t *e_state, int num_threads, vector<region_t> &fpga_regions)
{
	extern int nx, ny;

	int xsplit = split;
	int ysplit = (float)ny/nx*xsplit; // take aspect ratio into account
	int width = ceil((float)(nx+2)/xsplit);
	int height = ceil((float)(ny+2)/ysplit);

	int gid = 0;

	for (int x = 0; x < nx+2; x += width) {
		int xmax = std::min(x+width-1, nx+2);

		for (int y = 0; y < ny+2; y += height) {
			int ymax = std::min(y+height-1, ny+2);

			fpga_regions.emplace_back();
			auto &r = fpga_regions.back();

			r.gid = gid;
			bg::assign_values(r.bb, x, y, xmax, ymax);
			det_mutex_init(r.mutex, e_state, num_threads);
			r.mods.resize(num_threads);
			r.num_committed_mods.resize(num_threads);
			for (auto &n : r.num_committed_mods) {
				n.resize(num_threads, 0);
			}
			++gid;

			zlog_level(delta_log, ROUTER_V3, "FPGA region %d, BB %d,%d %d,%d\n", r.gid,
					bg::get<bg::min_corner, 0>(r.bb), bg::get<bg::min_corner, 1>(r.bb),
					bg::get<bg::max_corner, 0>(r.bb), bg::get<bg::max_corner, 1>(r.bb));
		}
	}

	for (int i = 0; i < fpga_regions.size(); ++i) {
		for (int j = i+1; j < fpga_regions.size(); ++j) {
			assert(!bg::intersects(fpga_regions[i].bb, fpga_regions[j].bb));
		}
	}
}

static void init_rr_node_regions(const RRGraph &g, vector<region_t> &fpga_regions, vector<region_t *> &rr_regions)
{
	rr_regions.resize(num_vertices(g), nullptr);	

	for (int i = 0; i < rr_regions.size(); ++i) {
		const auto &props = get_vertex_props(g, i);

		int num_regions = 0;
		for (auto &fpga_r : fpga_regions) {
			if (inside_region(props, fpga_r.bb)) {
				assert(!rr_regions[i]);
				rr_regions[i] = &fpga_r;
				zlog_level(delta_log, ROUTER_V3, "RR %d is inside region %d\n", i, fpga_r.gid);
				++num_regions;
			}
		}

		assert(num_regions == 1);
	}
}

bool speculative_deterministic_route_hb_fine(t_router_opts *opts)
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

	g_router = &router;

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

	for (int i = 0; i < nets.size(); ++i) {
		router.nets.push_back(&nets[i]);
	}
	std::sort(begin(router.nets), end(router.nets), [] (const net_t *a, const net_t *b) -> bool {
			return a->sinks.size() > b->sinks.size();
			});

	router.e_state.resize(opts->num_threads);
	for (auto &e : router.e_state) {
		e.logical_clock = new std::atomic<unsigned long>(0);
		e.inet = new std::atomic<int>(0);
	}

	router.window_size = router.num_threads * 2;

	vector<region_t> fpga_regions;
	init_fpga_regions(4, router.e_state.data(), opts->num_threads, fpga_regions);
	for (auto &r : fpga_regions) {
		router.fpga_regions.push_back(&r);	
	}
	init_rr_node_regions(router.g, fpga_regions, router.rr_regions);

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

	router.state.congestion = new congestion_t *[opts->num_threads];
	for (int i = 0; i < opts->num_threads; ++i) {
		router.state.congestion[i] = new congestion_t[num_vertices(router.g)];
		for (int j = 0; j < num_vertices(router.g); ++j) {
			router.state.congestion[i][j].acc_cost = 1;
			router.state.congestion[i][j].pres_cost = 1;
			router.state.congestion[i][j].occ = 0;
		}
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
	router.time.resize(opts->num_threads);

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
	SpeculativeDeterministicRouter *net_router = new SpeculativeDeterministicRouter(router.g, router.state.congestion, router.fpga_regions, router.rr_regions, router.e_state.data(), router.pres_fac);

	sprintf(buffer, "%d", net_router->get_instance());
	zlog_put_mdc("tid", buffer);

	bool routed = false;
	bool has_resched = false;
	for (; router.iter < opts->max_router_iterations && !routed; ++router.iter) {
		sprintf(buffer, "%d", router.iter.load());
		zlog_put_mdc("iter", buffer);

		printf("Routing iteration: %d\n", router.iter.load());

		//assert(router.sync.num_incoming_edges.empty());

		//router.sync.handles.resize(router.topo_nets.size());

		//for (const auto &vnet : router.topo_nets) {
			//router.sync.handles[vnet->global_index] = router.sync.num_incoming_edges.emplace(make_pair(router.num_incoming_edges[vnet->global_index], vnet));
		//}
		
		//assert(router.sync.fast_global_pending_nets.size_approx() > 0);

		router.sync.net_index = 0;

		for (int i = 0; i < opts->num_threads; ++i) {
			router.time[i].route_time = timer::duration::zero();
			router.time[i].wait_time = timer::duration::zero();
			router.time[i].pure_route_time = timer::duration::zero();
			router.time[i].last_barrier_wait_time = timer::duration::zero();
			router.time[i].last_sync_time = timer::duration::zero();

			router.stats[i].num_heap_pops = 0;
			router.stats[i].num_heap_pushes = 0;
			router.stats[i].num_neighbor_visits = 0;
		}

		for (int i = 0; i < nets.size(); ++i) {
			router.net_route_time[i] = timer::duration::zero();

			router.net_stats[i].num_heap_pops = std::numeric_limits<int>::max();
			router.net_stats[i].num_heap_pushes = std::numeric_limits<int>::max();
			router.net_stats[i].num_neighbor_visits = std::numeric_limits<int>::max();

			router.net_router[i] = -1;
		}

		for (auto &net : nets) {
			update_sink_criticalities(net, router.state.net_timing[net.vpr_id], router.params);
		}

		for (int i = 0; i < router.fpga_regions.size(); ++i) {
			region_clear_mods(*router.fpga_regions[i]);
		}

		for (int i = 0; i < router.e_state.size(); ++i) {
			*router.e_state[i].logical_clock = 0;
			*router.e_state[i].inet = i;
		}

		//do_virtual_work(&router, *net_router);
		do_work_2(&router, *net_router);

		/* checking */
		for (int i = 0; i < router.fpga_regions.size(); ++i) {
			for (int k = 0; k < opts->num_threads; ++k) {
				for (int l = 0; l < opts->num_threads; ++l) {
					if (l == k) {
						continue;
					}

					if (router.fpga_regions[i]->num_committed_mods[k][l] != router.fpga_regions[i]->mods[l].size()) {
						printf("%d %d %d %d %d %d\n", i, k, l, router.fpga_regions[i]->gid, router.fpga_regions[i]->num_committed_mods[k][l],
								router.fpga_regions[i]->mods[l].size());
						assert(false);
					}
				}
			}
		}

		for (int i = 0; i < opts->num_threads; ++i) {
			for (int j = 0; j < num_vertices(router.g); ++j) {
				router.state.congestion[i][j].recalc_occ = 0; 
			}

			for (const auto &net : nets) {
				check_route_tree(router.state.route_trees[net.local_id], net, router.g);
				recalculate_occ(router.state.route_trees[net.local_id], router.g, router.state.congestion[i]);
			}

			bool valid = true;
			for (int j = 0; j < num_vertices(router.g); ++j) {
				if (router.state.congestion[i][j].recalc_occ != router.state.congestion[i][j].occ) {
					sprintf_rr_node_impl(j, buffer);
					printf("Node %s occ mismatch, recalc: %d original: %d\n", buffer, router.state.congestion[i][j].recalc_occ, router.state.congestion[i][j].occ);
					valid = false;
				}
			}
			assert(valid);
		}

		//printf("Resched time [has_resched %d]: %g\n", has_resched ? 1 : 0, duration_cast<nanoseconds>(resched_time).count()/1e9);

		printf("Route time: ");
		for (int i = 0; i < opts->num_threads; ++i) {
			printf("%g ", duration_cast<nanoseconds>(router.time[i].route_time).count()/1e9);
		}
		printf("\n");

		//printf("Route time debug: ");
		//for (auto iter = begin(router.route_time); iter != end(router.route_time); ++iter) {
			//printf("%g ", duration_cast<nanoseconds>(*iter).count()/1e9);
		//}
		//printf("\n");
		timer::duration max_route_time = router.time[0].route_time;
		for (int i = 1;  i < opts->num_threads; ++i) {
			max_route_time = std::max(router.time[i].route_time, max_route_time);
		}
		total_route_time += max_route_time;

		printf("Max route time: %lu %g\n", max_route_time.count(), max_route_time.count() / 1e9);

		printf("Wait time: ");
		for (int i = 0; i < opts->num_threads; ++i) {
			printf("%g (%g) ", duration_cast<nanoseconds>(router.time[i].wait_time).count()/1e9, 100.0*router.time[i].wait_time.count()/router.time[i].route_time.count());
		}
		printf("\n");

		printf("Pure route time: ");
		for (int i = 0; i < opts->num_threads; ++i) {
			printf("%g (%g) ", duration_cast<nanoseconds>(router.time[i].pure_route_time).count()/1e9, 100.0*router.time[i].pure_route_time.count()/router.time[i].route_time.count());
		}
		printf("\n");

		printf("Last barrier wait time: ");
		for (int i = 0; i < opts->num_threads; ++i) {
			printf("%g (%g) ", duration_cast<nanoseconds>(router.time[i].last_barrier_wait_time).count()/1e9, 100.0*router.time[i].last_barrier_wait_time.count()/router.time[i].route_time.count());
		}
		printf("\n");

		printf("Last sync time: ");
		for (int i = 0; i < opts->num_threads; ++i) {
			printf("%g (%g) ", duration_cast<nanoseconds>(router.time[i].last_sync_time).count()/1e9, 100.0*router.time[i].last_sync_time.count()/router.time[i].route_time.count());
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

		if (feasible_routing(router.g, router.state.congestion[0])) {
			//dump_route(*current_traces_ptr, "route.txt");
			routed = true;
		} else {
			for (int i = 0; i < num_vertices(router.g); ++i) {
				if (router.state.congestion[0][i].occ > get_vertex_props(router.g, i).capacity) {
					++num_overused_nodes;
					const auto &v_p = get_vertex_props(router.g, i);
					++overused_nodes_by_type[v_p.type];
				}
			}

			//auto update_cost_start = timer::now();

			if (router.iter == 0) {
				router.pres_fac = opts->initial_pres_fac;
				update_costs(router.g, router.state.congestion[0], router.pres_fac, 0);
			} else {
				router.pres_fac *= opts->pres_fac_mult;

				/* Avoid overflow for high iteration counts, even if acc_cost is big */
				router.pres_fac = std::min(router.pres_fac, static_cast<float>(HUGE_POSITIVE_FLOAT / 1e5));

				update_costs(router.g, router.state.congestion[0], router.pres_fac, opts->acc_fac);
			}

			//update_cost_time = timer::now()-update_cost_start;
		}

		if (routed) {
			sprintf(buffer, "%s/routes_final.txt", dirname);
		} else {
			sprintf(buffer, "%s/routes_iter_%d.txt", dirname, router.iter.load());
		}
		write_routes(buffer, router.state.route_trees, nets.size());

		if (routed) {
			sprintf(buffer, "%s/net_dijkstra_stats_final.txt", dirname);
		} else {
			sprintf(buffer, "%s/net_dijkstra_stats_%d.txt", dirname, router.iter.load());
		}
		write_net_dijkstra_stats(buffer, router.net_stats);

		if (routed) {
			sprintf(buffer, "%s/congestion_state_final.txt", dirname);
		} else {
			sprintf(buffer, "%s/congestion_state_%d.txt", dirname, router.iter.load());
		}
		write_congestion_state(buffer, router.state.congestion[0], num_vertices(router.g));

		//auto analyze_timing_start = timer::now();

		float crit_path_delay = analyze_timing(router.state.net_timing);

		//analyze_timing_time = timer::now()-analyze_timing_start;
		
		printf("Overused: %lu/%lu (%g) Crit path delay: %g\n", num_overused_nodes, num_vertices(router.g), 100.0*num_overused_nodes/num_vertices(router.g), crit_path_delay);
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
		printf("Routed in %d iterations\n", router.iter.load());
		printf("Total route time: %g\n", duration_cast<nanoseconds>(total_route_time).count() / 1e9);
	}

	return false;
}
