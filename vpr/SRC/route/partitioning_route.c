#include <assert.h>
#include <zlog.h>
#include <boost/numeric/interval.hpp>
#include <boost/timer/timer.hpp>
#include "tbb/task_scheduler_init.h"
#include "tbb/tbb.h"
#include "util.h"
#include "vpr_types.h"
/*#include "globals.h"*/
//#include "route_export.h"
#include "route_common_local.h"
#include "route_tree_timing_local.h"
/*#include "route_timing.h"*/
#include "heapsort.h"
#include "path_delay.h"
#include "net_delay.h"
#include "stats.h"
#include "ReadOptions.h"
#include "rr_graph_util.h"
#include "rr_graph2.h"
#include "rr_graph.h"
#include "parallel_route_timing.h"
#include "barrier.h"
/*#include "bounding_box.h"*/
#include "graph.h"
#include "route.h"
#include "scheduler.h"
/*#include <boost/geometry.hpp>*/
/*#include <boost/geometry/geometries/point.hpp>*/
/*#include <boost/geometry/geometries/box.hpp>*/
/*#include <boost/geometry/index/rtree.hpp>*/

/* TODO: check whether nets are global before routing */
using namespace boost::timer;

enum {
	ROUTER_V1 = ZLOG_LEVEL_DEBUG+3,
	ROUTER_V2 = ZLOG_LEVEL_DEBUG+2, 
	ROUTER_V3 = ZLOG_LEVEL_DEBUG+1
};

#define zlog_level(cat, level, ...) \
	zlog(cat, __FILE__, sizeof(__FILE__)-1, __func__, sizeof(__func__)-1, __LINE__, \
	level, __VA_ARGS__)

#define zlog_level(cat, level, ...)
#define zlog_debug(cat, level, ...)

void sprintf_rr_node(int inode, char *buffer);

bool operator<(const sink_t &a, const sink_t &b)
{
	return a.criticality_fac > b.criticality_fac;
}

void load_overlapping_nets_rtree(vector<net_t *> &nets);

/*typedef struct heap_item_t {*/
	/*int node;*/
	/*float cost;*/
	/*float known_cost;*/
	/*float upstream_R;*/
/*} heap_item_t;*/

zlog_category_t *delta_log;
zlog_category_t *rr_log;
zlog_category_t *net_log;
zlog_category_t *schedule_log;
zlog_category_t *independent_log;

static void init_logging()
{
	delta_log = zlog_get_category("delta");
	rr_log = zlog_get_category("rr");
	net_log = zlog_get_category("net");
	schedule_log = zlog_get_category("schedule");
	independent_log = zlog_get_category("independent");
}

/*int get_num_nets()*/
/*{*/
	/*return num_nets;*/
/*}*/

/*int get_net_rr_terminals(int inet, int term)*/
/*{*/
	/*assert(inet < num_nets && term <= clb_net[inet].num_sinks);*/
	/*return net_rr_terminals[inet][term];*/
/*}*/

/*struct s_net *get_net(int inet)*/
/*{*/
	/*return &clb_net[inet];*/
/*}*/

bool operator<(const route_state_t &a, const route_state_t &b)
{
	return a.cost > b.cost;
}

float get_delay(const RREdge &e, const RRNode &v, float unbuffered_upstream_R)
{
	float upstream_R = e.properties.R;
	if (!e.properties.buffered) {
		upstream_R += unbuffered_upstream_R;
	}

	float delay = e.properties.switch_delay;
	delay += v.properties.C * (upstream_R + 0.5 * v.properties.R);

	zlog_level(delta_log, ROUTER_V3, " [edge_delay: %g edge_R: %g node_R: %g node_C: %g] ", e.properties.switch_delay, e.properties.R, v.properties.R, v.properties.C);

	return delay;
}

float get_congestion_cost(const RRNode &v)
{
	extern t_rr_indexed_data *rr_indexed_data;
	return rr_indexed_data[v.properties.cost_index].base_cost * v.properties.acc_cost * v.properties.pres_cost;
}

float get_known_cost(const RRGraph &g, const RREdge &e, float criticality_fac, float unbuffered_upstream_R)
{
	const auto &target = get_target(g, e);

	float delay = get_delay(e, target, unbuffered_upstream_R);
	float congestion = get_congestion_cost(target);

	zlog_level(delta_log, ROUTER_V3, " [delay: %g congestion %g crit_fac: %g] ", delay, congestion, criticality_fac);

	return criticality_fac * delay + (1 - criticality_fac) * congestion;
}

/* Macro used below to ensure that fractions are rounded up, but floating   *
 * point values very close to an integer are rounded to that integer.       */

#define ROUND_UP(x) (ceil (x - 0.001))

static int get_expected_segs_to_target(const RRNode &current, const RRNode &target,
		int *num_segs_ortho_dir_ptr) {

	/* Returns the number of segments the same type as inode that will be needed *
	 * to reach target_node (not including inode) in each direction (the same    *
	 * direction (horizontal or vertical) as inode and the orthogonal direction).*/

	t_rr_type rr_type;
	int target_x, target_y, num_segs_same_dir, cost_index, ortho_cost_index;
	int no_need_to_pass_by_clb;
	float inv_length, ortho_inv_length, ylow, yhigh, xlow, xhigh;

	extern t_rr_indexed_data *rr_indexed_data;

	target_x = target.properties.xlow;
	target_y = target.properties.ylow;
	cost_index = current.properties.cost_index;
	inv_length = rr_indexed_data[cost_index].inv_length;
	ortho_cost_index = rr_indexed_data[cost_index].ortho_cost_index;
	ortho_inv_length = rr_indexed_data[ortho_cost_index].inv_length;
	rr_type = current.properties.type;

	if (rr_type == CHANX) {
		ylow = current.properties.ylow;
		xhigh = current.properties.xhigh;
		xlow = current.properties.xlow;

		/* Count vertical (orthogonal to inode) segs first. */

		if (ylow > target_y) { /* Coming from a row above target? */
			*num_segs_ortho_dir_ptr =
					(int)(ROUND_UP((ylow - target_y + 1.) * ortho_inv_length));
			no_need_to_pass_by_clb = 1;
		} else if (ylow < target_y - 1) { /* Below the CLB bottom? */
			*num_segs_ortho_dir_ptr = (int)(ROUND_UP((target_y - ylow) *
					ortho_inv_length));
			no_need_to_pass_by_clb = 1;
		} else { /* In a row that passes by target CLB */
			*num_segs_ortho_dir_ptr = 0;
			no_need_to_pass_by_clb = 0;
		}

		/* Now count horizontal (same dir. as inode) segs. */

		if (xlow > target_x + no_need_to_pass_by_clb) {
			num_segs_same_dir = (int)(ROUND_UP((xlow - no_need_to_pass_by_clb -
							target_x) * inv_length));
		} else if (xhigh < target_x - no_need_to_pass_by_clb) {
			num_segs_same_dir = (int)(ROUND_UP((target_x - no_need_to_pass_by_clb -
							xhigh) * inv_length));
		} else {
			num_segs_same_dir = 0;
		}
	}

	else { /* inode is a CHANY */
		ylow = current.properties.ylow;
		yhigh = current.properties.yhigh;
		xlow = current.properties.xlow;

		/* Count horizontal (orthogonal to inode) segs first. */

		if (xlow > target_x) { /* Coming from a column right of target? */
			*num_segs_ortho_dir_ptr = (int)(
					ROUND_UP((xlow - target_x + 1.) * ortho_inv_length));
			no_need_to_pass_by_clb = 1;
		} else if (xlow < target_x - 1) { /* Left of and not adjacent to the CLB? */
			*num_segs_ortho_dir_ptr = (int)(ROUND_UP((target_x - xlow) *
					ortho_inv_length));
			no_need_to_pass_by_clb = 1;
		} else { /* In a column that passes by target CLB */
			*num_segs_ortho_dir_ptr = 0;
			no_need_to_pass_by_clb = 0;
		}

		/* Now count vertical (same dir. as inode) segs. */

		if (ylow > target_y + no_need_to_pass_by_clb) {
			num_segs_same_dir = (int)(ROUND_UP((ylow - no_need_to_pass_by_clb -
							target_y) * inv_length));
		} else if (yhigh < target_y - no_need_to_pass_by_clb) {
			num_segs_same_dir = (int)(ROUND_UP((target_y - no_need_to_pass_by_clb -
							yhigh) * inv_length));
		} else {
			num_segs_same_dir = 0;
		}
	}

	return (num_segs_same_dir);
}

static float get_timing_driven_expected_cost(const RRNode &current, const RRNode &target,
		float criticality_fac, float R_upstream) {

	/* Determines the expected cost (due to both delay and resouce cost) to reach *
	 * the target node from inode.  It doesn't include the cost of inode --       *
	 * that's already in the "known" path_cost.                                   */

	t_rr_type rr_type;
	int cost_index, ortho_cost_index, num_segs_same_dir, num_segs_ortho_dir;
	float expected_cost, cong_cost, Tdel;

	extern t_rr_indexed_data *rr_indexed_data;

	rr_type = current.properties.type;

	if (rr_type == CHANX || rr_type == CHANY) {
		num_segs_same_dir = get_expected_segs_to_target(current, target,
				&num_segs_ortho_dir);
		cost_index = current.properties.cost_index;
		ortho_cost_index = rr_indexed_data[cost_index].ortho_cost_index;

		cong_cost = num_segs_same_dir * rr_indexed_data[cost_index].saved_base_cost
				+ num_segs_ortho_dir
						* rr_indexed_data[ortho_cost_index].saved_base_cost;
		cong_cost += rr_indexed_data[IPIN_COST_INDEX].base_cost
				+ rr_indexed_data[SINK_COST_INDEX].base_cost;

		Tdel =
				num_segs_same_dir * rr_indexed_data[cost_index].T_linear
						+ num_segs_ortho_dir
								* rr_indexed_data[ortho_cost_index].T_linear
						+ num_segs_same_dir * num_segs_same_dir
								* rr_indexed_data[cost_index].T_quadratic
						+ num_segs_ortho_dir * num_segs_ortho_dir
								* rr_indexed_data[ortho_cost_index].T_quadratic
						+ R_upstream
								* (num_segs_same_dir
										* rr_indexed_data[cost_index].C_load
										+ num_segs_ortho_dir
												* rr_indexed_data[ortho_cost_index].C_load);

		Tdel += rr_indexed_data[IPIN_COST_INDEX].T_linear;

		expected_cost = criticality_fac * Tdel
				+ (1. - criticality_fac) * cong_cost;
		return (expected_cost);
	}

	else if (rr_type == IPIN) { /* Change if you're allowing route-throughs */
		return (rr_indexed_data[SINK_COST_INDEX].base_cost);
	}

	else { /* Change this if you want to investigate route-throughs */
		return (0.);
	}
}

template<typename ShouldExpandFunc>
void expand_neighbors(const RRGraph &g, const route_state_t &current, const RRNode &target, float criticality_fac, float astar_fac, std::priority_queue<route_state_t> &heap, const ShouldExpandFunc &should_expand, perf_t *perf)
{
	for_all_out_edges(g, get_vertex(g, current.rr_node), [&heap, &g, &current, &target, &criticality_fac, &astar_fac, &should_expand, &perf] (const RREdge &e) -> void {
			auto &neighbor = get_target(g, e);

			char buffer[256];
			sprintf_rr_node(id(neighbor), buffer);
			zlog_level(delta_log, ROUTER_V3, "\tNeighbor: %s ", buffer);
			
			if (!should_expand(neighbor)) {
				return;
			}

			route_state_t item;

			item.rr_node = id(neighbor);
			item.prev_edge = &e;

			float unbuffered_upstream_R = current.upstream_R;
			float upstream_R = e.properties.R + neighbor.properties.R;
			if (!e.properties.buffered) {
				upstream_R += unbuffered_upstream_R;
			}
			item.upstream_R = upstream_R;
			
			item.delay = current.delay + get_delay(e, neighbor, unbuffered_upstream_R);

			float known_cost = current.known_cost + get_known_cost(g, e, criticality_fac, unbuffered_upstream_R);
			item.known_cost = known_cost;

			float expected_cost = get_timing_driven_expected_cost(get_vertex(g, item.rr_node), target, criticality_fac, upstream_R);
			item.cost = known_cost + astar_fac * expected_cost;

			heap.push(item);

			if (perf) {
				++perf->num_heap_pushes;
			}
			
			zlog_level(delta_log, ROUTER_V3, "added with prev_edge: %X upstream_R: %g delay: %g known_cost: %g expected_cost: %g cost: %g\n", item.prev_edge, item.upstream_R, item.delay, item.known_cost, expected_cost, item.cost);
	});
}

void route_tree_init(route_tree_t &rt)
{
	rt.root = -1;
}

/*RouteTree &route_tree_get_sink_root(route_tree_t &rt, int sink_rr_node)*/
/*{*/
	/*assert(rt.sink_rr_node_to_start_rt_node.find(sink_rr_node) != rt.sink_rr_node_to_start_rt_node.end());*/
	/*return get_vertex(rt.graph, rt.sink_rr_node_to_start_rt_node[sink_rr_node]);*/
/*}*/

/*const vector<int> *route_tree_get_path_to_sink(const route_tree_t &rt, int sink_rr_node)*/
/*{*/
	/*const auto &iter = rt.sink_rr_node_to_path.find(sink_rr_node);*/
	/*const vector<int> *res = nullptr;*/
	/*if (iter != rt.sink_rr_node_to_path.end()) {*/
		/*res = &iter->second;*/
	/*}*/
	/*return res;*/
/*}*/

/*void route_tree_remove_path(route_tree_t &rt, int sink_rr_node)*/
/*{*/
	/*auto iter = rt.sink_rr_node_to_path.find(sink_rr_node);*/
	/*assert(iter != rt.sink_rr_node_to_path.end());*/
	/*vector<int> &path = iter->second;*/
	 /*[>important to sort because <]*/
	/*[>sort(path.begin(), path.end(), std::greater<int>());<]*/
	/*for (const auto &edge : path) {*/
		/*remove_edge(rt.graph, get_edge(rt.graph, edge));*/
		/*[> TODO: need to invalidate "path" after removing each edge <]*/
	/*}*/

	/*rt.sink_rr_node_to_path.erase(iter);*/
/*}*/

pair<RouteTreeNode *, bool> route_tree_add_rr_node(route_tree_t &rt, int rr_node)
{
	auto iter = rt.rr_node_to_rt_node.find(rr_node);

	RouteTreeNode *v = nullptr;
	bool added;
	if (iter == rt.rr_node_to_rt_node.end()) {
		add_vertex(rt.graph);
		int rt_node = num_vertices(rt.graph)-1;
		rt.rr_node_to_rt_node[rr_node] = rt_node;
		v = &get_vertex(rt.graph, rt_node);
		v->properties.valid = true;
		added = true;
	} else {
		v = &get_vertex(rt.graph, iter->second);
		added = !v->properties.valid;
		if (!v->properties.valid) {
			v->properties.valid = true;
		}
	}
	return make_pair(v, added);
}

/*RouteTreeEdge &route_tree_add_edge(route_tree_t &rt, int rt_node_a, int rt_node_b)*/
/*{*/
	/*[>int rt_node_a = route_tree_add_rr_node(rt, rr_node_a);<]*/
	/*[>int rt_node_b = route_tree_add_rr_node(rt, rr_node_b);<]*/
	/*[>assert(rt_node_a < num_vertices(rt.graph));<]*/
	/*[>assert(rt_node_b < num_vertices(rt.graph));<]*/
	/*return add_edge(rt.graph, rt_node_a, rt_node_b);*/
/*}*/

RouteTreeNode *route_tree_get_rt_node(route_tree_t &rt, int rr_node)
{
	auto iter = rt.rr_node_to_rt_node.find(rr_node);
	assert(iter != rt.rr_node_to_rt_node.end());
	RouteTreeNode *res;
	RouteTreeNode &v = get_vertex(rt.graph, iter->second);
	if (v.properties.valid) {
		res = &v;
	} else {
		res = nullptr;
	}
	return res;
}

const RouteTreeNode *route_tree_get_rt_node(const route_tree_t &rt, int rr_node)
{
	auto iter = rt.rr_node_to_rt_node.find(rr_node);
	assert(iter != rt.rr_node_to_rt_node.end());
	const RouteTreeNode *res;
	const RouteTreeNode &v = get_vertex(rt.graph, iter->second);
	if (v.properties.valid) {
		res = &v;
	} else {
		res = nullptr;
	}
	return res;
}

RouteTreeEdge &route_tree_add_edge_between_rr_node(route_tree_t &rt, int rr_node_a, int rr_node_b)
{
	/*int rt_node_a = route_tree_add_rr_node(rt, rr_node_a);*/
	/*int rt_node_b = route_tree_add_rr_node(rt, rr_node_b);*/
	/*assert(rt_node_a < num_vertices(rt.graph));*/
	/*assert(rt_node_b < num_vertices(rt.graph));*/
	RouteTreeNode *rt_node_a = route_tree_get_rt_node(rt, rr_node_a);
	RouteTreeNode *rt_node_b = route_tree_get_rt_node(rt, rr_node_b);

	if (get_edge(rt.graph, *rt_node_a, *rt_node_b)) {
		char buffer_a[256];
		char buffer_b[256];
		sprintf_rr_node(rt_node_a->properties.rr_node, buffer_a);
		sprintf_rr_node(rt_node_b->properties.rr_node, buffer_b);
		zlog_error(delta_log, "Existing edge between %s and %s\n", buffer_a, buffer_b);
		assert(false);
	}
	return add_edge(rt.graph, *rt_node_a, *rt_node_b);
}

void route_tree_merge(route_tree_t &merged_rt, const route_tree_t &other_rt)
{
	/*for_all_edges(other_rt.graph, [&merged_rt, &other_rt] (const RouteTreeEdge &e) -> void {*/
			/*const RouteTreeNode &other_source_rt_node = get_source(other_rt.graph, e);*/
			/*const RouteTreeNode &other_target_rt_node = get_target(other_rt.graph, e);*/
			/*int source_rr_node = other_source_rt_node.properties.rr_node;*/
			/*int target_rr_node = other_target_rt_node.properties.rr_node;*/

			/*RouteTreeNode &merged_source_rt_node = route_tree_add_rr_node(merged_rt, source_rr_node);*/
			/*merged_source_rt_node.properties = other_source_rt_node.properties;*/
			/*RouteTreeNode &merged_target_rt_node = route_tree_add_rr_node(merged_rt, target_rr_node);*/
			/*merged_target_rt_node.properties = other_target_rt_node.properties;*/

			/*route_tree_add_edge_between_rr_node(merged_rt, source_rr_node, target_rr_node);*/
			/*}*/
			/*);*/
}

void test_heap()
{
}

void route_tree_add_to_heap_internal(const route_tree_t &rt, const RouteTreeNode *rt_node, const RRGraph &g, const RRNode &target, float criticality_fac, float astar_fac, const bounding_box_t &current_bounding_box, std::priority_queue<route_state_t> &heap, perf_t *perf)
{
	assert(rt_node->properties.valid);

	const auto &prop = get_vertex(g, rt_node->properties.rr_node).properties;

	char buffer[256];
	sprintf_rr_node(rt_node->properties.rr_node, buffer);

	if (prop.xhigh < current_bounding_box.xmin
			|| prop.xlow > current_bounding_box.xmax
			|| prop.yhigh < current_bounding_box.ymin
			|| prop.ylow > current_bounding_box.ymax) {
		zlog_level(delta_log, ROUTER_V2, "Existing route tree node %s outside of current bounding box\n", buffer);
	} else if (rt_node->properties.reexpand) {
		route_state_t item;

		item.rr_node = rt_node->properties.rr_node;
		item.known_cost = criticality_fac * rt_node->properties.delay;
		float expected_cost = get_timing_driven_expected_cost(get_vertex(g, rt_node->properties.rr_node), target, criticality_fac, rt_node->properties.upstream_R);
		item.cost = item.known_cost + astar_fac * expected_cost;
		item.prev_edge = nullptr;
		item.upstream_R = rt_node->properties.upstream_R;
		item.delay = rt_node->properties.delay;

		zlog_level(delta_log, ROUTER_V2, "Adding route tree node %s delay: %g crit_fac: %g known_cost: %g expected_cost: %g astar_fac: %g cost: %g upstream_R: %g to heap\n", buffer, item.delay, criticality_fac, item.known_cost, expected_cost, astar_fac, item.cost, item.upstream_R);

		if (perf) {
			++perf->num_heap_pushes;
		}
		heap.push(item);
	}

	for_all_out_edges(rt.graph, *rt_node, [&rt, &heap, &criticality_fac, &astar_fac, &g, &target, &current_bounding_box, &perf] (const RouteTreeEdge &e) -> void {
		const auto &neighbor = get_target(rt.graph, e);

		route_tree_add_to_heap_internal(rt, &neighbor, g, target, criticality_fac, astar_fac, current_bounding_box, heap, perf);
	}
	);
}

void route_tree_add_to_heap(const route_tree_t &rt, const RRGraph &g, const RRNode &target, float criticality_fac, float astar_fac, const bounding_box_t &current_bounding_box, std::priority_queue<route_state_t> &heap, perf_t *perf)
{
	route_tree_add_to_heap_internal(rt, &get_vertex(rt.graph, rt.root), g, target, criticality_fac, astar_fac, current_bounding_box, heap, perf);
}

void trace_init(trace_t &trace)
{
	/*trace.num_sources = 0;*/
	/*trace.first_sink_rr_node = -1;*/
}

bool trace_empty(const trace_t &trace)
{
	return trace.existing_nodes.empty() && trace.paths_starting_with_source.empty() && trace.segments.empty();
}

const Segment &trace_add_path(trace_t &trace, const RRGraph &g, const route_state_t *state, int sink_rr_node, int vpr_net_id)
{
	assert(trace.segments.find(sink_rr_node) == trace.segments.end());

	auto &new_segment = trace.segments.insert(make_pair(sink_rr_node, Segment())).first->second;

	char buffer[256];
	int prev_rr_node = -1;
	int current_rr_node = sink_rr_node;
	while (state[current_rr_node].prev_edge) {
		if (trace.existing_nodes.find(current_rr_node) != trace.existing_nodes.end()) {
			sprintf_rr_node(current_rr_node, buffer);
			char prev_s[256];
			sprintf_rr_node(prev_rr_node, prev_s);
			zlog_warn(delta_log, "Warning: Trying to add existing node %s to trace. Stopping traceback\n", buffer);
			/*assert(false);*/
			/*current_rr_node = prev_rr_node;*/
			break;
		}

		new_segment.push_back(current_rr_node);
		trace.existing_nodes.insert(current_rr_node);

		sprintf_rr_node(current_rr_node, buffer);
		zlog_level(delta_log, ROUTER_V2, "Net %d trace: %s\n", vpr_net_id, buffer);

		int parent_rr_node = id(get_source(g, *state[current_rr_node].prev_edge));
		prev_rr_node = current_rr_node;
		current_rr_node = parent_rr_node;
	}

	/* we're checking for various conditions: */
	const RRNode &first_node = get_vertex(g, current_rr_node);
	sprintf_rr_node(current_rr_node, buffer);
	if (first_node.properties.type == SOURCE) {
		if (trace.paths_starting_with_source.size() == 0) {
			/* this is a first segment */
			/*assert(trace.segments.size() == 1);*/ /* this is not true after the first iteration because we're ripping up segment by segment */
			assert(trace.existing_nodes.find(current_rr_node) == trace.existing_nodes.end());

			zlog_level(delta_log, ROUTER_V2, "Net %d first segment trace: %s\n", vpr_net_id, buffer);

			/* shouldn't reach this case actually because paths_starting_with_source is empty */
			/*if (trace.existing_nodes.find(current_rr_node) != trace.existing_nodes.end()) {*/
				/*zlog_error(delta_log, "Error: Trying to add existing node %s prev_edge %d to trace\n", buffer, state[current_rr_node].prev_edge ? 1 : 0);*/
				/*assert(false);*/
			/*}*/

			new_segment.push_back(current_rr_node);
			trace.existing_nodes.insert(current_rr_node);
			trace.paths_starting_with_source.push_back(&new_segment);
		} else { /* we have existing paths that start with SOURCE */
			if (trace.existing_nodes.find(current_rr_node) != trace.existing_nodes.end()) {
				/* this is an existing SOURCE
				 * we need to check whether it has capacity of > 1
				 * >1 : no error 
				 * <=1 : error */
				assert(any_of(trace.paths_starting_with_source.begin(), trace.paths_starting_with_source.end(),
							[&current_rr_node] (const Segment *path) -> bool {
							return path->back() == current_rr_node;
							}));

				if (first_node.properties.capacity > 1) {
					zlog_level(delta_log, ROUTER_V2, "Reconnecting to existing %s\n", buffer);
					new_segment.push_back(current_rr_node);
					trace.paths_starting_with_source.push_back(&new_segment);
				} else {
					zlog_error(delta_log, "Error: Trying to reconnect to %s with only capacity of %d\n", buffer, first_node.properties.capacity);
					assert(false);
				}
			} else {
				/* this is a new SOURCE 
				 * error because we now have mutliple sources for the net */
				zlog_error(delta_log, "Error: Trying to start the trace with another source %s when existing sources ", buffer);
				for (const auto &p : trace.paths_starting_with_source) {
					char old_source[256];
					sprintf_rr_node(p->back(), old_source);
					zlog_error(delta_log, "%s ", buffer);
				}
				zlog_error(delta_log, "exist\n");
				assert(false);
			}
		}
	} else {
		assert(first_node.properties.type == CHANX || first_node.properties.type == CHANY || first_node.properties.type == OPIN);

		if (trace.segments.size() == 1) {
			zlog_error(delta_log, "Error: First path started with %s instead of a SOURCE\n", buffer);
			assert(false);
		}
	}

	/*if (get_vertex(g, current_rr_node).properties.type == SOURCE) {*/
		/*if (trace.num_sources != 0) {*/
			/*zlog_error(delta_log, "Error: Trying to add another path that starts from SOURCE [num_sources: %d]\n", trace.num_sources);*/
			/*assert(false);*/
		/*}*/
		/*new_segment.push_back(current_rr_node);*/
		/*++trace.num_sources;*/
	/*}*/

	return new_segment;
}

void adjust_bb_factor(trace_t &trace, net_t &net, const RRGraph &g, int threshold)
{
	assert(net.sinks.size() == trace.segments.size());

	for (const auto &item : trace.segments) {
		const Segment &segment = item.second;
		int sink_rr_node = item.first;

		int num_overused_nodes = 0;
		for (const auto &rr_node : segment) {
			auto &v = get_vertex(g, rr_node);
			if (v.properties.occ > v.properties.capacity) {
				++num_overused_nodes;
			}
		}

		int isink = -1;
		int num_matches = 0;
		for (int i = 0; i < net.sinks.size(); ++i) {
			if (net.sinks[i].rr_node == sink_rr_node) {
				isink = i;
				++num_matches;
			}
		}
		assert(num_matches == 1 && isink != -1);
		auto &sink = net.sinks[isink];

		char buffer[256];
		sprintf_rr_node(sink_rr_node, buffer);

		if (num_overused_nodes > 0) {
			const bounding_box_t &bb = sink.current_bounding_box;
			zlog_level(delta_log, ROUTER_V2, "Net %d %s (distance rank %d/%d BB %d-%d %d-%d) segment has %d/%d (%g) overused nodes\n", net.vpr_id, buffer, sink.distance_to_source_rank, net.sinks.size(), bb.xmin, bb.xmax, bb.ymin, bb.ymax, num_overused_nodes, segment.size(), (float)num_overused_nodes/segment.size()*100);

			++sink.congested_iterations;
			if (sink.congested_iterations >= threshold) {
				sink.congested_iterations = 0;
				zlog_level(delta_log, ROUTER_V2, "Increasing net %d %s bb_factor by 1\n", net.vpr_id, buffer);
				++sink.bb_factor;
			}
		} else {
			/*if (sink.congested_iterations > 0) {*/
				/*zlog_level(delta_log, ROUTER_V2, "Decreasing net %d %s bb_factor by 1\n", net.vpr_id, buffer);*/
				/*--sink.congested_iterations;*/
			/*}*/
		}

		assert(sink.congested_iterations >= 0);
	}
}

/*const vector<int> &route_tree_add_path(route_tree_t &rt, const RRGraph &g, const route_state_t *state, pair<vector<int>::const_iterator, vector<int>::const_iterator> &path)*/
/*{*/
	/*int new_branch = -1;*/

	/*auto iter = path.first*/

	/*for (; iter != path.second; ++iter) {*/
		/*int current_rr_node = *iter;*/
		/*int prev_rr_node = *(iter + 1);*/

		/*auto &prev_rt_node = route_tree_add_rr_node(rt, prev_rr_node);*/
		/*prev_rt_node.properties.reexpand = get_vertex(g, prev_rr_node).properties.type != IPIN;*/
		/*prev_rt_node.properties.rr_node = prev_rr_node;*/
		/*prev_rt_node.properties.prev_edge = nullptr;*/
		/*prev_rt_node.properties.upstream_R = state[prev_rr_node].upstream_R;*/
		/*prev_rt_node.properties.delay = state[prev_rr_node].delay;*/

		/*auto &current_rt_node = route_tree_add_rr_node(rt, current_rr_node);*/
		/*current_rt_node.properties.reexpand = get_vertex(g, current_rr_node).properties.type != IPIN;*/
		/*current_rt_node.properties.rr_node = current_rr_node;*/
		/*current_rt_node.properties.prev_edge = state[current_rr_node].prev_edge;*/
		/*current_rt_node.properties.upstream_R = state[current_rr_node].upstream_R;*/
		/*current_rt_node.properties.delay = state[current_rr_node].delay;*/

		/*char p[256];*/
		/*char c[256];*/
		/*sprintf_rr_node(prev_rr_node, p);*/
		/*sprintf_rr_node(current_rr_node, c);*/
		/*zlog_level(delta_log, ROUTER_V2, "Route tree add path: %s -> %s\n", p, c);*/

		/*auto &e = route_tree_add_edge_between_rr_node(rt, prev_rr_node, current_rr_node);*/
		/*e.properties.rr_edge = state[current_rr_node].prev_edge;*/

		/*path.push_back(current_rr_node);*/

		/*new_branch = current_rr_node;*/
	/*}*/

	/*[>assert(rt.sink_rr_node_to_start_rt_node.find(sink_rr_node) != rt.sink_rr_node_to_start_rt_node.end());<]*/
	/*RouteTreeNode &new_branch_rt_node = route_tree_get_rt_node(rt, new_branch);*/
	/*[>rt.sink_rr_node_to_start_rt_node[sink_rr_node] = id(new_branch_rt_node);<]*/

	/*assert(rt.sink_rr_node_to_path.find(sink_rr_node) == rt.sink_rr_node_to_path.end());*/
	/*const auto &iter = rt.sink_rr_node_to_path.insert(make_pair(sink_rr_node, path)).first;*/

	/*return iter->second;*/
/*}*/

void route_tree_add_path(route_tree_t &rt, const RRGraph &g, const route_state_t *state, const vector<int> &rr_nodes, int vpr_net_id)
{
	char p[256];
	char c[256];

	for (int i = 0; i < rr_nodes.size()-1; ++i) {
		int current_rr_node = rr_nodes[i];
		RouteTreeNode *current_rt_node;
		bool added;

		std::tie(current_rt_node, added) = route_tree_add_rr_node(rt, current_rr_node);
		if (!added) {
			if (current_rt_node->properties.upstream_R != state[current_rr_node].upstream_R ||
					current_rt_node->properties.delay != state[current_rr_node].delay) {
				zlog_error(delta_log, "Error: Adding a new RT node with different properties [old_upstream_R: %g new_upstream_R: %g] [old_delay: %g new_delay: %g]\n", current_rt_node->properties.upstream_R, state[current_rr_node].upstream_R, current_rt_node->properties.delay, state[current_rr_node].delay); 
				assert(false);
			}
		}
		current_rt_node->properties.reexpand = get_vertex(g, current_rr_node).properties.type != IPIN;
		current_rt_node->properties.rr_node = current_rr_node;
		current_rt_node->properties.prev_edge = state[current_rr_node].prev_edge;
		current_rt_node->properties.upstream_R = state[current_rr_node].upstream_R;
		current_rt_node->properties.delay = state[current_rr_node].delay;

		const RRNode &parent = get_source(g, *state[current_rr_node].prev_edge);
		int parent_rr_node = id(parent);
		assert(parent_rr_node == rr_nodes[i+1]);

		RouteTreeNode *parent_rt_node;
		std::tie(parent_rt_node, added) = route_tree_add_rr_node(rt, parent_rr_node);
		if (!added) {
			if (parent_rt_node->properties.upstream_R != state[parent_rr_node].upstream_R ||
					parent_rt_node->properties.delay != state[parent_rr_node].delay) {
				zlog_error(delta_log, "Error: Adding a new RT node with different properties [old_upstream_R: %g new_upstream_R: %g] [old_delay: %g new_delay: %g]\n", parent_rt_node->properties.upstream_R, state[parent_rr_node].upstream_R, parent_rt_node->properties.delay, state[parent_rr_node].delay); 
				assert(false);
			}
		}

		auto &e = route_tree_add_edge_between_rr_node(rt, parent_rr_node, current_rr_node);
		e.properties.rr_edge = state[current_rr_node].prev_edge;

		sprintf_rr_node(parent_rr_node, p);
		sprintf_rr_node(current_rr_node, c);
		zlog_level(delta_log, ROUTER_V2, "Net %d route tree edge: %s -> %s\n", vpr_net_id, p, c);

		parent_rt_node->properties.reexpand = parent.properties.type != IPIN;
		parent_rt_node->properties.rr_node = parent_rr_node;
		if (i+1 == rr_nodes.size()-1) {
			if (parent.properties.type == SOURCE) {
				parent_rt_node->properties.prev_edge = nullptr;
			} else {
				const RREdge *edge_to_connection = state[parent_rr_node].prev_edge;
				int connection_rr_node = id(get_source(g, *edge_to_connection));
				if (state[connection_rr_node].prev_edge != nullptr) {
					sprintf_rr_node(connection_rr_node, c);
					sprintf_rr_node(id(get_source(g, *state[connection_rr_node].prev_edge)), p);
					zlog_warn(delta_log, "Warning: Connection %s is not from existing route tree because we found a shorter path to it via %s\n", c, p);
				}

				auto &e = route_tree_add_edge_between_rr_node(rt, connection_rr_node, parent_rr_node);
				e.properties.rr_edge = edge_to_connection;

				parent_rt_node->properties.prev_edge = edge_to_connection;

				sprintf_rr_node(connection_rr_node, p);
				sprintf_rr_node(parent_rr_node, c);

				zlog_level(delta_log, ROUTER_V2, "Net %d connection route tree edge: %s -> %s\n", vpr_net_id, p, c);
			}
		} else {
			/* will be updated next iteration */
			parent_rt_node->properties.prev_edge = nullptr;
		}
		parent_rt_node->properties.upstream_R = state[parent_rr_node].upstream_R;
		parent_rt_node->properties.delay = state[parent_rr_node].delay;
	}
}

void route_tree_add_path(route_tree_t &rt, const RRGraph &g, const route_state_t *state, int sink_rr_node, int vpr_net_id)
{
	int current_rr_node = sink_rr_node;
	int new_branch = -1;

	/*vector<int> path;*/

	while (state[current_rr_node].prev_edge) {
		int prev_rr_node = id(get_source(g, *state[current_rr_node].prev_edge));

		char p[256];
		char c[256];
		sprintf_rr_node(prev_rr_node, p);
		sprintf_rr_node(current_rr_node, c);
		zlog_level(delta_log, ROUTER_V2, "Net %d route tree edge: %s -> %s\n", vpr_net_id, p, c);

		RouteTreeNode *prev_rt_node;
	   	bool added;
		std::tie(prev_rt_node, added) = route_tree_add_rr_node(rt, prev_rr_node);
		if (!added) {
			if (prev_rt_node->properties.upstream_R != state[prev_rr_node].upstream_R ||
					prev_rt_node->properties.delay != state[prev_rr_node].delay) {
				zlog_error(delta_log, "Error: Adding a new RT node with different properties [old_upstream_R: %g new_upstream_R: %g] [old_delay: %g new_delay: %g]\n", prev_rt_node->properties.upstream_R, state[prev_rr_node].upstream_R, prev_rt_node->properties.delay, state[prev_rr_node].delay); 
			}
		}
		prev_rt_node->properties.reexpand = get_vertex(g, prev_rr_node).properties.type != IPIN;
		prev_rt_node->properties.rr_node = prev_rr_node;
		prev_rt_node->properties.prev_edge = nullptr;
		prev_rt_node->properties.upstream_R = state[prev_rr_node].upstream_R;
		prev_rt_node->properties.delay = state[prev_rr_node].delay;

		RouteTreeNode *current_rt_node;
		std::tie(current_rt_node, added) = route_tree_add_rr_node(rt, current_rr_node);
		if (!added) {
			if (current_rt_node->properties.upstream_R != state[current_rr_node].upstream_R ||
					current_rt_node->properties.delay != state[current_rr_node].delay) {
				zlog_error(delta_log, "Error: Adding a new RT node with different properties [old_upstream_R: %g new_upstream_R: %g] [old_delay: %g new_delay: %g]\n", current_rt_node->properties.upstream_R, state[current_rr_node].upstream_R, current_rt_node->properties.delay, state[current_rr_node].delay); 
			}
		}
		current_rt_node->properties.reexpand = get_vertex(g, current_rr_node).properties.type != IPIN;
		current_rt_node->properties.rr_node = current_rr_node;
		current_rt_node->properties.prev_edge = state[current_rr_node].prev_edge;
		current_rt_node->properties.upstream_R = state[current_rr_node].upstream_R;
		current_rt_node->properties.delay = state[current_rr_node].delay;

		auto &e = route_tree_add_edge_between_rr_node(rt, prev_rr_node, current_rr_node);
		e.properties.rr_edge = state[current_rr_node].prev_edge;

		/*path.push_back(current_rr_node);*/

		new_branch = current_rr_node;
		current_rr_node = prev_rr_node;
	}

	/*assert(rt.sink_rr_node_to_start_rt_node.find(sink_rr_node) != rt.sink_rr_node_to_start_rt_node.end());*/
	/*RouteTreeNode &new_branch_rt_node = route_tree_get_rt_node(rt, new_branch);*/
	/*rt.sink_rr_node_to_start_rt_node[sink_rr_node] = id(new_branch_rt_node);*/

	/*assert(rt.sink_rr_node_to_path.find(sink_rr_node) == rt.sink_rr_node_to_path.end());*/
	/*const auto &iter = rt.sink_rr_node_to_path.insert(make_pair(sink_rr_node, path)).first;*/

	/*return iter->second;*/
}

void route_tree_clear(route_tree_t &rt)
{
	rt.root = -1;
	clear_edges(rt.graph);
	for (const auto &item : rt.rr_node_to_rt_node) {
		RouteTreeNode &rt_node = get_vertex(rt.graph, item.second);
		rt_node.properties.valid = false;
	}
}

void route_tree_set_source(route_tree_t &rt, const RRNode &source)
{
	const auto &rt_node = route_tree_add_rr_node(rt, id(source)); 

	rt_node.first->properties.upstream_R = source.properties.R;
	rt_node.first->properties.delay = 0.5 * source.properties.R * source.properties.C;
	rt_node.first->properties.rr_node = id(source);
	rt_node.first->properties.reexpand = true;
	rt_node.first->properties.prev_edge = nullptr;

	rt.root = id(*rt_node.first);
}

void update_costs(RRGraph &g, float pres_fac, float acc_fac)
{
	for_all_vertices(g, [&pres_fac, &acc_fac] (RRNode &v) -> void {
			int occ = v.properties.occ;
			int capacity = v.properties.capacity;
			if (occ > capacity) {
			v.properties.acc_cost += (occ - capacity) * acc_fac;
			v.properties.pres_cost = 1.
			+ (occ + 1 - capacity) * pres_fac;
			}

			/* If occ == capacity, we don't need to increase acc_cost, but a change    *
			 * in pres_fac could have made it necessary to recompute the cost anyway.  */

			else if (occ == capacity) {
			v.properties.pres_cost = 1. + pres_fac;
			}
			});
}

void update_one_cost_internal(RRNode &rr_node, int delta, float pres_fac)
{
	rr_node.properties.occ += delta;

	assert(rr_node.properties.occ >= 0);

	if (rr_node.properties.occ < rr_node.properties.capacity) {
		rr_node.properties.pres_cost = 1;
	} else {
		rr_node.properties.pres_cost = 1 + (rr_node.properties.occ + 1 - rr_node.properties.capacity) * pres_fac;
	}
		
	char buffer[256];
	sprintf_rr_node(id(rr_node), buffer);
	zlog_level(delta_log, ROUTER_V2, "Update cost of %s delta: %d pres_fac: %g\n", buffer, delta, pres_fac);
}

void update_one_cost(RRGraph &g, const vector<int>::const_iterator &rr_nodes_begin, const vector<int>::const_iterator &rr_nodes_end, int delta, float pres_fac)
{
	/*const RouteTreeNode *last = nullptr;*/
	for (auto iter = rr_nodes_begin; iter != rr_nodes_end; ++iter) {
		/* we don't update get_source because that's the link to existing route tree and the cost is handled by update_one_cost_internal */
		/*const RouteTreeNode &rt_node = get_target(rt.graph, get_edge(rt.graph, rt_edge_id));*/
		/*last = &get_target(rt.graph, get_edge(rt.graph, rt_edge_id));*/
		/*RRNode &rr_node = get_vertex(g, rt_node.properties.rr_node);*/
		RRNode &rr_node = get_vertex(g, *iter);
		update_one_cost_internal(rr_node, delta, pres_fac);
	}
	/*RRNode &rr_node = get_vertex(g, last->properties.rr_node);*/
	/*update_one_cost_internal(rr_node, delta, pres_fac);*/
}

void update_one_cost(RRGraph &g, route_tree_t &rt, const RouteTreeNode &node, int delta, float pres_fac)
{
	assert(node.properties.valid);

	RRNode &rr_node = get_vertex(g, node.properties.rr_node);

	update_one_cost_internal(rr_node, delta, pres_fac);

	for_all_out_edges(rt.graph, node, [&g, &rt, &delta, &pres_fac] (const RouteTreeEdge &e) -> void {
			update_one_cost(g, rt, get_target(rt.graph, e), delta, pres_fac);
			});
}

void update_R()
{
}

/*void route_tree_update_timing(route_tree_t &rt, const RouteTreeNode &node)*/
/*{*/
	/*for_all_out_edges(rt.graph, node, [] (const RouteTreeEdge &e) -> void {*/
			/*e.properties.rr_edge;*/
	/*});*/
/*}*/

void update_sink_criticalities(net_t &net, const t_net_timing &net_timing, const route_parameters_t &params)
{
	for (int ipin = 1; ipin <= net.sinks.size(); ipin++) { 
		float pin_criticality;
		/*if (!net_timing) {*/
			/* Use criticality of 1. This makes all nets critical.  Note: There is a big difference between setting pin criticality to 0*/
			/*compared to 1.  If pin criticality is set to 0, then the current path delay is completely ignored during routing.  By setting*/
			/*pin criticality to 1, the current path delay to the pin will always be considered and optimized for */
			/*pin_criticality = 1.0;*/
		/*} else { */
#ifdef PATH_COUNTING
			/* Pin criticality is based on a weighted sum of timing and path criticalities. */	
			pin_criticality =		 ROUTE_PATH_WEIGHT	* net_timing.path_criticality[ipin]
								  + (1 - ROUTE_PATH_WEIGHT) * net_timing.timing_criticality[ipin]; 
#else
			/* Pin criticality is based on only timing criticality. */
			pin_criticality = net_timing.timing_criticality[ipin];
#endif
			/* Currently, pin criticality is between 0 and 1. Now shift it downwards 
			by 1 - max_criticality (max_criticality is 0.99 by default, so shift down
			by 0.01) and cut off at 0.  This means that all pins with small criticalities 
			(<0.01) get criticality 0 and are ignored entirely, and everything
			else becomes a bit less critical. This effect becomes more pronounced if
			max_criticality is set lower. */
			assert(pin_criticality > -0.01 && pin_criticality < 1.01);
			pin_criticality = std::max(pin_criticality - (1.0 - params.max_criticality), 0.0);

			/* Take pin criticality to some power (1 by default). */
			pin_criticality = pow(pin_criticality, params.criticality_exp);
			
			/* Cut off pin criticality at max_criticality. */
			pin_criticality = std::min(pin_criticality, params.max_criticality);
		/*}*/

		net.sinks[ipin-1].criticality_fac = pin_criticality;
	}
}

bool route_tree_empty(const route_tree_t &rt)
{
	return rt.root == -1;
}

void check_net_route(const route_tree_t &rt, const RouteTreeNode &node, RRGraph &g, vector<int> &overused_rr_node)
{
	assert(node.properties.valid);

	const auto &rr_node = get_vertex(g, node.properties.rr_node);

	if (rr_node.properties.occ > rr_node.properties.capacity) {
		overused_rr_node.push_back(id(rr_node));
	}
	
	for_all_out_edges(rt.graph, node, [&rt, &g, &overused_rr_node] (const RouteTreeEdge &e) -> void {
			const auto &neighbor = get_target(rt.graph, e);
			check_net_route(rt, neighbor, g, overused_rr_node);
			});
}

void check_route_tree_internal(const route_tree_t &rt, const RouteTreeNode &node, RRGraph &g, vector<int> &visited_sinks, vector<int> &visited_nodes)
{
	assert(node.properties.valid);

	auto &rr_node = get_vertex(g, node.properties.rr_node);
	char buffer[256];
	if (rr_node.properties.type == SINK) {
		visited_sinks.push_back(id(rr_node));
		sprintf_rr_node(id(rr_node), buffer);
		/*zlog_level(delta_log, ROUTER_V2, "route_tree_check: %s\n", buffer);*/
	}

	visited_nodes.push_back(node.properties.rr_node);

	if (rr_node.properties.type == SOURCE) {
		rr_node.properties.recalc_occ += num_out_edges(rt.graph, node);
	} else {
		++rr_node.properties.recalc_occ;
	}

	for_all_out_edges(rt.graph, node, [&rt, &g, &visited_sinks, &visited_nodes] (const RouteTreeEdge &e) -> void {
			const auto &neighbor = get_target(rt.graph, e);
			check_route_tree_internal(rt, neighbor, g, visited_sinks, visited_nodes);
			});
}

void check_route_tree(const route_tree_t &rt, const net_t &net, RRGraph &g)
{
	const auto &rt_root = *route_tree_get_rt_node(rt, net.source.rr_node);
	vector<int> sinks;
	for (const auto &s : net.sinks) {
		sinks.push_back(s.rr_node);
	}
	vector<int> visited_sinks;

	vector<int> visited_nodes;
	check_route_tree_internal(rt, rt_root, g, visited_sinks, visited_nodes);
	set<int> visited_nodes_set(visited_nodes.begin(), visited_nodes.end());
	assert(visited_nodes_set.size() == visited_nodes.size());

	sort(visited_sinks.begin(), visited_sinks.end());
	sort(sinks.begin(), sinks.end());
	
	if (visited_sinks != sinks) {
		zlog_error(delta_log, "Error: Visited %d sinks out of %d sinks of net %d\n", visited_sinks.size(), sinks.size(), net.vpr_id);
		vector<int> only_in_required;
		vector<int> only_in_visited;
		vector<int> sym_difference;
		set_difference(sinks.begin(), sinks.end(), visited_sinks.begin(), visited_sinks.end(), back_inserter(only_in_required));
		set_difference(visited_sinks.begin(), visited_sinks.end(), sinks.begin(), sinks.end(), back_inserter(only_in_visited));
		/*set_symmetric_difference(sinks.begin(), sinks.end(), visited_sinks.begin(), visited_sinks.end(), back_inserter(sym_difference));*/
		/*assert(difference == sym_difference);*/

		write_graph(rt.graph, "error.dot",
				[] (const RouteTreeNode &rt_node) -> string {
					char buffer[256];
					sprintf_rr_node(rt_node.properties.rr_node, buffer);
					return buffer;
				},
				[] (const RouteTreeEdge &rt_edge) -> string {
					return "";
				},
				[] (const RouteTreeNode &rt_node) -> bool {
					return !rt_node.properties.valid;
				});


		zlog_error(delta_log, "Only in required: ");
		for (const int d : only_in_required) {
			zlog_error(delta_log, "%d ", d);
		}
		zlog_error(delta_log, "\n");

		zlog_error(delta_log, "Only in visited: ");
		for (const int d : only_in_visited) {
			zlog_error(delta_log, "%d ", d);
		}
		zlog_error(delta_log, "\n");

		zlog_error(delta_log, "Visited sinks: ");
		for (const int d : visited_sinks) {
			zlog_error(delta_log, "%d ", d);
		}
		zlog_error(delta_log, "\n");

		assert(false);
	}
}

void route_net_fast(RRGraph &g, const net_t &net, const route_parameters_t &params, route_state_t *state, route_tree_t &rt, trace_t &trace, t_net_timing &net_timing, perf_t *perf)
{
	std::priority_queue<route_state_t> heap;

	vector<int> modified;

	zlog_level(delta_log, ROUTER_V1, "Routing VPR net id %d\n", net.vpr_id);

	if (route_tree_empty(rt)) {
		/* special case */
		/*update_one_cost_internal(get_vertex(g, net.current_source.rr_node), 1, params.pres_fac);*/

		route_tree_set_source(rt, get_vertex(g, net.source.rr_node));
	} else {
		RouteTreeNode rt_root = get_vertex(rt.graph, rt.root);
		if (rt_root.properties.rr_node != net.source.rr_node) {
			char root[256];
			char source[256];
			sprintf_rr_node(rt_root.properties.rr_node, root);
			sprintf_rr_node(net.source.rr_node, source);
			zlog_error(delta_log, "Error: Non empty route tree node has root %s that is different from net source %s\n",
					root, source);
			assert(false);
		}
	}

	vector<sink_t> sorted_sinks = net.sinks;
	std::sort(sorted_sinks.begin(), sorted_sinks.end());

	char buffer[256];

	sprintf_rr_node(net.source.rr_node, buffer);
	zlog_level(delta_log, ROUTER_V1, "Source: %s\n", buffer);

	int isink = 0;
	for (const auto &sink : sorted_sinks) {
		sprintf_rr_node(sink.rr_node, buffer);
		zlog_level(delta_log, ROUTER_V1, "Sink %d: %s criticality: %g BB: %d-%d %d-%d\n", sink.id, buffer, sink.criticality_fac, sink.current_bounding_box.xmin, sink.current_bounding_box.xmax, sink.current_bounding_box.ymin, sink.current_bounding_box.ymax);
		route_tree_add_to_heap(rt, g, get_vertex(g, sink.rr_node), sink.criticality_fac, params.astar_fac, sink.current_bounding_box, heap, perf);

		bool found_sink = false;
		while (!heap.empty() && !found_sink) {
			auto item = heap.top();
			heap.pop();

			sprintf_rr_node(item.rr_node, buffer);
			zlog_level(delta_log, ROUTER_V2, "Current: %s prev: %d old_delay: %g new_delay: %g old_known: %g new_known: %g old_cost: %g new_cost: %g\n", buffer, item.prev_edge ? id(get_source(g, *item.prev_edge)) : -1, state[item.rr_node].delay, item.delay, state[item.rr_node].known_cost, item.known_cost, state[item.rr_node].cost, item.cost);

			if (item.rr_node == sink.rr_node) {
				state[item.rr_node] = item;
				modified.push_back(item.rr_node);
				found_sink = true;
			} else if (item.known_cost < state[item.rr_node].known_cost && item.cost < state[item.rr_node].cost) {
				/*if (route_tree_get_rt_node(rt, item.rr_node)) {*/
					/*sprintf_rr_node(item.rr_node, buffer);*/
					/*zlog_warn(delta_log, "Warning: Found a lower cost path to existing route tree node %s\n", buffer);*/
					/*assert(false);*/
				/*}*/
				const auto &sink_vertex = get_vertex(g, sink.rr_node);

				state[item.rr_node] = item;
				modified.push_back(item.rr_node);

				for (const auto &e_i : get_vertex(g, item.rr_node).edges) {
					auto &e = get_edge(g, e_i);
					auto &neighbor = get_target(g, e);

					char buffer[256];
					sprintf_rr_node(id(neighbor), buffer);
					zlog_level(delta_log, ROUTER_V3, "\tNeighbor: %s ", buffer);

					const auto &prop = neighbor.properties;

					if (prop.xhigh < sink.current_bounding_box.xmin
							|| prop.xlow > sink.current_bounding_box.xmax
							|| prop.yhigh < sink.current_bounding_box.ymin
							|| prop.ylow > sink.current_bounding_box.ymax) {
						zlog_level(delta_log, ROUTER_V3, "outside of bounding box\n");
						continue;
					}

					if (prop.type == IPIN
							&& (prop.xhigh != sink_vertex.properties.xhigh 
								|| prop.yhigh != sink_vertex.properties.yhigh)) {
						zlog_level(delta_log, ROUTER_V3, "not target IPIN\n");
						continue;
					}

					route_state_t new_item;

					new_item.rr_node = id(neighbor);
					new_item.prev_edge = &e;

					float unbuffered_upstream_R = item.upstream_R;
					float upstream_R = e.properties.R + neighbor.properties.R;
					if (!e.properties.buffered) {
						upstream_R += unbuffered_upstream_R;
					}
					new_item.upstream_R = upstream_R;

					new_item.delay = item.delay + get_delay(e, neighbor, unbuffered_upstream_R);

					float known_cost = item.known_cost + get_known_cost(g, e, sink.criticality_fac, unbuffered_upstream_R);
					new_item.known_cost = known_cost;

					float expected_cost = get_timing_driven_expected_cost(get_vertex(g, new_item.rr_node), sink_vertex, sink.criticality_fac, upstream_R);
					new_item.cost = known_cost + params.astar_fac * expected_cost;

					heap.push(new_item);

					if (perf) {
						++perf->num_heap_pushes;
					}

					zlog_level(delta_log, ROUTER_V3, "added with prev_edge: %X upstream_R: %g delay: %g known_cost: %g expected_cost: %g cost: %g\n", new_item.prev_edge, new_item.upstream_R, new_item.delay, new_item.known_cost, expected_cost, new_item.cost);
				}
			} else {
				zlog_level(delta_log, ROUTER_V3, "\tNot expanding neighbor because known_cost %g > %g or cost %g > %g\n", item.known_cost, state[item.rr_node].known_cost, item.cost, state[item.rr_node].cost);
			}
		}

		assert(found_sink);

		const auto &new_path = trace_add_path(trace, g, state, sink.rr_node, net.vpr_id);
		/* new_path.second-1 because we do not update the last node because it's the node in existing route tree (which was already updated previously */
		/*int new_path_first_rr_node = *(new_path.end()-1);*/
		/* if the first node is a connection to existing route, don't update its cost again */
		/*if (route_tree_get_rt_node(rt, new_path_first_rr_node)) {*/
			/*update_one_cost(g, new_path.begin(), new_path.end()-1, 1, params.pres_fac);*/
		/*} else {*/
			/*update_one_cost(g, new_path.begin(), new_path.end(), 1, params.pres_fac);*/
		/*}*/
		update_one_cost(g, new_path.begin(), new_path.end(), 1, params.pres_fac);

		/* important to update route tree only after updating cost because the previous check relies on whether the
		 * first node is already in the route tree */
		route_tree_add_path(rt, g, state, new_path, net.vpr_id);

		net_timing.delay[sink.id+1] = route_tree_get_rt_node(rt, sink.rr_node)->properties.delay;

		for (const auto &m : modified)  {
			state[m].known_cost = std::numeric_limits<float>::max();
			state[m].cost = std::numeric_limits<float>::max();
		}

		heap = std::priority_queue<route_state_t>();
	}

	/*delete [] state;*/

	/*check_route_tree(rt, net, g);*/
}

void route_net(RRGraph &g, const net_t &net, const route_parameters_t &params, route_state_t *state, route_tree_t &rt, trace_t &trace, t_net_timing &net_timing, perf_t *perf)
{
	std::priority_queue<route_state_t> heap;

	vector<int> modified;

	zlog_level(delta_log, ROUTER_V1, "Routing VPR net id %d\n", net.vpr_id);

	if (route_tree_empty(rt)) {
		/* special case */
		/*update_one_cost_internal(get_vertex(g, net.current_source.rr_node), 1, params.pres_fac);*/

		route_tree_set_source(rt, get_vertex(g, net.source.rr_node));
	} else {
		RouteTreeNode rt_root = get_vertex(rt.graph, rt.root);
		if (rt_root.properties.rr_node != net.source.rr_node) {
			char root[256];
			char source[256];
			sprintf_rr_node(rt_root.properties.rr_node, root);
			sprintf_rr_node(net.source.rr_node, source);
			zlog_error(delta_log, "Error: Non empty route tree node has root %s that is different from net source %s\n",
					root, source);
			assert(false);
		}
	}

	vector<sink_t> sorted_sinks = net.sinks;
	std::sort(sorted_sinks.begin(), sorted_sinks.end());

	char buffer[256];

	sprintf_rr_node(net.source.rr_node, buffer);
	zlog_level(delta_log, ROUTER_V1, "Source: %s\n", buffer);

	int isink = 0;
	for (const auto &sink : sorted_sinks) {
		sprintf_rr_node(sink.rr_node, buffer);
		zlog_level(delta_log, ROUTER_V1, "Sink %d: %s criticality: %g BB: %d-%d %d-%d\n", sink.id, buffer, sink.criticality_fac, sink.current_bounding_box.xmin, sink.current_bounding_box.xmax, sink.current_bounding_box.ymin, sink.current_bounding_box.ymax);
		route_tree_add_to_heap(rt, g, get_vertex(g, sink.rr_node), sink.criticality_fac, params.astar_fac, sink.current_bounding_box, heap, perf);

		bool found_sink = false;
		while (!heap.empty() && !found_sink) {
			auto item = heap.top();
			heap.pop();

			sprintf_rr_node(item.rr_node, buffer);
			zlog_level(delta_log, ROUTER_V2, "Current: %s prev: %d old_delay: %g new_delay: %g old_known: %g new_known: %g old_cost: %g new_cost: %g\n", buffer, item.prev_edge ? id(get_source(g, *item.prev_edge)) : -1, state[item.rr_node].delay, item.delay, state[item.rr_node].known_cost, item.known_cost, state[item.rr_node].cost, item.cost);

			if (item.rr_node == sink.rr_node) {
				state[item.rr_node] = item;
				modified.push_back(item.rr_node);
				found_sink = true;
			} else if (item.known_cost < state[item.rr_node].known_cost && item.cost < state[item.rr_node].cost) {
				/*if (route_tree_get_rt_node(rt, item.rr_node)) {*/
					/*sprintf_rr_node(item.rr_node, buffer);*/
					/*zlog_warn(delta_log, "Warning: Found a lower cost path to existing route tree node %s\n", buffer);*/
					/*assert(false);*/
				/*}*/
				const auto &sink_vertex = get_vertex(g, sink.rr_node);

				state[item.rr_node] = item;
				modified.push_back(item.rr_node);

				expand_neighbors(g, item, sink_vertex, sink.criticality_fac, params.astar_fac, heap, [&net, &sink, &sink_vertex] (const RRNode &v) -> bool {
					const auto &prop = v.properties;

					if (prop.xhigh < sink.current_bounding_box.xmin
							|| prop.xlow > sink.current_bounding_box.xmax
							|| prop.yhigh < sink.current_bounding_box.ymin
							|| prop.ylow > sink.current_bounding_box.ymax) {
					zlog_level(delta_log, ROUTER_V3, "outside of bounding box\n");
					return false;
					}

					if (prop.type == IPIN
							&& (prop.xhigh != sink_vertex.properties.xhigh 
								|| prop.yhigh != sink_vertex.properties.yhigh)) {
					zlog_level(delta_log, ROUTER_V3, "not target IPIN\n");
					return false;
					}

					/*if (net.sinks.size() >= HIGH_FANOUT_NET_LIM) {*/
					/*if (prop.xhigh < target_x - highfanout_rlim*/
							/*|| prop.xlow > target_x + highfanout_rlim*/
							/*|| prop.yhigh < target_y - highfanout_rlim*/
							/*|| prop.ylow > target_y + highfanout_rlim) {*/
						/*return false;*/
					/*}*/
					/*}*/
					return true;
				}, perf);
			} else {
				zlog_level(delta_log, ROUTER_V3, "\tNot expanding neighbor because known_cost %g > %g or cost %g > %g\n", item.known_cost, state[item.rr_node].known_cost, item.cost, state[item.rr_node].cost);
			}
		}

		assert(found_sink);

		const auto &new_path = trace_add_path(trace, g, state, sink.rr_node, net.vpr_id);
		/* new_path.second-1 because we do not update the last node because it's the node in existing route tree (which was already updated previously */
		/*int new_path_first_rr_node = *(new_path.end()-1);*/
		/* if the first node is a connection to existing route, don't update its cost again */
		/*if (route_tree_get_rt_node(rt, new_path_first_rr_node)) {*/
			/*update_one_cost(g, new_path.begin(), new_path.end()-1, 1, params.pres_fac);*/
		/*} else {*/
			/*update_one_cost(g, new_path.begin(), new_path.end(), 1, params.pres_fac);*/
		/*}*/
		update_one_cost(g, new_path.begin(), new_path.end(), 1, params.pres_fac);

		/* important to update route tree only after updating cost because the previous check relies on whether the
		 * first node is already in the route tree */
		route_tree_add_path(rt, g, state, new_path, net.vpr_id);

		net_timing.delay[sink.id+1] = route_tree_get_rt_node(rt, sink.rr_node)->properties.delay;

		for (const auto &m : modified)  {
			state[m].known_cost = std::numeric_limits<float>::max();
			state[m].cost = std::numeric_limits<float>::max();
		}

		heap = std::priority_queue<route_state_t>();
	}

	/*delete [] state;*/

	/*check_route_tree(rt, net, g);*/
}

void init_net_timing(vector<net_t> &nets, vector<net_t> &global_nets, t_net_timing *net_timing)
{
	for (const auto &net : nets) {
	   	net_timing[net.vpr_id].delay = new float[net.sinks.size() + 1];
	   	net_timing[net.vpr_id].timing_criticality = new float[net.sinks.size() + 1];
	   	net_timing[net.vpr_id].slack = new float[net.sinks.size() + 1];
	}
	for (const auto &net : global_nets) {
	   	net_timing[net.vpr_id].delay = new float[net.sinks.size() + 1];
	   	net_timing[net.vpr_id].timing_criticality = new float[net.sinks.size() + 1];
	   	net_timing[net.vpr_id].slack = new float[net.sinks.size() + 1];
	}

	float init_timing_criticality_val = 1;

	for (const auto &net : nets) {
		for (int ipin = 1; ipin <= net.sinks.size(); ipin++) {
			net_timing[net.vpr_id].timing_criticality[ipin] = init_timing_criticality_val;
#ifdef PATH_COUNTING
			net_timing[net.vpr_id].path_criticality[ipin] = init_timing_criticality_val;
#endif		
		}
	} 
	
	for (const auto &net : global_nets) {
		/* Set delay of global signals to zero. Non-global net 
		 * 			delays are set by update_net_delays_from_route_tree() 
		 * 						inside timing_driven_route_net(), which is only called
		 * 									for non-global nets. */
		for (int ipin = 1; ipin <= net.sinks.size(); ipin++) {
			net_timing[net.vpr_id].delay[ipin] = 0.;
		}
	}
}

void fine_grained_route_nets(RRGraph &g, vector<net_t> &nets, const route_parameters_t &params, route_tree_t *route_trees, t_net_timing *net_timing)
{
	int inet = 0;
	for (auto &net : nets) {
		/*route_net(g, net, params, route_trees[net.local_id], net_timing[net.vpr_id]);*/
		++inet;
	}
}

bool use_net_level_parallelism()
{
	return true;
}

void init_graph(RRGraph &g)
{
	extern s_rr_node *rr_node;
	extern int num_rr_nodes;
	extern t_rr_indexed_data *rr_indexed_data;
	extern struct s_switch_inf *switch_inf;

	add_vertex(g, num_rr_nodes);
	for (int i = 0; i < num_vertices(g); ++i) {
		auto &v = get_vertex(g, i);
		v.properties.type = rr_node[i].type;
		v.properties.inc_direction = rr_node[i].direction == INC_DIRECTION;
		v.properties.xlow = rr_node[i].xlow;
		v.properties.ylow = rr_node[i].ylow;
		v.properties.xhigh = rr_node[i].xhigh;
		v.properties.yhigh = rr_node[i].yhigh;

		extern struct s_grid_tile **grid;
		auto type = &grid[v.properties.xlow][v.properties.ylow];

		v.properties.real_xlow = rr_node[i].xlow;
		v.properties.real_ylow = rr_node[i].ylow;
		v.properties.real_xhigh = rr_node[i].xhigh;
		v.properties.real_yhigh = rr_node[i].ylow + type->offset;
		v.properties.R = rr_node[i].R;
		v.properties.C = rr_node[i].C;
		v.properties.cost_index = rr_node[i].cost_index;
		v.properties.capacity = rr_node[i].capacity;

		char buffer[256];
		sprintf_rr_node(i, buffer);
		zlog_debug(rr_log, "%s: real_xlow: %d real_xhigh: %d real_ylow: %d real_yhigh: %d\n", buffer, v.properties.real_xlow, v.properties.real_xhigh, v.properties.real_ylow, v.properties.real_yhigh);

		for (int j = 0; j < rr_node[i].num_edges; ++j) {
			auto &e = add_edge(g, v, get_vertex(g, rr_node[i].edges[j]));

			int si = rr_node[i].switches[j];
			
			e.properties.buffered = switch_inf[si].buffered; 
			e.properties.switch_delay = switch_inf[si].Tdel; 
			e.properties.R = switch_inf[si].R; 
		}
	}
	zlog_info(delta_log, "RR graph num vertices: %d\n", num_vertices(g));
	zlog_info(delta_log, "RR graph num edges: %d\n", num_edges(g));
}

template<typename PointA, typename PointB>
bounding_box_t get_bounding_box(const PointA &a, const PointB &b, int bb_factor)
{
	bounding_box_t bb;

	bb.xmin = std::min(a.x, b.x);
	bb.xmax = std::max(a.x, b.x);
	bb.ymin = std::min(a.y, b.y);
	bb.ymax = std::max(a.y, b.y);

	bb.xmin -= 1;
	bb.ymin -= 1;

	extern int nx;
	extern int ny;

	bb.xmin = std::max(bb.xmin - bb_factor, 0);
	bb.xmax = std::min(bb.xmax + bb_factor, nx + 1);
	bb.ymin = std::max(bb.ymin - bb_factor, 0);
	bb.ymax = std::min(bb.ymax + bb_factor, ny + 1);

	return bb;
}

void test_interval()
{
	using namespace boost::numeric;
	interval<int> i1(1,2);
	interval<int> i2(3,3);

	if (overlap(i1, i2)) {
		printf("1 overlap\n");
	}
}

void test_scheduler()
{
	using namespace boost::numeric;
	vector<pair<interval<int>, int>> intervals = {
		{ { 1, 5 }, 0 },
		{ { 6, 7 }, 1 },
		{ { 3, 10 }, 2 },
		{ { 8, 9 }, 3 }
	};
	vector<pair<interval<int>, int> *> chosen;
	max_independent_intervals(intervals, chosen);
	for (const auto &c : chosen) {
		printf("chosen: [%d, %d]\n", c->first.lower(), c->first.upper());
	}

	vector<bounding_box_t> bbs = {
		{ 1, 5, 3, 4 },
		{ 9, 10, 3, 4 },
		{ 1, 5, 3, 4 },
		{ 9, 10, 3, 4 }
	};
	vector<pair<const bounding_box_t *, int> *> bbs_ptr;

	for (const auto &b : bbs) {
		bbs_ptr.push_back(new pair<const bounding_box_t *, int>(&b, 0));
	}
	vector<pair<const bounding_box_t *, int> *> chosen_bbs;
	max_independent_rectangles(bbs_ptr, chosen_bbs);
	for (const auto &c : chosen_bbs) {
		printf("chosen: x %d-%d y %d-%d]\n", c->first->xmin, c->first->xmax, c->first->ymin, c->first->ymax);
	}
}

void init_nets(vector<net_t> &nets, vector<net_t> &global_nets, int bb_factor)
{
	extern struct s_net *clb_net;
	extern int num_nets;
	extern int **net_rr_terminals; /* [0..num_nets-1][0..num_pins-1] */
	extern struct s_rr_node *rr_node;
	extern struct s_bb *route_bb;
	extern struct s_block *block;

	int local_id = 0;
	for (int i = 0; i < num_nets; ++i) {
		net_t net;

		net.vpr_id = i;
		net.local_id = local_id;

		net.source.rr_node = net_rr_terminals[i][0];
		int b = clb_net[i].node_block[0];
		int p = clb_net[i].node_block_pin[0];
		net.source.x = block[b].x;
		net.source.y = block[b].y + block[b].type->pin_height[p];

		net.current_source = net.source;
		/*net.previous_source.rr_node = -1;*/
		/*net.previous_source.x = -1;*/
		/*net.previous_source.y = -1;*/
		/*net.previous_source_valid = false;*/
		 
		for (int j = 1; j <= clb_net[i].num_sinks; j++) {
			sink_t sink;

			sink.rr_node = net_rr_terminals[i][j];
			sink.id = j-1;
			int b = clb_net[i].node_block[j];
			int p = clb_net[i].node_block_pin[j];
			sink.x = block[b].x;
			sink.y = block[b].y + block[b].type->pin_height[p];
			sink.criticality_fac = std::numeric_limits<float>::max();
			/* initially all sink will be reached from the net's source */
			/* this is later updated during scheduling */
			sink.source = net.source;
			sink.bb_factor = bb_factor;
			sink.current_bounding_box = get_bounding_box(sink.source, sink, sink.bb_factor); 
			sink.congested_iterations = 0;

			net.sinks.push_back(sink);

			/*int inode = net_rr_terminals[i][j];*/
			/*extern struct s_block *block;*/
			/*assert(clb_net[i].node_block_pin[j] == rr_node[inode].ptc_num);*/
		}

		net.sink_routed.resize(net.sinks.size(), false);

		assert(net.sinks.size() > 0);
		net.has_sink = true;

		int rank = 0;
		for (auto &sink : net.sinks) {
			sink.distance_to_source_rank = rank++;
		}

		net.current_sink_index = 0;

		/*net.current_bounding_box = get_bounding_box(net.current_source, net.sinks[net.current_sink_index]);*/
		/*net.previous_sink_index = -1;*/
		/*net.previous_sink_valid = false;*/

		/*net.previous_bounding_box_valid = false;*/

		char buffer[256];
		zlog_debug(net_log, "Net %d sorted\n", i);
		sprintf_rr_node(net.source.rr_node, buffer);
		zlog_debug(net_log, "%s x: %d y: %d\n", buffer, net.source.x, net.source.y);
		zlog_debug(net_log, "Sorted sinks:\n");
		for (const auto &s : net.sinks) {
			sprintf_rr_node(s.rr_node, buffer);
			zlog_debug(net_log, "%s x: %d y: %d\n", buffer, s.x, s.y);
		}

		/*net.box.xmin = route_bb[i].xmin;*/
		/*net.box.ymin = route_bb[i].ymin;*/
		/*net.box.xmax = route_bb[i].xmax;*/
		/*net.box.ymax = route_bb[i].ymax;*/

		if (clb_net[i].is_global) {
			global_nets.push_back(net);

			zlog_info(net_log, "Global net %d\n", i);
		} else {
			nets.push_back(net);
			++local_id;
		}
	}

	/* update pointers */
	for (auto &net : nets) {
		net.source.net = &net;
		for (auto &sink : net.sinks) {
			sink.net = &net;
		}
	}
	/*int num_local_nets = local_id;*/
	/*for (auto &net : nets) {*/
		/*net.num_local_nets = num_local_nets;*/
		/*net.overlapping_nets = new bool[num_local_nets];*/
		/*net.non_overlapping_nets = new bool[num_local_nets];*/
	/*}*/
}

void write_metis_file(const graph_t<int, int> &dep_g, const char *filename)
{
	FILE *file = fopen("metis.graph", "w");
	assert(file);
	assert(num_edges(dep_g) % 2 == 0);
	fprintf(file, "%d %d\n", num_vertices(dep_g), num_edges(dep_g)/2);
	for_all_vertices(dep_g, [&dep_g, &file] (const vertex_t<int, int> &v) -> void {
			for_all_out_edges(dep_g, v, [&dep_g, &file] (const edge_t<int> &e) -> void {
					fprintf(file, "%d ", id(get_target(dep_g, e))+1);
					});
			fprintf(file, "\n");
			});
	fclose(file);
}

void load_overlapping_nets_vec(vector<net_t *> &nets)
{
	for (int i = 0; i < nets.size(); ++i) {
		nets[i]->overlapping_nets_vec.clear();
		nets[i]->non_overlapping_nets_vec.clear();
	}

	for (int i = 0; i < nets.size(); ++i) {
		for (int j = i+1; j < nets.size(); ++j) {
			if (box_overlap(nets[i]->sinks[nets[i]->current_sink_index].current_bounding_box, nets[j]->sinks[nets[j]->current_sink_index].current_bounding_box)) {
				nets[i]->overlapping_nets_vec.push_back(nets[j]);
				nets[j]->overlapping_nets_vec.push_back(nets[i]);
			} else {
				nets[i]->non_overlapping_nets_vec.push_back(nets[j]);
				nets[j]->non_overlapping_nets_vec.push_back(nets[i]);
			}
		}
		/*zlog_debug(net_log, "Net %i overlaps with %d nets\n", i, num_overlaps);*/
	}
}

void load_overlapping_nets(vector<net_t *> &nets)
{
	for (auto &net : nets) {
		/*for (int i = 0; i < n.num_local_nets; ++i) {*/
			/*n.overlapping_nets[i] = false;*/
			/*n.non_overlapping_nets[i] = false;*/
			/*n.num_overlapping_nets = 0;*/
			/*n.num_non_overlapping_nets = 0;*/
		/*}*/
		net->overlapping_nets = vector<bool>(nets.size(), false);
		net->non_overlapping_nets = vector<bool>(nets.size(), false);
		net->num_overlapping_nets = 0;
		net->num_non_overlapping_nets = 0;
	}

	for (int i = 0; i < nets.size(); ++i) {
		for (int j = i+1; j < nets.size(); ++j) {
			if (box_overlap(nets[i]->sinks[nets[i]->current_sink_index].current_bounding_box, nets[j]->sinks[nets[j]->current_sink_index].current_bounding_box)) {
				assert(nets[j]->current_local_id < nets[i]->overlapping_nets.size());
				nets[i]->overlapping_nets[nets[j]->current_local_id] = true;
				assert(nets[i]->current_local_id < nets[j]->overlapping_nets.size());
				nets[j]->overlapping_nets[nets[i]->current_local_id] = true;

				++nets[i]->num_overlapping_nets;
				++nets[j]->num_overlapping_nets;
			} else {
				assert(nets[j]->current_local_id < nets[i]->non_overlapping_nets.size());
				nets[i]->non_overlapping_nets[nets[j]->current_local_id] = true;
				assert(nets[i]->current_local_id < nets[j]->non_overlapping_nets.size());
				nets[j]->non_overlapping_nets[nets[i]->current_local_id] = true;

				++nets[i]->num_non_overlapping_nets;
				++nets[j]->num_non_overlapping_nets;
			}
		}
		/*zlog_debug(net_log, "Net %i overlaps with %d nets\n", i, num_overlaps);*/
	}
}

void init_overlapping_nets(vector<net_t> &nets, graph_t<int, int> &dep_g)
{
	/*add_vertex(dep_g, nets.size());*/
	/*for (auto &n : nets) {*/
		/*n.overlapping_nets = vector<bool>(nets.size(), false);*/
		/*n.non_overlapping_nets = vector<bool>(nets.size(), false);*/
	/*}*/

	/*for (int i = 0; i < nets.size(); ++i) {*/
		/*int num_overlaps = 0;*/
		/*for (int j = i+1; j < nets.size(); ++j) {*/
			/*if (box_overlap(bb_a, bb_b)) {*/
				/*++num_overlaps;*/
				/*add_edge(dep_g, i, j);*/
				/*add_edge(dep_g, j, i);*/

				/*nets[i].overlapping_nets[nets[i].id] = true;*/
				/*nets[j].overlapping_nets[nets[i].id] = true;*/
			/*} else {*/
				/*nets[i].non_overlapping_nets[nets[j].id] = true;*/
				/*nets[j].non_overlapping_nets[nets[i].id] = true;*/
			/*}*/
		/*}*/
		/*zlog_debug(net_log, "Net %i overlaps with %d nets\n", i, num_overlaps);*/
	/*}*/
}

typedef struct partition_t {
	vector<const net_t *> nets;
} partition_t;

void init_nets_partition(vector<net_t> &nets, const char *part_file, int num_parts, vector<partition_t> &partitions)
{
	FILE *file = fopen(part_file, "r");
	int i = 0;
	printf("Partitions\n");
	while (!feof(file)) {
		fscanf(file, "%d\n", &nets[i].pid);
		printf("Net %d pid %d\n", i, nets[i].pid);	
		assert(nets[i].pid < num_parts && nets[i].pid >= 0);
		partitions[nets[i].pid].nets.push_back(&nets[i]);
		++i;
	}
	
	assert(i == nets.size());
}

template<typename T>
inline bool in_set(const set<T> &s, const T &item)
{
	return s.find(item) != s.end();
}

net_t *schedule_net_greedy(const vector<net_t *> &nets, vector<bool> &scheduled_nets, vector<const net_t *> &current_time_scheduled_nets, int time)
{
	int max_num_non_overlapping_nets = std::numeric_limits<int>::min();	
	net_t *res = nullptr;

	/*zlog_debug(schedule_log, "-- Start scheduling for time %d\n", time);*/

	zlog_debug(schedule_log, "Scheduled nets: ");
	for (const auto &n : nets) {
		if (scheduled_nets[n->current_local_id]) {
			zlog_debug(schedule_log, "%d ", n->current_local_id);
		}
	}
	zlog_debug(schedule_log, "\n");

	zlog_debug(schedule_log, "Current time scheduled nets: ");
	for (const auto &n : current_time_scheduled_nets) {
		zlog_debug(schedule_log, "[%d BB: %d] ", n->current_local_id, get_bounding_box_area(n->current_bounding_box));
	}
	zlog_debug(schedule_log, "\n");

	for (const auto &net : nets) {
		bool net_scheduled = scheduled_nets[net->current_local_id];
		/*zlog_debug(schedule_log, "Checking net %d non_overlapping_nets %d net_scheduled %d ", net->id, num_non_overlapping_nets, net_scheduled ? 1 : 0);*/

		if (!net_scheduled && net->num_non_overlapping_nets > max_num_non_overlapping_nets) {
			bool overlap_scheduled_nets;

			/*if (current_time_scheduled_nets.empty()) {*/
				/*overlap_scheduled_nets = false;*/
			/*} else {*/
				/*overlap_scheduled_nets = any_of(current_time_scheduled_nets.begin(), current_time_scheduled_nets.end(), [&net] (const net_t *current_time_net) -> bool {*/
					/*return net->overlapping_nets[current_time_net->id];*/
					/*});*/
			/*}*/
			overlap_scheduled_nets = any_of(current_time_scheduled_nets.begin(), current_time_scheduled_nets.end(), [&net] (const net_t *current_time_net) -> bool {
					return net->overlapping_nets[current_time_net->current_local_id];
					});

			zlog_debug(schedule_log, "Net %d non_overlapping_nets %d not scheduled and overlap_scheduled_nets %d\n", net->current_local_id, net->num_non_overlapping_nets, overlap_scheduled_nets ? 1 : 0);

			if (!overlap_scheduled_nets) {
				max_num_non_overlapping_nets = net->num_non_overlapping_nets;
				res = net;
				zlog_debug(schedule_log, "\tRecording net %d\n", net->current_local_id);
			}
		}

		/*zlog_debug(schedule_log, "\n");*/
	}

	if (res) {
		res->schedule = time;
		scheduled_nets[res->current_local_id] = true;
		current_time_scheduled_nets.push_back(res);
		zlog_debug(schedule_log, "Scheduled net %d\n", res->current_local_id);
	}

	/*zlog_debug(schedule_log, "-- End scheduling for time %d\n", time);*/

	return res;
}

void schedule_nets_num_sinks(vector<net_t *> &nets, vector<vector<const net_t *>> &net_scheduled_at_time)
{
	sort(nets.begin(), nets.end(), [] (const net_t *a, const net_t *b) -> bool {
			return a->sinks.size() > b->sinks.size();
			});

	if (nets.size() > 0) {
		net_scheduled_at_time.push_back(vector<const net_t *>());
		for (const auto &net : nets) {
			net_scheduled_at_time[0].push_back(net);
		}
	}
}

void schedule_nets_greedy(vector<net_t *> &nets, vector<vector<const net_t *>> &net_scheduled_at_time)
{
	cpu_timer timer;
	timer.start();
	load_overlapping_nets(nets);
	timer.stop();
	cpu_times elapsed = timer.elapsed();
	/*zlog_info(delta_log, "Scheduling took %s\n", format(elapsed).c_str());*/
	zlog_info(delta_log, "Load overlapping net took %g\n", elapsed.wall / 1e9);

	/* we only schedule nets that has sinks to be routed */
	vector<bool> scheduled_nets(nets.size(), false);
	
	vector<const net_t *> current_time_scheduled_nets;

	timer.start();

	int time = 0;

	zlog_debug(schedule_log, "-- Start scheduling nets for time %d\n", time);
	net_t *scheduled_net = schedule_net_greedy(nets, scheduled_nets, current_time_scheduled_nets, time);
	if (!scheduled_net) {
		/* no more nets to schedule */
		return;
	}
	net_scheduled_at_time.push_back(vector<const net_t *>());
	net_scheduled_at_time[time].push_back(scheduled_net);
	int num_scheduled_nets = 1;

	while (num_scheduled_nets < nets.size()) {
		vector<net_t *> temp;
		for (const auto &net : nets) {
			if (scheduled_net->non_overlapping_nets[net->current_local_id]) {
				temp.push_back(net);
			}
		}
		scheduled_net = schedule_net_greedy(temp, scheduled_nets, current_time_scheduled_nets, time);
		if (!scheduled_net) {
			current_time_scheduled_nets.clear();
			++time;

			zlog_debug(schedule_log, "-- Start scheduling nets for time %d\n", time);
			scheduled_net = schedule_net_greedy(nets, scheduled_nets, current_time_scheduled_nets, time);
			assert(scheduled_net);
			net_scheduled_at_time.push_back(vector<const net_t *>());
			net_scheduled_at_time[time].push_back(scheduled_net);
		} else {
			net_scheduled_at_time[time].push_back(scheduled_net);
		}
		++num_scheduled_nets;
	}
	for (int i = 0; i < nets.size(); ++i) {
		/*zlog_info(schedule_log, "Net %d scheduled %d\n", i, scheduled_nets[i] ? 1 : 0);*/
	}
	assert(all_of(scheduled_nets.begin(), scheduled_nets.end(), [] (const bool &item) -> bool { return item; }));

	timer.stop();
	elapsed = timer.elapsed();
	/*zlog_info(delta_log, "Scheduling took %s\n", format(elapsed).c_str());*/
	zlog_info(delta_log, "Scheduling took %g\n", elapsed.wall / 1e9);

	/*delete [] scheduled_nets;*/
}

void verify_scheduling()
{
}

template<typename Nets>
net_t *schedule_net(const Nets &nets, set<int> &scheduled_net, set<int> &scheduled_partition, int time)
{
	int max_num_overlapping_nets = std::numeric_limits<int>::min();	
	net_t *res = nullptr;

	zlog_debug(schedule_log, "-- Start scheduling for time %d\n", time);

	zlog_debug(schedule_log, "Scheduled partitions: ");
	for (const auto &p : scheduled_partition) {
		zlog_debug(schedule_log, "%d ", p);
	}
	zlog_debug(schedule_log, "\n");

	zlog_debug(schedule_log, "Scheduled nets: ");
	for (const auto &n : scheduled_net) {
		zlog_debug(schedule_log, "%d ", n);
	}
	zlog_debug(schedule_log, "\n");

	for (const auto &net : nets) {
		bool net_scheduled = in_set(scheduled_net, net->local_id);
		bool part_scheduled = in_set(scheduled_partition, net->pid);

		zlog_debug(schedule_log, "Checking net %d pid %d non_overlapping_nets %d net_scheduled %d part_scheduled %d\n", net->local_id, net->pid, net->num_non_overlapping_nets, net_scheduled ? 1 : 0, part_scheduled ? 1 : 0);
		
		if (!net_scheduled && !part_scheduled
				&& net->num_non_overlapping_nets > max_num_overlapping_nets) {
			max_num_overlapping_nets = net->num_non_overlapping_nets ;
			res = net;
			zlog_debug(schedule_log, "\tRecording net %d\n", net->local_id);
		}
	}

	if (res) {
		res->schedule = time;
		scheduled_net.insert(res->local_id);
		scheduled_partition.insert(res->pid);
		zlog_debug(schedule_log, "Scheduled net %d pid %d\n", res->local_id, res->pid);
	}

	zlog_debug(schedule_log, "-- End scheduling for time %d\n", time);

	return res;
}

void schedule_nets(vector<net_t> &nets, int num_parts)
{
	for (auto &n : nets) {
		n.schedule = -1;
	}

	vector<net_t *> bootstrap;
	for (auto &n : nets) {
		bootstrap.push_back(&n);
	}
	set<int> scheduled_nets;
	set<int> scheduled_partitions;

	int time = 0;
	net_t *scheduled_net = schedule_net(bootstrap, scheduled_nets, scheduled_partitions, time);
	while (scheduled_nets.size() < nets.size()) {
		/*scheduled_net = schedule_net(scheduled_net->non_overlapping_nets, scheduled_nets, scheduled_partitions, time);*/
		if (!scheduled_net) {

			scheduled_partitions.clear();
			++time;
			scheduled_net = schedule_net(bootstrap, scheduled_nets, scheduled_partitions, time);
			assert(scheduled_net);
		} else if (scheduled_partitions.size() == num_parts) {
			scheduled_partitions.clear();
			++time;
		}
	}
}

/*void schedule_nets(vector<net_t> &nets, int num_parts)*/
/*{*/
	/*for (auto &n : nets) {*/
		/*n.schedule = -1;*/
	/*}*/

	/*int max_time;*/
	/*set<int> scheduled;*/
	/*for (int time = 0; time < max_time; ++time) {*/
		/*const net_t *cur_part_min = &nets[0];*/

		/*for (int i = 1; i < nets.size(); ++i) {*/
			/*auto iter = scheduled.find(nets[i].id);*/
			/*if (nets[i].non_overlapping_nets.size() < cur_part_min->non_overlapping_nets.size()*/
					/*&& iter == scheduled.end()) {*/
				/*cur_part_min = &nets[i];*/
			/*}*/
		/*}*/
		/*cur_part_min->schedule = time;*/
		/*scheduled.insert(cur_part_min->id);*/

		/*for (int i = 0; i < num_parts; ++i) { */
			/*const net_t *other_part_min = nullptr;*/
			/*for (int j = 0; j < cur_part_min->non_overlapping_nets.size(); ++j) {*/
				/*if (other_part_min == nullptr) {*/
					/*if (cur_part_min->non_overlapping_nets[j]->pid != cur_part_min->pid) {*/
						/*other_part_min = &cur_part_min->non_overlapping_nets[j];*/
					/*}*/
				/*} else if (cur_part_min->non_overlapping_nets[j]->pid != cur_part_min->pid && cur_part_min->non_overlapping_nets[j].size() < other_part_min->non_overlapping_nets.size()) {*/
					/*other_part_min = &cur_part_min->non_overlapping_nets[j];*/
				/*}*/
			/*}*/
		/*}*/
	/*}*/
/*}*/

void schedule_nets_fast(vector<net_t *> &nets, vector<vector<const net_t *>> &net_scheduled_at_time)
{
	/*load_overlapping_nets_vec(nets);*/

	cpu_timer timer;
	timer.start();

	/*sort(nets.begin(), nets.end(), [] (const net_t *a, const net_t *b) { return a->non_overlapping_nets_vec.size() < b->non_overlapping_nets_vec.size(); });*/
	sort(nets.begin(), nets.end(), [] (const net_t *a, const net_t *b) { return get_bounding_box_area(a->sinks[a->current_sink_index].current_bounding_box) > get_bounding_box_area(b->sinks[b->current_sink_index].current_bounding_box); });

	/*printf("Num overlapping_nets:\n");*/
	/*for (const auto &net : nets) {*/
		/*printf("Net %d: %d\n", net->current_local_id, net->overlapping_nets_vec.size());*/
	/*}*/

	vector<bool> scheduled_nets(nets.size(), false);
	net_scheduled_at_time.reserve(nets.size()); //worst case
	int num_scheduled_nets = 0;
	while (num_scheduled_nets < nets.size()) {
		net_scheduled_at_time.push_back(vector<const net_t *>());

		for (int i = 0; i < nets.size(); ++i) {
			if (!scheduled_nets[i]) {
				bool overlap_scheduled_nets = any_of(net_scheduled_at_time.back().begin(), net_scheduled_at_time.back().end(), [&i, &nets] (const net_t *scheduled_net) -> bool {
						return box_overlap(nets[i]->sinks[nets[i]->current_sink_index].current_bounding_box,
								scheduled_net->sinks[scheduled_net->current_sink_index].current_bounding_box);
					});

				if (!overlap_scheduled_nets) {
					scheduled_nets[i] = true;
					net_scheduled_at_time.back().push_back(nets[i]);
					++num_scheduled_nets;
				}
			}
		}
	}
	timer.stop();
	cpu_times elapsed = timer.elapsed();
	zlog_info(delta_log, "Scheduling took %g\n", elapsed.wall / 1e9);
	assert(all_of(scheduled_nets.begin(), scheduled_nets.end(), [] (const bool &item) -> bool { return item; }));

	for (const auto &at_time : net_scheduled_at_time) {
		verify_ind(at_time);
	}
}

void analyze_timing(t_net_timing *net_timing) 
{
	load_timing_graph_net_delays_new(net_timing); 
#ifdef HACK_LUT_PIN_SWAPPING
	do_timing_analysis_new(net_timing, FALSE, TRUE, FALSE);
#else
	do_timing_analysis_new(net_timing, FALSE, FALSE, FALSE);
#endif

	float critical_path_delay = get_critical_path_delay();
	zlog_info(delta_log, "Critical path: %g ns\n", critical_path_delay);
}

bool feasible_routing(const RRGraph &g)
{
	bool feasible = true;
	for_all_vertices_breakable(g, [&feasible] (const RRNode &v) -> bool {
			if (v.properties.occ > v.properties.capacity) {
				feasible = false;
			}
			return feasible;
			});

	return feasible;
}

void run_metis(const char *graph_filename, int num_partitions)
{
	char buffer[256];
	sprintf(buffer, "./gpmetis %s %d", graph_filename, num_partitions);
	FILE *pipe = popen(buffer, "r");
	while (!feof(pipe)) {
		fgets(buffer, 256, pipe);
		printf("%s", buffer);
	}
	pclose(pipe);
}

void set_current_source(net_t &net, route_tree_t &rt, const RRGraph &g, const sink_t &target, float astar_fac)
{
	/*net.previous_source = net.current_source;*/
	/*net.previous_source_valid = true;*/
	/*assert(net.previous_sink_valid);*/
	/*net.previous_bounding_box = get_bounding_box(net.previous_source, net.sinks[net.previous_sink_index]);*/
	/*net.previous_bounding_box_valid = true;*/

	/*if (route_tree_empty(rt)) {*/
		/*net.current_source = net.source;*/
	/*} else {*/
		/*std::priority_queue<route_state_t> heap;*/
		/*route_tree_add_to_heap(rt, g, get_vertex(g, target.rr_node), target.criticality_fac, astar_fac, net.previous_bounding_box, heap);*/
		/*net.current_source.rr_node = heap.top().rr_node;*/
		/*extern struct s_net *clb_net;*/
		/*extern struct s_block *block;*/
		/*int b = clb_net[net.id].node_block[0];*/
		/*int p = clb_net[net.id].node_block_pin[0];*/
		/*net.current_source.x = block[b].x;*/
		/*net.current_source.y = block[b].y + block[b].type->pin_height[p];*/
	/*}*/
	/*net.current_bounding_box = get_bounding_box(net.current_source, net.sinks[net.current_sink_index]);*/
}

void adjust_bounding_box(net_t &net)
{
	assert(net.current_sink_index >= 0 && net.current_sink_index < net.sinks.size());
	auto &current_sink = net.sinks[net.current_sink_index];

	if (current_sink.current_bounding_box.xmin != current_sink.previous_bounding_box.xmin ||
			current_sink.current_bounding_box.ymin != current_sink.previous_bounding_box.ymin || 
			current_sink.current_bounding_box.xmax != current_sink.previous_bounding_box.xmax ||
			current_sink.current_bounding_box.ymax != current_sink.previous_bounding_box.ymax) {
		char buffer[256];
		sprintf_rr_node(current_sink.rr_node, buffer);
		zlog_level(delta_log, ROUTER_V2, "Net %d %s bounding box differs. Prev: %d-%d %d-%d Cur: %d-%d %d-%d\n", net.vpr_id, buffer,
				current_sink.current_bounding_box.xmin, current_sink.current_bounding_box.ymin,
				current_sink.current_bounding_box.xmax, current_sink.current_bounding_box.ymax,
				current_sink.previous_bounding_box.xmin, current_sink.previous_bounding_box.ymin,
				current_sink.previous_bounding_box.xmax, current_sink.previous_bounding_box.ymax);
	}

	current_sink.current_bounding_box.xmin = std::min(current_sink.current_bounding_box.xmin, current_sink.previous_bounding_box.xmin);
	current_sink.current_bounding_box.ymin = std::min(current_sink.current_bounding_box.ymin, current_sink.previous_bounding_box.ymin);
	current_sink.current_bounding_box.xmax = std::max(current_sink.current_bounding_box.xmax, current_sink.previous_bounding_box.xmax);
	current_sink.current_bounding_box.ymax = std::max(current_sink.current_bounding_box.ymax, current_sink.previous_bounding_box.ymax);
}

int adjust_bounding_box_2(net_t &net)
{
	int num_changed = 0;
	for (auto &sink : net.sinks) {
		if (!net.sink_routed[sink.id]) {
			zlog_level(delta_log, ROUTER_V2, "Adjust bounding box for net %d sink %d. Current x %d-%d y %d-%d Previous x %d-%d y %d-%d\n", net.vpr_id, sink.id, sink.current_bounding_box.xmin, sink.current_bounding_box.xmax, sink.current_bounding_box.ymin, sink.current_bounding_box.ymax, sink.previous_bounding_box.xmin, sink.previous_bounding_box.xmax, sink.previous_bounding_box.ymin, sink.previous_bounding_box.ymax);

			bool changed = false;
			if (sink.current_bounding_box.xmin != sink.previous_bounding_box.xmin) {
				changed = true;
			}
			if (sink.current_bounding_box.xmax != sink.previous_bounding_box.xmax) {
				changed = true;
			}
			if (sink.current_bounding_box.ymin != sink.previous_bounding_box.ymin) {
				changed = true;
			}
			if (sink.current_bounding_box.ymax != sink.previous_bounding_box.ymax) {
				changed = true;
			}

			if (changed) {
				++num_changed;
			}

			sink.current_bounding_box.xmin = std::min(sink.current_bounding_box.xmin, sink.previous_bounding_box.xmin);
			sink.current_bounding_box.ymin = std::min(sink.current_bounding_box.ymin, sink.previous_bounding_box.ymin);
			sink.current_bounding_box.xmax = std::max(sink.current_bounding_box.xmax, sink.previous_bounding_box.xmax);
			sink.current_bounding_box.ymax = std::max(sink.current_bounding_box.ymax, sink.previous_bounding_box.ymax);
		}
	}
	return num_changed;
}

void update_sink_bounding_boxes(net_t &net, route_tree_t &rt, const RRGraph &g, float astar_fac, perf_t &perf)
{
	/*net.previous_sink_index = net.current_sink_index;*/
	/*net.previous_source = net.current_source;*/
	/*net.previous_bounding_box = get_bounding_box(net.previous_source, net.sinks[net.previous_sink_index]);*/

	/*current_sink.previous_bounding_box = current_sink.current_bounding_box;*/

	assert(!route_tree_empty(rt));

	net.has_sink = false;
	for (auto &sink : net.sinks) {
		if (!net.sink_routed[sink.id]) {
			std::priority_queue<route_state_t> heap;
			extern int nx;
			extern int ny;
			bounding_box_t fpga_bounding_box;
			fpga_bounding_box.xmin = 0;
			fpga_bounding_box.xmax = nx+2;
			fpga_bounding_box.ymin = 0;
			fpga_bounding_box.ymax = ny+2;
			route_tree_add_to_heap(rt, g, get_vertex(g, sink.rr_node), sink.criticality_fac, astar_fac, fpga_bounding_box, heap, &perf);
			const RRNode *current_rr_node = &get_vertex(g, heap.top().rr_node);
			while (current_rr_node->properties.type != CHANX && current_rr_node->properties.type != CHANY) {
				heap.pop();
				current_rr_node = &get_vertex(g, heap.top().rr_node);
			}
			sink.source.rr_node = id(*current_rr_node);
			if (current_rr_node->properties.inc_direction) {
				sink.source.x = current_rr_node->properties.xlow; 
				sink.source.y = current_rr_node->properties.ylow;  
			} else {
				sink.source.x = current_rr_node->properties.xhigh; 
				sink.source.y = current_rr_node->properties.yhigh;  
			}
			/*sink.previous_bounding_box = sink.current_bounding_box;*/
			sink.current_bounding_box = get_bounding_box(sink.source, sink, sink.bb_factor);
			char source_s[256];
			char sink_s[256];
			sprintf_rr_node(sink.rr_node, sink_s);
			sprintf_rr_node(sink.source.rr_node, source_s);
			zlog_level(delta_log, ROUTER_V2, "Updating source of sink %s to %s. New bouding box %d-%d %d-%d\n", sink_s, source_s,
					sink.current_bounding_box.xmin, sink.current_bounding_box.xmax, sink.current_bounding_box.ymin, sink.current_bounding_box.ymax);
			net.has_sink = true;
		}
	}
}

void set_to_next_sink(net_t &net, route_tree_t &rt, const RRGraph &g, float astar_fac, perf_t &perf)
{
	/*net.previous_sink_index = net.current_sink_index;*/
	/*net.previous_source = net.current_source;*/
	/*net.previous_bounding_box = get_bounding_box(net.previous_source, net.sinks[net.previous_sink_index]);*/

	assert(net.current_sink_index >= 0 && net.current_sink_index < net.sinks.size());

	auto &old_sink = net.sinks[net.current_sink_index];
	old_sink.previous_bounding_box = old_sink.current_bounding_box;

	net.current_sink_index++;
	if (net.current_sink_index >= net.sinks.size()) {
		net.has_sink = false;
	} else {
		net.has_sink = true;
	}

	if (!net.has_sink) {
		return;
	}

	zlog_level(delta_log, ROUTER_V2, "Set to next sink for net %d\n", net.vpr_id);

	assert(!route_tree_empty(rt));
	
	std::priority_queue<route_state_t> heap;
	extern int nx;
	extern int ny;
	bounding_box_t fpga_bounding_box;
	fpga_bounding_box.xmin = 0;
	fpga_bounding_box.xmax = nx+2;
	fpga_bounding_box.ymin = 0;
	fpga_bounding_box.ymax = ny+2;
	auto &new_sink = net.sinks[net.current_sink_index];
	/*route_tree_add_to_heap(rt, g, get_vertex(g, new_sink.rr_node), new_sink.criticality_fac, astar_fac, fpga_bounding_box, heap, &perf);*/
	/*const RRNode *current_rr_node = &get_vertex(g, heap.top().rr_node);*/
	/*zlog_info(delta_log, "Set to next sink for net %d\n", net.vpr_id);*/
	/*while (current_rr_node->properties.type != CHANX && current_rr_node->properties.type != CHANY) {*/
		/*heap.pop();*/
		/*current_rr_node = &get_vertex(g, heap.top().rr_node);*/
		/*zlog_info(delta_log, "Current: %d", net.vpr_id);*/
	/*}*/

	auto &sink_vertex = get_vertex(g, new_sink.rr_node);
	/*if (sink_vertex.properties.xlow != sink_vertex.properties.xhigh || sink_vertex.properties.ylow != sink_vertex.properties.yhigh) {*/

		/*assert(false);*/
	/*}*/

	int min_area = std::numeric_limits<int>::max();
	const RRNode *new_from = nullptr;
	bounding_box_t new_bounding_box;
	for_all_vertices(rt.graph, [&min_area, &new_from, &new_bounding_box, &new_sink, &g] (const RouteTreeNode &rt_node) -> void {
			if (!rt_node.properties.valid) {
				return;
			}

			const auto &from_node = get_vertex(g, rt_node.properties.rr_node);

			if (from_node.properties.type == IPIN || from_node.properties.type == SINK) {
				return;
			}

			struct point {
				int x;
				int y;
			} bottom_left, top_right;
			
			bottom_left.x = std::min(new_sink.x, from_node.properties.xlow);
			bottom_left.y = std::min(new_sink.y, from_node.properties.ylow);
			top_right.x = std::max(new_sink.x, from_node.properties.xhigh);
			top_right.y = std::max(new_sink.y, from_node.properties.yhigh);

			new_bounding_box = get_bounding_box(bottom_left, top_right, new_sink.bb_factor);

			int area = get_bounding_box_area(new_bounding_box);
			if (area < min_area) {
				new_from = &from_node;
				min_area = area;
			}
			});

	char buffer[256];

	if (new_from == nullptr) {
		sprintf_rr_node(new_sink.rr_node, buffer);
		zlog_error(delta_log, "Error: Unable to find node in route tree that is closest to %s\n", buffer);
		assert(false);
	}

	sprintf_rr_node(id(*new_from), buffer);
	zlog_level(delta_log, ROUTER_V2, "New source %s\n", buffer);
	/*source_t new_source;*/
	/*int min_distance = std::numeric_limits<int>::max();*/
	/*for_all_vertices(rt.graph, [] (const RouteTreeNode &v) -> void {*/
			/*if (v.properties.valid) {*/
			/*v.properties.xlow*/
			/*}*/
			/*});*/
	/*const sink_t &target = net.sinks[net.current_sink_index];*/
	/*route_tree_add_to_heap(rt, g, get_vertex(g, target.rr_node), target.criticality_fac, astar_fac, fpga_bounding_box, heap, NULL);*/
	/*net.current_source.rr_node = heap.top().rr_node;*/
	/*assert(get_vertex(g, net.current_source.rr_node).properties.type == SOURCE);*/
	/*extern struct s_net *clb_net;*/
	/*extern struct s_block *block;*/
	/*int b = clb_net[net.vpr_id].node_block[0];*/
	/*int p = clb_net[net.vpr_id].node_block_pin[0];*/
	/* TODO: fix the code below, the calculation of coordinates is wrong */
	/*net.current_source.x = block[b].x;*/
	/*net.current_source.y = block[b].y + block[b].type->pin_height[p];*/

	new_sink.current_bounding_box = new_bounding_box;

	/*adjust_bounding_box(net);*/
}

void reset_current_source_sink_2(net_t &net)
{
	net.current_sink = nullptr;

	for (auto &sink : net.sinks) {
		sink.source = net.source;
		sink.current_bounding_box = get_bounding_box(sink.source, sink, sink.bb_factor);
	}

	fill(net.sink_routed.begin(), net.sink_routed.end(), false);
	net.has_sink = true;
}

void reset_current_source_sink(net_t &net)
{
	/*net.previous_source_valid = false;*/
	/*net.previous_sink_valid = false;*/
	/*net.previous_bounding_box_valid = false;*/

	net.current_source = net.source;
	/*net.previous_source.rr_node = -1;*/
	/*net.previous_source.x = -1;*/
	/*net.previous_source.y = -1;*/

	net.current_sink_index = 0;
	/*net.previous_sink_index = -1;*/
	assert(net.sinks.size() > 0);
	net.has_sink = true;

	auto &sink = net.sinks[net.current_sink_index];
	char buffer[256];
	sprintf_rr_node(sink.rr_node, buffer);
	zlog_level(delta_log, ROUTER_V2, "Resetting %s bounding box with bb_factor %d\n", buffer, sink.bb_factor);

	sink.current_bounding_box = get_bounding_box(net.current_source, sink, sink.bb_factor);
}

/*void temp_route_nets(RRGraph &g, vector<net_t *> &nets, const route_parameters_t &params, route_tree_t *route_trees, t_net_timing *net_timing)*/
/*{*/
	/*[> rip up previous routing <]*/
	/*if (route_tree_initiated(rt)) {*/
		/*update_one_cost(g, rt, route_tree_get_rt_node(rt, net.source.rr_node), -1, params.pres_fac);*/
	/*} */

	/*route_tree_clear(rt);*/
/*}*/

void trace_rip_up_net(trace_t &trace, RRGraph &g, float pres_fac)
{
	for (const auto &segment : trace.segments) {
		update_one_cost(g, segment.second.begin(), segment.second.end(), -1, pres_fac);
	}
	trace.segments.clear();
	/*trace.first_sink_rr_node = -1;*/
	trace.paths_starting_with_source.clear();
	trace.existing_nodes.clear();
}

void trace_rip_up_segment(trace_t &trace, RRGraph &g, int sink_rr_node, float pres_fac)
{
	/*const vector<int> *path = route_tree_get_path_to_sink(rt, sink_rr_node);*/
	/*if (path) {*/
		/*update_one_cost(g, path->begin(), path->end()-1, -1, pres_fac);*/
		/*route_tree_remove_path(rt, sink_rr_node);*/
	/*}*/
	auto iter = trace.segments.find(sink_rr_node);
	if (iter != trace.segments.end()) {
		update_one_cost(g, iter->second.begin(), iter->second.end(), -1, pres_fac);
		/*if (get_vertex(g, iter->second.back()).properties.type == SOURCE) {*/
			/*--trace.num_sources;*/
			/*if (trace.num_sources != 0) {*/
				/*char buffer[256];*/
				/*sprintf_rr_node(sink_rr_node, buffer);*/
				/*zlog_error(delta_log, "Error: Ripped up a segment for %s but num_sources (%d) != 0\n", buffer, trace.num_sources);*/
				/*assert(false);*/
			/*}*/
		/*}*/

		for (const auto &node : iter->second) {
			char buffer[256];
			sprintf_rr_node(node, buffer);
			if (trace.existing_nodes.find(node) == trace.existing_nodes.end()) {
				/* this is possible because SOURCE with capacity > 1 only appears
				 * once in the set */
				const auto &prop = get_vertex(g, node).properties;
				if (prop.type == SOURCE && prop.capacity > 1 && !trace.paths_starting_with_source.empty()) {
					zlog_level(delta_log, ROUTER_V2, "Ripping up %s with capacity > 1\n", buffer);
				} else {
					zlog_error(delta_log, "Error: Ripping up non-existing node %s\n", buffer);
					assert(false);
				}
			} else {
				zlog_level(delta_log, ROUTER_V2, "Ripping up node %s from the trace\n", buffer);
				trace.existing_nodes.erase(node);
			}
		}

		if (get_vertex(g, iter->second.back()).properties.type == SOURCE) {
			auto ptr_iter = find_if(trace.paths_starting_with_source.begin(), trace.paths_starting_with_source.end(), [&sink_rr_node] (const Segment *path) {
					return path->front() == sink_rr_node;
					});
			assert(ptr_iter != trace.paths_starting_with_source.end());
			trace.paths_starting_with_source.erase(ptr_iter);
		}

		trace.segments.erase(iter);
	}
}

/*bool route_tree_has_path(const route_tree_t &rt, int sink_rr_node)*/
/*{*/
	/*return route_tree_get_path_to_sink(rt, sink_rr_node) ? true : false;*/
/*}*/

void test_graph()
{
	graph_t<int, int> g;
	add_vertex(g, 5);
	add_edge(g, get_vertex(g, 0), get_vertex(g, 1));
	add_edge(g, get_vertex(g, 0), get_vertex(g, 2));
	add_edge(g, get_vertex(g, 1), get_vertex(g, 3));
	add_edge(g, get_vertex(g, 1), get_vertex(g, 4));
	remove_edge(g, get_edge(g, 0));
}

void test_rtree();

void sort_sinks(vector<net_t> &nets)
{
	for (auto net : nets) {
		std::sort(net.sinks.begin(), net.sinks.end(), [&net] (const sink_t &a, const sink_t &b) -> bool {
				int a_to_src_dist = abs(a.x - net.source.x) + abs(a.y - net.source.y);
				int b_to_src_dist = abs(b.x - net.source.x) + abs(b.y - net.source.y);
				return a_to_src_dist > b_to_src_dist;
				});
	}
}

bool partitioning_route_bounding_box(t_router_opts *opts)
{
	const int max_iter = 50;

	init_logging();

	/*test_scheduler();*/

	/*test_rtree();*/

	/*test_graph();*/

	RRGraph g;
	init_graph(g);

	vector<net_t> nets;
	vector<net_t> global_nets;
	init_nets(nets, global_nets, opts->bb_factor);

	vector<net_t *> nets_ptr(nets.size());
	for (int i = 0; i < nets.size(); ++i) {
		nets_ptr[i] = &nets[i];
	}
	extern s_bb *route_bb;
	sort(nets_ptr.begin(), nets_ptr.end(), [] (const net_t *a, const net_t *b) -> bool {
			return get_bounding_box_area(route_bb[a->vpr_id]) > get_bounding_box_area(route_bb[b->vpr_id]);
			});
	int rank = 0;
	for (auto &net : nets_ptr) {
		net->bb_area_rank = rank++;
	}

	extern int nx;
	extern int ny;
	zlog_info(delta_log, "Num nets: %d Num global nets: %d nx: %d ny: %d\n", nets.size(), global_nets.size(), nx, ny);
	
	t_net_timing *net_timing = new t_net_timing[nets.size()+global_nets.size()];
	init_net_timing(nets, global_nets, net_timing);

	vector<route_tree_t> route_trees(nets.size());
	for (int i = 0; i < nets.size(); ++i) {
		route_tree_init(route_trees[i]);
	}

	vector<trace_t> traces(nets.size());
	for (int i = 0; i < nets.size(); ++i) {
		trace_init(traces[i]);
	}

	route_state_t *state = new route_state_t[num_vertices(g)];

	for (int i = 0; i < num_vertices(g); ++i) {
		state[i].rr_node = -1;
		state[i].known_cost = std::numeric_limits<float>::max();
		state[i].cost = std::numeric_limits<float>::max();
		state[i].prev_edge = nullptr;
		state[i].upstream_R = -1;
		state[i].delay = std::numeric_limits<float>::max();
	}


	/*route_tree_t **route_trees = new route_tree_t *[nets.size()+global_nets.size()];*/
	/*for (int i = 0; i < nets.size(); ++i) {*/
		/*route_trees[nets[i].id] = new route_tree_t[nets[i].sinks.size()+1];*/
		/*for (int j = 0; j <= nets[i].sinks.size(); ++j) {*/
			/*route_tree_init(route_trees[nets[i].id][j]);*/
		/*}*/
	/*}*/
	/*for (int i = 0; i < global_nets.size(); ++i) {*/
		/*route_trees[nets[i].id] = new route_tree_t[nets[i].sinks.size()+1];*/
		/*for (int j = 0; j <= nets[i].sinks.size(); ++j) {*/
			/*route_tree_init(route_trees[nets[i].id][j]);*/
		/*}*/
	/*}*/

	route_parameters_t params;
	params.criticality_exp = opts->criticality_exp;
	params.astar_fac = opts->astar_fac;
	params.max_criticality = opts->max_criticality;
	params.bend_cost = opts->bend_cost;
	params.pres_fac = opts->first_iter_pres_fac; /* Typically 0 -> ignore cong. */

	for_all_vertices(g, [] (RRNode &v) -> void {
			v.properties.acc_cost = 1;
			v.properties.pres_cost = 1;
			v.properties.occ = 0;
			});

	bool routed = false;
	for (int iter = 0; iter < max_iter && !routed; ++iter) {
		char iter_s[64];
		sprintf(iter_s, "%d", iter);
		assert(!zlog_put_mdc("iter", iter_s));

		zlog_info(delta_log, "Routing iteration: %d\n", iter);

		cpu_timer iter_timer;
		iter_timer.start();

		for (int i = 0; i < nets.size(); ++i) {
			route_tree_clear(route_trees[i]);
		}

		for (auto &net : nets) {
			update_sink_criticalities(net, net_timing[net.vpr_id], params);
		}

		cpu_timer schedule_timer;
		schedule_timer.start();

		vector<net_t *> nets_to_route;
		int current_local_id = 0;
		for (auto &net : nets) {
			if (net.has_sink) {
				nets_to_route.push_back(&net);
				net.current_local_id = current_local_id++;
			}
		}
		if (iter > 0) {
			for (auto &net : nets_to_route) {
				adjust_bounding_box_2(*net);
			}
		}
		vector<pair<sink_t *, net_t *>> scheduled_sinks;
		schedule_nets_bounding_box(nets_to_route, scheduled_sinks);
		zlog_info(delta_log, "Num scheduled sink: %d\n", scheduled_sinks.size());

		schedule_timer.stop();
		nanosecond_type schedule_time = schedule_timer.elapsed().wall;

		int time = 0;

		perf_t perf;
		perf.num_heap_pushes = 0;

		nanosecond_type route_time = 0;

		int total_num_sinks = 0;

		while (scheduled_sinks.size() > 0) {
			zlog_info(delta_log, "Time: %d\n", time);
			cpu_timer route_timer;
			route_timer.start();

			/* this for loop can be parallelized with parallel_for */
			for (const auto &sink_pair : scheduled_sinks) {
				net_t *net;
				sink_t *sink;
				std::tie(sink, net) = sink_pair;
				net_t virtual_net = *net;
				virtual_net.sinks.clear();
				virtual_net.sinks.push_back(*sink);
				virtual_net.current_source = sink->source;
				virtual_net.current_bounding_box = sink->current_bounding_box;
				/*if (time == 0) {*/
				/*trace_rip_up_net(traces[net->local_id], g, params.pres_fac);*/
				/*}*/
				trace_rip_up_segment(traces[net->local_id], g, sink->rr_node, params.pres_fac);
				/*printf("Net: %d\n", virtual_net.vpr_id);*/
				route_net(g, virtual_net, params, state, route_trees[net->local_id], traces[net->local_id], net_timing[net->vpr_id], &perf);

				net->sink_routed[sink->id] = true;
				/* save this so that we can use it in next iteration */
				sink->previous_bounding_box = sink->current_bounding_box;
				/*route_tree_merge(route_trees[net->id][0], route_trees[net->id][virtual_net.sinks[0].id]);*/
				++total_num_sinks;
			}

			/* TODO: fix the loop range and update_sink_bounding_boxes */
			for (auto &sink_pair : scheduled_sinks) {
				net_t *net = sink_pair.second;
				update_sink_bounding_boxes(*net, route_trees[net->local_id], g, params.astar_fac, perf);
			}

			route_timer.stop();
			cpu_times elapsed = route_timer.elapsed();
			route_time += elapsed.wall;
			/*zlog_info(delta_log, "Scheduling took %s\n", format(elapsed).c_str());*/
			zlog_info(delta_log, "Routing took %g\n", elapsed.wall / 1e9);

			schedule_timer.start();
		
			nets_to_route.clear();
			int current_local_id = 0;
			for (auto &net : nets) {
				if (net.has_sink) {
					nets_to_route.push_back(&net);
					net.current_local_id = current_local_id++;
				}
			}
			scheduled_sinks.clear();
			if (iter > 0) {
				for (auto &net : nets_to_route) {
					adjust_bounding_box_2(*net);
				}
			}
			schedule_nets_bounding_box(nets_to_route, scheduled_sinks);

			schedule_timer.stop();
			schedule_time += schedule_timer.elapsed().wall;

			zlog_info(delta_log, "Num scheduled sink: %d\n", scheduled_sinks.size());

			++time;
		}

		zlog_info(delta_log, "Average concurrency: %g Num heap pushes: %d\n", (float)total_num_sinks/time, perf.num_heap_pushes);

		bool fail = false;
		for (const auto &net : nets) {
			for (int i = 0; i < net.sinks.size(); ++i) {
				if (!net.sink_routed[net.sinks[i].id]) {
					char buffer[256];
					sprintf_rr_node(net.sinks[i].rr_node, buffer);
					zlog_error(delta_log, "Error: %s of net %d is not routed\n", buffer, net.vpr_id);
					fail = true;
				}
			}
		}
		assert(!fail);

		for_all_vertices(g, [] (RRNode &v) -> void { v.properties.recalc_occ = 0; });

		for (const auto &net : nets) {
			check_route_tree(route_trees[net.local_id], net, g);
		}

		bool valid = true;
		for_all_vertices(g, [&valid] (RRNode &v) -> void {
				char buffer[256];
				sprintf_rr_node(id(v), buffer);
				if (v.properties.recalc_occ != v.properties.occ) {
					zlog_error(delta_log, "Node %s occ mismatch, recalc: %d original: %d\n", buffer, v.properties.recalc_occ, v.properties.occ);
					valid = false;
				}
				});
		assert(valid);

		for (const auto &net : nets) {
			vector<int> overused_rr_node;
			check_net_route(route_trees[net.local_id], get_vertex(route_trees[net.local_id].graph, route_trees[net.local_id].root), g, overused_rr_node);
			if (!overused_rr_node.empty()) {
				zlog_level(delta_log, ROUTER_V3, "Net %d bb_rank %d overused nodes:\n", net.vpr_id, net.bb_area_rank);
				for (const auto &item : overused_rr_node) {
					zlog_level(delta_log, ROUTER_V3, "%d ", item);
				}
				zlog_level(delta_log, ROUTER_V3, "\n");
			}
		}

		if (feasible_routing(g)) {
			zlog_info(delta_log, "Routed in %d iterations\n", iter+1);
			routed = true;
		} else {
			int num_overused_nodes = 0;
			for_all_vertices(g, [&num_overused_nodes] (const RRNode &v) -> void {
					if (v.properties.occ > v.properties.capacity) {
					++num_overused_nodes;
					}
					});
			zlog_info(delta_log, "Num overused nodes: %d/%d (%.2f)\n", num_overused_nodes, num_vertices(g), num_overused_nodes*100.0/num_vertices(g));

			if (iter == 0) {
				params.pres_fac = opts->initial_pres_fac;
				update_costs(g, params.pres_fac, 0);
			} else {
				params.pres_fac *= opts->pres_fac_mult;

				/* Avoid overflow for high iteration counts, even if acc_cost is big */
				params.pres_fac = std::min(params.pres_fac, static_cast<float>(HUGE_POSITIVE_FLOAT / 1e5));

				update_costs(g, params.pres_fac, opts->acc_fac);
			}

			for (auto &net : nets) { 
				reset_current_source_sink_2(net);
			}
		}

		analyze_timing(net_timing);

		iter_timer.stop();
		cpu_times elapsed = iter_timer.elapsed();

		zlog_info(delta_log, "Iteration time: %g s. Route time: %g s. Schedule time: %g s.\n", elapsed.wall / 1e9, route_time / 1e9, schedule_time / 1e9);
	}
	return routed;
}

void dump_congestion_map(const RRGraph &g, const char *filename)
{
	FILE *file = fopen(filename, "w");

	for_all_vertices(g, [&file] (const RRNode &v) -> void {
			auto &prop = v.properties;
			if (prop.occ > prop.capacity) {
			for (int x = prop.xlow; x <= prop.xhigh; ++x) {
			for (int y = prop.ylow; y <= prop.yhigh; ++y) {
				for (int i = 0; i < prop.occ-prop.capacity; ++i) {
					fprintf(file, "%d %d\n", x, y);
				}
			}
			}
			}
			});

	fclose(file);
}

void print_parameters(t_router_opts *opts)
{

}	

void dump_route(const vector<trace_t> &traces, const char *filename)
{
	FILE *routes = fopen(filename, "w");
	assert(routes);

	int inet = 0;
	char buffer[256];
	for (const auto &trace : traces) {
		fprintf(routes, "Net %d\n", inet);
		for (const auto &item : trace.segments) {
			int sink_rr_node = item.first;
			const Segment &segment = item.second;

			sprintf_rr_node(sink_rr_node, buffer);
			fprintf(routes, "Sink: %s\n", buffer);

			for (const auto &node : segment) {
				sprintf_rr_node(node, buffer);
				fprintf(routes, "\t%s\n", buffer);
			}
		}
		++inet;
	}

	fclose(routes);
}

bool partitioning_route(t_router_opts *opts)
{
	tbb::task_scheduler_init init(opts->num_threads);
	
	init_logging();

	/*test_scheduler();*/

	/*test_rtree();*/

	/*test_graph();*/

	RRGraph g;
	init_graph(g);

	vector<net_t> nets;
	vector<net_t> global_nets;
	init_nets(nets, global_nets, opts->bb_factor);

	sort_sinks(nets);

	vector<net_t *> nets_ptr(nets.size());
	for (int i = 0; i < nets.size(); ++i) {
		nets_ptr[i] = &nets[i];
	}
	extern s_bb *route_bb;
	sort(nets_ptr.begin(), nets_ptr.end(), [] (const net_t *a, const net_t *b) -> bool {
			return get_bounding_box_area(route_bb[a->vpr_id]) > get_bounding_box_area(route_bb[b->vpr_id]);
			});
	int rank = 0;
	for (auto &net : nets_ptr) {
		net->bb_area_rank = rank++;
	}

	extern int nx;
	extern int ny;
	zlog_info(delta_log, "Num nets: %d Num global nets: %d nx: %d ny: %d\n", nets.size(), global_nets.size(), nx, ny);
	
	t_net_timing *net_timing = new t_net_timing[nets.size()+global_nets.size()];
	init_net_timing(nets, global_nets, net_timing);

	vector<route_tree_t> route_trees(nets.size());
	for (int i = 0; i < nets.size(); ++i) {
		route_tree_init(route_trees[i]);
	}

	vector<trace_t> prev_traces(nets.size());
	vector<trace_t> traces(nets.size());
	for (int i = 0; i < nets.size(); ++i) {
		trace_init(traces[i]);
		trace_init(prev_traces[i]);
	}

	vector<trace_t> *prev_traces_ptr = &prev_traces;
	vector<trace_t> *current_traces_ptr = &traces;

	route_state_t *state = new route_state_t[num_vertices(g)];

	for (int i = 0; i < num_vertices(g); ++i) {
		state[i].rr_node = -1;
		state[i].known_cost = std::numeric_limits<float>::max();
		state[i].cost = std::numeric_limits<float>::max();
		state[i].prev_edge = nullptr;
		state[i].upstream_R = -1;
		state[i].delay = std::numeric_limits<float>::max();
	}


	/*route_tree_t **route_trees = new route_tree_t *[nets.size()+global_nets.size()];*/
	/*for (int i = 0; i < nets.size(); ++i) {*/
		/*route_trees[nets[i].id] = new route_tree_t[nets[i].sinks.size()+1];*/
		/*for (int j = 0; j <= nets[i].sinks.size(); ++j) {*/
			/*route_tree_init(route_trees[nets[i].id][j]);*/
		/*}*/
	/*}*/
	/*for (int i = 0; i < global_nets.size(); ++i) {*/
		/*route_trees[nets[i].id] = new route_tree_t[nets[i].sinks.size()+1];*/
		/*for (int j = 0; j <= nets[i].sinks.size(); ++j) {*/
			/*route_tree_init(route_trees[nets[i].id][j]);*/
		/*}*/
	/*}*/

	route_parameters_t params;
	params.criticality_exp = opts->criticality_exp;
	params.astar_fac = opts->astar_fac;
	params.max_criticality = opts->max_criticality;
	params.bend_cost = opts->bend_cost;
	params.pres_fac = opts->first_iter_pres_fac; /* Typically 0 -> ignore cong. */

	for_all_vertices(g, [] (RRNode &v) -> void {
			v.properties.acc_cost = 1;
			v.properties.pres_cost = 1;
			v.properties.occ = 0;
			});

	char buffer[256];

	bool routed = false;
	for (int iter = 0; iter < opts->max_router_iterations && !routed; ++iter) {
		zlog_info(delta_log, "Routing iteration: %d\n", iter);

		sprintf(buffer, "concurrency_dump_%d.txt", iter);
		FILE *concurrency_dump = fopen(buffer, "w");
		assert(concurrency_dump);

		cpu_timer iter_timer;
		iter_timer.start();

		for (int i = 0; i < nets.size(); ++i) {
			route_tree_clear(route_trees[i]);
		}

		for (auto &net : nets) {
			update_sink_criticalities(net, net_timing[net.vpr_id], params);
		}

		cpu_timer schedule_timer;
		schedule_timer.start();

		vector<net_t *> nets_to_route;
		int current_local_id = 0;
		for (auto &net : nets) {
			if (net.has_sink) {
				nets_to_route.push_back(&net);
				net.current_local_id = current_local_id++;
			}
		}
		vector<vector<const net_t *>> scheduled_nets_at_time;
		if (iter > 0) {
			for (auto &net : nets_to_route) {
				adjust_bounding_box(*net);
			}
		}
		switch (opts->scheduler) {
			case SchedulerType::IND:
				schedule_nets_ind(nets_to_route, scheduled_nets_at_time);
				break;
			case SchedulerType::FAST:
				schedule_nets_fast(nets_to_route, scheduled_nets_at_time);
				break;
			default:
				zlog_error(delta_log, "Error: Unknown scheduler\n");
				assert(false);
				break;
		}

		schedule_timer.stop();
		nanosecond_type schedule_time = schedule_timer.elapsed().wall;

		zlog_info(delta_log, "Average concurrency: %g\n", (float)nets_to_route.size()/scheduled_nets_at_time.size());

		/*zlog_level(delta_log, ROUTER_V3, "-- Start concurrency stats --\n");*/
		/*int temp = 0;*/
		/*for (const auto &nets : scheduled_nets_at_time) {*/
			/*zlog_level(delta_log, ROUTER_V3, "Concurrency at time %d: %d\n", temp, nets.size());*/
			/*++temp;*/
		/*}*/
		/*zlog_level(delta_log, ROUTER_V3, "Average concurrency: %g\n", (float)nets_to_route.size()/scheduled_nets_at_time.size());*/
		/*zlog_level(delta_log, ROUTER_V3, "-- End concurrency stats --\n");*/

		int isink = 0;

		perf_t perf;
		perf.num_heap_pushes = 0;

		nanosecond_type route_time = 0;

		while (scheduled_nets_at_time.size() > 0) {
			zlog_info(delta_log, "isink: %d\n", isink);
			tbb::atomic<int> num_nets_routed = 0;
			cpu_timer route_timer;
			route_timer.start();

			for (const auto &nets_at_time : scheduled_nets_at_time) {
				/* this for loop can be parallelized with parallel_for */
				if (false) {
					tbb::parallel_for(tbb::blocked_range<size_t>(0,nets_at_time.size()), [&] (const tbb::blocked_range<size_t> &r) -> void {
						for (auto inet = r.begin(); inet != r.end(); ++inet) {
							const net_t *net = nets_at_time[inet];
							net_t virtual_net = *net;
							virtual_net.sinks.clear();
							virtual_net.sinks.push_back(net->sinks[net->current_sink_index]);
							/*if (time == 0) {*/
							/*trace_rip_up_net(traces[net->local_id], g, params.pres_fac);*/
							/*}*/
							trace_rip_up_segment((*prev_traces_ptr)[net->local_id], g, net->sinks[net->current_sink_index].rr_node, params.pres_fac);
							/*printf("Net: %d\n", virtual_net.vpr_id);*/
							route_net(g, virtual_net, params, state, route_trees[net->local_id], (*current_traces_ptr)[net->local_id], net_timing[net->vpr_id], &perf);
							/*route_tree_merge(route_trees[net->id][0], route_trees[net->id][virtual_net.sinks[0].id]);*/
							++num_nets_routed;
						}
							});
				} else {
					for (const auto &net : nets_at_time) {
						net_t virtual_net = *net;
						virtual_net.sinks.clear();
						virtual_net.sinks.push_back(net->sinks[net->current_sink_index]);
						/*if (time == 0) {*/
						/*trace_rip_up_net(traces[net->local_id], g, params.pres_fac);*/
						/*}*/
						trace_rip_up_segment((*prev_traces_ptr)[net->local_id], g, net->sinks[net->current_sink_index].rr_node, params.pres_fac);
						/*printf("Net: %d\n", virtual_net.vpr_id);*/
						route_net(g, virtual_net, params, state, route_trees[net->local_id], (*current_traces_ptr)[net->local_id], net_timing[net->vpr_id], &perf);
						/*route_tree_merge(route_trees[net->id][0], route_trees[net->id][virtual_net.sinks[0].id]);*/
						++num_nets_routed;
					}
				}

				fprintf(concurrency_dump, "%d\n", nets_at_time.size());
			}

			if (isink == 0) {
				assert(num_nets_routed == nets.size());
			}
			assert(num_nets_routed == nets_to_route.size());

			for (auto &net : nets_to_route) {
				set_to_next_sink(*net, route_trees[net->local_id], g, params.astar_fac, perf);
			}

			route_timer.stop();
			cpu_times elapsed = route_timer.elapsed();
			route_time += elapsed.wall;
			/*zlog_info(delta_log, "Scheduling took %s\n", format(elapsed).c_str());*/
			zlog_info(delta_log, "Routing took %g\n", elapsed.wall / 1e9);

			schedule_timer.start();
		
			nets_to_route.clear();
			int current_local_id = 0;
			for (auto &net : nets) {
				if (net.has_sink) {
					nets_to_route.push_back(&net);
					net.current_local_id = current_local_id++;
				}
			}
			scheduled_nets_at_time.clear();
			if (iter > 0) {
				for (auto &net : nets_to_route) {
					adjust_bounding_box(*net);
				}
			}
			switch (opts->scheduler) {
				case SchedulerType::IND:
					schedule_nets_ind(nets_to_route, scheduled_nets_at_time);
					break;
				case SchedulerType::FAST:
					schedule_nets_fast(nets_to_route, scheduled_nets_at_time);
					break;
				default:
					zlog_error(delta_log, "Error: Unknown scheduler\n");
					assert(false);
					break;
			}

			schedule_timer.stop();
			schedule_time += schedule_timer.elapsed().wall;
			zlog_info(delta_log, "Average concurrency: %g\n", (float)nets_to_route.size()/scheduled_nets_at_time.size());

			++isink;
		}

		fclose(concurrency_dump);

		zlog_info(delta_log, "Num heap pushes: %d\n", perf.num_heap_pushes);

		for_all_vertices(g, [] (RRNode &v) -> void { v.properties.recalc_occ = 0; });

		for (const auto &net : nets) {
			check_route_tree(route_trees[net.local_id], net, g);
		}

		bool valid = true;
		for_all_vertices(g, [&valid, &buffer] (RRNode &v) -> void {
				sprintf_rr_node(id(v), buffer);
				if (v.properties.recalc_occ != v.properties.occ) {
					zlog_error(delta_log, "Node %s occ mismatch, recalc: %d original: %d\n", buffer, v.properties.recalc_occ, v.properties.occ);
					valid = false;
				}
				});
		assert(valid);

		for (const auto &net : nets) {
			vector<int> overused_rr_node;
			check_net_route(route_trees[net.local_id], get_vertex(route_trees[net.local_id].graph, route_trees[net.local_id].root), g, overused_rr_node);
			if (!overused_rr_node.empty()) {
				zlog_level(delta_log, ROUTER_V1, "Net %d bb_rank %d overused nodes:\n", net.vpr_id, net.bb_area_rank);
				for (const auto &item : overused_rr_node) {
					zlog_level(delta_log, ROUTER_V1, "%d ", item);
				}
				zlog_level(delta_log, ROUTER_V1, "\n");
			}
		}

		char filename[256];
		sprintf(filename, "congestion_%d.txt", iter);
		dump_congestion_map(g, filename);

		for (auto &net : nets) {
			/*zlog_level(delta_log, ROUTER_V1, "Net %d trace stats:\n", net.vpr_id);*/
			adjust_bb_factor((*current_traces_ptr)[net.local_id], net, g, 3);
		}

		if (feasible_routing(g)) {
			zlog_info(delta_log, "Routed in %d iterations\n", iter+1);
			dump_route(*current_traces_ptr, "route.txt");
			routed = true;
		} else {
			int num_overused_nodes = 0;
			for_all_vertices(g, [&num_overused_nodes] (const RRNode &v) -> void {
					if (v.properties.occ > v.properties.capacity) {
					++num_overused_nodes;
					}
					});
			zlog_info(delta_log, "Num overused nodes: %d/%d (%.2f)\n", num_overused_nodes, num_vertices(g), num_overused_nodes*100.0/num_vertices(g));

			if (iter == 0) {
				params.pres_fac = opts->initial_pres_fac;
				update_costs(g, params.pres_fac, 0);
			} else {
				params.pres_fac *= opts->pres_fac_mult;

				/* Avoid overflow for high iteration counts, even if acc_cost is big */
				params.pres_fac = std::min(params.pres_fac, static_cast<float>(HUGE_POSITIVE_FLOAT / 1e5));

				update_costs(g, params.pres_fac, opts->acc_fac);
			}

			for (auto &net : nets) { 
				reset_current_source_sink(net);
			}

			/* swap the pointers */
			std::swap(current_traces_ptr, prev_traces_ptr);
			for (const auto &trace : *current_traces_ptr) {
				assert(trace_empty(trace));
			}
		}

		analyze_timing(net_timing);

		iter_timer.stop();
		cpu_times elapsed = iter_timer.elapsed();

		zlog_info(delta_log, "Iteration time: %g s. Route time: %g s. Schedule time: %g s.\n", elapsed.wall / 1e9, route_time / 1e9, schedule_time / 1e9);
	}
	return routed;
}

bool hybrid_route(t_router_opts *opts)
{
	const int max_iter = 50;

	test_interval();

	init_logging();

	RRGraph g;
	vector<net_t> nets;
	vector<net_t> global_nets;
	init_graph(g);
	init_nets(nets, global_nets, opts->bb_factor);
	
	const int num_parts = 4;

	graph_t<int, int> dep_g;
	/*load_overlapping_nets(nets);*/

	/*const char *metis_filename = "metis.graph";*/
	/*write_metis_file(dep_g, metis_filename);*/

	/*run_metis(metis_filename, num_parts);*/

	/*char part_filename[256];*/
	/*sprintf(part_filename, "%s.part.%d", metis_filename, num_parts);*/
	/*vector<partition_t> partitions(num_parts);*/
	/*init_nets_partition(nets, part_filename, num_parts, partitions);*/

	/*schedule_nets_greedy(nets);*/

	std::sort(nets.begin(), nets.end(), [] (const net_t &a, const net_t &b) -> bool {
		return a.sinks.size() > b.sinks.size(); }
		);

	t_net_timing *net_timing = new t_net_timing[nets.size()];
	init_net_timing(nets, global_nets, net_timing);

	exit(-1);

	route_tree_t *route_trees = new route_tree_t[nets.size()];

	route_parameters_t params;
	params.criticality_exp = opts->criticality_exp;
	params.astar_fac = opts->astar_fac;
	params.max_criticality = opts->max_criticality;
	params.bend_cost = opts->bend_cost;

	for_all_vertices(g, [] (RRNode &v) -> void {
			v.properties.acc_cost = 1;
			v.properties.pres_cost = 1;
			v.properties.occ = 0;
			});

	params.pres_fac = opts->first_iter_pres_fac; /* Typically 0 -> ignore cong. */

	bool routed = false;
	for (int iter = 0; iter < max_iter && !routed; ++iter) {
		zlog_info(delta_log, "Routing iteration: %d\n", iter);
		fine_grained_route_nets(g, nets, params, route_trees, net_timing);
		/*if (use_net_level_parallelism()) {*/
			/*[>delta_stepping_route_nets();<]*/
		/*} else {*/
			/*fine_grained_route_nets();*/
		/*}*/

		if (feasible_routing(g)) {
			zlog_info(delta_log, "Routed in %d iterations\n", iter);
			routed = true;
		} else {
			int num_overused_nodes = 0;
			for_all_vertices(g, [&num_overused_nodes] (const RRNode &v) -> void {
					if (v.properties.occ > v.properties.capacity) {
					++num_overused_nodes;
					}
					});
			zlog_info(delta_log, "Num overused nodes: %d/%d (%.2f)\n", num_overused_nodes, num_vertices(g), num_overused_nodes*100.0/num_vertices(g));

			if (iter == 0) {
				params.pres_fac = opts->initial_pres_fac;
				update_costs(g, params.pres_fac, 0);
			} else {
				params.pres_fac *= opts->pres_fac_mult;

				/* Avoid overflow for high iteration counts, even if acc_cost is big */
				params.pres_fac = std::min(params.pres_fac, static_cast<float>(HUGE_POSITIVE_FLOAT / 1e5));

				update_costs(g, params.pres_fac, opts->acc_fac);
			}
		}

		analyze_timing(net_timing);
	}
	return routed;
}
