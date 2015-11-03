#include <assert.h>
#include <zlog.h>
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
#include "bounding_box.h"
#include "graph.h"

void sprintf_rr_node(int inode, char *buffer);

/*typedef struct heap_item_t {*/
	/*int node;*/
	/*float cost;*/
	/*float known_cost;*/
	/*float upstream_R;*/
/*} heap_item_t;*/

typedef struct route_parameters_t {
	float pres_fac;
	float astar_fac;
	float bend_cost;
	float max_criticality;
	float criticality_exp;
} route_parameters_t;

typedef struct bounding_box_t {
	int xmin;
	int xmax;
	int ymin;
	int ymax;
} bounding_box_t;

typedef struct net_t {
	int id;
	bool global;
	int source;
	vector<int> sinks;
	bounding_box_t box;
} net_t;

typedef struct rr_node_property_t {
	t_rr_type type;
	int xlow;
	int ylow;
	int xhigh;
	int yhigh;
	float R;
	float C;
	int cost_index;
	int capacity;
	float base_cost;
	int occ;
	float pres_cost;
	float acc_cost;
} rr_node_property_t;

typedef struct rr_edge_property_t {
	bool buffered;
	float switch_delay;
	float R;
} rr_edge_property_t;

typedef struct route_state_t {
	int rr_node;
	const edge_t<rr_edge_property_t> *prev_edge;
	float upstream_R;
	float delay;
	float known_cost;
	float cost;
} route_state_t;

typedef struct rt_node_property_t {
	bool reexpand;
	int rr_node;
	const edge_t<rr_edge_property_t> *prev_edge;
	float upstream_R;	
	/*float upstream_R_from_route_state;*/
	float downstream_C;
	float delay;
} rt_node_property_t;

typedef struct rt_edge_property_t {
	const edge_t<rr_edge_property_t> *rr_edge;
} rt_edge_property_t;

typedef graph_t<rr_node_property_t, rr_edge_property_t> RRGraph;
typedef vertex_t<rr_node_property_t, rr_edge_property_t> RRNode;
typedef edge_t<rr_edge_property_t> RREdge;

typedef graph_t<rt_node_property_t, rt_edge_property_t> RouteTree;
typedef vertex_t<rt_node_property_t, rt_edge_property_t> RouteTreeNode;
typedef edge_t<rt_edge_property_t> RouteTreeEdge;

typedef struct route_tree_t {
	route_tree_t() : root(-1) {}
	RouteTree graph;
	map<int, int> rr_node_to_rt_node;
	int root;
} route_tree_t;

zlog_category_t *delta_log;

static void init_logging()
{
	delta_log = zlog_get_category("delta");
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

float get_delay(const edge_t<rr_edge_property_t> &e, const RRNode &v, float unbuffered_upstream_R)
{
	float upstream_R = e.properties.R;
	if (!e.properties.buffered) {
		upstream_R += unbuffered_upstream_R;
	}

	float delay = e.properties.switch_delay;
	delay += v.properties.C * (upstream_R + 0.5 * v.properties.R);

	zlog_debug(delta_log, " [edge_delay: %g edge_R: %g node_R: %g node_C: %g] ", e.properties.switch_delay, e.properties.R, v.properties.R, v.properties.C);

	return delay;
}

float get_congestion_cost(const RRNode &v)
{
	return v.properties.base_cost * v.properties.acc_cost * v.properties.pres_cost;
}

float get_known_cost(const RRGraph &g, const edge_t<rr_edge_property_t> &e, float criticality_fac, float unbuffered_upstream_R)
{
	const auto &target = get_target(g, e);

	float delay = get_delay(e, target, unbuffered_upstream_R);
	float congestion = get_congestion_cost(target);

	zlog_debug(delta_log, " [delay: %g congestion %g crit_fac: %g] ", delay, congestion, criticality_fac);

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
void expand_neighbors(const RRGraph &g, const route_state_t &current, const RRNode &target, float criticality_fac, float astar_fac, std::priority_queue<route_state_t> &heap, ShouldExpandFunc should_expand)
{
	for_all_out_edges(g, get_vertex(g, current.rr_node), [&heap, &g, &current, &target, &criticality_fac, &astar_fac, &should_expand] (const edge_t<rr_edge_property_t> &e) -> void {
			auto &neighbor = get_target(g, e);

			char buffer[256];
			sprintf_rr_node(id(neighbor), buffer);
			zlog_debug(delta_log, "\tNeighbor: %s ", buffer);
			
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
			
			zlog_debug(delta_log, "added with prev_edge: %X upstream_R: %g delay: %g known_cost: %g expected_cost: %g cost: %g\n", item.prev_edge, item.upstream_R, item.delay, item.known_cost, expected_cost, item.cost);
	});
}

RouteTreeNode &route_tree_add_rr_node(route_tree_t &rt, int rr_node)
{
	auto iter = rt.rr_node_to_rt_node.find(rr_node);

	RouteTreeNode *v = nullptr;
	if (iter == rt.rr_node_to_rt_node.end()) {
		add_vertex(rt.graph);
		int rt_node = num_vertices(rt.graph)-1;
		rt.rr_node_to_rt_node[rr_node] = rt_node;
		v = &get_vertex(rt.graph, rt_node);
	} else {
		v = &get_vertex(rt.graph, iter->second);
	}
	return *v;
}

/*RouteTreeEdge &route_tree_add_edge(route_tree_t &rt, int rt_node_a, int rt_node_b)*/
/*{*/
	/*[>int rt_node_a = route_tree_add_rr_node(rt, rr_node_a);<]*/
	/*[>int rt_node_b = route_tree_add_rr_node(rt, rr_node_b);<]*/
	/*[>assert(rt_node_a < num_vertices(rt.graph));<]*/
	/*[>assert(rt_node_b < num_vertices(rt.graph));<]*/
	/*return add_edge(rt.graph, rt_node_a, rt_node_b);*/
/*}*/

RouteTreeEdge &route_tree_add_edge_between_rt_node(route_tree_t &rt, const RouteTreeNode &a, const RouteTreeNode &b)
{
	/*int rt_node_a = route_tree_add_rr_node(rt, rr_node_a);*/
	/*int rt_node_b = route_tree_add_rr_node(rt, rr_node_b);*/
	/*assert(rt_node_a < num_vertices(rt.graph));*/
	/*assert(rt_node_b < num_vertices(rt.graph));*/
	return add_edge(rt.graph, id(a), id(b));
}

void route_tree_add_to_heap(const route_tree_t &rt, const RouteTreeNode *rt_node, const RRGraph &g, const RRNode &target, float criticality_fac, float astar_fac, std::priority_queue<route_state_t> &heap)
{
	if (!rt_node) {
		rt_node = &get_vertex(rt.graph, rt.root);
	}

	if (rt_node->properties.reexpand) {
		route_state_t item;

		item.rr_node = rt_node->properties.rr_node;
		item.known_cost = criticality_fac * rt_node->properties.delay;
		float expected_cost = get_timing_driven_expected_cost(get_vertex(g, rt_node->properties.rr_node), target, criticality_fac, rt_node->properties.upstream_R);
		item.cost = item.known_cost + astar_fac * expected_cost;
		item.prev_edge = nullptr;
		item.upstream_R = rt_node->properties.upstream_R;
		item.delay = rt_node->properties.delay;

		char buffer[256];

		sprintf_rr_node(item.rr_node, buffer);
		zlog_debug(delta_log, "Adding route tree node %s delay: %g known_cost: %g expected_cost: %g astar_fac: %g cost: %g upstream_R: %g to heap\n", buffer, item.delay, item.known_cost, expected_cost, astar_fac, item.cost, item.upstream_R);

		heap.push(item);
	}

	for_all_out_edges(rt.graph, *rt_node, [&rt, &heap, &criticality_fac, &astar_fac, &g, &target] (const RouteTreeEdge &e) -> void {
		const auto &neighbor = get_target(rt.graph, e);

		route_tree_add_to_heap(rt, &neighbor, g, target, criticality_fac, astar_fac, heap);
	}
	);
}

RouteTreeNode &route_tree_get_rt_node(route_tree_t &rt, int rr_node)
{
	auto iter = rt.rr_node_to_rt_node.find(rr_node);
	assert(iter != rt.rr_node_to_rt_node.end());
	return get_vertex(rt.graph, iter->second);
}

const RouteTreeNode &route_tree_get_rt_node(const route_tree_t &rt, int rr_node)
{
	auto iter = rt.rr_node_to_rt_node.find(rr_node);
	assert(iter != rt.rr_node_to_rt_node.end());
	return get_vertex(rt.graph, iter->second);
}

const RouteTreeNode &route_tree_update(route_tree_t &rt, const RRGraph &g, const route_state_t *state, int sink)
{
	int current = sink;
	int new_branch = -1;

	while (state[current].prev_edge) {
		int prev = id(get_source(g, *state[current].prev_edge));

		auto &prev_rt_node = route_tree_add_rr_node(rt, prev);
		prev_rt_node.properties.reexpand = get_vertex(g, prev).properties.type != IPIN;
		prev_rt_node.properties.rr_node = prev;
		prev_rt_node.properties.prev_edge = nullptr;
		prev_rt_node.properties.upstream_R = state[prev].upstream_R;
		prev_rt_node.properties.delay = state[prev].delay;

		auto &current_rt_node = route_tree_add_rr_node(rt, current);
		current_rt_node.properties.reexpand = get_vertex(g, current).properties.type != IPIN;
		current_rt_node.properties.rr_node = current;
		current_rt_node.properties.prev_edge = state[current].prev_edge;
		current_rt_node.properties.upstream_R = state[current].upstream_R;
		current_rt_node.properties.delay = state[current].delay;

		char p[256];
		char c[256];
		sprintf_rr_node(prev, p);
		sprintf_rr_node(current, c);
		zlog_debug(delta_log, "Traceback: %s -> %s\n", p, c);

		auto &e = route_tree_add_edge_between_rt_node(rt, route_tree_get_rt_node(rt, prev), route_tree_get_rt_node(rt, current));
		e.properties.rr_edge = state[current].prev_edge;

		new_branch = current;
		current = prev;
	}

	return route_tree_get_rt_node(rt, new_branch);
}

void route_tree_clear(route_tree_t &rt)
{
	rt.root = -1;
	clear_edges(rt.graph);
}

void route_tree_init(route_tree_t &rt, const RRNode &source)
{
	auto &rt_node = route_tree_add_rr_node(rt, id(source)); 

	rt_node.properties.upstream_R = source.properties.R;
	rt_node.properties.delay = 0.5 * source.properties.R * source.properties.C;
	rt_node.properties.rr_node = id(source);
	rt_node.properties.reexpand = true;
	rt_node.properties.prev_edge = nullptr;

	rt.root = id(rt_node);
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

void update_one_cost(RRGraph &g, route_tree_t &rt, const RouteTreeNode &node, int delta, float pres_fac)
{
	RRNode &rr_node = get_vertex(g, node.properties.rr_node);

	rr_node.properties.occ += delta;

	assert(rr_node.properties.occ >= 0);

	if (rr_node.properties.occ < rr_node.properties.capacity) {
		rr_node.properties.pres_cost = 1;
	} else {
		rr_node.properties.pres_cost = 1 + (rr_node.properties.occ + 1 - rr_node.properties.capacity) * pres_fac;
	}

	for_all_out_edges(rt.graph, node, [&g, &rt, &delta, &pres_fac] (const RouteTreeEdge &e) -> void {
			update_one_cost(g, rt, get_target(rt.graph, e), delta, pres_fac);
			});
}

void update_R()
{
}

void route_tree_update_timing(route_tree_t &rt, const RouteTreeNode &node)
{
	for_all_out_edges(rt.graph, node, [] (const RouteTreeEdge &e) -> void {
			e.properties.rr_edge;
	});
}

typedef struct sink_t {
	int id;
	float criticality_fac;
	int rr_node;
} sink_t;

bool operator<(const sink_t &a, const sink_t &b) {
	return a.criticality_fac > b.criticality_fac;
}

void sort_sinks(const net_t &net, const t_net_timing *net_timing, const route_parameters_t &params, vector<sink_t> &sorted_sinks)
{
	for (int ipin = 1; ipin <= net.sinks.size(); ipin++) { 
		float pin_criticality;
		if (!net_timing) {
			/* Use criticality of 1. This makes all nets critical.  Note: There is a big difference between setting pin criticality to 0
			compared to 1.  If pin criticality is set to 0, then the current path delay is completely ignored during routing.  By setting
			pin criticality to 1, the current path delay to the pin will always be considered and optimized for */
			pin_criticality = 1.0;
		} else { 
#ifdef PATH_COUNTING
			/* Pin criticality is based on a weighted sum of timing and path criticalities. */	
			pin_criticality =		 ROUTE_PATH_WEIGHT	* net_timing->path_criticality[ipin]
								  + (1 - ROUTE_PATH_WEIGHT) * net_timing->timing_criticality[ipin]; 
#else
			/* Pin criticality is based on only timing criticality. */
			pin_criticality = net_timing->timing_criticality[ipin];
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
		}
		sorted_sinks.push_back({ ipin, pin_criticality, net.sinks[ipin-1] });
	}

	std::sort(sorted_sinks.begin(), sorted_sinks.end());
}

bool route_tree_initiated(const route_tree_t &rt)
{
	return rt.root != -1;
}

void check_route_tree_internal(const route_tree_t &rt, const RouteTreeNode &node, const RRGraph &g, vector<int> &visited_sinks)
{
	const auto &rr_node = get_vertex(g, node.properties.rr_node);
	char buffer[256];
	if (rr_node.properties.type == SINK) {
		visited_sinks.push_back(id(rr_node));
		sprintf_rr_node(id(rr_node), buffer);
		zlog_debug(delta_log, "route_tree_check: %s\n", buffer);
	}

	for_all_out_edges(rt.graph, node, [&rt, &g, &visited_sinks] (const RouteTreeEdge &e) -> void {
			const auto &neighbor = get_target(rt.graph, e);
			check_route_tree_internal(rt, neighbor, g, visited_sinks);
			});
}

void check_route_tree(const route_tree_t &rt, const net_t &net, const RRGraph &g)
{
	const auto &rt_root = route_tree_get_rt_node(rt, net.source);
	vector<int> sinks = net.sinks;
	vector<int> visited_sinks;

	check_route_tree_internal(rt, rt_root, g, visited_sinks);

	sort(visited_sinks.begin(), visited_sinks.end());
	sort(sinks.begin(), sinks.end());
	
	assert(visited_sinks == sinks);
}

void route_net(RRGraph &g, const net_t &net, const route_parameters_t &params, route_tree_t &rt, t_net_timing *net_timing)

{
	std::priority_queue<route_state_t> heap;
	route_state_t *state = new route_state_t[num_vertices(g)];

	for (int i = 0; i < num_vertices(g); ++i) {
		state[i].rr_node = -1;
		state[i].known_cost = std::numeric_limits<float>::max();
		state[i].cost = std::numeric_limits<float>::max();
		state[i].prev_edge = nullptr;
		state[i].upstream_R = -1;
		state[i].delay = std::numeric_limits<float>::max();
	}

	vector<int> modified;

	/* rip up previous routing */
	if (route_tree_initiated(rt)) {
		update_one_cost(g, rt, route_tree_get_rt_node(rt, net.source), -1, params.pres_fac);
	} 

	route_tree_clear(rt);
	route_tree_init(rt, get_vertex(g, net.source));
	update_one_cost(g, rt, route_tree_get_rt_node(rt, net.source), 1, params.pres_fac);

	vector<sink_t> sorted_sinks;
	sort_sinks(net, net_timing, params, sorted_sinks);

	zlog_debug(delta_log, "Routing net %d\n", net.id);

	char buffer[256];

	sprintf_rr_node(net.source, buffer);
	zlog_debug(delta_log, "Source: %s\n", buffer);

	int isink = 0;
	for (const auto &sink : sorted_sinks) {
		sprintf_rr_node(sink.rr_node, buffer);
		zlog_debug(delta_log, "%d sink: %s criticality: %g\n", isink++, buffer, sink.criticality_fac);
		route_tree_add_to_heap(rt, nullptr, g, get_vertex(g, sink.rr_node), sink.criticality_fac, params.astar_fac, heap);

		bool found_sink = false;
		while (!heap.empty() && !found_sink) {
			auto item = heap.top();
			heap.pop();

			sprintf_rr_node(item.rr_node, buffer);
			zlog_debug(delta_log, "Current: %s prev: %d old_delay: %g new_delay: %g old_known: %g new_known: %g old_cost: %g new_cost: %g\n", buffer, item.prev_edge ? id(get_source(g, *item.prev_edge)) : -1, state[item.rr_node].delay, item.delay, state[item.rr_node].known_cost, item.known_cost, state[item.rr_node].cost, item.cost);

			if (item.rr_node == sink.rr_node) {
				state[item.rr_node] = item;
				modified.push_back(item.rr_node);
				found_sink = true;
			} else if (item.known_cost < state[item.rr_node].known_cost && item.cost < state[item.rr_node].cost) {
				const auto &sink_vertex = get_vertex(g, sink.rr_node);

				state[item.rr_node] = item;
				modified.push_back(item.rr_node);

				expand_neighbors(g, item, sink_vertex, sink.criticality_fac, params.astar_fac, heap, [&net, &sink_vertex] (const RRNode &v) -> bool {
					const auto &prop = v.properties;

					if (prop.xhigh < net.box.xmin
							|| prop.xlow > net.box.xmax
							|| prop.yhigh < net.box.ymin
							|| prop.ylow > net.box.ymax) {
					zlog_debug(delta_log, "outside of bounding box\n");
					return false;
					}

					if (prop.type == IPIN
							&& (prop.xhigh != sink_vertex.properties.xhigh 
								|| prop.yhigh != sink_vertex.properties.yhigh)) {
					zlog_debug(delta_log, "not target IPIN\n");
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
				});
			} else {
				zlog_debug(delta_log, "\tNot expanding neighbor because known_cost %g > %g or cost %g > %g\n", item.known_cost, state[item.rr_node].known_cost, item.cost, state[item.rr_node].cost);
			}
		}

		assert(found_sink);

		const RouteTreeNode &new_branch = route_tree_update(rt, g, state, sink.rr_node);
		update_one_cost(g, rt, new_branch, 1, params.pres_fac);

		net_timing->delay[sink.id] = route_tree_get_rt_node(rt, sink.rr_node).properties.delay;

		for (const auto &m : modified)  {
			state[m].known_cost = std::numeric_limits<float>::max();
			state[m].cost = std::numeric_limits<float>::max();
		}

		heap = std::priority_queue<route_state_t>();
	}

	check_route_tree(rt, net, g);
}

void init_net_timing(const vector<net_t> &nets, t_net_timing *net_timing)
{
	for (const auto &net : nets) {
	   	net_timing[net.id].delay = new float[net.sinks.size() + 1];
	   	net_timing[net.id].timing_criticality = new float[net.sinks.size() + 1];
	   	net_timing[net.id].slack = new float[net.sinks.size() + 1];
	}

	float init_timing_criticality_val = 1;

	for (const auto &net : nets) {
		if (!net.global) {
			for (int ipin = 1; ipin <= net.sinks.size(); ipin++)
				net_timing[net.id].timing_criticality[ipin] = init_timing_criticality_val;
#ifdef PATH_COUNTING
			net_timing[net.id].path_criticality[ipin] = init_timing_criticality_val;
#endif		
		} else { 
			/* Set delay of global signals to zero. Non-global net 
			 * 			delays are set by update_net_delays_from_route_tree() 
			 * 						inside timing_driven_route_net(), which is only called
			 * 									for non-global nets. */
			for (int ipin = 1; ipin <= net.sinks.size(); ipin++) {
				net_timing[net.id].delay[ipin] = 0.;
			}
		}
	}
}

void fine_grained_route_nets(RRGraph &g, const vector<net_t> &nets, const route_parameters_t &params, route_tree_t *route_trees, t_net_timing *net_timing)
{
	int inet = 0;
	for (const auto &net : nets) {
		route_net(g, net, params, route_trees[net.id], &net_timing[net.id]);
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
		v.properties.xlow = rr_node[i].xlow;
		v.properties.ylow = rr_node[i].ylow;
		v.properties.xhigh = rr_node[i].xhigh;
		v.properties.yhigh = rr_node[i].yhigh;
		v.properties.R = rr_node[i].R;
		v.properties.C = rr_node[i].C;
		v.properties.cost_index = rr_node[i].cost_index;
		v.properties.capacity = rr_node[i].capacity;

		int ci = rr_node[i].cost_index;
		v.properties.base_cost = rr_indexed_data[ci].saved_base_cost;

		for (int j = 0; j < rr_node[i].num_edges; ++j) {
			auto &e = add_edge(g, i, rr_node[i].edges[j]);

			int si = rr_node[i].switches[j];
			
			e.properties.buffered = switch_inf[si].buffered; 
			e.properties.switch_delay = switch_inf[si].Tdel; 
			e.properties.R = switch_inf[si].R; 
		}
	}
}

void init_nets(vector<net_t> &nets)
{
	extern struct s_net *clb_net;
	extern int num_nets;
	extern int **net_rr_terminals; /* [0..num_nets-1][0..num_pins-1] */
	extern struct s_bb *route_bb;

	nets.resize(num_nets);
	for (int i = 0; i < num_nets; ++i) {
		nets[i].id = i;
		nets[i].global = clb_net[i].is_global ? true : false;
		nets[i].source = net_rr_terminals[i][0];
		for (int j = 1; j <= clb_net[i].num_sinks; j++) {
			nets[i].sinks.push_back(net_rr_terminals[i][j]);
		}
		nets[i].box.xmin = route_bb[i].xmin;
		nets[i].box.ymin = route_bb[i].ymin;
		nets[i].box.xmax = route_bb[i].xmax;
		nets[i].box.ymax = route_bb[i].ymax;
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

bool hybrid_route(t_router_opts *opts)
{
	const int max_iter = 50;

	init_logging();

	RRGraph g;
	vector<net_t> nets;
	init_graph(g);
	init_nets(nets);

	std::sort(nets.begin(), nets.end(), [] (const net_t &a, const net_t &b) -> bool {
		return a.sinks.size() > b.sinks.size(); }
		);

	t_net_timing *net_timing = new t_net_timing[nets.size()];
	init_net_timing(nets, net_timing);

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
