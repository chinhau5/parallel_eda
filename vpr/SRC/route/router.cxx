#include "pch.h"

#include "vpr_types.h"
#include "path_delay.h"

#include "route.h"
#include "route_tree.h"
#include "congestion.h"
#include "trace.h"

bool operator<(const sink_t &a, const sink_t &b)
{
	return a.criticality_fac > b.criticality_fac;
}

bool operator<(const route_state_t &a, const route_state_t &b)
{
	return a.cost > b.cost;
}

bool feasible_routing(const RRGraph &g, const congestion_t *congestion)
{
	bool feasible = true;

	for (int i = 0; i < num_vertices(g) && feasible; ++i) {
		if (congestion[i].occ > get_vertex(g, i).properties.capacity) {
			feasible = false;
		}
	}

	return feasible;
}

float analyze_timing(t_net_timing *net_timing) 
{
	load_timing_graph_net_delays_new(net_timing); 
#ifdef HACK_LUT_PIN_SWAPPING
	do_timing_analysis_new(net_timing, FALSE, TRUE, FALSE);
#else
	do_timing_analysis_new(net_timing, FALSE, FALSE, FALSE);
#endif

	return get_critical_path_delay();
	/*zlog_info(delta_log, "Critical path: %g ns\n", critical_path_delay);*/
	/*printf("Critical path: %g ns\n", critical_path_delay);*/
}

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

void get_overused_nodes(const route_tree_t &rt, const RouteTreeNode &node, const RRGraph &g, const congestion_t *congestion, vector<int> &overused_rr_node)
{
	assert(node.properties.valid);

	const auto &rr_node = get_vertex(g, node.properties.rr_node);

	int rr_node_id = id(rr_node);
	if (congestion[rr_node_id].occ > rr_node.properties.capacity) {
		overused_rr_node.push_back(rr_node_id);
	}
	
	for_all_out_edges(rt.graph, node, [&rt, &g, &congestion, &overused_rr_node] (const RouteTreeEdge &e) -> void {
			const auto &neighbor = get_vertex(rt.graph, get_target(rt.graph, id(e)));
			get_overused_nodes(rt, neighbor, g, congestion, overused_rr_node);
			});
}

void recalculate_occ_internal(const route_tree_t &rt, const RouteTreeNode &node, const RRGraph &g, congestion_t *congestion)
{
	assert(node.properties.valid);

	auto &rr_node = get_vertex(g, node.properties.rr_node);
	int rr_node_id = id(rr_node);
	if (rr_node.properties.type == SOURCE) {
		congestion[rr_node_id].recalc_occ += num_out_edges(rt.graph, id(node));
	} else {
		++congestion[rr_node_id].recalc_occ;
	}

	for (const auto &branch : route_tree_get_branches(rt, id(node))) {
		/*for_all_out_edges(rt.graph, node, [&rt, &g, &visited_sinks, &visited_nodes] (const RouteTreeEdge &e) -> void {*/
		const auto &child = get_vertex(rt.graph, get_target(rt.graph, branch));
		recalculate_occ_internal(rt, child, g, congestion);
	}
}

void recalculate_occ(const route_tree_t &rt, const RRGraph &g, congestion_t *congestion)
{
	const auto &rt_root = get_vertex(rt.graph, rt.root_rt_node_id);

	recalculate_occ_internal(rt, rt_root, g, congestion);
}

void check_route_tree_internal(const route_tree_t &rt, const RouteTreeNode &node, RRGraph &g, vector<int> &visited_sinks, vector<int> &visited_nodes)
{
	assert(node.properties.valid);

	auto &rr_node = get_vertex(g, node.properties.rr_node);
	if (rr_node.properties.type == SINK) {
		visited_sinks.push_back(id(rr_node));
		char buffer[256];
		sprintf_rr_node(id(rr_node), buffer);
		/*zlog_level(delta_log, ROUTER_V2, "route_tree_check: %s\n", buffer);*/
	}

	visited_nodes.push_back(node.properties.rr_node);

	for (const auto &branch : route_tree_get_branches(rt, id(node))) {
		/*for_all_out_edges(rt.graph, node, [&rt, &g, &visited_sinks, &visited_nodes] (const RouteTreeEdge &e) -> void {*/
		const auto &child = get_vertex(rt.graph, get_target(rt.graph, branch));
		check_route_tree_internal(rt, child, g, visited_sinks, visited_nodes);
	}
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
	vector<int> duplicated_nodes;
	sort(begin(visited_nodes), end(visited_nodes));
	int current = visited_nodes[0];
	bool added_duplicated = false;
	for (int i = 1; i < visited_nodes.size(); ++i) {
		/*zlog_info(delta_log, "visited nodes %d\n", visited_nodes[i]);*/
		if (visited_nodes[i] == current) {
			if (!added_duplicated) {
				duplicated_nodes.push_back(visited_nodes[i]);
				added_duplicated = true;
			}
		} else {
			current = visited_nodes[i];
			added_duplicated = false;
		}
	}
	char buffer[256];
	if (!duplicated_nodes.empty()) {
		zlog_error(delta_log, "Error: net %d visited_nodes has duplicates: \n", net.vpr_id);
		for (const auto &d : duplicated_nodes) {
			sprintf_rr_node(d, buffer);
			zlog_error(delta_log, "%s\n", buffer);
		}
		write_graph(rt.graph, "duplicate.dot",
				[&duplicated_nodes] (const RouteTreeNode &rt_node) -> string {
				char s_rr_node[256];
				char buffer[256];
				if (find(begin(duplicated_nodes), end(duplicated_nodes), rt_node.properties.rr_node) !=
						end(duplicated_nodes)) {
				sprintf_rr_node(rt_node.properties.rr_node, s_rr_node);
				sprintf(buffer, "label=\"%s\" style=filled fillcolor=red", s_rr_node);
				} else {
				sprintf_rr_node(rt_node.properties.rr_node, s_rr_node);
				sprintf(buffer, "label=\"%s\"", s_rr_node);
				}
				return buffer;
				},
				[] (const RouteTreeEdge &rt_edge) -> string {
				return "";
				},
				[] (const RouteTreeNode &rt_node) -> bool {
				return false;
				});
		assert(false);
	}

	sort(visited_sinks.begin(), visited_sinks.end());
	sort(sinks.begin(), sinks.end());
	
	if (visited_sinks != sinks) {
		zlog_error(delta_log, "Error: Visited %lu sinks out of %lu sinks of net %d\n", visited_sinks.size(), sinks.size(), net.vpr_id);
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

float get_delay(const RREdge &e, const RRNode &v, float unbuffered_upstream_R)
{
	float upstream_R = e.properties.R;
	if (!e.properties.buffered) {
		upstream_R += unbuffered_upstream_R;
	}

	float delay = e.properties.switch_delay;
	delay += v.properties.C * (upstream_R + 0.5 * v.properties.R);

	/*zlog_level(delta_log, ROUTER_V3, " [edge_delay: %g edge_R: %g node_R: %g node_C: %g] ", e.properties.switch_delay, e.properties.R, v.properties.R, v.properties.C);*/

	return delay;
}

float get_congestion_cost(const congestion_t &congestion, int cost_index)
{
	extern t_rr_indexed_data *rr_indexed_data;
	/*zlog_level(delta_log, ROUTER_V3, " [pres: %g acc: %g] ", v.properties.pres_cost, v.properties.acc_cost);*/
	return rr_indexed_data[cost_index].base_cost * congestion.acc_cost * congestion.pres_cost;
}

float get_known_cost(const RRGraph &g, const RREdge &e, float criticality_fac, float unbuffered_upstream_R)
{
	/*const auto &target = get_target(g, e);*/

	/*float delay = get_delay(e, target, unbuffered_upstream_R);*/
	/*float congestion = get_congestion_cost(target);*/

	/*zlog_level(delta_log, ROUTER_V3, " [delay: %g congestion %g crit_fac: %g] ", delay, congestion, criticality_fac);*/

	/*return criticality_fac * delay + (1 - criticality_fac) * congestion;*/
	assert(false);
	return 0;
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

float get_timing_driven_expected_cost(const RRNode &current, const RRNode &target,
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
void expand_neighbors_fast(const RRGraph &g, const route_state_t &current, const RRNode &target, float criticality_fac, float astar_fac, std::priority_queue<route_state_t> &heap, const ShouldExpandFunc &should_expand, perf_t *perf)
{
	for_all_out_edges(g, get_vertex(g, current.rr_node), [&heap, &g, &current, &target, &criticality_fac, &astar_fac, &should_expand, &perf] (const RREdge &e) -> void {
			int neighbor_id = get_target(g, id(e));
			auto &neighbor = get_vertex(g, neighbor_id);

			char buffer[256];
			sprintf_rr_node(neighbor_id, buffer);
			zlog_level(delta_log, ROUTER_V3, "\tNeighbor: %s ", buffer);
			
			if (!should_expand(neighbor)) {
				return;
			}

			route_state_t item;

			item.rr_node = id(neighbor);
			item.prev_edge = id(e);

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

template<typename ShouldExpandFunc>
void expand_neighbors(const RRGraph &g, const RRNode &current, const route_state_t *state, const congestion_t *congestion, const RRNode &target, float criticality_fac, float astar_fac, std::priority_queue<route_state_t> &heap, const ShouldExpandFunc &should_expand, bool lock, perf_t *perf, lock_perf_t *lock_perf)
{
	for_all_out_edges(g, current, [&heap, &g, &current, &state, &congestion, &target, &criticality_fac, &astar_fac, &should_expand, &perf, &lock_perf, &lock] (const RREdge &e) -> void {
			int neighbor_id = get_target(g, id(e));
			const auto &neighbor = get_vertex(g, neighbor_id);

			char buffer[256];
			sprintf_rr_node(neighbor_id, buffer);
			zlog_level(delta_log, ROUTER_V3, "\tNeighbor: %s ", buffer);

			if (perf) {
				++perf->num_neighbor_visits;
			}
			
			if (!should_expand(neighbor)) {
				return;
			}

			route_state_t item;

			item.rr_node = neighbor_id;
			item.prev_edge = id(e);

			const route_state_t *current_state = &state[id(current)];

			float unbuffered_upstream_R = current_state->upstream_R;
			float upstream_R = e.properties.R + neighbor.properties.R;
			if (!e.properties.buffered) {
				upstream_R += unbuffered_upstream_R;
			}
			item.upstream_R = upstream_R;
			
			/*if (lock) {*/
				/*[>if (lock_perf) {<]*/
					/*[>++lock_perf->num_lock_tries;<]*/
				/*[>}<]*/
				/*[>if (!neighbor.properties.lock->try_lock()) {<]*/
					/*[>if (lock_perf) {<]*/
						/*[>++lock_perf->num_lock_waits;<]*/
					/*[>}<]*/
					/*[>neighbor.properties.lock->lock();<]*/
				/*[>} <]*/
				/*neighbor.properties.lock->lock();*/
			/*}*/
			float congestion_cost = get_congestion_cost(congestion[item.rr_node], neighbor.properties.cost_index);
			/*if (lock) {*/
				/*neighbor.properties.lock->unlock();*/
			/*}*/
			float delay = get_delay(e, neighbor, unbuffered_upstream_R);

			item.delay = current_state->delay + delay;

			float known_cost = criticality_fac * delay + (1 - criticality_fac) * congestion_cost;
			item.known_cost = current_state->known_cost + known_cost;

			float expected_cost = get_timing_driven_expected_cost(get_vertex(g, item.rr_node), target, criticality_fac, upstream_R);
			item.cost = item.known_cost + astar_fac * expected_cost;

			heap.push(item);

			if (perf) {
				++perf->num_heap_pushes;
			}

			zlog_level(delta_log, ROUTER_V3, " [cost: %g known_cost: %g][occ/cap: %d/%d pres: %g acc: %g][edge_delay: %g edge_R: %g node_R: %g node_C: %g] \n",
					item.cost, item.known_cost, 
					congestion[item.rr_node].occ, neighbor.properties.capacity, congestion[item.rr_node].pres_cost, congestion[item.rr_node].acc_cost,
					e.properties.switch_delay, e.properties.R, neighbor.properties.R, neighbor.properties.C);
	});
}


void route_net_fast(RRGraph &g, const net_t &net, const route_parameters_t &params, route_state_t *state, route_tree_t &rt, trace_t &trace, t_net_timing &net_timing, perf_t *perf)
{
	std::priority_queue<route_state_t> heap;

	vector<int> modified;

	zlog_level(delta_log, ROUTER_V1, "Routing VPR net id %d\n", net.vpr_id);

	if (route_tree_empty(rt)) {
		/* special case */
		/*update_one_cost_internal(get_vertex(g, net.current_source.rr_node), 1, params.pres_fac);*/

		/*route_tree_set_source(rt, get_vertex(g, net.source.rr_node));*/
	} else {
		RouteTreeNode rt_root = get_vertex(rt.graph, rt.root_rt_node_id);
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
			zlog_level(delta_log, ROUTER_V2, "Current: %s prev: %d old_delay: %g new_delay: %g old_known: %g new_known: %g old_cost: %g new_cost: %g\n", buffer, item.prev_edge != -1 ? get_source(g, item.prev_edge) : -1, state[item.rr_node].delay, item.delay, state[item.rr_node].known_cost, item.known_cost, state[item.rr_node].cost, item.cost);

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
					auto &neighbor = get_vertex(g, get_target(g, e_i));

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
					new_item.prev_edge = id(e);

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
		/*update_one_cost(g, new_path.begin(), new_path.end(), 1, params.pres_fac, false, nullptr);*/

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

int get_previous_edge(int rr_node_id, const route_state_t *state, const route_tree_t &rt, const RRGraph &g)
{
	int previous_edge;

	if (state[rr_node_id].prev_edge == -1) {
		previous_edge = -1;
	} else {
		auto rt_node = route_tree_get_rt_node(rt, rr_node_id);

		char buffer[256];
		sprintf_rr_node(rr_node_id, buffer);

		if (rt_node) {
			/* already reach an existing route tree node but the current state suggests there's a possible path of 
			 * lower cost */
			/*if (rt_node->properties.rr_edge_to_parent != -1) {*/
				/*char parent[256];*/
				/*sprintf_rr_node(get_source(g, rt_node->properties.rr_edge_to_parent), parent);*/
				/*zlog_error(delta_log, "Error: Existing route tree node %s has non-null rr_edge_to_parent that connects to %s\n", buffer, parent);*/
				/*assert(false);*/
			/*}*/

			char s_state[256];
			char s_rt[256];
			sprintf_rr_node(get_source(g, state[rr_node_id].prev_edge), s_state);
			sprintf_rr_node(get_source(g, rt_node->properties.rr_edge_to_parent), s_rt);
			zlog_warn(delta_log, "Warning: Existing route tree node %s does not have a matching route state. (state.prev_edge: %s rt_node.properties.rr_edge_to_parent: %s) because we have found a shorter path to that node\n", buffer, s_state, s_rt);

			previous_edge = -1;
		} else {
			previous_edge = state[rr_node_id].prev_edge;
		}
	} 

	return previous_edge;
}

vector<path_node_t> get_path(int sink_rr_node_id, const route_state_t *state, const route_tree_t &rt, const RRGraph &g, int vpr_net_id)
{
	assert(get_vertex(g, sink_rr_node_id).properties.type == SINK);

	int current_rr_node_id = sink_rr_node_id;
	int previous_edge;
	vector<path_node_t> path;

	while ((previous_edge = get_previous_edge(current_rr_node_id, state, rt, g)) != -1) {
		/* parent */
		int parent_rr_node_id = get_source(g, previous_edge);

		path_node_t node;
		node.rr_node_id = current_rr_node_id;
		node.prev_edge = previous_edge;
		node.update_cost = true;

		path.emplace_back(node);

		char p[256];
		char c[256];
		/* printing */
		sprintf_rr_node(parent_rr_node_id, p);
		sprintf_rr_node(current_rr_node_id, c);
		zlog_level(delta_log, ROUTER_V2, "Net %d path: %s -> %s\n", vpr_net_id, p, c);

		current_rr_node_id = parent_rr_node_id;
	}

	path_node_t node;
	node.rr_node_id = current_rr_node_id;
	node.prev_edge = -1;
	node.update_cost = get_vertex(g, current_rr_node_id).properties.type == SOURCE;

	path.emplace_back(node);

	return path;
}

vector<sink_t *> route_net_3(const RRGraph &g, const vector<int> &pid, int this_pid, int vpr_id, const source_t *source, const vector<sink_t *> &sinks, const route_parameters_t &params, route_state_t *state, congestion_t *congestion, route_tree_t &rt, t_net_timing &net_timing, vector<vector<int>> &unrouted_sinks_boundary_nodes, int num_partitions, bool lock, perf_t *perf, lock_perf_t *lock_perf)
{
	std::priority_queue<route_state_t> heap;

	vector<int> modified;

	zlog_level(delta_log, ROUTER_V1, "Routing net %d (%lu sinks)\n", vpr_id, sinks.size());

	vector<sink_t *> sorted_sinks = sinks;
	std::sort(begin(sorted_sinks), end(sorted_sinks), [] (const sink_t *a, const sink_t *b) -> bool {
			return a->criticality_fac > b->criticality_fac;
			});

	char buffer[256];

	if (route_tree_empty(rt)) {
		/* special case */
		sprintf_rr_node(source->rr_node, buffer);
		zlog_level(delta_log, ROUTER_V2, "Empty route tree. Setting root to %s\n", buffer);

		const auto &source_rr_node = get_vertex(g, source->rr_node);
		RouteTreeNode *root_rt_node = route_tree_add_rr_node(rt, source_rr_node);
		route_tree_set_node_properties(*root_rt_node, true, -1, source_rr_node.properties.R, 0.5 * source_rr_node.properties.R * source_rr_node.properties.C);
		route_tree_add_root(rt, source->rr_node);

		/*update_one_cost_internal(source_rr_node, 1, params.pres_fac);*/
	} else {
		if (rt.root_rt_node_id != -1) {
			const RouteTreeNode &rt_root = get_vertex(rt.graph, rt.root_rt_node_id);
			if (source && rt_root.properties.rr_node != source->rr_node) {
				char root[256];
				char source_str[256];
				sprintf_rr_node(rt_root.properties.rr_node, root);
				sprintf_rr_node(source->rr_node, source_str);
				zlog_error(delta_log, "Error: Non empty route tree node has root %s that is different from net source %s\n",
						root, source_str);
				assert(false);
			}
		}
	}

	int isink = 0;
	int num_routed_sinks = 0;
	vector<sink_t *> unrouted_sinks;
	for (const auto &sink : sorted_sinks) {
		sprintf_rr_node(sink->rr_node, buffer);
		zlog_level(delta_log, ROUTER_V1, "Net %d Sink %d: %s criticality: %g BB: %d-%d %d-%d Prev BB: %d-%d %d-%d\n", vpr_id, sink->id, buffer, sink->criticality_fac, sink->current_bounding_box.xmin, sink->current_bounding_box.xmax, sink->current_bounding_box.ymin, sink->current_bounding_box.ymax, sink->previous_bounding_box.xmin, sink->previous_bounding_box.xmax, sink->previous_bounding_box.ymin, sink->previous_bounding_box.ymax);

		const auto &sink_rr_node = get_vertex(g, sink->rr_node);

		route_tree_multi_root_add_to_heap(rt, g, sink_rr_node, sink->criticality_fac, params.astar_fac, sink->current_bounding_box, heap, perf);

		priority_queue<pair<float, int>, vector<pair<float, int>>, std::greater<pair<float, int>>> boundary_nodes;

		bool found_sink = false;
		while (!heap.empty() && !found_sink) {
			auto item = heap.top();
			heap.pop();

			assert(pid[item.rr_node] == -1 || pid[item.rr_node] == this_pid);

			if (perf) {
				++perf->num_heap_pops;
			}

			sprintf_rr_node(item.rr_node, buffer);
			const auto &v = get_vertex(g, item.rr_node);
			zlog_level(delta_log, ROUTER_V3, "Current: %s occ/cap: %d/%d prev: %d new_cost: %g new_known: %g new_delay: %g old_cost: %g old_known: %g old_delay: %g\n", buffer, congestion[item.rr_node].occ, v.properties.capacity, item.prev_edge != -1 ? get_source(g, item.prev_edge) : -1, item.cost, item.known_cost, item.delay, state[item.rr_node].cost, state[item.rr_node].known_cost, state[item.rr_node].delay);

			if (item.rr_node == sink->rr_node) {
				state[item.rr_node] = item;
				modified.push_back(item.rr_node);
				found_sink = true;
			} else if (item.known_cost < state[item.rr_node].known_cost && item.cost < state[item.rr_node].cost) {
				/*if (route_tree_get_rt_node(rt, item.rr_node)) {*/
					/*sprintf_rr_node(item.rr_node, buffer);*/
					/*zlog_warn(delta_log, "Warning: Found a lower cost path to existing route tree node %s\n", buffer);*/
					/*assert(false);*/
				/*}*/

				state[item.rr_node] = item;
				modified.push_back(item.rr_node);

				expand_neighbors(g, v, state, congestion, sink_rr_node, sink->criticality_fac, params.astar_fac, heap, [&sink, &sink_rr_node, &v, &pid, &boundary_nodes, &item, &this_pid] (const RRNode &n) -> bool {

					/*if (trace_has_node(prev_trace, id(n))) {*/
						/*zlog_level(delta_log, ROUTER_V3, " existing node route tree ");*/
					/*}*/
					const auto &prop = n.properties;

					if (prop.xhigh < sink->current_bounding_box.xmin
							|| prop.xlow > sink->current_bounding_box.xmax
							|| prop.yhigh < sink->current_bounding_box.ymin
							|| prop.ylow > sink->current_bounding_box.ymax) {
					zlog_level(delta_log, ROUTER_V3, "outside of bounding box\n");
					return false;
					}

					if (prop.type == IPIN
							&& (prop.xhigh != sink_rr_node.properties.xhigh 
								|| prop.yhigh != sink_rr_node.properties.yhigh)) {
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
					if (pid[id(n)] != -1 && pid[id(n)] != this_pid) {
						zlog_level(delta_log, ROUTER_V3, " boundary node %d with cost %g\n", id(v), item.cost);
						boundary_nodes.push(make_pair(item.cost, id(v)));
						return false;
					}

					return true;
				}, lock, perf, lock_perf);
			} else {
				zlog_level(delta_log, ROUTER_V3, "\tNot expanding neighbor because known_cost %g > %g or cost %g > %g\n", item.known_cost, state[item.rr_node].known_cost, item.cost, state[item.rr_node].cost);
			}
		}

		if (!found_sink) {
			assert(heap.empty());
			zlog_error(delta_log, "Error: Failed to find sink %d\n", sink->rr_node);
			vector<int> lowest_cost_boundary_nodes;
			int num_visited_partitions = 0;
			vector<bool> visited_partitions(num_partitions, false);
			set<int> added_boundary_nodes;
			while (!boundary_nodes.empty() && num_visited_partitions < num_partitions-1) {
				auto bn = boundary_nodes.top(); boundary_nodes.pop();
				if (added_boundary_nodes.find(bn.second) == added_boundary_nodes.end()) {
					lowest_cost_boundary_nodes.push_back(bn.second);
					added_boundary_nodes.insert(bn.second);
					for (const auto &e : get_out_edges(g, bn.second)) {
						int to = get_target(g, e);
						int to_pid = pid[to];
						if (to_pid != -1 && to_pid != this_pid && !visited_partitions[to_pid]) {
							visited_partitions[to_pid] = true;
							++num_visited_partitions;
						}
					}
				}
			}
			if (num_visited_partitions < num_partitions-1) {
				zlog_level(delta_log, ROUTER_V3, "Net %d sink %d does't have boundaries nodes that covers all partitions\n");
			}
			unrouted_sinks_boundary_nodes.emplace_back(lowest_cost_boundary_nodes);
			unrouted_sinks.push_back(sink);
		} else {
			const auto &path = get_path(sink->rr_node, state, rt, g, vpr_id);

			route_tree_add_path(rt, path, g, state);

			vector<int> added_rr_nodes;
			for (const auto &node : path) {
				if (node.update_cost) {
					added_rr_nodes.push_back(node.rr_node_id);
				}
			}

			update_one_cost(g, congestion, added_rr_nodes.begin(), added_rr_nodes.end(), 1, params.pres_fac, lock, lock_perf);

			net_timing.delay[sink->id+1] = route_tree_get_rt_node(rt, sink->rr_node)->properties.delay;

			++num_routed_sinks;
		}

		for (const auto &m : modified)  {
			state[m].known_cost = std::numeric_limits<float>::max();
			state[m].cost = std::numeric_limits<float>::max();
		}

		heap = std::priority_queue<route_state_t>();
	}

	zlog_level(delta_log, ROUTER_V1, "\n");

	return unrouted_sinks;

	/*delete [] state;*/

	/*check_route_tree(rt, net, g);*/
}

vector<sink_t *> route_net_2(const RRGraph &g, int vpr_id, const source_t *source, const vector<sink_t *> &sinks, const route_parameters_t &params, route_state_t *state, congestion_t *congestion, route_tree_t &rt, t_net_timing &net_timing, bool lock, perf_t *perf, lock_perf_t *lock_perf)
{
	std::priority_queue<route_state_t> heap;

	vector<int> modified;

	zlog_level(delta_log, ROUTER_V1, "Routing net %d (%lu sinks)\n", vpr_id, sinks.size());

	vector<sink_t *> sorted_sinks = sinks;
	std::sort(begin(sorted_sinks), end(sorted_sinks), [] (const sink_t *a, const sink_t *b) -> bool {
			return a->criticality_fac > b->criticality_fac;
			});

	char buffer[256];

	sprintf_rr_node(source->rr_node, buffer);
	zlog_level(delta_log, ROUTER_V1, "Source: %s\n", buffer);

	if (route_tree_empty(rt)) {
		/* special case */
		zlog_level(delta_log, ROUTER_V2, "Empty route tree. Setting root to %s\n", buffer);

		const auto &source_rr_node = get_vertex(g, source->rr_node);
		RouteTreeNode *root_rt_node = route_tree_add_rr_node(rt, source_rr_node);
		route_tree_set_node_properties(*root_rt_node, true, -1, source_rr_node.properties.R, 0.5 * source_rr_node.properties.R * source_rr_node.properties.C);
		route_tree_set_root(rt, source->rr_node);

		/*update_one_cost_internal(source_rr_node, 1, params.pres_fac);*/
	} else {
		const RouteTreeNode &rt_root = get_vertex(rt.graph, rt.root_rt_node_id);
		if (source && rt_root.properties.rr_node != source->rr_node) {
			char root[256];
			char source_str[256];
			sprintf_rr_node(rt_root.properties.rr_node, root);
			sprintf_rr_node(source->rr_node, source_str);
			zlog_error(delta_log, "Error: Non empty route tree node has root %s that is different from net source %s\n",
					root, source_str);
			assert(false);
		}
	}

	int isink = 0;
	int num_routed_sinks = 0;
	vector<sink_t *> unrouted_sinks;
	for (const auto &sink : sorted_sinks) {
		sprintf_rr_node(sink->rr_node, buffer);
		zlog_level(delta_log, ROUTER_V1, "Net %d Sink %d: %s criticality: %g BB: %d-%d %d-%d Prev BB: %d-%d %d-%d\n", vpr_id, sink->id, buffer, sink->criticality_fac, sink->current_bounding_box.xmin, sink->current_bounding_box.xmax, sink->current_bounding_box.ymin, sink->current_bounding_box.ymax, sink->previous_bounding_box.xmin, sink->previous_bounding_box.xmax, sink->previous_bounding_box.ymin, sink->previous_bounding_box.ymax);

		const auto &sink_rr_node = get_vertex(g, sink->rr_node);

		route_tree_add_to_heap(rt, g, sink_rr_node, sink->criticality_fac, params.astar_fac, sink->current_bounding_box, heap, perf);

		bool found_sink = false;
		while (!heap.empty() && !found_sink) {
			auto item = heap.top();
			heap.pop();

			if (perf) {
				++perf->num_heap_pops;
			}

			sprintf_rr_node(item.rr_node, buffer);
			const auto &v = get_vertex(g, item.rr_node);
			zlog_level(delta_log, ROUTER_V3, "Current: %s occ/cap: %d/%d prev: %d new_cost: %g new_known: %g new_delay: %g old_cost: %g old_known: %g old_delay: %g\n", buffer, congestion[item.rr_node].occ, v.properties.capacity, item.prev_edge != -1 ? get_source(g, item.prev_edge) : -1, item.cost, item.known_cost, item.delay, state[item.rr_node].cost, state[item.rr_node].known_cost, state[item.rr_node].delay);

			if (item.rr_node == sink->rr_node) {
				state[item.rr_node] = item;
				modified.push_back(item.rr_node);
				found_sink = true;
			} else if (item.known_cost < state[item.rr_node].known_cost && item.cost < state[item.rr_node].cost) {
				/*if (route_tree_get_rt_node(rt, item.rr_node)) {*/
					/*sprintf_rr_node(item.rr_node, buffer);*/
					/*zlog_warn(delta_log, "Warning: Found a lower cost path to existing route tree node %s\n", buffer);*/
					/*assert(false);*/
				/*}*/

				state[item.rr_node] = item;
				modified.push_back(item.rr_node);

				expand_neighbors(g, v, state, congestion, sink_rr_node, sink->criticality_fac, params.astar_fac, heap, [&sink, &sink_rr_node] (const RRNode &v) -> bool {

					/*if (trace_has_node(prev_trace, id(v))) {*/
						/*zlog_level(delta_log, ROUTER_V3, " existing node route tree ");*/
					/*}*/
					const auto &prop = v.properties;

					if (prop.xhigh < sink->current_bounding_box.xmin
							|| prop.xlow > sink->current_bounding_box.xmax
							|| prop.yhigh < sink->current_bounding_box.ymin
							|| prop.ylow > sink->current_bounding_box.ymax) {
					zlog_level(delta_log, ROUTER_V3, "outside of bounding box\n");
					return false;
					}

					if (prop.type == IPIN
							&& (prop.xhigh != sink_rr_node.properties.xhigh 
								|| prop.yhigh != sink_rr_node.properties.yhigh)) {
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
				}, lock, perf, lock_perf);
			} else {
				zlog_level(delta_log, ROUTER_V3, "\tNot expanding neighbor because known_cost %g > %g or cost %g > %g\n", item.known_cost, state[item.rr_node].known_cost, item.cost, state[item.rr_node].cost);
			}
		}

		if (!found_sink) {
			zlog_error(delta_log, "Error: Failed to find sink\n");
			unrouted_sinks.push_back(sink);
			assert(heap.empty());
		} else {
			vector<int> added_rr_nodes;
			route_tree_add_path(rt, g, state, sink->rr_node, vpr_id, added_rr_nodes);

			update_one_cost(g, congestion, added_rr_nodes.begin(), added_rr_nodes.end(), 1, params.pres_fac, lock, lock_perf);

			net_timing.delay[sink->id+1] = route_tree_get_rt_node(rt, sink->rr_node)->properties.delay;

			++num_routed_sinks;
		}

		for (const auto &m : modified)  {
			state[m].known_cost = std::numeric_limits<float>::max();
			state[m].cost = std::numeric_limits<float>::max();
		}

		heap = std::priority_queue<route_state_t>();
	}

	zlog_level(delta_log, ROUTER_V1, "\n");

	return unrouted_sinks;

	/*delete [] state;*/

	/*check_route_tree(rt, net, g);*/
}

void route_net(RRGraph &g, int vpr_id, const source_t *source, const vector<sink_t *> &sinks, const route_parameters_t &params, route_state_t *state, route_tree_t &rt, trace_t &trace, const trace_t &prev_trace, t_net_timing &net_timing, perf_t *perf)
{
	std::priority_queue<route_state_t> heap;

	vector<int> modified;

	zlog_level(delta_log, ROUTER_V1, "Routing VPR net id %d\n", vpr_id);

	if (route_tree_empty(rt)) {
		/* special case */
		/*update_one_cost_internal(get_vertex(g, net.current_source.rr_node), 1, params.pres_fac);*/

		/*route_tree_set_source(rt, get_vertex(g, source->rr_node));*/
		auto &source_rr_node = get_vertex(g, source->rr_node);
		RouteTreeNode *root_rt_node = route_tree_add_rr_node(rt, source_rr_node);
		route_tree_set_node_properties(*root_rt_node, true, -1, source_rr_node.properties.R, 0.5 * source_rr_node.properties.R * source_rr_node.properties.C);
		route_tree_set_root(rt, source->rr_node);
	} else {
		RouteTreeNode rt_root = get_vertex(rt.graph, rt.root_rt_node_id);
		if (rt_root.properties.rr_node != source->rr_node) {
			char root[256];
			char source_str[256];
			sprintf_rr_node(rt_root.properties.rr_node, root);
			sprintf_rr_node(source->rr_node, source_str);
			zlog_error(delta_log, "Error: Non empty route tree node has root %s that is different from net source %s\n",
					root, source_str);
			assert(false);
		}
	}

	vector<sink_t *> sorted_sinks = sinks;
	std::sort(sorted_sinks.begin(), sorted_sinks.end());

	char buffer[256];

	sprintf_rr_node(source->rr_node, buffer);
	zlog_level(delta_log, ROUTER_V1, "Source: %s\n", buffer);

	int isink = 0;
	for (const auto &sink : sorted_sinks) {
		sprintf_rr_node(sink->rr_node, buffer);
		zlog_level(delta_log, ROUTER_V1, "Net %d Sink %d: %s criticality: %g BB: %d-%d %d-%d Prev BB: %d-%d %d-%d\n", vpr_id, sink->id, buffer, sink->criticality_fac, sink->current_bounding_box.xmin, sink->current_bounding_box.xmax, sink->current_bounding_box.ymin, sink->current_bounding_box.ymax, sink->previous_bounding_box.xmin, sink->previous_bounding_box.xmax, sink->previous_bounding_box.ymin, sink->previous_bounding_box.ymax);
		route_tree_add_to_heap(rt, g, get_vertex(g, sink->rr_node), sink->criticality_fac, params.astar_fac, sink->current_bounding_box, heap, perf);

		bool found_sink = false;
		while (!heap.empty() && !found_sink) {
			auto item = heap.top();
			heap.pop();

			sprintf_rr_node(item.rr_node, buffer);
			const auto &v = get_vertex(g, item.rr_node);
			/*zlog_level(delta_log, ROUTER_V2, "Current: %s occ/cap: %d/%d prev: %d old_cost: %g new_cost: %g old_delay: %g new_delay: %g old_known: %g new_known: %g \n", buffer, congestion[item.rr_node].occ, v.properties.capacity, item.prev_edge ? id(get_source(g, *item.prev_edge)) : -1, state[item.rr_node].cost, item.cost, state[item.rr_node].delay, item.delay, state[item.rr_node].known_cost, item.known_cost);*/

			if (item.rr_node == sink->rr_node) {
				state[item.rr_node] = item;
				modified.push_back(item.rr_node);
				found_sink = true;
			} else if (item.known_cost < state[item.rr_node].known_cost && item.cost < state[item.rr_node].cost) {
				/*if (route_tree_get_rt_node(rt, item.rr_node)) {*/
					/*sprintf_rr_node(item.rr_node, buffer);*/
					/*zlog_warn(delta_log, "Warning: Found a lower cost path to existing route tree node %s\n", buffer);*/
					/*assert(false);*/
				/*}*/
				const auto &sink_vertex = get_vertex(g, sink->rr_node);

				state[item.rr_node] = item;
				modified.push_back(item.rr_node);

				expand_neighbors(g, v, nullptr, nullptr, sink_vertex, sink->criticality_fac, params.astar_fac, heap, [&sink, &prev_trace, &sink_vertex] (const RRNode &v) -> bool {

					if (trace_has_node(prev_trace, id(v))) {
						zlog_level(delta_log, ROUTER_V3, " existing node route tree ");
					}
					const auto &prop = v.properties;

					if (prop.xhigh < sink->current_bounding_box.xmin
							|| prop.xlow > sink->current_bounding_box.xmax
							|| prop.yhigh < sink->current_bounding_box.ymin
							|| prop.ylow > sink->current_bounding_box.ymax) {
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
				}, false, perf, nullptr);
			} else {
				zlog_level(delta_log, ROUTER_V3, "\tNot expanding neighbor because known_cost %g > %g or cost %g > %g\n", item.known_cost, state[item.rr_node].known_cost, item.cost, state[item.rr_node].cost);
			}
		}

		if (!found_sink) {
			zlog_error(delta_log, "Failed to find sink\n");
			assert(false);
		}

		const auto &new_path = trace_add_path(trace, g, state, sink->rr_node, vpr_id);
		/* new_path.second-1 because we do not update the last node because it's the node in existing route tree (which was already updated previously */
		/*int new_path_first_rr_node = *(new_path.end()-1);*/
		/* if the first node is a connection to existing route, don't update its cost again */
		/*if (route_tree_get_rt_node(rt, new_path_first_rr_node)) {*/
			/*update_one_cost(g, new_path.begin(), new_path.end()-1, 1, params.pres_fac);*/
		/*} else {*/
			/*update_one_cost(g, new_path.begin(), new_path.end(), 1, params.pres_fac);*/
		/*}*/
		/*update_one_cost(g, new_path.begin(), new_path.end(), 1, params.pres_fac, false, nullptr);*/

		/* important to update route tree only after updating cost because the previous check relies on whether the
		 * first node is already in the route tree */
		route_tree_add_path(rt, g, state, new_path, vpr_id);

		net_timing.delay[sink->id+1] = route_tree_get_rt_node(rt, sink->rr_node)->properties.delay;

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

		/*route_tree_set_source(rt, get_vertex(g, net.source.rr_node));*/
	} else {
		RouteTreeNode rt_root = get_vertex(rt.graph, rt.root_rt_node_id);
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
			zlog_level(delta_log, ROUTER_V2, "Current: %s prev: %d old_delay: %g new_delay: %g old_known: %g new_known: %g old_cost: %g new_cost: %g\n", buffer, item.prev_edge != -1 ? get_source(g, item.prev_edge) : -1, state[item.rr_node].delay, item.delay, state[item.rr_node].known_cost, item.known_cost, state[item.rr_node].cost, item.cost);

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

				expand_neighbors(g, get_vertex(g, item.rr_node), nullptr, nullptr, sink_vertex, sink.criticality_fac, params.astar_fac, heap, [&net, &sink, &sink_vertex] (const RRNode &v) -> bool {
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
				}, false, perf, nullptr);
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
		/*update_one_cost(g, new_path.begin(), new_path.end(), 1, params.pres_fac, false, nullptr);*/

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

