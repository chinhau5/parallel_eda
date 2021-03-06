#include "pch.h"

#include "vpr_types.h"
#include "path_delay.h"

#include "route.h"
#include "route_tree.h"
#include "congestion.h"
#include "trace.h"
#include "router.h"
#include "filtered_graph.h"
#include "log.h"
#include "expand.h"
#include "router.h"

using namespace std;

bool operator<(const sink_t &a, const sink_t &b)
{
	return a.criticality_fac > b.criticality_fac;
}

bool operator<(const route_state_t &a, const route_state_t &b)
{
	return a.cost > b.cost;
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

void check_route_tree_internal(const route_tree_t &rt, RouteTreeNode rt_node, RRGraph &g, vector<int> &visited_sinks, vector<int> &visited_nodes)
{
	const auto &rt_node_p = get_vertex_props(rt.graph, rt_node);

	assert(rt_node_p.valid);

	int rr_node = rt_node_p.rr_node;
	auto &rr_node_p = get_vertex_props(g, rr_node);
	if (rr_node_p.type == SINK) {
		visited_sinks.push_back(rr_node);
		char buffer[256];
		sprintf_rr_node(rr_node, buffer);
		/*zlog_level(delta_log, ROUTER_V2, "route_tree_check: %s\n", buffer);*/
	}

	visited_nodes.push_back(rr_node);

	for (const auto &branch : route_tree_get_branches(rt, rt_node)) {
		/*for_all_out_edges(rt.graph, node, [&rt, &g, &visited_sinks, &visited_nodes] (const RouteTreeEdge &e) -> void {*/
		const auto &child = get_target(rt.graph, branch);
		check_route_tree_internal(rt, child, g, visited_sinks, visited_nodes);
	}
}

void check_route_tree(const route_tree_t &rt, const net_t &net, const vector<sink_t *> &routed_sinks, RRGraph &g)
{
	RouteTreeNode rt_root = route_tree_get_rt_node(rt, net.source.rr_node);

	vector<int> sinks;
	for (const auto &sink : routed_sinks) {
		sinks.push_back(sink->rr_node);
	}

	vector<int> visited_sinks;
	vector<int> visited_nodes;

	check_route_tree_internal(rt, rt_root, g, visited_sinks, visited_nodes);

	int num_rt_nodes = 0;
	for (const auto &rt_node : route_tree_get_nodes(rt)) {
		++num_rt_nodes;
	}
	assert(num_rt_nodes == num_edges(rt.graph)+1);
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
				[&duplicated_nodes, &rt] (RouteTreeNode rt_node) -> string {
				char s_rr_node[256];
				char buffer[256];
				const auto &rt_node_p = get_vertex_props(rt.graph, rt_node);
				if (find(begin(duplicated_nodes), end(duplicated_nodes), rt_node_p.rr_node) !=
						end(duplicated_nodes)) {
				sprintf_rr_node(rt_node_p.rr_node, s_rr_node);
				sprintf(buffer, "label=\"%s\" style=filled fillcolor=red", s_rr_node);
				} else {
				sprintf_rr_node(rt_node_p.rr_node, s_rr_node);
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
				[&rt] (const RouteTreeNode &rt_node) -> string {
					char buffer[256];
					sprintf_rr_node(get_vertex_props(rt.graph, rt_node).rr_node, buffer);
					return buffer;
				},
				[] (const RouteTreeEdge &rt_edge) -> string {
					return "";
				},
				[&rt] (const RouteTreeNode &rt_node) -> bool {
					return !get_vertex_props(rt.graph, rt_node).valid;
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

void check_route_tree(const route_tree_t &rt, const net_t &net, RRGraph &g)
{
	RouteTreeNode rt_root = route_tree_get_rt_node(rt, net.source.rr_node);
	vector<int> sinks;
	for (const auto &s : net.sinks) {
		sinks.push_back(s.rr_node);
	}
	vector<int> visited_sinks;

	vector<int> visited_nodes;
	check_route_tree_internal(rt, rt_root, g, visited_sinks, visited_nodes);
	int num_rt_nodes = 0;
	for (const auto &rt_node : route_tree_get_nodes(rt)) {
		++num_rt_nodes;
	}
	assert(num_rt_nodes == num_edges(rt.graph)+1);
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
				[&duplicated_nodes, &rt] (RouteTreeNode rt_node) -> string {
				char s_rr_node[256];
				char buffer[256];
				const auto &rt_node_p = get_vertex_props(rt.graph, rt_node);
				if (find(begin(duplicated_nodes), end(duplicated_nodes), rt_node_p.rr_node) !=
						end(duplicated_nodes)) {
				sprintf_rr_node(rt_node_p.rr_node, s_rr_node);
				sprintf(buffer, "label=\"%s\" style=filled fillcolor=red", s_rr_node);
				} else {
				sprintf_rr_node(rt_node_p.rr_node, s_rr_node);
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

	assert(!sinks.empty());
	
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
				[&rt] (const RouteTreeNode &rt_node) -> string {
					char buffer[256];
					sprintf_rr_node(get_vertex_props(rt.graph, rt_node).rr_node, buffer);
					return buffer;
				},
				[] (const RouteTreeEdge &rt_edge) -> string {
					return "";
				},
				[&rt] (const RouteTreeNode &rt_node) -> bool {
					return !get_vertex_props(rt.graph, rt_node).valid;
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

float get_delay(const rr_edge_property_t &e, const rr_node_property_t &v, float unbuffered_upstream_R)
{
	extern struct s_switch_inf *switch_inf;
	const struct s_switch_inf *sw = &switch_inf[e.switch_index];

	float upstream_R = sw->R;

	if (!sw->buffered) {
		upstream_R += unbuffered_upstream_R;
	}

	/*zlog_level(delta_log, ROUTER_V3, " [edge_delay: %g edge_R: %g node_R: %g node_C: %g] ", sw->Tdel, sw->R, v.R, v.C);*/

	return sw->Tdel + v.C * (upstream_R + 0.5 * v.R);
}

float get_congestion_cost_mpi_recv(RRNode rr_node, const RRGraph &g, congestion_t *congestion, MPI_Comm comm, int this_pid, int num_procs, float pres_fac)
{
	extern t_rr_indexed_data *rr_indexed_data;
	/*zlog_level(delta_log, ROUTER_V3, " [pres: %g acc: %g] ", v.pres_cost, v.acc_cost);*/
	const auto &rr_node_p = get_vertex_props(g, rr_node);

	//for (int i = 0; i < num_procs; ++i) {
		//if (i != this_pid) {
			//int flag;
			//MPI_Status status;
			//assert(MPI_Iprobe(i, rr_node, comm, &flag, &status) == MPI_SUCCESS);

			//if (flag) {
				//int delta;

				//assert(MPI_Recv(&delta, 1, MPI_INT, status.MPI_SOURCE, status.MPI_TAG, comm, MPI_STATUS_IGNORE) == MPI_SUCCESS);

				//congestion[rr_node].occ += delta;

				////assert(MPI_Iprobe(i, rr_node, comm, &flag, &status) == MPI_SUCCESS);
			//}
		//}
	//}

	//assert(MPI_Iprobe(MPI_ANY_SOURCE, rr_node, comm, &flag, &status) == MPI_SUCCESS);

	//while (flag) {
		//int delta;

		//assert(MPI_Recv(&delta, 1, MPI_INT, status.MPI_SOURCE, status.MPI_TAG, comm, MPI_STATUS_IGNORE) == MPI_SUCCESS);

		//congestion[rr_node].occ += delta;

		//assert(MPI_Iprobe(MPI_ANY_SOURCE, rr_node, comm, &flag, &status) == MPI_SUCCESS);
	//}

	float pres_cost;
	if (congestion[rr_node].occ < rr_node_p.capacity) {
		pres_cost = 1;
	} else {
		pres_cost = 1 + (congestion[rr_node].occ + 1 - rr_node_p.capacity) * pres_fac;
	}

	float cong_cost = rr_indexed_data[rr_node_p.cost_index].base_cost * congestion[rr_node].acc_cost * pres_cost;

	return cong_cost;
}

float get_congestion_cost_mpi_rma(RRNode rr_node, const RRGraph &g, const vector<int> &pid, int this_pid, congestion_t *congestion, MPI_Win win, int cost_index, float pres_fac)
{
	extern t_rr_indexed_data *rr_indexed_data;
	/*zlog_level(delta_log, ROUTER_V3, " [pres: %g acc: %g] ", v.pres_cost, v.acc_cost);*/
	const auto &rr_node_p = get_vertex_props(g, rr_node);

	int rr_node_pid = pid[rr_node];
	int from_pid = rr_node_pid == -1 ? 0 : rr_node_pid;

	assert(MPI_Win_lock(MPI_LOCK_SHARED, from_pid, 0, win) == MPI_SUCCESS);

	if (from_pid != this_pid) {
		congestion[rr_node].occ = std::numeric_limits<int>::min();
		assert(MPI_Get(&congestion[rr_node].occ, 1, get_occ_dt(),
					from_pid,
					get_occ_disp(rr_node), 1, get_occ_dt(),
					win) == MPI_SUCCESS);
		assert(MPI_Win_flush(from_pid, win) == MPI_SUCCESS);
		assert(congestion[rr_node].occ != std::numeric_limits<int>::min());
	} 

	float pres_cost = 1 + (congestion[rr_node].occ + 1 - rr_node_p.capacity) * pres_fac;
	float cong_cost = rr_indexed_data[cost_index].base_cost * congestion[rr_node].acc_cost * pres_cost;

	assert(MPI_Win_unlock(from_pid, win) == MPI_SUCCESS);

	return cong_cost;
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

static int get_expected_segs_to_target(const rr_node_property_t &current, const rr_node_property_t &target,
		int *num_segs_ortho_dir_ptr) {

	/* Returns the number of segments the same type as inode that will be needed *
	 * to reach target_node (not including inode) in each direction (the same    *
	 * direction (horizontal or vertical) as inode and the orthogonal direction).*/

	t_rr_type rr_type;
	int target_x, target_y, num_segs_same_dir, cost_index, ortho_cost_index;
	int no_need_to_pass_by_clb;
	float inv_length, ortho_inv_length, ylow, yhigh, xlow, xhigh;

	extern t_rr_indexed_data *rr_indexed_data;

	target_x = target.xlow;
	target_y = target.ylow;
	cost_index = current.cost_index;
	inv_length = rr_indexed_data[cost_index].inv_length;
	ortho_cost_index = rr_indexed_data[cost_index].ortho_cost_index;
	ortho_inv_length = rr_indexed_data[ortho_cost_index].inv_length;
	rr_type = current.type;

	zlog_level(delta_log, ROUTER_V3, "\t[inv %g ortho_inv %g] [Target xlow %d ylow %d real_xlow %d real_ylow %d] ", inv_length, ortho_inv_length, target.xlow, target.ylow,
			target.real_xlow, target.real_ylow);

	if (rr_type == CHANX) {
		ylow = current.ylow;
		xhigh = current.xhigh;
		xlow = current.xlow;

		/* Count vertical (orthogonal to inode) segs first. */

		if (ylow > target_y) { /* Coming from a row above target? */
			*num_segs_ortho_dir_ptr =
					(int)(ROUND_UP((ylow - target_y + 1.) * ortho_inv_length));
			no_need_to_pass_by_clb = 1;
			zlog_level(delta_log, ROUTER_V3, "[ylow > target_y] ");
		} else if (ylow < target_y - 1) { /* Below the CLB bottom? */
			*num_segs_ortho_dir_ptr = (int)(ROUND_UP((target_y - ylow) *
					ortho_inv_length));
			no_need_to_pass_by_clb = 1;
			zlog_level(delta_log, ROUTER_V3, "[ylow < target_y-1] ");
		} else { /* In a row that passes by target CLB */
			*num_segs_ortho_dir_ptr = 0;
			no_need_to_pass_by_clb = 0;
			zlog_level(delta_log, ROUTER_V3, "[same row] ");
		}

		/* Now count horizontal (same dir. as inode) segs. */

		if (xlow > target_x + no_need_to_pass_by_clb) {
			num_segs_same_dir = (int)(ROUND_UP((xlow - no_need_to_pass_by_clb -
							target_x) * inv_length));
			zlog_level(delta_log, ROUTER_V3, "[xlow > target_x] ");
		} else if (xhigh < target_x - no_need_to_pass_by_clb) {
			num_segs_same_dir = (int)(ROUND_UP((target_x - no_need_to_pass_by_clb -
							xhigh) * inv_length));
			zlog_level(delta_log, ROUTER_V3, "[xhigh < target_x] ");
		} else {
			num_segs_same_dir = 0;
			zlog_level(delta_log, ROUTER_V3, "[same column] ");
		}
	}

	else { /* inode is a CHANY */
		ylow = current.ylow;
		yhigh = current.yhigh;
		xlow = current.xlow;

		/* Count horizontal (orthogonal to inode) segs first. */

		if (xlow > target_x) { /* Coming from a column right of target? */
			*num_segs_ortho_dir_ptr = (int)(
					ROUND_UP((xlow - target_x + 1.) * ortho_inv_length));
			no_need_to_pass_by_clb = 1;
			zlog_level(delta_log, ROUTER_V3, "[xlow > target_x] ");
		} else if (xlow < target_x - 1) { /* Left of and not adjacent to the CLB? */
			*num_segs_ortho_dir_ptr = (int)(ROUND_UP((target_x - xlow) *
					ortho_inv_length));
			no_need_to_pass_by_clb = 1;
			zlog_level(delta_log, ROUTER_V3, "[xlow < target_x-1] ");
		} else { /* In a column that passes by target CLB */
			*num_segs_ortho_dir_ptr = 0;
			no_need_to_pass_by_clb = 0;
			zlog_level(delta_log, ROUTER_V3, "[same column] ");
		}

		/* Now count vertical (same dir. as inode) segs. */

		if (ylow > target_y + no_need_to_pass_by_clb) {
			num_segs_same_dir = (int)(ROUND_UP((ylow - no_need_to_pass_by_clb -
							target_y) * inv_length));
			zlog_level(delta_log, ROUTER_V3, "[ylow > target_y] ");
		} else if (yhigh < target_y - no_need_to_pass_by_clb) {
			num_segs_same_dir = (int)(ROUND_UP((target_y - no_need_to_pass_by_clb -
							yhigh) * inv_length));
			zlog_level(delta_log, ROUTER_V3, "[yhigh < target_y] ");
		} else {
			num_segs_same_dir = 0;
			zlog_level(delta_log, ROUTER_V3, "[same row] ");
		}
	}
	
	zlog_level(delta_log, ROUTER_V3, "\n");

	return (num_segs_same_dir);
}

float get_timing_driven_expected_cost_new(const rr_node_property_t &current, const rr_node_property_t &target,
		float criticality_fac, float R_upstream) {

	/* Determines the expected cost (due to both delay and resouce cost) to reach *
	 * the target node from inode.  It doesn't include the cost of inode --       *
	 * that's already in the "known" path_cost.                                   */

	t_rr_type rr_type;
	int cost_index, ortho_cost_index, num_segs_same_dir, num_segs_ortho_dir;
	float expected_cost, cong_cost, Tdel;

	extern t_rr_indexed_data *rr_indexed_data;

	rr_type = current.type;

	if (rr_type == CHANX || rr_type == CHANY) {
		num_segs_same_dir = get_expected_segs_to_target(current, target,
				&num_segs_ortho_dir);
		cost_index = current.cost_index;
		ortho_cost_index = rr_indexed_data[cost_index].ortho_cost_index;

		cong_cost = num_segs_same_dir * rr_indexed_data[cost_index].saved_base_cost
				+ num_segs_ortho_dir
						* rr_indexed_data[ortho_cost_index].saved_base_cost;
		//cong_cost += rr_indexed_data[IPIN_COST_INDEX].base_cost;
				//+ rr_indexed_data[SINK_COST_INDEX].base_cost;

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

		//Tdel += rr_indexed_data[IPIN_COST_INDEX].T_linear;

		expected_cost = criticality_fac * Tdel
				+ (1. - criticality_fac) * cong_cost;
		return (expected_cost);
	}

	//else if (rr_type == IPIN) { [> Change if you're allowing route-throughs <]
		//return (rr_indexed_data[SINK_COST_INDEX].base_cost);
	//}

	else { /* Change this if you want to investigate route-throughs */
		return (0.);
	}
}

float get_timing_driven_expected_cost(const rr_node_property_t &current, const rr_node_property_t &target, float criticality_fac, float R_upstream, int *same, int *ortho)
{

	/* Determines the expected cost (due to both delay and resouce cost) to reach *
	 * the target node from inode.  It doesn't include the cost of inode --       *
	 * that's already in the "known" path_cost.                                   */

	t_rr_type rr_type;
	int cost_index, ortho_cost_index, num_segs_same_dir, num_segs_ortho_dir;
	float expected_cost, cong_cost, Tdel;

	extern t_rr_indexed_data *rr_indexed_data;

	rr_type = current.type;

	if (rr_type == CHANX || rr_type == CHANY) {
		num_segs_same_dir = get_expected_segs_to_target(current, target,
				&num_segs_ortho_dir);
		cost_index = current.cost_index;
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

		if (same) {
			*same = num_segs_same_dir;
		}
		if (ortho) {
			*ortho = num_segs_ortho_dir;
		}
		return (expected_cost);
	}

	else if (rr_type == IPIN) { /* Change if you're allowing route-throughs */
		return (rr_indexed_data[SINK_COST_INDEX].base_cost);
	}

	else { /* Change this if you want to investigate route-throughs */
		return (0.);
	}
}

//template<typename ShouldExpandFunc>
//void expand_neighbors_fast(const RRGraph &g, const route_state_t &current, const RRNode &target, float criticality_fac, float astar_fac, std::priority_queue<route_state_t> &heap, const ShouldExpandFunc &should_expand, perf_t *perf)
//{
	//for_all_out_edges(g, get_vertex_props(g, current.rr_node), [&heap, &g, &current, &target, &criticality_fac, &astar_fac, &should_expand, &perf] (const RREdge &e) -> void {
			//int neighbor_id = get_target(g, id(e));
			//auto &neighbor = get_vertex_props(g, neighbor_id);

			//char buffer[256];
			//sprintf_rr_node(neighbor_id, buffer);
			//zlog_level(delta_log, ROUTER_V3, "\tNeighbor: %s ", buffer);
			
			//if (!should_expand(neighbor)) {
				//return;
			//}

			//route_state_t item;

			//item.rr_node = id(neighbor);
			//item.prev_edge = id(e);

			//float unbuffered_upstream_R = current.upstream_R;
			//float upstream_R = e.R + neighbor.R;
			//if (!e.buffered) {
				//upstream_R += unbuffered_upstream_R;
			//}
			//item.upstream_R = upstream_R;
			
			//item.delay = current.delay + get_delay(e, neighbor, unbuffered_upstream_R);

			//float known_cost = current.known_cost + get_known_cost(g, e, criticality_fac, unbuffered_upstream_R);
			//item.known_cost = known_cost;

			//float expected_cost = get_timing_driven_expected_cost(get_vertex_props(g, item.rr_node), target, criticality_fac, upstream_R);
			//item.cost = known_cost + astar_fac * expected_cost;

			//heap.push(item);

			//if (perf) {
				//++perf->num_heap_pushes;
			//}
			
			//zlog_level(delta_log, ROUTER_V3, "added with prev_edge: %X upstream_R: %g delay: %g known_cost: %g expected_cost: %g cost: %g\n", item.prev_edge, item.upstream_R, item.delay, item.known_cost, expected_cost, item.cost);
	//});
//}
//
template<typename ShouldExpandFunc>
void expand_neighbors_one_pass(const RRGraph &g, RRNode current, const route_state_t *state, const congestion_t *congestion, const rr_node_property_t &target, float criticality_fac, float astar_fac, std::priority_queue<route_state_t> &heap, const ShouldExpandFunc &should_expand, perf_t *perf)
{
	for (const auto &e : get_out_edges(g, current)) {
		int neighbor = get_target(g, e);
		const auto &neighbor_p = get_vertex_props(g, neighbor);

		char buffer[256];
		sprintf_rr_node(neighbor, buffer);
		zlog_level(delta_log, ROUTER_V3, "\tNeighbor: %s ", buffer);

		if (perf) {
			++perf->num_neighbor_visits;
		}

		if (!should_expand(neighbor)) {
			continue;
		}

		route_state_t item;

		item.rr_node = neighbor;
		item.prev_edge = e;

		const auto &e_p = get_edge_props(g, e);

		extern struct s_switch_inf *switch_inf;
		const struct s_switch_inf *sw = &switch_inf[e_p.switch_index];

		const route_state_t *current_state = &state[current];

		float unbuffered_upstream_R = current_state->upstream_R;
		float upstream_R = sw->R + neighbor_p.R;
		if (!sw->buffered) {
			upstream_R += unbuffered_upstream_R;
		}
		item.upstream_R = upstream_R;

		float congestion_cost = get_congestion_cost(congestion[neighbor], neighbor_p.cost_index);

		float delay = get_delay(e_p, neighbor_p, unbuffered_upstream_R);

		item.delay = current_state->delay + delay;

		float known_cost = criticality_fac * delay + (1 - criticality_fac) * congestion_cost;
		item.known_cost = current_state->known_cost + known_cost;

		float expected_cost = get_timing_driven_expected_cost(neighbor_p, target, criticality_fac, upstream_R);
		item.cost = item.known_cost + astar_fac * expected_cost;

		heap.push(item);

		if (perf) {
			++perf->num_heap_pushes;
		}

		zlog_level(delta_log, ROUTER_V3, " [cost: %g known_cost: %g][occ/cap: %d/%d pres: %g acc: %g][edge_delay: %g edge_R: %g node_R: %g node_C: %g] \n",
				item.cost, item.known_cost, 
				congestion[item.rr_node].occ, neighbor_p.capacity, congestion[item.rr_node].pres_cost, congestion[item.rr_node].acc_cost,
				sw->Tdel, sw->R, neighbor_p.R, neighbor_p.C);
	}
}

template<typename ShouldExpandFunc>
void expand_neighbors_lockless(const RRGraph &g, RRNode current, const route_state_t *state, const congestion_t *congestion, const rr_node_property_t &target, float criticality_fac, float astar_fac, std::priority_queue<route_state_t> &heap, const ShouldExpandFunc &should_expand, perf_t *perf)
{
	for (const auto &e : get_out_edges(g, current)) {
		int neighbor = get_target(g, e);
		const auto &neighbor_p = get_vertex_props(g, neighbor);

		char buffer[256];
		sprintf_rr_node(neighbor, buffer);
		zlog_level(delta_log, ROUTER_V3, "\tNeighbor: %s ", buffer);

		if (perf) {
			++perf->num_neighbor_visits;
		}

		if (!should_expand(neighbor)) {
			continue;
		}

		route_state_t item;

		item.rr_node = neighbor;
		item.prev_edge = e;

		const auto &e_p = get_edge_props(g, e);

		extern struct s_switch_inf *switch_inf;
		const struct s_switch_inf *sw = &switch_inf[e_p.switch_index];

		const route_state_t *current_state = &state[current];

		float unbuffered_upstream_R = current_state->upstream_R;
		float upstream_R = sw->R + neighbor_p.R;
		if (!sw->buffered) {
			upstream_R += unbuffered_upstream_R;
		}
		item.upstream_R = upstream_R;

		float congestion_cost = get_congestion_cost(congestion[neighbor], neighbor_p.cost_index);

		float delay = get_delay(e_p, neighbor_p, unbuffered_upstream_R);

		item.delay = current_state->delay + delay;

		float known_cost = criticality_fac * delay + (1 - criticality_fac) * congestion_cost;
		item.known_cost = current_state->known_cost + known_cost;

		float expected_cost = get_timing_driven_expected_cost(neighbor_p, target, criticality_fac, upstream_R);
		item.cost = item.known_cost + astar_fac * expected_cost;

		heap.push(item);

		if (perf) {
			++perf->num_heap_pushes;
		}

		zlog_level(delta_log, ROUTER_V3, " [cost: %g known_cost: %g][occ/cap: %d/%d pres: %g acc: %g][edge_delay: %g edge_R: %g node_R: %g node_C: %g] \n",
				item.cost, item.known_cost, 
				congestion[item.rr_node].occ, neighbor_p.capacity, congestion[item.rr_node].pres_cost, congestion[item.rr_node].acc_cost,
				sw->Tdel, sw->R, neighbor_p.R, neighbor_p.C);
	}
}

template<typename ShouldExpandFunc>
void expand_neighbors_mpi_rma(const RRGraph &g, RRNode current, const route_state_t *state, congestion_t *congestion, MPI_Win win, const rr_node_property_t &target, float criticality_fac, float astar_fac, float pres_fac, const ShouldExpandFunc &should_expand, const vector<int> &pid, int this_pid, std::priority_queue<route_state_t> &heap, bool &explored_interpartition_neighbor, perf_t *perf)
{
	for (const auto &e : get_out_edges(g, current)) {
		int neighbor = get_target(g, e);
		const auto &neighbor_p = get_vertex_props(g, neighbor);
		const auto &current_p = get_vertex_props(g, current);

		char buffer[256];
		sprintf_rr_node(neighbor, buffer);
		zlog_level(delta_log, ROUTER_V3, "\tNeighbor: %s pid: %d ", buffer, pid[neighbor]);

		if (perf) {
			++perf->num_neighbor_visits;
		}

		if (!should_expand(neighbor)) {
			continue;
		}

		const route_state_t *current_state = &state[current];

		float inter_partition_factor;
		int pc = pid[current];
		int pn = pid[neighbor];
		bool inter_channel = pc != -1 && pn != -1;
		bool opin_to_chan = pc == -1 && pn != -1;
		if (inter_channel) {
			assert((current_p.type == CHANX || current_p.type == CHANY) && (neighbor_p.type == CHANX || neighbor_p.type == CHANY));
		}
		if (opin_to_chan) {
			assert(current_p.type == OPIN && (neighbor_p.type == CHANX || neighbor_p.type == CHANY));
		}
		if ((inter_channel && pc == this_pid && pn != this_pid)
				|| (opin_to_chan && pn != this_pid)) {
			inter_partition_factor = 10;
			explored_interpartition_neighbor = true;
			zlog_level(delta_log, ROUTER_V3, " INTERPART ");
		} else {
			inter_partition_factor = 1;
		}

		route_state_t item;

		item.rr_node = neighbor;
		item.prev_edge = e;

		const auto &e_p = get_edge_props(g, e);

		extern struct s_switch_inf *switch_inf;
		const struct s_switch_inf *sw = &switch_inf[e_p.switch_index];

		float unbuffered_upstream_R = current_state->upstream_R;
		float upstream_R = sw->R + neighbor_p.R;
		if (!sw->buffered) {
			upstream_R += unbuffered_upstream_R;
		}
		item.upstream_R = upstream_R;

		float congestion_cost = get_congestion_cost_mpi_rma(item.rr_node, g, pid, this_pid, congestion, win, neighbor_p.cost_index, pres_fac);
		float delay = get_delay(e_p, neighbor_p, unbuffered_upstream_R);

		item.delay = current_state->delay + delay;

		float known_cost = criticality_fac * inter_partition_factor * delay + (1 - criticality_fac) * congestion_cost;
		item.known_cost = current_state->known_cost + known_cost;

		float expected_cost = get_timing_driven_expected_cost(get_vertex_props(g, item.rr_node), target, criticality_fac, upstream_R);
		item.cost = item.known_cost + astar_fac * expected_cost;

		heap.push(item);

		if (perf) {
			++perf->num_heap_pushes;
		}

		zlog_level(delta_log, ROUTER_V3, " [cost: %g known_cost: %g][occ/cap: %d/%d pres: %g acc: %g][edge_delay: %g edge_R: %g node_R: %g node_C: %g] \n",
				item.cost, item.known_cost, 
				congestion[item.rr_node].occ, neighbor_p.capacity, congestion[item.rr_node].pres_cost, congestion[item.rr_node].acc_cost,
				sw->Tdel, sw->R, neighbor_p.R, neighbor_p.C);
	}
}

using bpqueue = priority_queue<pair<float, int>, vector<pair<float, int>>, std::greater<pair<float, int>>>;

/* high interpartition edge cost */
template<typename ShouldExpandFunc>
void expand_neighbors_with_high_interpartition_cost(const RRGraph &g, RRNode current, const route_state_t *state, const congestion_locked_t *congestion, const rr_node_property_t &target, float criticality_fac, float astar_fac, const ShouldExpandFunc &should_expand, const vector<int> &pid, int this_pid, std::priority_queue<route_state_t> &heap, bool &explored_interpartition_neighbor, perf_t *perf)
{
	for (const auto &e : get_out_edges(g, current)) {
		int neighbor = get_target(g, e);
		const auto &neighbor_p = get_vertex_props(g, neighbor);
		const auto &current_p = get_vertex_props(g, current);

		char buffer[256];
		sprintf_rr_node(neighbor, buffer);
		zlog_level(delta_log, ROUTER_V3, "\tNeighbor: %s pid: %d ", buffer, pid[neighbor]);

		if (perf) {
			++perf->num_neighbor_visits;
		}

		if (!should_expand(neighbor)) {
			continue;
		}

		const route_state_t *current_state = &state[current];

		float inter_partition_factor;
		int pc = pid[current];
		int pn = pid[neighbor];
		bool inter_channel = pc != -1 && pn != -1;
		bool opin_to_chan = pc == -1 && pn != -1;
		if (inter_channel) {
			assert((current_p.type == CHANX || current_p.type == CHANY) && (neighbor_p.type == CHANX || neighbor_p.type == CHANY));
		}
		if (opin_to_chan) {
			assert(current_p.type == OPIN && (neighbor_p.type == CHANX || neighbor_p.type == CHANY));
		}
		if ((inter_channel && pc == this_pid && pn != this_pid)
				|| (opin_to_chan && pn != this_pid)) {
			inter_partition_factor = 10;
			explored_interpartition_neighbor = true;
			zlog_level(delta_log, ROUTER_V3, " INTERPART ");
		} else {
			inter_partition_factor = 1;
		}

		route_state_t item;

		item.rr_node = neighbor;
		item.prev_edge = e;

		const auto &e_p = get_edge_props(g, e);

		extern struct s_switch_inf *switch_inf;
		const struct s_switch_inf *sw = &switch_inf[e_p.switch_index];

		float unbuffered_upstream_R = current_state->upstream_R;
		float upstream_R = sw->R + neighbor_p.R;
		if (!sw->buffered) {
			upstream_R += unbuffered_upstream_R;
		}
		item.upstream_R = upstream_R;

		float congestion_cost = get_congestion_cost(congestion[item.rr_node].cong, neighbor_p.cost_index);
		float delay = get_delay(e_p, neighbor_p, unbuffered_upstream_R);

		item.delay = current_state->delay + delay;

		float known_cost = criticality_fac * inter_partition_factor * delay + (1 - criticality_fac) * congestion_cost;
		item.known_cost = current_state->known_cost + known_cost;

		float expected_cost = get_timing_driven_expected_cost(get_vertex_props(g, item.rr_node), target, criticality_fac, upstream_R);
		item.cost = item.known_cost + astar_fac * expected_cost;

		heap.push(item);

		if (perf) {
			++perf->num_heap_pushes;
		}

		zlog_level(delta_log, ROUTER_V3, " [cost: %g known_cost: %g][occ/cap: %d/%d pres: %g acc: %g][edge_delay: %g edge_R: %g node_R: %g node_C: %g] \n",
				item.cost, item.known_cost, 
				congestion[item.rr_node].cong.occ, neighbor_p.capacity, congestion[item.rr_node].cong.pres_cost, congestion[item.rr_node].cong.acc_cost,
				sw->Tdel, sw->R, neighbor_p.R, neighbor_p.C);
	}
}

/* with tracking of bounding nodes for partition hopping */
template<typename ShouldExpandFunc>
void expand_neighbors_with_partition_hopping(const RRGraph &g, int current, const route_state_t *state, const congestion_locked_t *congestion, const rr_node_property_t &target, float criticality_fac, float astar_fac, std::priority_queue<route_state_t> &heap, const ShouldExpandFunc &should_expand, const vector<int> &pid, int this_pid, bpqueue &boundary_nodes, bool lock, perf_t *perf, lock_perf_t *lock_perf)
{
	bool added_boundary_node = false; 
	for (const auto &e : get_out_edges(g, current)) {
		int neighbor_id = get_target(g, e);
		const auto &neighbor_p = get_vertex_props(g, neighbor_id);

		char buffer[256];
		sprintf_rr_node(neighbor_id, buffer);
		zlog_level(delta_log, ROUTER_V3, "\tNeighbor: %s ", buffer);

		if (perf) {
			++perf->num_neighbor_visits;
		}

		const route_state_t *current_state = &state[current];

		if (pid[neighbor_id] != -1 && pid[neighbor_id] != this_pid && !added_boundary_node) {
			boundary_nodes.push(make_pair(current_state->cost, current));
			added_boundary_node = true;
		}

		if (!should_expand(neighbor_id)) {
			continue;
		}

		route_state_t item;

		item.rr_node = neighbor_id;
		item.prev_edge = e;

		const auto &e_p = get_edge_props(g, e);

		extern struct s_switch_inf *switch_inf;
		const struct s_switch_inf *sw = &switch_inf[e_p.switch_index];

		float unbuffered_upstream_R = current_state->upstream_R;
		float upstream_R = sw->R + neighbor_p.R;
		if (!sw->buffered) {
			upstream_R += unbuffered_upstream_R;
		}
		item.upstream_R = upstream_R;

		/*if (lock) {*/
		/*[>if (lock_perf) {<]*/
		/*[>++lock_perf->num_lock_tries;<]*/
		/*[>}<]*/
		/*[>if (!neighbor_p.lock->try_lock()) {<]*/
		/*[>if (lock_perf) {<]*/
		/*[>++lock_perf->num_lock_waits;<]*/
		/*[>}<]*/
		/*[>neighbor_p.lock->lock();<]*/
		/*[>} <]*/
		/*neighbor_p.lock->lock();*/
		/*}*/
		float congestion_cost = get_congestion_cost(congestion[item.rr_node].cong, neighbor_p.cost_index);
		/*if (lock) {*/
		/*neighbor_p.lock->unlock();*/
		/*}*/
		float delay = get_delay(e_p, neighbor_p, unbuffered_upstream_R);

		item.delay = current_state->delay + delay;

		float known_cost = criticality_fac * delay + (1 - criticality_fac) * congestion_cost;
		item.known_cost = current_state->known_cost + known_cost;

		float expected_cost = get_timing_driven_expected_cost(get_vertex_props(g, item.rr_node), target, criticality_fac, upstream_R);
		item.cost = item.known_cost + astar_fac * expected_cost;

		heap.push(item);

		if (perf) {
			++perf->num_heap_pushes;
		}

		zlog_level(delta_log, ROUTER_V3, " [cost: %g known_cost: %g][occ/cap: %d/%d pres: %g acc: %g][edge_delay: %g edge_R: %g node_R: %g node_C: %g] \n",
				item.cost, item.known_cost, 
				congestion[item.rr_node].cong.occ, neighbor_p.capacity, congestion[item.rr_node].cong.pres_cost, congestion[item.rr_node].cong.acc_cost,
				sw->Tdel, sw->R, neighbor_p.R, neighbor_p.C);
	}
}

RREdge get_previous_edge(int rr_node_id, const route_state_t *state, const route_tree_t &rt, const RRGraph &g)
{
	RREdge previous_edge;

	if (!valid(state[rr_node_id].prev_edge)) {
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
			sprintf_rr_node(get_source(g, state[rr_node_id].prev_edge), s_state);
			sprintf_rr_node(get_vertex_props(rt.graph, get_source(rt.graph, get_vertex_props(rt.graph, rt_node).rt_edge_to_parent)).rr_node, s_rt);
			zlog_warn(delta_log, "Warning: Existing route tree node %s does not have a matching route state. (state.prev_edge: %s rt_node.rr_edge_to_parent: %s) because we have found a shorter path to that node\n", buffer, s_state, s_rt);

			previous_edge = RRGraph::null_edge();
		} else {
			previous_edge = state[rr_node_id].prev_edge;
		}
	} 

	return previous_edge;
}

std::shared_ptr<vector<path_node_t>> get_path(int sink_rr_node_id, const route_state_t *state, const route_tree_t &rt, const RRGraph &g, int vpr_net_id)
{
	int current_rr_node_id = sink_rr_node_id;
	RREdge previous_edge;
	std::shared_ptr<vector<path_node_t>> path = make_shared<vector<path_node_t>>();

	char p[256];
	char c[256];

	while (valid((previous_edge = get_previous_edge(current_rr_node_id, state, rt, g)))) {
		/* parent */
		int parent_rr_node_id = get_source(g, previous_edge);

		path_node_t node;
		node.rr_node_id = current_rr_node_id;
		node.prev_edge = previous_edge;
		//node.update_cost = true;

		path->emplace_back(node);

		/* printing */
		sprintf_rr_node(parent_rr_node_id, p);
		sprintf_rr_node(current_rr_node_id, c);
		zlog_level(delta_log, ROUTER_V2, "Net %d get path: %s\n", vpr_net_id, c);

		current_rr_node_id = parent_rr_node_id;
	}

	path_node_t node;
	node.rr_node_id = current_rr_node_id;
	node.prev_edge = RRGraph::null_edge();
	//node.update_cost = get_vertex_props(g, current_rr_node_id).type == SOURCE;
	//node.update_cost = false;
	sprintf_rr_node(current_rr_node_id, c);
	zlog_level(delta_log, ROUTER_V2, "Net %d get path: %s\n", vpr_net_id, c);

	path->emplace_back(node);

	return path;
}

//vector<sink_t *> fg_route_net_3(const FastRRGraph &g, const vector<int> &pid, int this_pid, int vpr_id, const source_t *source, const vector<sink_t *> &sinks, const route_parameters_t &params, route_state_t *state, congestion_t *congestion, route_tree_t &rt, t_net_timing &net_timing, vector<vector<int>> &unrouted_sinks_boundary_nodes, int num_partitions, bool lock, perf_t *perf, lock_perf_t *lock_perf)
//{
	//std::priority_queue<route_state_t> heap;

	//vector<int> modified;

	//zlog_level(delta_log, ROUTER_V1, "Routing net %d (%lu sinks)\n", vpr_id, sinks.size());

	//vector<sink_t *> sorted_sinks = sinks;
	//std::sort(begin(sorted_sinks), end(sorted_sinks), [] (const sink_t *a, const sink_t *b) -> bool {
			//return a->criticality_fac > b->criticality_fac;
			//});

	//char buffer[256];

	//if (route_tree_empty(rt)) {
		//[> special case <]
		//sprintf_rr_node(source->rr_node, buffer);
		//zlog_level(delta_log, ROUTER_V2, "Empty route tree. Setting root to %s\n", buffer);

		//const auto &source_rr_node = get_vertex_props(g, source->rr_node);
		//RouteTreeNode *root_rt_node = route_tree_add_rr_node(rt, source_rr_node);
		//route_tree_set_node_properties(*root_rt_node, true, -1, source_rr_node.R, 0.5 * source_rr_node.R * source_rr_node.C);
		//route_tree_add_root(rt, source->rr_node);

		//[>update_one_cost_internal(source_rr_node, 1, params.pres_fac);<]
	//} else {
		//if (rt.root_rt_node_id != -1) {
			//const RouteTreeNode &rt_root = get_vertex_props(rt.graph, rt.root_rt_node_id);
			//if (source && rt_root.rr_node != source->rr_node) {
				//char root[256];
				//char source_str[256];
				//sprintf_rr_node(rt_root.rr_node, root);
				//sprintf_rr_node(source->rr_node, source_str);
				//zlog_error(delta_log, "Error: Non empty route tree node has root %s that is different from net source %s\n",
						//root, source_str);
				//assert(false);
			//}
		//}
	//}

	//int isink = 0;
	//int num_routed_sinks = 0;
	//vector<sink_t *> unrouted_sinks;
	//for (const auto &sink : sorted_sinks) {
		//sprintf_rr_node(sink->rr_node, buffer);
		//zlog_level(delta_log, ROUTER_V1, "Net %d Sink %d: %s criticality: %g BB: %d-%d %d-%d Prev BB: %d-%d %d-%d\n", vpr_id, sink->id, buffer, sink->criticality_fac, sink->current_bounding_box.xmin, sink->current_bounding_box.xmax, sink->current_bounding_box.ymin, sink->current_bounding_box.ymax, sink->previous_bounding_box.xmin, sink->previous_bounding_box.xmax, sink->previous_bounding_box.ymin, sink->previous_bounding_box.ymax);

		//const auto &sink_rr_node = get_vertex_props(g, sink->rr_node);

		//route_tree_multi_root_add_to_heap(rt, g, sink_rr_node, sink->criticality_fac, params.astar_fac, sink->current_bounding_box, heap, perf);

		//priority_queue<pair<float, int>, vector<pair<float, int>>, std::greater<pair<float, int>>> boundary_nodes;

		//bool found_sink = false;
		//while (!heap.empty() && !found_sink) {
			//auto item = heap.top();
			//heap.pop();

			//assert(pid[item.rr_node] == -1 || pid[item.rr_node] == this_pid);

			//if (perf) {
				//++perf->num_heap_pops;
			//}

			//sprintf_rr_node(item.rr_node, buffer);
			//const auto &v = get_vertex_props(g, item.rr_node);
			//zlog_level(delta_log, ROUTER_V3, "Current: %s occ/cap: %d/%d prev: %d new_cost: %g new_known: %g new_delay: %g old_cost: %g old_known: %g old_delay: %g\n", buffer, congestion[item.rr_node].occ, v.capacity, item.prev_edge != -1 ? get_source(g, item.prev_edge) : -1, item.cost, item.known_cost, item.delay, state[item.rr_node].cost, state[item.rr_node].known_cost, state[item.rr_node].delay);

			//if (item.rr_node == sink->rr_node) {
				//state[item.rr_node] = item;
				//modified.push_back(item.rr_node);
				//found_sink = true;
			//} else if (item.known_cost < state[item.rr_node].known_cost && item.cost < state[item.rr_node].cost) {
				//[>if (route_tree_get_rt_node(rt, item.rr_node)) {<]
					//[>sprintf_rr_node(item.rr_node, buffer);<]
					//[>zlog_warn(delta_log, "Warning: Found a lower cost path to existing route tree node %s\n", buffer);<]
					//[>assert(false);<]
				//[>}<]

				//state[item.rr_node] = item;
				//modified.push_back(item.rr_node);

				//expand_neighbors(g, v, state, congestion, sink_rr_node, sink->criticality_fac, params.astar_fac, heap, [&sink, &sink_rr_node, &v, &pid, &boundary_nodes, &item, &this_pid] (const RRNode &n) -> bool {

					//[>if (trace_has_node(prev_trace, id(n))) {<]
						//[>zlog_level(delta_log, ROUTER_V3, " existing node route tree ");<]
					//[>}<]
					//const auto &prop = n.properties;

					//if (prop.xhigh < sink->current_bounding_box.xmin
							//|| prop.xlow > sink->current_bounding_box.xmax
							//|| prop.yhigh < sink->current_bounding_box.ymin
							//|| prop.ylow > sink->current_bounding_box.ymax) {
					//zlog_level(delta_log, ROUTER_V3, "outside of bounding box\n");
					//return false;
					//}

					//if (prop.type == IPIN
							//&& (prop.xhigh != sink_rr_node.xhigh 
								//|| prop.yhigh != sink_rr_node.yhigh)) {
					//zlog_level(delta_log, ROUTER_V3, "not target IPIN\n");
					//return false;
					//}

					//[>if (net.sinks.size() >= HIGH_FANOUT_NET_LIM) {<]
					//[>if (prop.xhigh < target_x - highfanout_rlim<]
							//[>|| prop.xlow > target_x + highfanout_rlim<]
							//[>|| prop.yhigh < target_y - highfanout_rlim<]
							//[>|| prop.ylow > target_y + highfanout_rlim) {<]
						//[>return false;<]
					//[>}<]
					//[>}<]
					//if (pid[id(n)] != -1 && pid[id(n)] != this_pid) {
						//zlog_level(delta_log, ROUTER_V3, " boundary node %d with cost %g\n", id(v), item.cost);
						//boundary_nodes.push(make_pair(item.cost, id(v)));
						//return false;
					//}

					//return true;
				//}, lock, perf, lock_perf);
			//} else {
				//zlog_level(delta_log, ROUTER_V3, "\tNot expanding neighbor because known_cost %g > %g or cost %g > %g\n", item.known_cost, state[item.rr_node].known_cost, item.cost, state[item.rr_node].cost);
			//}
		//}

		//if (!found_sink) {
			//assert(heap.empty());
			//zlog_error(delta_log, "Error: Failed to find sink %d\n", sink->rr_node);
			//vector<int> lowest_cost_boundary_nodes;
			//int num_visited_partitions = 0;
			//vector<bool> visited_partitions(num_partitions, false);
			//set<int> added_boundary_nodes;
			//while (!boundary_nodes.empty() && num_visited_partitions < num_partitions-1) {
				//auto bn = boundary_nodes.top(); boundary_nodes.pop();
				//if (added_boundary_nodes.find(bn.second) == added_boundary_nodes.end()) {
					//lowest_cost_boundary_nodes.push_back(bn.second);
					//added_boundary_nodes.insert(bn.second);
					//for (const auto &e : get_out_edges(g, bn.second)) {
						//int to = get_target(g, e);
						//int to_pid = pid[to];
						//if (to_pid != -1 && to_pid != this_pid && !visited_partitions[to_pid]) {
							//visited_partitions[to_pid] = true;
							//++num_visited_partitions;
						//}
					//}
				//}
			//}
			//if (num_visited_partitions < num_partitions-1) {
				//zlog_level(delta_log, ROUTER_V3, "Net %d sink %d does't have boundaries nodes that covers all partitions\n");
			//}
			//unrouted_sinks_boundary_nodes.emplace_back(lowest_cost_boundary_nodes);
			//unrouted_sinks.push_back(sink);
		//} else {
			//const auto &path = get_path(sink->rr_node, state, rt, g, vpr_id);

			//route_tree_add_path(rt, path, g, state);

			//vector<int> added_rr_nodes;
			//for (const auto &node : path) {
				//if (node.update_cost) {
					//added_rr_nodes.push_back(node.rr_node_id);
				//}
			//}

			//update_one_cost(g, congestion, added_rr_nodes.begin(), added_rr_nodes.end(), 1, params.pres_fac, lock, lock_perf);

			//net_timing.delay[sink->id+1] = route_tree_get_rt_node(rt, sink->rr_node)->delay;

			//++num_routed_sinks;
		//}

		//for (const auto &m : modified)  {
			//state[m].known_cost = std::numeric_limits<float>::max();
			//state[m].cost = std::numeric_limits<float>::max();
		//}

		//heap = std::priority_queue<route_state_t>();
	//}

	//zlog_level(delta_log, ROUTER_V1, "\n");

	//return unrouted_sinks;

	//[>delete [] state;<]

	//[>check_route_tree(rt, net, g);<]
//}
//
void route_net_one_pass(const RRGraph &g, int vpr_id, const source_t *source, const vector<sink_t *> &sinks, const route_parameters_t &params, float pres_fac, route_state_t *state, congestion_t *congestion, route_tree_t &rt, t_net_timing &net_timing, vector<sink_t *> &routed_sinks, vector<sink_t *> &unrouted_sinks, perf_t *perf)
{
	std::priority_queue<route_state_t> heap;

	vector<int> modified;
	vector<bool> visited(num_vertices(g), false);

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

		RouteTreeNode root_rt_node = route_tree_add_rr_node(rt, source->rr_node);
		const auto &source_rr_node_p = get_vertex_props(g, source->rr_node);
		route_tree_set_node_properties(rt, root_rt_node, true, source_rr_node_p.R, 0.5 * source_rr_node_p.R * source_rr_node_p.C);
		route_tree_add_root(rt, source->rr_node);

		/*update_one_cost_internal(source_rr_node_p, 1, params.pres_fac);*/
	} else {
		if (rt.root_rt_node_id != -1) {
			const auto &rt_root_p = get_vertex_props(rt.graph, rt.root_rt_node_id);
			if (source && rt_root_p.rr_node != source->rr_node) {
				char root[256];
				char source_str[256];
				sprintf_rr_node(rt_root_p.rr_node, root);
				sprintf_rr_node(source->rr_node, source_str);
				zlog_error(delta_log, "Error: Non empty route tree node has root %s that is different from net source %s\n",
						root, source_str);
				assert(false);
			}
		}
	}

	int isink = 0;
	for (const auto &sink : sorted_sinks) {
		sprintf_rr_node(sink->rr_node, buffer);
		//zlog_level(delta_log, ROUTER_V1, "Net %d current sink %d: %s criticality: %g BB: %d-%d %d-%d Prev BB: %d-%d %d-%d\n", vpr_id, sink->id, buffer, sink->criticality_fac, sink->current_bounding_box.xmin, sink->current_bounding_box.xmax, sink->current_bounding_box.ymin, sink->current_bounding_box.ymax, sink->previous_bounding_box.xmin, sink->previous_bounding_box.xmax, sink->previous_bounding_box.ymin, sink->previous_bounding_box.ymax);

		const auto &sink_rr_node = get_vertex_props(g, sink->rr_node);

		route_tree_multi_root_add_to_heap(rt, g, sink->rr_node, sink->criticality_fac, params.astar_fac, sink->current_bounding_box, heap, perf);

		bool found_sink = false;
		while (!heap.empty() && !found_sink) {
			auto item = heap.top();
			heap.pop();

			sprintf_rr_node(item.rr_node, buffer);
			const auto &v = get_vertex_props(g, item.rr_node);
			zlog_level(delta_log, ROUTER_V3, "Current: %s occ/cap: %d/%d prev: %d new_cost: %g new_known: %g new_delay: %g old_cost: %g old_known: %g old_delay: %g\n", buffer, congestion[item.rr_node].occ, v.capacity, valid(item.prev_edge) ? get_source(g, item.prev_edge) : -1, item.cost, item.known_cost, item.delay, state[item.rr_node].cost, state[item.rr_node].known_cost, state[item.rr_node].delay);

			if (route_tree_get_rt_node(rt, item.rr_node) == RouteTree::null_vertex()) {
				assert(!visited[item.rr_node]);
			}
			visited[item.rr_node] = true;

			if (perf) {
				++perf->num_heap_pops;
			}
			//assert(pid[item.rr_node] == -1 || pid[item.rr_node] == this_pid);

			if (item.rr_node == sink->rr_node) {
				state[item.rr_node] = item;
				modified.push_back(item.rr_node);
				found_sink = true;
			} else {
				if (route_tree_get_rt_node(rt, item.rr_node) == RouteTree::null_vertex()) {
					assert(item.known_cost < state[item.rr_node].known_cost);
				}
				state[item.rr_node] = item;
				modified.push_back(item.rr_node);

				expand_neighbors_one_pass(g, item.rr_node, state, congestion, sink_rr_node, sink->criticality_fac, params.astar_fac, heap, [&g, &sink, &sink_rr_node, &v, &item, &visited] (const RRNode &n) -> bool {

					/*if (trace_has_node(prev_trace, id(n))) {*/
						/*zlog_level(delta_log, ROUTER_V3, " existing node route tree ");*/
					/*}*/

					if (visited[n]) {
					return false;
					}
					
					const auto &prop = get_vertex_props(g, n);

					if (prop.xhigh < sink->current_bounding_box.xmin
							|| prop.xlow > sink->current_bounding_box.xmax
							|| prop.yhigh < sink->current_bounding_box.ymin
							|| prop.ylow > sink->current_bounding_box.ymax) {
					zlog_level(delta_log, ROUTER_V3, "outside of bounding box\n");
					return false;
					}

					if (prop.type == IPIN
							&& (prop.xhigh != sink_rr_node.xhigh 
								|| prop.yhigh != sink_rr_node.yhigh)) {
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
					//if (pid[n] != -1 && pid[n] != this_pid) {
						//zlog_level(delta_log, ROUTER_V3, " pid %d not in current partition %d\n", pid[n], this_pid);
						//return false;
					//}

					return true;
				}, perf);
			} 
			//else {
				//zlog_level(delta_log, ROUTER_V3, "\tNot expanding neighbor because known_cost %g > %g or cost %g > %g\n", item.known_cost, state[item.rr_node].known_cost, item.cost, state[item.rr_node].cost);
			//}
		}

		if (!found_sink) {
			assert(heap.empty());
			zlog_error(delta_log, "Error: Failed to find sink %d\n", sink->rr_node);

			assert(find(begin(unrouted_sinks), end(unrouted_sinks), sink) == end(unrouted_sinks));
			unrouted_sinks.push_back(sink);

			if (sink->id != -1) {
				net_timing.delay[sink->id+1] = 0;
			}
		} else {
			assert(find(begin(routed_sinks), end(routed_sinks), sink) == end(routed_sinks));

			routed_sinks.push_back(sink);

			const auto &path = get_path(sink->rr_node, state, rt, g, vpr_id);

			route_tree_add_path(rt, path, g, state);

			vector<RRNode> added_nodes;
			for (const auto &n : *path) {
				if (valid(n.prev_edge) || get_vertex_props(g, n.rr_node_id).type == SOURCE) {
					added_nodes.push_back(n.rr_node_id);
				} 
			}

			update_one_cost(g, congestion, added_nodes.begin(), added_nodes.end(), 1, pres_fac);

			if (sink->id != -1) {
				net_timing.delay[sink->id+1] = get_vertex_props(rt.graph, route_tree_get_rt_node(rt, sink->rr_node)).delay;
			}

			zlog_level(delta_log, ROUTER_V2, "Routed sink %d\n", sink->rr_node);
		}

		for (const auto &m : modified)  {
			state[m].known_cost = std::numeric_limits<float>::max();
			state[m].cost = std::numeric_limits<float>::max();
			visited[m] = false;
		}

		modified.clear();

		heap = std::priority_queue<route_state_t>();
	}
	//zlog_level(delta_log, ROUTER_V1, "\n");

	//return unrouted_sinks;

	/*delete [] state;*/

	/*check_route_tree(rt, net, g);*/
}

void route_net_with_partitioned_fine_grain_lock(const RRGraph &g, const vector<int> &pid, int this_pid, int vpr_id, const source_t *source, const vector<sink_t *> &sinks, const route_parameters_t &params, float pres_fac, route_state_t *state, congestion_locked_t *congestion, route_tree_t &rt, t_net_timing &net_timing, vector<sink_t *> &routed_sinks, vector<sink_t *> &unrouted_sinks, bool lock, perf_t *perf, lock_perf_t *lock_perf)
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

		const auto &source_rr_node = get_vertex_props(g, source->rr_node);
		RouteTreeNode root_rt_node = route_tree_add_rr_node(rt, source->rr_node);
		route_tree_set_node_properties(rt, root_rt_node, true, source_rr_node.R, 0.5 * source_rr_node.R * source_rr_node.C);
		route_tree_add_root(rt, source->rr_node);

		/*update_one_cost_internal(source_rr_node, 1, params.pres_fac);*/
	} else {
		assert(rt.root_rt_nodes.size() == 1);
		const auto &rt_root_p = get_vertex_props(rt.graph, rt.root_rt_nodes[0]);
		if (source && rt_root_p.rr_node != source->rr_node) {
			char root[256];
			char source_str[256];
			sprintf_rr_node(rt_root_p.rr_node, root);
			sprintf_rr_node(source->rr_node, source_str);
			zlog_error(delta_log, "Error: Non empty route tree node has root %s that is different from net source %s\n",
					root, source_str);
			assert(false);
		}
	}

	int isink = 0;
	int num_routed_sinks = 0;
	for (const auto &sink : sorted_sinks) {
		sprintf_rr_node(sink->rr_node, buffer);
		//zlog_level(delta_log, ROUTER_V1, "Net %d Current sink %d: %s criticality: %g BB: %d-%d %d-%d Prev BB: %d-%d %d-%d\n", vpr_id, sink->id, buffer, sink->criticality_fac, sink->current_bounding_box.xmin, sink->current_bounding_box.xmax, sink->current_bounding_box.ymin, sink->current_bounding_box.ymax, sink->previous_bounding_box.xmin, sink->previous_bounding_box.xmax, sink->previous_bounding_box.ymin, sink->previous_bounding_box.ymax);

		const auto &sink_rr_node = get_vertex_props(g, sink->rr_node);

		route_tree_multi_root_add_to_heap(rt, g, sink->rr_node, sink->criticality_fac, params.astar_fac, sink->current_bounding_box, heap, perf);

		bool found_sink = false;
		while (!heap.empty() && !found_sink) {
			auto item = heap.top();
			heap.pop();

			if (perf) {
				++perf->num_heap_pops;
			}

			sprintf_rr_node(item.rr_node, buffer);
			const auto &v = get_vertex_props(g, item.rr_node);
			zlog_level(delta_log, ROUTER_V3, "Current: %s pid: %d occ/cap: %d/%d prev: %d new_cost: %g new_known: %g new_delay: %g old_cost: %g old_known: %g old_delay: %g\n", buffer, pid[item.rr_node], congestion[item.rr_node].cong.occ, v.capacity, valid(item.prev_edge) ? get_source(g, item.prev_edge) : -1, item.cost, item.known_cost, item.delay, state[item.rr_node].cost, state[item.rr_node].known_cost, state[item.rr_node].delay);

			assert(pid[item.rr_node] == -1 || pid[item.rr_node] == this_pid);

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

				expand_neighbors_with_fine_grain_lock(g, item.rr_node, state, congestion, sink_rr_node, sink->criticality_fac, params.astar_fac, heap, [&g, &sink, &sink_rr_node, &pid, &this_pid] (const RRNode &n) -> bool {

					/*if (trace_has_node(prev_trace, id(v))) {*/
						/*zlog_level(delta_log, ROUTER_V3, " existing node route tree ");*/
					/*}*/
					const auto &prop = get_vertex_props(g, n);

					if (prop.xhigh < sink->current_bounding_box.xmin
							|| prop.xlow > sink->current_bounding_box.xmax
							|| prop.yhigh < sink->current_bounding_box.ymin
							|| prop.ylow > sink->current_bounding_box.ymax) {
					zlog_level(delta_log, ROUTER_V3, "outside of bounding box\n");
					return false;
					}

					if (prop.type == IPIN
							&& (prop.xhigh != sink_rr_node.xhigh 
								|| prop.yhigh != sink_rr_node.yhigh)) {
					zlog_level(delta_log, ROUTER_V3, "not target IPIN\n");
					return false;
					}

					if (pid[n] != -1 && pid[n] != this_pid) {
						zlog_level(delta_log, ROUTER_V3, " pid %d not in current partition %d\n", pid[n], this_pid);
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
			const auto &path = get_path(sink->rr_node, state, rt, g, vpr_id);

			route_tree_add_path(rt, path, g, state);

			vector<RRNode> added_nodes;
			for (const auto &n : *path) {
				const auto &prop = get_vertex_props(g, n.rr_node_id);

				if (valid(n.prev_edge) || prop.type == SOURCE) {
					added_nodes.push_back(n.rr_node_id);
				} 
			}

			update_one_cost(g, congestion, added_nodes.begin(), added_nodes.end(), 1, pres_fac, lock, lock_perf);

			net_timing.delay[sink->id+1] = get_vertex_props(rt.graph, route_tree_get_rt_node(rt, sink->rr_node)).delay;

			routed_sinks.push_back(sink);

			++num_routed_sinks;
		}

		for (const auto &m : modified)  {
			state[m].known_cost = std::numeric_limits<float>::max();
			state[m].cost = std::numeric_limits<float>::max();
		}

		heap = std::priority_queue<route_state_t>();
	}

	zlog_level(delta_log, ROUTER_V1, "\n");

	/*delete [] state;*/

	/*check_route_tree(rt, net, g);*/
}

void broadcast_costs(const vector<RRNode>::const_iterator &rr_nodes_begin, const vector<RRNode>::const_iterator &rr_nodes_end, int delta, int this_pid, int num_procs, MPI_Comm comm, vector<ongoing_transaction_t> &transactions)
{
	ongoing_transaction_t trans;

	trans.data = make_shared<vector<node_update_t>>();

	for (auto iter = rr_nodes_begin; iter != rr_nodes_end; ++iter) {
		node_update_t d;
		d.rr_node = *iter;
		d.delta = delta;

		trans.data->push_back(d);
	}

	assert(!trans.data->empty());

	for (int i = 0; i < num_procs; ++i) {
		if (i != this_pid) {
			zlog_level(delta_log, ROUTER_V3, "MPI update from %d to %d\n", this_pid, i);
			assert(MPI_Isend(trans.data->data(), trans.data->size()*2, MPI_INT, i, 0, comm, &trans.req) == MPI_SUCCESS);
			//assert(MPI_ISend(trans.data->data(), trans.data->size()*2, MPI_INT, i, 0, comm) == MPI_SUCCESS);
		}
	}

	transactions.push_back(trans);
}

void route_net_lockless(const RRGraph &g, const vector<int> &pid, int this_pid, int vpr_id, const source_t *source, const vector<sink_t *> &sinks, const route_parameters_t &params, float pres_fac, route_state_t *state, congestion_t *congestion, route_tree_t &rt, t_net_timing &net_timing, vector<sink_t *> &routed_sinks, vector<sink_t *> &unrouted_sinks, perf_t *perf)
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

		RouteTreeNode root_rt_node = route_tree_add_rr_node(rt, source->rr_node);
		const auto &source_rr_node_p = get_vertex_props(g, source->rr_node);
		route_tree_set_node_properties(rt, root_rt_node, true, source_rr_node_p.R, 0.5 * source_rr_node_p.R * source_rr_node_p.C);
		route_tree_add_root(rt, source->rr_node);

		/*update_one_cost_internal(source_rr_node_p, 1, params.pres_fac);*/
	} else {
		if (rt.root_rt_node_id != -1) {
			const auto &rt_root_p = get_vertex_props(rt.graph, rt.root_rt_node_id);
			if (source && rt_root_p.rr_node != source->rr_node) {
				char root[256];
				char source_str[256];
				sprintf_rr_node(rt_root_p.rr_node, root);
				sprintf_rr_node(source->rr_node, source_str);
				zlog_error(delta_log, "Error: Non empty route tree node has root %s that is different from net source %s\n",
						root, source_str);
				assert(false);
			}
		}
	}

	int isink = 0;
	for (const auto &sink : sorted_sinks) {
		sprintf_rr_node(sink->rr_node, buffer);
		//zlog_level(delta_log, ROUTER_V1, "Net %d current sink %d: %s criticality: %g BB: %d-%d %d-%d Prev BB: %d-%d %d-%d\n", vpr_id, sink->id, buffer, sink->criticality_fac, sink->current_bounding_box.xmin, sink->current_bounding_box.xmax, sink->current_bounding_box.ymin, sink->current_bounding_box.ymax, sink->previous_bounding_box.xmin, sink->previous_bounding_box.xmax, sink->previous_bounding_box.ymin, sink->previous_bounding_box.ymax);

		const auto &sink_rr_node = get_vertex_props(g, sink->rr_node);

		route_tree_multi_root_add_to_heap(rt, g, sink->rr_node, sink->criticality_fac, params.astar_fac, sink->current_bounding_box, heap, perf);

		bool found_sink = false;
		while (!heap.empty() && !found_sink) {
			auto item = heap.top();
			heap.pop();

			if (perf) {
				++perf->num_heap_pops;
			}

			sprintf_rr_node(item.rr_node, buffer);
			const auto &v = get_vertex_props(g, item.rr_node);
			zlog_level(delta_log, ROUTER_V3, "Current: %s pid: %d occ/cap: %d/%d prev: %d new_cost: %g new_known: %g new_delay: %g old_cost: %g old_known: %g old_delay: %g\n", buffer, pid[item.rr_node], congestion[item.rr_node].occ, v.capacity, valid(item.prev_edge) ? get_source(g, item.prev_edge) : -1, item.cost, item.known_cost, item.delay, state[item.rr_node].cost, state[item.rr_node].known_cost, state[item.rr_node].delay);

			assert(pid[item.rr_node] == -1 || pid[item.rr_node] == this_pid);

			if (item.rr_node == sink->rr_node) {
				state[item.rr_node] = item;
				modified.push_back(item.rr_node);
				found_sink = true;
			} else if (item.cost < state[item.rr_node].cost) {
				assert(item.known_cost < state[item.rr_node].known_cost);
				state[item.rr_node] = item;
				modified.push_back(item.rr_node);

				expand_neighbors_lockless(g, item.rr_node, state, congestion, sink_rr_node, sink->criticality_fac, params.astar_fac, heap, [&g, &sink, &sink_rr_node, &v, &pid, &item, &this_pid] (const RRNode &n) -> bool {

					/*if (trace_has_node(prev_trace, id(n))) {*/
						/*zlog_level(delta_log, ROUTER_V3, " existing node route tree ");*/
					/*}*/
					const auto &prop = get_vertex_props(g, n);

					if (prop.xhigh < sink->current_bounding_box.xmin
							|| prop.xlow > sink->current_bounding_box.xmax
							|| prop.yhigh < sink->current_bounding_box.ymin
							|| prop.ylow > sink->current_bounding_box.ymax) {
					zlog_level(delta_log, ROUTER_V3, "outside of bounding box\n");
					return false;
					}

					if (prop.type == IPIN
							&& (prop.xhigh != sink_rr_node.xhigh 
								|| prop.yhigh != sink_rr_node.yhigh)) {
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
					if (pid[n] != -1 && pid[n] != this_pid) {
						zlog_level(delta_log, ROUTER_V3, " pid %d not in current partition %d\n", pid[n], this_pid);
						return false;
					}

					return true;
				}, perf);
			} else {
				zlog_level(delta_log, ROUTER_V3, "\tNot expanding neighbor because known_cost %g > %g or cost %g > %g\n", item.known_cost, state[item.rr_node].known_cost, item.cost, state[item.rr_node].cost);
			}
		}

		if (!found_sink) {
			assert(heap.empty());
			zlog_error(delta_log, "Error: Failed to find sink %d\n", sink->rr_node);

			assert(find(begin(unrouted_sinks), end(unrouted_sinks), sink) == end(unrouted_sinks));
			unrouted_sinks.push_back(sink);

			if (sink->id != -1) {
				net_timing.delay[sink->id+1] = 0;
			}
		} else {
			assert(find(begin(routed_sinks), end(routed_sinks), sink) == end(routed_sinks));

			routed_sinks.push_back(sink);

			const auto &path = get_path(sink->rr_node, state, rt, g, vpr_id);

			route_tree_add_path(rt, path, g, state);

			vector<RRNode> added_nodes;
			for (const auto &n : *path) {
				if (valid(n.prev_edge) || get_vertex_props(g, n.rr_node_id).type == SOURCE) {
					added_nodes.push_back(n.rr_node_id);
				} 
			}

			update_one_cost(g, congestion, added_nodes.begin(), added_nodes.end(), 1, pres_fac);

			if (sink->id != -1) {
				net_timing.delay[sink->id+1] = get_vertex_props(rt.graph, route_tree_get_rt_node(rt, sink->rr_node)).delay;
			}

			zlog_level(delta_log, ROUTER_V2, "Routed sink %d\n", sink->rr_node);
		}

		for (const auto &m : modified)  {
			state[m].known_cost = std::numeric_limits<float>::max();
			state[m].cost = std::numeric_limits<float>::max();
		}

		modified.clear();

		heap = std::priority_queue<route_state_t>();
	}
	//zlog_level(delta_log, ROUTER_V1, "\n");

	//return unrouted_sinks;

	/*delete [] state;*/

	/*check_route_tree(rt, net, g);*/
}

void route_net_mpi_rma(const RRGraph &g, const vector<int> &pid, int this_pid, MPI_Win win, int vpr_id, const source_t *source, const vector<sink_t *> &sinks, const route_parameters_t &params, float pres_fac, route_state_t *state, congestion_t *congestion, route_tree_t &rt, t_net_timing &net_timing, vector<interpartition_sink_t> &interpartition_sinks, perf_t *perf)
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

		RouteTreeNode root_rt_node = route_tree_add_rr_node(rt, source->rr_node);
		const auto &source_rr_node_p = get_vertex_props(g, source->rr_node);
		route_tree_set_node_properties(rt, root_rt_node, true, source_rr_node_p.R, 0.5 * source_rr_node_p.R * source_rr_node_p.C);
		route_tree_add_root(rt, source->rr_node);

		/*update_one_cost_internal(source_rr_node_p, 1, params.pres_fac);*/
	} else {
		if (rt.root_rt_node_id != -1) {
			const auto &rt_root_p = get_vertex_props(rt.graph, rt.root_rt_node_id);
			if (source && rt_root_p.rr_node != source->rr_node) {
				char root[256];
				char source_str[256];
				sprintf_rr_node(rt_root_p.rr_node, root);
				sprintf_rr_node(source->rr_node, source_str);
				zlog_error(delta_log, "Error: Non empty route tree node has root %s that is different from net source %s\n",
						root, source_str);
				assert(false);
			}
		}
	}
	
	vector<bool> touched(num_vertices(g), false);

	int isink = 0;
	int num_routed_sinks = 0;
	bpqueue boundary_node_heap;
	for (const auto &sink : sorted_sinks) {
		sprintf_rr_node(sink->rr_node, buffer);
		//zlog_level(delta_log, ROUTER_V1, "Net %d Sink %d: %s criticality: %g BB: %d-%d %d-%d Prev BB: %d-%d %d-%d\n", vpr_id, sink->id, buffer, sink->criticality_fac, sink->current_bounding_box.xmin, sink->current_bounding_box.xmax, sink->current_bounding_box.ymin, sink->current_bounding_box.ymax, sink->previous_bounding_box.xmin, sink->previous_bounding_box.xmax, sink->previous_bounding_box.ymin, sink->previous_bounding_box.ymax);

		const auto &sink_rr_node = get_vertex_props(g, sink->rr_node);

		route_tree_multi_root_add_to_heap(rt, g, sink->rr_node, sink->criticality_fac, params.astar_fac, sink->current_bounding_box, heap, perf);

		bool expanded_interpartition_node = false;
		bool explored_interpartition_neighbor = false;

		bool found_sink = false;
		while (!heap.empty() && !found_sink) {
			auto item = heap.top();
			heap.pop();

			if (perf) {
				++perf->num_heap_pops;
			}

			sprintf_rr_node(item.rr_node, buffer);
			const auto &v = get_vertex_props(g, item.rr_node);
			zlog_level(delta_log, ROUTER_V3, "Current: %s %s pid: %d occ/cap: %d/%d prev: %d new_cost: %g new_known: %g new_delay: %g old_cost: %g old_known: %g old_delay: %g\n", buffer, (pid[item.rr_node] != -1 && pid[item.rr_node] != this_pid) ? "OTHER": "", pid[item.rr_node], congestion[item.rr_node].occ, v.capacity, valid(item.prev_edge) ? get_source(g, item.prev_edge) : -1, item.cost, item.known_cost, item.delay, state[item.rr_node].cost, state[item.rr_node].known_cost, state[item.rr_node].delay);

			if (pid[item.rr_node] != -1 && pid[item.rr_node] != this_pid) {
				expanded_interpartition_node = true;
			}

			if (item.rr_node == sink->rr_node) {
				state[item.rr_node] = item;
				modified.push_back(item.rr_node);
				found_sink = true;
			} else if (item.cost < state[item.rr_node].cost) {
				//if (touched[item.rr_node]) {
					//RRGraph temp;
					//add_vertex(temp, num_vertices(g));
					//for (const auto &m : modified) {
						//if (valid(state[m].prev_edge)) {
							//add_edge(temp, get_source(g, state[m].prev_edge), get_target(g, state[m].prev_edge));
						//}
					//}
					//auto ft = make_filtered_graph(temp, [&modified] (RRNode v) -> bool {
							//return find(begin(modified), end(modified), v) != end(modified);
							//},
							//[] (const RREdge &e) -> bool {
							//return true;
							//});

					//write_graph(ft, "weird.dot", [&state, &item] (RRNode n) -> string {
						//char buffer[256];
						//if (n == item.rr_node) {
						//sprintf(buffer, "label=\"%d %g %g %g\" style=filled fillcolor=red", n, state[n].cost, state[n].known_cost, state[n].cost-state[n].known_cost);
						//} else {
						//sprintf(buffer, "label=\"%d %g %g %g\"", n, state[n].cost, state[n].known_cost, state[n].cost-state[n].known_cost);
						//}
						//return string(buffer);
						//}, [] (const RREdge &e) -> string {
						//return string();
						//});
					//assert(false);
				//}
				assert(item.known_cost < state[item.rr_node].known_cost);
				//touched[item.rr_node] = true;
				state[item.rr_node] = item;
				modified.push_back(item.rr_node);

				expand_neighbors_mpi_rma(g, item.rr_node, state, congestion, win, sink_rr_node, sink->criticality_fac, params.astar_fac, pres_fac,
						[&g, &sink, &sink_rr_node, &v, &pid, &boundary_node_heap, &item, &this_pid, &congestion] (const RRNode &n) -> bool {
					const auto &prop = get_vertex_props(g, n);

					//if (congestion[n].occ >= prop.capacity) {
					//zlog_level(delta_log, ROUTER_V3, "congested\n");
					//return false;
					//}
					/*if (trace_has_node(prev_trace, id(n))) {*/
						/*zlog_level(delta_log, ROUTER_V3, " existing node route tree ");*/
					/*}*/

					if (prop.xhigh < sink->current_bounding_box.xmin
							|| prop.xlow > sink->current_bounding_box.xmax
							|| prop.yhigh < sink->current_bounding_box.ymin
							|| prop.ylow > sink->current_bounding_box.ymax) {
					zlog_level(delta_log, ROUTER_V3, "outside of bounding box\n");
					return false;
					}

					if (prop.type == IPIN
							&& (prop.xhigh != sink_rr_node.xhigh 
								|| prop.yhigh != sink_rr_node.yhigh)) {
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
				}, pid, this_pid, heap, explored_interpartition_neighbor, perf);
			} else {
				zlog_level(delta_log, ROUTER_V3, "\tNot expanding neighbor because known_cost %g > %g or cost %g > %g\n", item.known_cost, state[item.rr_node].known_cost, item.cost, state[item.rr_node].cost);
			}
		}

		if (!found_sink) {
			assert(heap.empty());
			zlog_error(delta_log, "Error: Failed to find sink %d\n", sink->rr_node);
			assert(false);
		} else {
			const auto &path = get_path(sink->rr_node, state, rt, g, vpr_id);

			route_tree_add_path(rt, path, g, state);

			vector<RRNode> added_nodes;
			for (const auto &n : *path) {
				if (valid(n.prev_edge) || get_vertex_props(g, n.rr_node_id).type == SOURCE) {
					added_nodes.push_back(n.rr_node_id);
				} 
			}

			update_one_cost_mpi_rma(added_nodes.begin(), added_nodes.end(), g, pid, this_pid, congestion, win, 1, pres_fac);

			if (!expanded_interpartition_node) {
			} else {
				interpartition_sink_t is;
				is.sink = sink;
				is.path = *path;
				interpartition_sinks.emplace_back(is);
			}

			if (sink->id != -1) {
				net_timing.delay[sink->id+1] = get_vertex_props(rt.graph, route_tree_get_rt_node(rt, sink->rr_node)).delay;
			}

			++num_routed_sinks;
		}

		for (const auto &m : modified)  {
			state[m].known_cost = std::numeric_limits<float>::max();
			state[m].cost = std::numeric_limits<float>::max();
			touched[m] = false;
		}

		modified.clear();

		heap = std::priority_queue<route_state_t>();
	}

	//zlog_level(delta_log, ROUTER_V1, "\n");

	//return unrouted_sinks;

	/*delete [] state;*/

	/*check_route_tree(rt, net, g);*/
}

void route_net_with_high_interpartition_cost(const RRGraph &g, const vector<int> &pid, int this_pid, int vpr_id, const source_t *source, const vector<sink_t *> &sinks, const route_parameters_t &params, float pres_fac, route_state_t *state, congestion_locked_t *congestion, route_tree_t &rt, t_net_timing &net_timing, vector<interpartition_sink_t> &interpartition_sinks, bool lock, perf_t *perf, lock_perf_t *lock_perf)
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

		RouteTreeNode root_rt_node = route_tree_add_rr_node(rt, source->rr_node);
		const auto &source_rr_node_p = get_vertex_props(g, source->rr_node);
		route_tree_set_node_properties(rt, root_rt_node, true, source_rr_node_p.R, 0.5 * source_rr_node_p.R * source_rr_node_p.C);
		route_tree_add_root(rt, source->rr_node);

		/*update_one_cost_internal(source_rr_node_p, 1, params.pres_fac);*/
	} else {
		if (rt.root_rt_node_id != -1) {
			const auto &rt_root_p = get_vertex_props(rt.graph, rt.root_rt_node_id);
			if (source && rt_root_p.rr_node != source->rr_node) {
				char root[256];
				char source_str[256];
				sprintf_rr_node(rt_root_p.rr_node, root);
				sprintf_rr_node(source->rr_node, source_str);
				zlog_error(delta_log, "Error: Non empty route tree node has root %s that is different from net source %s\n",
						root, source_str);
				assert(false);
			}
		}
	}
	
	vector<bool> touched(num_vertices(g), false);

	int isink = 0;
	int num_routed_sinks = 0;
	bpqueue boundary_node_heap;
	for (const auto &sink : sorted_sinks) {
		sprintf_rr_node(sink->rr_node, buffer);
		//zlog_level(delta_log, ROUTER_V1, "Net %d Sink %d: %s criticality: %g BB: %d-%d %d-%d Prev BB: %d-%d %d-%d\n", vpr_id, sink->id, buffer, sink->criticality_fac, sink->current_bounding_box.xmin, sink->current_bounding_box.xmax, sink->current_bounding_box.ymin, sink->current_bounding_box.ymax, sink->previous_bounding_box.xmin, sink->previous_bounding_box.xmax, sink->previous_bounding_box.ymin, sink->previous_bounding_box.ymax);

		const auto &sink_rr_node = get_vertex_props(g, sink->rr_node);

		route_tree_multi_root_add_to_heap(rt, g, sink->rr_node, sink->criticality_fac, params.astar_fac, sink->current_bounding_box, heap, perf);

		bool expanded_interpartition_node = false;
		bool explored_interpartition_neighbor = false;

		bool found_sink = false;
		while (!heap.empty() && !found_sink) {
			auto item = heap.top();
			heap.pop();

			if (perf) {
				++perf->num_heap_pops;
			}

			sprintf_rr_node(item.rr_node, buffer);
			const auto &v = get_vertex_props(g, item.rr_node);
			zlog_level(delta_log, ROUTER_V3, "Current: %s %s pid: %d occ/cap: %d/%d prev: %d new_cost: %g new_known: %g new_delay: %g old_cost: %g old_known: %g old_delay: %g\n", buffer, (pid[item.rr_node] != -1 && pid[item.rr_node] != this_pid) ? "OTHER": "", pid[item.rr_node], congestion[item.rr_node].cong.occ, v.capacity, valid(item.prev_edge) ? get_source(g, item.prev_edge) : -1, item.cost, item.known_cost, item.delay, state[item.rr_node].cost, state[item.rr_node].known_cost, state[item.rr_node].delay);

			if (pid[item.rr_node] != -1 && pid[item.rr_node] != this_pid) {
				expanded_interpartition_node = true;
			}

			if (item.rr_node == sink->rr_node) {
				state[item.rr_node] = item;
				modified.push_back(item.rr_node);
				found_sink = true;
			} else if (item.cost < state[item.rr_node].cost) {
				//if (touched[item.rr_node]) {
					//RRGraph temp;
					//add_vertex(temp, num_vertices(g));
					//for (const auto &m : modified) {
						//if (valid(state[m].prev_edge)) {
							//add_edge(temp, get_source(g, state[m].prev_edge), get_target(g, state[m].prev_edge));
						//}
					//}
					//auto ft = make_filtered_graph(temp, [&modified] (RRNode v) -> bool {
							//return find(begin(modified), end(modified), v) != end(modified);
							//},
							//[] (const RREdge &e) -> bool {
							//return true;
							//});

					//write_graph(ft, "weird.dot", [&state, &item] (RRNode n) -> string {
						//char buffer[256];
						//if (n == item.rr_node) {
						//sprintf(buffer, "label=\"%d %g %g %g\" style=filled fillcolor=red", n, state[n].cost, state[n].known_cost, state[n].cost-state[n].known_cost);
						//} else {
						//sprintf(buffer, "label=\"%d %g %g %g\"", n, state[n].cost, state[n].known_cost, state[n].cost-state[n].known_cost);
						//}
						//return string(buffer);
						//}, [] (const RREdge &e) -> string {
						//return string();
						//});
					//assert(false);
				//}
				assert(item.known_cost < state[item.rr_node].known_cost);
				//touched[item.rr_node] = true;
				state[item.rr_node] = item;
				modified.push_back(item.rr_node);

				expand_neighbors_with_high_interpartition_cost(g, item.rr_node, state, congestion, sink_rr_node, sink->criticality_fac, params.astar_fac, 
						[&g, &sink, &sink_rr_node, &v, &pid, &boundary_node_heap, &item, &this_pid, &congestion] (const RRNode &n) -> bool {
					const auto &prop = get_vertex_props(g, n);

					//if (congestion[n].occ >= prop.capacity) {
					//zlog_level(delta_log, ROUTER_V3, "congested\n");
					//return false;
					//}
					/*if (trace_has_node(prev_trace, id(n))) {*/
						/*zlog_level(delta_log, ROUTER_V3, " existing node route tree ");*/
					/*}*/

					if (prop.xhigh < sink->current_bounding_box.xmin
							|| prop.xlow > sink->current_bounding_box.xmax
							|| prop.yhigh < sink->current_bounding_box.ymin
							|| prop.ylow > sink->current_bounding_box.ymax) {
					zlog_level(delta_log, ROUTER_V3, "outside of bounding box\n");
					return false;
					}

					if (prop.type == IPIN
							&& (prop.xhigh != sink_rr_node.xhigh 
								|| prop.yhigh != sink_rr_node.yhigh)) {
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
				}, pid, this_pid, heap, explored_interpartition_neighbor, perf);
			} else {
				zlog_level(delta_log, ROUTER_V3, "\tNot expanding neighbor because known_cost %g > %g or cost %g > %g\n", item.known_cost, state[item.rr_node].known_cost, item.cost, state[item.rr_node].cost);
			}
		}

		if (!found_sink) {
			assert(heap.empty());
			zlog_error(delta_log, "Error: Failed to find sink %d\n", sink->rr_node);
			assert(false);
		} else {
			const auto &path = get_path(sink->rr_node, state, rt, g, vpr_id);

			route_tree_add_path(rt, path, g, state);

			vector<RRNode> added_nodes;
			for (const auto &n : *path) {
				if (valid(n.prev_edge) || get_vertex_props(g, n.rr_node_id).type == SOURCE) {
					added_nodes.push_back(n.rr_node_id);
				} 
			}

			update_one_cost(g, congestion, added_nodes.begin(), added_nodes.end(), 1, pres_fac, lock, lock_perf);

			if (!expanded_interpartition_node) {
			} else {
				interpartition_sink_t is;
				is.sink = sink;
				is.path = *path;
				interpartition_sinks.emplace_back(is);
			}

			if (sink->id != -1) {
				net_timing.delay[sink->id+1] = get_vertex_props(rt.graph, route_tree_get_rt_node(rt, sink->rr_node)).delay;
			}

			++num_routed_sinks;
		}

		for (const auto &m : modified)  {
			state[m].known_cost = std::numeric_limits<float>::max();
			state[m].cost = std::numeric_limits<float>::max();
			touched[m] = false;
		}

		modified.clear();

		heap = std::priority_queue<route_state_t>();
	}

	//zlog_level(delta_log, ROUTER_V1, "\n");

	//return unrouted_sinks;

	/*delete [] state;*/

	/*check_route_tree(rt, net, g);*/
}

void route_net_with_partition_hopping(const RRGraph &g, const vector<int> &pid, int this_pid, int vpr_id, const source_t *source, const vector<sink_t *> &sinks, const route_parameters_t &params, float pres_fac, route_state_t *state, congestion_locked_t *congestion, route_tree_t &rt, t_net_timing &net_timing, unrouted_t &unrouted, int num_partitions, bool lock, perf_t *perf, lock_perf_t *lock_perf)
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

		RouteTreeNode root_rt_node = route_tree_add_rr_node(rt, source->rr_node);
		const auto &source_rr_node_p = get_vertex_props(g, source->rr_node);
		route_tree_set_node_properties(rt, root_rt_node, true, source_rr_node_p.R, 0.5 * source_rr_node_p.R * source_rr_node_p.C);
		route_tree_add_root(rt, source->rr_node);

		/*update_one_cost_internal(source_rr_node_p, 1, params.pres_fac);*/
	} else {
		if (rt.root_rt_node_id != -1) {
			const auto &rt_root_p = get_vertex_props(rt.graph, rt.root_rt_node_id);
			if (source && rt_root_p.rr_node != source->rr_node) {
				char root[256];
				char source_str[256];
				sprintf_rr_node(rt_root_p.rr_node, root);
				sprintf_rr_node(source->rr_node, source_str);
				zlog_error(delta_log, "Error: Non empty route tree node has root %s that is different from net source %s\n",
						root, source_str);
				assert(false);
			}
		}
	}

	int isink = 0;
	int num_routed_sinks = 0;
	bpqueue boundary_node_heap;
	for (const auto &sink : sorted_sinks) {
		sprintf_rr_node(sink->rr_node, buffer);
		//zlog_level(delta_log, ROUTER_V1, "Net %d Sink %d: %s criticality: %g BB: %d-%d %d-%d Prev BB: %d-%d %d-%d\n", vpr_id, sink->id, buffer, sink->criticality_fac, sink->current_bounding_box.xmin, sink->current_bounding_box.xmax, sink->current_bounding_box.ymin, sink->current_bounding_box.ymax, sink->previous_bounding_box.xmin, sink->previous_bounding_box.xmax, sink->previous_bounding_box.ymin, sink->previous_bounding_box.ymax);

		const auto &sink_rr_node = get_vertex_props(g, sink->rr_node);

		route_tree_multi_root_add_to_heap(rt, g, sink->rr_node, sink->criticality_fac, params.astar_fac, sink->current_bounding_box, heap, perf);

		bool found_sink = false;
		while (!heap.empty() && !found_sink) {
			auto item = heap.top();
			heap.pop();

			if (perf) {
				++perf->num_heap_pops;
			}

			sprintf_rr_node(item.rr_node, buffer);
			const auto &v = get_vertex_props(g, item.rr_node);
			zlog_level(delta_log, ROUTER_V3, "Current: %s pid: %d occ/cap: %d/%d prev: %d new_cost: %g new_known: %g new_delay: %g old_cost: %g old_known: %g old_delay: %g\n", buffer, pid[item.rr_node], congestion[item.rr_node].cong.occ, v.capacity, valid(item.prev_edge) ? get_source(g, item.prev_edge) : -1, item.cost, item.known_cost, item.delay, state[item.rr_node].cost, state[item.rr_node].known_cost, state[item.rr_node].delay);

			assert(pid[item.rr_node] == -1 || pid[item.rr_node] == this_pid);

			if (item.rr_node == sink->rr_node) {
				state[item.rr_node] = item;
				modified.push_back(item.rr_node);
				found_sink = true;
			} else if (item.cost < state[item.rr_node].cost) {
				assert(item.known_cost < state[item.rr_node].known_cost);
				state[item.rr_node] = item;
				modified.push_back(item.rr_node);

				expand_neighbors_with_partition_hopping(g, item.rr_node, state, congestion, sink_rr_node, sink->criticality_fac, params.astar_fac, heap, [&g, &sink, &sink_rr_node, &v, &pid, &boundary_node_heap, &item, &this_pid] (const RRNode &n) -> bool {

					/*if (trace_has_node(prev_trace, id(n))) {*/
						/*zlog_level(delta_log, ROUTER_V3, " existing node route tree ");*/
					/*}*/
					const auto &prop = get_vertex_props(g, n);

					if (prop.xhigh < sink->current_bounding_box.xmin
							|| prop.xlow > sink->current_bounding_box.xmax
							|| prop.yhigh < sink->current_bounding_box.ymin
							|| prop.ylow > sink->current_bounding_box.ymax) {
					zlog_level(delta_log, ROUTER_V3, "outside of bounding box\n");
					return false;
					}

					if (prop.type == IPIN
							&& (prop.xhigh != sink_rr_node.xhigh 
								|| prop.yhigh != sink_rr_node.yhigh)) {
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
					if (pid[n] != -1 && pid[n] != this_pid) {
						zlog_level(delta_log, ROUTER_V3, " pid %d not in current partition %d\n", pid[n], this_pid);
						return false;
					}

					return true;
				}, pid, this_pid, boundary_node_heap, lock, perf, lock_perf);
			} else {
				zlog_level(delta_log, ROUTER_V3, "\tNot expanding neighbor because known_cost %g > %g or cost %g > %g\n", item.known_cost, state[item.rr_node].known_cost, item.cost, state[item.rr_node].cost);
			}
		}

		if (!found_sink) {
			assert(heap.empty());
			zlog_error(delta_log, "Error: Failed to find sink %d\n", sink->rr_node);

			unrouted.unrouted_sinks.push_back(sink);
		} else {
			const auto &path = get_path(sink->rr_node, state, rt, g, vpr_id);

			route_tree_add_path(rt, path, g, state);

			vector<RRNode> added_nodes;
			for (const auto &n : *path) {
				if (valid(n.prev_edge) || get_vertex_props(g, n.rr_node_id).type == SOURCE) {
					added_nodes.push_back(n.rr_node_id);
				} 
			}

			update_one_cost(g, congestion, added_nodes.begin(), added_nodes.end(), 1, pres_fac, lock, lock_perf);

			if (sink->id != -1) {
				net_timing.delay[sink->id+1] = get_vertex_props(rt.graph, route_tree_get_rt_node(rt, sink->rr_node)).delay;
			}

			++num_routed_sinks;
		}

		for (const auto &m : modified)  {
			state[m].known_cost = std::numeric_limits<float>::max();
			state[m].cost = std::numeric_limits<float>::max();
		}

		modified.clear();

		heap = std::priority_queue<route_state_t>();
	}

	int num_visited_partitions = 0;
	vector<bool> visited_partitions(num_partitions, false);

	set<int> added_boundary_nodes;

	while (!boundary_node_heap.empty() && num_visited_partitions < num_partitions-1) {
		RRNode bn; float cost;
		std::tie(cost, bn) = boundary_node_heap.top(); boundary_node_heap.pop();

		zlog_level(delta_log, ROUTER_V3, "Current boundary node %d with cost %g", bn, cost);

		if (added_boundary_nodes.find(bn) == added_boundary_nodes.end()) {
			boundary_node_t new_bn;
			new_bn.rr_node = bn;
			//new_bn.path = get_path(new_bn.rr_node, state, rt, g, vpr_id);

			added_boundary_nodes.insert(bn);
			for (const auto &e : get_out_edges(g, bn)) {
				int to = get_target(g, e);
				int to_pid = pid[to];
				if (to_pid != -1 && to_pid != this_pid && !visited_partitions[to_pid]) {
					visited_partitions[to_pid] = true;
					++num_visited_partitions;
				}
			}
			zlog_level(delta_log, ROUTER_V3, " ADDED");

			unrouted.boundary_nodes.emplace_back(new_bn);
		}
		zlog_level(delta_log, ROUTER_V3, "\n");
	}

	if (num_visited_partitions < num_partitions-1) {
		zlog_level(delta_log, ROUTER_V3, "Net %d does't have boundaries nodes that covers all partitions\n", vpr_id);
	}

	//zlog_level(delta_log, ROUTER_V1, "\n");

	//return unrouted_sinks;

	/*delete [] state;*/

	/*check_route_tree(rt, net, g);*/
}
