#ifndef ROUTER_H
#define ROUTER_H

#include "route.h"
#include "route_tree.h"
#include "congestion.h"

bool operator<(const route_state_t &a, const route_state_t &b);

float get_timing_driven_expected_cost(const rr_node_property_t &current, const rr_node_property_t &target,
		float criticality_fac, float R_upstream);

float analyze_timing(t_net_timing *net_timing);

void update_sink_criticalities(net_t &net, const t_net_timing &net_timing, const route_parameters_t &params);


template<typename Congestion>
void recalculate_occ_internal(const route_tree_t &rt, RouteTreeNode rt_node, const RRGraph &g, Congestion *congestion)
{
	const auto &rt_node_p = get_vertex_props(rt.graph, rt_node);

	assert(rt_node_p.valid);

	int rr_node = rt_node_p.rr_node;
	auto &rr_node_p = get_vertex_props(g, rr_node);
	//if (rr_node_p.type == SOURCE) {
		//get_recalc_occ(congestion[rr_node]) += num_out_edges(rt.graph, rt_node);
	//} else {
		++get_recalc_occ(congestion[rr_node]);
	//}

	for (const auto &branch : route_tree_get_branches(rt, rt_node)) {
		/*for_all_out_edges(rt.graph, node, [&rt, &g, &visited_sinks, &visited_nodes] (const RouteTreeEdge &e) -> void {*/
		const auto &child = get_target(rt.graph, branch);
		recalculate_occ_internal(rt, child, g, congestion);
	}
}

template<typename Congestion>
void recalculate_occ(const route_tree_t &rt, const RRGraph &g, Congestion *congestion)
{
	for (const auto &root : rt.root_rt_nodes) {
		recalculate_occ_internal(rt, root, g, congestion);
	}
}

template<typename Congestion>
void recalculate_occ_internal_locking_route(const route_tree_t &rt, RouteTreeNode rt_node, const RRGraph &g, Congestion *congestion)
{
	const auto &rt_node_p = get_vertex_props(rt.graph, rt_node);

	assert(rt_node_p.valid);

	int rr_node = rt_node_p.rr_node;
	auto &rr_node_p = get_vertex_props(g, rr_node);
	if (rr_node_p.type == SOURCE) {
		get_recalc_occ(congestion[rr_node]) += num_out_edges(rt.graph, rt_node);
	} else {
		++get_recalc_occ(congestion[rr_node]);
	}

	for (const auto &branch : route_tree_get_branches(rt, rt_node)) {
		/*for_all_out_edges(rt.graph, node, [&rt, &g, &visited_sinks, &visited_nodes] (const RouteTreeEdge &e) -> void {*/
		const auto &child = get_target(rt.graph, branch);
		recalculate_occ_internal_locking_route(rt, child, g, congestion);
	}
}

template<typename Congestion>
void recalculate_occ_locking_route(const route_tree_t &rt, const RRGraph &g, Congestion *congestion)
{
	for (const auto &root : rt.root_rt_nodes) {
		recalculate_occ_internal_locking_route(rt, root, g, congestion);
	}
}

void check_route_tree(const route_tree_t &rt, const net_t &net, RRGraph &g);
void check_route_tree(const route_tree_t &rt, const net_t &net, const vector<sink_t *> &routed_sinks, RRGraph &g);

template<typename Congestion>
void get_overused_nodes(const route_tree_t &rt, RouteTreeNode rt_node, const RRGraph &g, const Congestion *congestion, vector<int> &overused_rr_node)
{
	const auto &rt_node_p = get_vertex_props(rt.graph, rt_node);

	assert(rt_node_p.valid);

	int rr_node = rt_node_p.rr_node;
	const auto &rr_node_p = get_vertex_props(g, rr_node);

	if (get_occ(congestion[rr_node]) > rr_node_p.capacity) {
		overused_rr_node.push_back(rr_node);
	}
	
	for (const auto &e : get_out_edges(rt.graph, rt_node)) {
		const auto &neighbor = get_target(rt.graph, e);
		get_overused_nodes(rt, neighbor, g, congestion, overused_rr_node);
	}
}


template<typename Congestion>
bool feasible_routing(const RRGraph &g, const Congestion *congestion)
{
	bool feasible = true;

	for (int i = 0; i < num_vertices(g) && feasible; ++i) {
		if (get_occ(congestion[i]) > get_vertex_props(g, i).capacity) {
			feasible = false;
		}
	}

	return feasible;
}

template<typename Congestion>
float get_congestion_cost(const Congestion &congestion, int cost_index)
{
	extern t_rr_indexed_data *rr_indexed_data;
	/*zlog_level(delta_log, ROUTER_V3, " [pres: %g acc: %g] ", v.pres_cost, v.acc_cost);*/
	return rr_indexed_data[cost_index].base_cost * congestion.acc_cost * congestion.pres_cost;
}

//void route_net(RRGraph &g, const net_t &net, const route_parameters_t &params, route_state_t *state, route_tree_t &rt, trace_t &trace, t_net_timing &net_timing, perf_t *perf);

//void route_net(RRGraph &g, int vpr_id, const source_t *source, const vector<sink_t *> &sinks, const route_parameters_t &params, route_state_t *state, route_tree_t &rt, trace_t &trace, const trace_t &prev_trace, t_net_timing &net_timing, perf_t *perf);
//
//
void route_net_one_pass(const RRGraph &g, int vpr_id, const source_t *source, const vector<sink_t *> &sinks, const route_parameters_t &params, route_state_t *state, congestion_t *congestion, route_tree_t &rt, t_net_timing &net_timing, vector<sink_t *> &routed_sinks, vector<sink_t *> &unrouted_sinks, perf_t *perf);

void route_net_with_partitioned_fine_grain_lock(const RRGraph &g, const vector<int> &pid, int this_pid, int vpr_id, const source_t *source, const vector<sink_t *> &sinks, const route_parameters_t &params, route_state_t *state, congestion_locked_t *congestion, route_tree_t &rt, t_net_timing &net_timing, vector<sink_t *> &routed_sinks, vector<sink_t *> &unrouted_sinks, bool lock, perf_t *perf, lock_perf_t *lock_perf);

void route_net_mpi_send_recv_improved(const RRGraph &g, int vpr_id, const source_t *source, const vector<sink_t *> &sinks, const route_parameters_t &params, route_state_t *state, congestion_t *congestion, route_tree_t &rt, t_net_timing &net_timing, vector<sink_t *> &routed_sinks, vector<sink_t *> &unrouted_sinks, mpi_context_t *mpi, perf_t *perf, mpi_perf_t *mpi_perf);

void route_net_mpi_send_recv(const RRGraph &g, int this_pid, int vpr_id, const source_t *source, const vector<sink_t *> &sinks, const route_parameters_t &params, route_state_t *state, congestion_t *congestion, int num_procs, MPI_Comm comm, vector<ongoing_transaction_t> &transactions, route_tree_t &rt, t_net_timing &net_timing, vector<sink_t *> &routed_sinks, vector<sink_t *> &unrouted_sinks, perf_t *perf, mpi_perf_t *mpi_perf);

void route_net_lockless(const RRGraph &g, const vector<int> &pid, int this_pid, int vpr_id, const source_t *source, const vector<sink_t *> &sinks, const route_parameters_t &params, route_state_t *state, congestion_t *congestion, route_tree_t &rt, t_net_timing &net_timing, vector<sink_t *> &routed_sinks, vector<sink_t *> &unrouted_sinks, perf_t *perf);

void route_net_mpi_rma(const RRGraph &g, const vector<int> &pid, int this_pid, MPI_Win win, int vpr_id, const source_t *source, const vector<sink_t *> &sinks, const route_parameters_t &params, route_state_t *state, congestion_t *congestion, route_tree_t &rt, t_net_timing &net_timing, vector<interpartition_sink_t> &interpartition_sinks, perf_t *perf);

void route_net_with_high_interpartition_cost(const RRGraph &g, const vector<int> &pid, int this_pid, int vpr_id, const source_t *source, const vector<sink_t *> &sinks, const route_parameters_t &params, route_state_t *state, congestion_locked_t *congestion, route_tree_t &rt, t_net_timing &net_timing, vector<interpartition_sink_t> &interpartition_sinks, bool lock, perf_t *perf, lock_perf_t *lock_perf);

void route_net_with_partition_hopping(const RRGraph &g, const vector<int> &pid, int this_pid, int vpr_id, const source_t *source, const vector<sink_t *> &sinks, const route_parameters_t &params, route_state_t *state, congestion_locked_t *congestion, route_tree_t &rt, t_net_timing &net_timing, unrouted_t &unrouted, int num_partitions, bool lock, perf_t *perf, lock_perf_t *lock_perf);

vector<const sink_t *> route_net_with_fine_grain_lock(const RRGraph &g, int vpr_id, const source_t *source, const vector<const sink_t *> &sinks, const route_parameters_t &params, route_state_t *state, congestion_locked_t *congestion, route_tree_t &rt, t_net_timing &net_timing, bool lock, perf_t *perf, lock_perf_t *lock_perf);


#endif
