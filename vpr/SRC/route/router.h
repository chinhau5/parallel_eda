#ifndef ROUTER_H
#define ROUTER_H

#include "route.h"
#include "route_tree.h"
#include "trace.h"

bool operator<(const route_state_t &a, const route_state_t &b);

float get_timing_driven_expected_cost(const rr_node_property_t &current, const rr_node_property_t &target,
		float criticality_fac, float R_upstream);

float analyze_timing(t_net_timing *net_timing);

void update_sink_criticalities(net_t &net, const t_net_timing &net_timing, const route_parameters_t &params);

void recalculate_occ(const route_tree_t &rt, const RRGraph &g, congestion_t *congestion);

void check_route_tree(const route_tree_t &rt, const net_t &net, RRGraph &g);

void get_overused_nodes(const route_tree_t &rt, RouteTreeNode rt_node, const RRGraph &g, const congestion_t *congestion, vector<int> &overused_rr_node);

bool feasible_routing(const RRGraph &g, const congestion_t *congestion);

void route_net(RRGraph &g, const net_t &net, const route_parameters_t &params, route_state_t *state, route_tree_t &rt, trace_t &trace, t_net_timing &net_timing, perf_t *perf);

void route_net(RRGraph &g, int vpr_id, const source_t *source, const vector<sink_t *> &sinks, const route_parameters_t &params, route_state_t *state, route_tree_t &rt, trace_t &trace, const trace_t &prev_trace, t_net_timing &net_timing, perf_t *perf);

vector<sink_t *> route_net_2(const RRGraph &g, int vpr_id, const source_t *source, const vector<sink_t *> &sinks, const route_parameters_t &params, route_state_t *state, congestion_t *congestion, route_tree_t &rt, t_net_timing &net_timing, bool lock, perf_t *perf, lock_perf_t *lock_perf);

vector<sink_t *> route_net_3(const RRGraph &g, const vector<int> &pid, int this_pid, int vpr_id, const source_t *source, const vector<sink_t *> &sinks, const route_parameters_t &params, route_state_t *state, congestion_t *congestion, route_tree_t &rt, t_net_timing &net_timing, vector<vector<int>> &unrouted_sinks_boundary_nodes, int next_pid, bool lock, perf_t *perf, lock_perf_t *lock_perf);

#endif
