#ifndef FUNC_H
#define FUNC_H

#include "route.h"
#include "route_tree.h"
#include "new_rr_graph.h"

void update_virtual_net_bounding_box(virtual_net_t &virtual_net, route_tree_t &rt, const RRGraph &g, float astar_fac, perf_t *perf);

void update_virtual_net_current_sinks(virtual_net_t &virtual_net, const route_tree_t &rt);

void update_virtual_net_scheduler_bounding_box(virtual_net_t &virtual_net, const box &initial_scheduler_bounding_box);

void update_sink_criticalities(net_t &net, const t_net_timing &net_timing, const route_parameters_t &params);

void check_route_tree(const route_tree_t &rt, const net_t &net, RRGraph &g);

void get_overused_nodes(const route_tree_t &rt, const RouteTreeNode &node, RRGraph &g, vector<int> &overused_rr_node);

bool feasible_routing(const RRGraph &g);

void route_net_2(RRGraph &g, int vpr_id, const source_t *source, const vector<sink_t *> &sinks, const route_parameters_t &params, route_state_t *state, route_tree_t &rt, t_net_timing &net_timing, bool lock, perf_t *perf, lock_perf_t *lock_perf);

#endif
