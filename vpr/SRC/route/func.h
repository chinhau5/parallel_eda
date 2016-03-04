#ifndef FUNC_H
#define FUNC_H

#include "route.h"
#include "route_tree.h"
#include "new_rr_graph.h"

void update_virtual_net_bounding_box(virtual_net_t &virtual_net, route_tree_t &rt, const RRGraph &g, float astar_fac, perf_t *perf);

void update_virtual_net_current_sinks(virtual_net_t &virtual_net, const route_tree_t &rt);

void update_virtual_net_scheduler_bounding_box(virtual_net_t &virtual_net, const box &initial_scheduler_bounding_box);

#endif
