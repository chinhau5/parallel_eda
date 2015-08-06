#include "route_common_types.h"
#include "route_tree_timing_types.h"

t_rt_node *
thread_safe_update_route_tree(const struct s_heap * hptr, const t_rr_node_route_inf *l_rr_node_route_inf, t_rt_node **l_rr_node_to_rt_node);

t_rt_node *
init_route_tree_to_source(t_rt_node **l_rr_node_to_rt_node, int inet);
