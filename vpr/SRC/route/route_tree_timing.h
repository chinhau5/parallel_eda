#include "route_tree_timing_types.h"

/**************** Subroutines exported by route_tree_timing.c ***************/

void alloc_route_tree_timing_structs(void);

void free_route_tree_timing_structs(void);

t_rt_node *init_route_tree_to_source(int inet);
t_rt_node *init_route_tree_to_source(t_rt_node **l_rr_node_to_rt_node, int inet);

void free_route_tree(t_rt_node * rt_node);

t_rt_node *update_route_tree(struct s_heap *hptr);

void update_net_delays_from_route_tree(float *net_delay,
		t_rt_node ** rt_node_of_sink, int inet);
