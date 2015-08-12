#include <queue>
#include "route_common_types.h"

void print_route(char *route_file, int **sink_order);

float get_rr_cong_cost(int inode, int num_sinks);

void node_to_heap(int inode, float cost, int prev_node, int prev_edge,
		float backward_path_cost, float R_upstream,
		t_rr_node_route_inf *l_rr_node_route_inf,
		std::priority_queue<struct s_heap> &heap);

void thread_safe_free_traceback(int inet);

void thread_safe_pathfinder_update_one_cost(int inet, struct s_trace *route_segment_start,
		int add_or_sub, float pres_fac);

void thread_safe_pathfinder_update_cost(float pres_fac, float acc_fac);

struct s_trace *
thread_safe_update_traceback(const struct s_heap *hptr, t_rr_node_route_inf *l_rr_node_route_inf, int inet);

void thread_safe_mark_ends(int inet, t_rr_node_route_inf *l_rr_node_route_inf);

void thread_safe_reserve_locally_used_opins(float pres_fac, boolean rip_up_local_opins,
		t_ivec ** clb_opins_used_locally, int inet);
