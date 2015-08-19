#include "route_common_types.h"
/**************** Variables shared by all route_files ***********************/

extern t_rr_node_route_inf *rr_node_route_inf; /* [0..num_rr_nodes-1] */
extern struct s_bb *route_bb; /* [0..num_nets-1]     */

/******* Subroutines in route_common used only by other router modules ******/

void pathfinder_update_one_cost(struct s_trace *route_segment_start,
		int add_or_sub, float pres_fac);

void pathfinder_update_cost(float pres_fac, float acc_fac);

struct s_trace *update_traceback(struct s_heap *hptr, int inet);

void reset_path_costs(void);

float get_rr_cong_cost(int inode);

void mark_ends(int inet);

void node_to_heap(int inode, float cost, int prev_node, int prev_edge,
		float backward_path_cost, float R_upstream);

boolean is_empty_heap(void);

void free_traceback(int inet);

void add_to_mod_list(float *fptr);

struct s_heap *get_heap_head(void);

void empty_heap(void);

void free_heap_data(struct s_heap *hptr);

void invalidate_heap_entries(int sink_node, int ipin_node);

void init_route_structs(int bb_factor);

void free_rr_node_route_structs(void);

void alloc_and_load_rr_node_route_structs(void);

void reset_rr_node_route_structs(void);

void alloc_route_static_structs(void);

void free_trace_structs(void);

void reserve_locally_used_opins(float pres_fac, boolean rip_up_local_opins,
		t_ivec ** clb_opins_used_locally);

void free_chunk_memory_trace(void);

void load_route_bb(int bb_factor);
