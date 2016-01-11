#ifndef ROUTE_H
#define ROUTE_H

#include <chrono>
#include "geometry.h"
#include "quadtree.h"
#include "new_rr_graph.h"

typedef struct lock_perf_t {
	unsigned long num_lock_waits;
	unsigned long num_lock_tries;
} lock_perf_t;

typedef struct perf_t {
	unsigned long num_heap_pushes;
	unsigned long num_heap_pops;
	unsigned long num_neighbor_visits;
	std::chrono::high_resolution_clock::duration total_wait_time;
	std::chrono::high_resolution_clock::duration total_rip_up_time;
	std::chrono::high_resolution_clock::duration total_route_time;
	std::chrono::high_resolution_clock::duration total_update_time;
	std::chrono::high_resolution_clock::duration total_push_time;
	std::chrono::high_resolution_clock::duration total_centroid_time;
	std::chrono::high_resolution_clock::duration total_get_nearest_time;
	std::chrono::high_resolution_clock::duration total_verificaion_time;
	std::chrono::high_resolution_clock::duration total_expansion_time;
	std::chrono::high_resolution_clock::duration total_scheduler_box_time;
} perf_t;

typedef struct sched_perf_t {
	unsigned long num_updates;
	unsigned long num_leaf_node_pred_calls;
	unsigned long num_internal_node_pred_calls;
	std::chrono::high_resolution_clock::duration total_build_time;
	std::chrono::high_resolution_clock::duration total_partitioning_time;
	std::chrono::high_resolution_clock::duration total_dispatch_time;
	std::chrono::high_resolution_clock::duration total_rtree_update_time;
	std::chrono::high_resolution_clock::duration total_wait_time;
} sched_perf_t;

typedef struct route_parameters_t {
	float pres_fac;
	float astar_fac;
	float bend_cost;
	float max_criticality;
	float criticality_exp;
} route_parameters_t;

typedef struct source_t {
	struct net_t *net;
	int rr_node;
	int x;
	int y;
} source_t;

typedef struct sink_t {
	struct net_t *net;
	int id;
	float criticality_fac;
	int rr_node;
	int x;
	int y;
	source_t source;
	bounding_box_t current_bounding_box;
	bounding_box_t previous_bounding_box;
	bounding_box_t scheduler_bounding_box;
	int bb_factor;
	int distance_to_source_rank;
	int congested_iterations;
} sink_t;

bool operator<(const sink_t &a, const sink_t &b);

typedef struct sink_schedule_t {
	bounding_box_t current_bounding_box;
	bounding_box_t previous_bounding_box;
} sink_schedule_t;

typedef struct net_t {
	//net_t(const net_t &other) = delete;
	//net_t &operator=(const net_t &other) = delete;
	//std::mutex lock;
	int vpr_id;
	int local_id;
	/*bool global;*/
	source_t source;
	//source_t previous_source;
	/*bool previous_source_valid;*/
	std::vector<sink_t> sinks;

	vector<struct virtual_net_t *> virtual_nets;

	int bb_area_rank;
	unsigned long num_bounding_box_updates;
	unsigned long num_nearest_iters;
	unsigned long total_point_tree_size;

	//int previous_sink_index;
	/*bool previous_sink_valid;*/
	//bounding_box_t current_bounding_box;
	/*bool previous_bounding_box_valid;*/
	int current_local_id;
	bool has_sink;
	int num_sinks_routed;
	int current_sink_index;
	sink_t *current_sink;
	vector<bool> sink_routed;
	source_t current_source;

	bounding_box_t bounding_box;
	std::vector<bool> overlapping_nets;
	std::vector<bool> non_overlapping_nets;
	int num_overlapping_nets;
	int num_non_overlapping_nets;
	std::vector<const net_t *> overlapping_nets_vec;
	std::vector<const net_t *> non_overlapping_nets_vec;
	//int num_local_nets;
	int pid;
	int schedule;
} net_t;

typedef struct virtual_net_t {
	//tbb::spin_mutex lock;
	int id;
	bool routed;
	bool dispatched;
	net_t *net;
	source_t *source;
	vector<sink_t *> sinks;
	vector<sink_t *> current_sinks;
	int nearest_rr_node;
	box sink_bounding_box;
	point centroid;
	box current_bounding_box;
	box scheduler_bounding_box;
	box saved_scheduler_bounding_box;
} virtual_net_t;

//typedef struct simple_net_t {
	//const net_t *parent_net;
	//source_t source
//} 

typedef struct route_state_t {
	int rr_node;
	const RREdge *prev_edge;
	float upstream_R;
	float delay;
	float known_cost;
	float cost;
} route_state_t;

bool hybrid_route(t_router_opts *opts);
bool partitioning_route(t_router_opts *opts);
bool greedy_route(t_router_opts *opts);
bool partitioning_route_bounding_box(t_router_opts *opts);

float get_timing_driven_expected_cost(const RRNode &current, const RRNode &target,
		float criticality_fac, float R_upstream);

bool operator<(const route_state_t &a, const route_state_t &b);

void update_costs(RRGraph &g, float pres_fac, float acc_fac);

float analyze_timing(t_net_timing *net_timing);

int zlog_custom_output(zlog_msg_t *msg);

int zlog_sched_custom_output(zlog_msg_t *msg);

#endif
