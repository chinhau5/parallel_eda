#ifndef ROUTE_H
#define ROUTE_H

#include <chrono>
#include <mpi.h>
#include "queue.h"
#include "geometry.h"
#include "quadtree.h"
#include "new_rr_graph.h"

typedef struct lock_perf_t {
	unsigned long num_lock_waits;
	unsigned long num_lock_tries;
	std::chrono::high_resolution_clock::duration total_wait_time;
} lock_perf_t;

typedef struct perf_t {
	unsigned long num_heap_pushes;
	unsigned long num_heap_pops;
	unsigned long num_neighbor_visits;
	unsigned long max_num_heap_pops_per_sink;
	unsigned long min_num_heap_pops_per_sink;
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

typedef struct mpi_perf_t {
	std::chrono::high_resolution_clock::duration total_sync_time;
	std::chrono::high_resolution_clock::duration total_broadcast_time;
	int num_syncs_while_expanding;

	std::chrono::high_resolution_clock::duration total_send_testsome_time;
	std::chrono::high_resolution_clock::duration total_testsome_time;
	std::chrono::high_resolution_clock::duration total_iprobe_time;
	std::chrono::high_resolution_clock::duration total_probe_time;
	std::chrono::high_resolution_clock::duration total_irecv_time;
	std::chrono::high_resolution_clock::duration total_recv_time;
	std::chrono::high_resolution_clock::duration total_update_time;
	int total_calls;
} mpi_perf_t;

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
	//int xhigh;
	//int yhigh;
} source_t;

typedef struct sink_t {
	struct net_t *net;
	int id;
	float criticality_fac;
	int rr_node;
	int x;
	int y;
	//int xhigh;
	//int yhigh;
	//source_t source;
	struct new_virtual_net_t *vnet;
	bounding_box_t current_bounding_box;
	//bounding_box_t previous_bounding_box;
	//bounding_box_t scheduler_bounding_box;
	//int bb_factor;
	//int distance_to_source_rank;
	//int congested_iterations;
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
	bounding_box_t bounding_box;

	//vector<struct virtual_net_t *> virtual_nets;

	//int bb_area_rank;
	//unsigned long num_bounding_box_updates;
	//unsigned long num_nearest_iters;
	//unsigned long total_point_tree_size;

	//int previous_sink_index;
	/*bool previous_sink_valid;*/
	//bounding_box_t current_bounding_box;
	/*bool previous_bounding_box_valid;*/
	//int current_local_id;
	//bool has_sink;
	//int num_sinks_routed;
	//int current_sink_index;
	//sink_t *current_sink;
	//vector<bool> sink_routed;
	//source_t current_source;

	//std::vector<bool> overlapping_nets;
	//std::vector<bool> non_overlapping_nets;
	//int num_overlapping_nets;
	//int num_non_overlapping_nets;
	//std::vector<const net_t *> overlapping_nets_vec;
	//std::vector<const net_t *> non_overlapping_nets_vec;
	//int num_local_nets;
	//int pid;
	//int schedule;
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
//

typedef struct congestion_t {
	//tbb::spin_mutex lock;
	int occ;
	float pres_cost;
	float acc_cost;
	int recalc_occ;
} congestion_t;

typedef struct alignas(64) congestion_locked_t {
	tbb::spin_mutex lock;
	congestion_t cong;
} congestion_locked_t;

typedef struct local_cost_t {
	int occ_delta;
	float pres_cost;
} local_cost_t;

typedef struct congestion_local_t {
	int tid;
	congestion_t *global;
	vector<local_cost_t> local;
	set<int> dirty_nodes;
} congestion_local_t;

typedef struct route_state_t {
	int rr_node;
	RREdge prev_edge;
	float upstream_R;
	float delay;
	float known_cost;
	float cost;
} route_state_t;

typedef struct path_node_t {
	RRNode rr_node_id;
	RREdge prev_edge;
} path_node_t;

typedef struct path_node_send_t {
	int rr_node_id;
	int parent_rr_node;
	float delay;
	float upstream_R;
} path_node_send_t;

typedef struct boundary_node_t {
	RRNode rr_node;
	vector<path_node_t> path;
} boundary_node_t;

typedef	struct node_update_t {
	int rr_node;
	int delta;
} node_update_t;

typedef struct ongoing_transaction_t {
	std::shared_ptr<vector<node_update_t>> data;
	MPI_Request req;
} ongoing_transaction_t;

typedef struct broadcast_data_t {
	bool is_meta;
	int meta_ref;
	int data_ref;
} broadcast_data_t;

typedef struct request_meta_t {
	int data_index;
	int data_offset;
	int data_count;
} request_meta_t;

typedef struct large_t {
	int data_index;
	int data_offset;
} large_t;

typedef struct mpi_context_t {
	MPI_Comm comm;
	int rank;
	int comm_size;

	int buffer_size;
	int bit_width;

	vector<vector<int> *> pending_send_data_nbc;
	vector<int> pending_send_data_ref_count;

	vector<MPI_Request> pending_send_req;
	vector<int> pending_send_req_data_ref;
	vector<int> pending_send_req_dst;
	vector<int> completed_indices;
	vector<MPI_Status> completed_recv_statuses;
	int num_pending_reqs;
	vector<int> num_pending_reqs_by_rank;
	vector<vector<int>> num_pending_reqs_by_time;

	vector<bool> received_last_update; 
	vector<large_t> large_packets;

	queue<int> free_send_data_index;
	queue<int> free_send_req_index;

	int max_send_data_size;
	unsigned long total_send_data_size;
	unsigned long total_send_count;
	int max_req_buffer_size;
	int max_active_reqs;
	vector<int> max_pending_send_reqs_by_rank;

	// -----------

	vector<MPI_Comm> ibcast_comm;
	vector<int> max_ibcast_count;
	vector<request_meta_t> pending_req_meta;

	vector<int> all_data_sizes;
	vector<int> packet_id;

	int num_broadcasts_required;

	vector<std::shared_ptr<vector<node_update_t>>> pending_send_data;
	int send_req_queue_size;
	int send_req_queue_head;
	int send_req_queue_tail;

	MPI_Request active_req;
	broadcast_data_t active_bcast_data;
	queue_t<broadcast_data_t> pending_bcast_data_q;

	vector<vector<int>> pending_send_req_refs;

	vector<pair<vector<int> *, int>> pending_send_data_raw;
	vector<int> pending_send_req_meta_data_ref;

	vector<int> completed_send_indices;

	vector<vector<std::shared_ptr<vector<node_update_t>>>> pending_recv_data; 
	vector<vector<MPI_Request>> pending_recv_req; 

	vector<vector<int>> pending_recv_data_flat; 
	vector<MPI_Request> pending_recv_req_flat; 

	vector<int> completed_recv_indices;
} mpi_context_t;

//typedef struct unrouted_sink_t {
	//sink_t *sink;
	//vector<boundary_node_t> boundary_nodes;
//} unrouted_sink_t;
//
typedef struct interpartition_sink_t {
	sink_t *sink;
	vector<path_node_t> path;
} interpartition_sink_t;

typedef struct unrouted_t {
	vector<sink_t *> unrouted_sinks;
	vector<boundary_node_t> boundary_nodes;
} unrouted_t;

typedef struct pseudo_source_t {
	RRNode node;
	RREdge prev_edge;
	boundary_node_t *bnode;
} pseudo_source_t;

typedef struct pseudo_sink_t {
	sink_t *sink;
	vector<path_node_t> path;
	pseudo_source_t *pseudo_source;
} pseudo_sink_t;

typedef struct pseudo_net_t {
	net_t *net;
	vector<pseudo_source_t> pseudo_sources;
	vector<pseudo_sink_t> pseudo_sinks;
} pseudo_net_t;

bool hybrid_route(t_router_opts *opts);
bool partitioning_route(t_router_opts *opts);
bool greedy_route(t_router_opts *opts);
bool partitioning_route_bounding_box(t_router_opts *opts);

#endif
