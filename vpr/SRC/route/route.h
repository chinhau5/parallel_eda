#ifndef ROUTE_H
#define ROUTE_H

#include <chrono>
#include <assert.h>
#include <boost/numeric/interval.hpp>
#include <boost/geometry.hpp>
#include <tbb/tbb.h>
#include "vpr_types.h"
#include "graph.h"
#include "geometry.h"

typedef struct perf_t {
	unsigned long num_heap_pushes;
	std::chrono::high_resolution_clock::duration total_route_time;
	std::chrono::high_resolution_clock::duration total_wait_time;
} perf_t;

typedef struct sched_perf_t {
	std::chrono::high_resolution_clock::duration total_rtree_build_time;
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

typedef struct bounding_box_t {
	int xmin;
	int xmax;
	int ymin;
	int ymax;
	bool operator==(const bounding_box_t &other) const {
		return other.xmin == xmin && 
			other.xmax == xmax && 
			other.ymin == ymin && 
			other.ymax == ymax;
	}
} bounding_box_t;

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
	//std::mutex lock;
	int vpr_id;
	int local_id;
	int current_local_id;
	int bb_area_rank;
	/*bool global;*/
	bool has_sink;
	source_t source;
	source_t current_source;
	//source_t previous_source;
	/*bool previous_source_valid;*/
	std::vector<sink_t> sinks;
	int num_sinks_routed;
	int current_sink_index;
	sink_t *current_sink;
	vector<bool> sink_routed;

	vector<struct virtual_net_t *> virtual_nets;
	//int previous_sink_index;
	/*bool previous_sink_valid;*/
	//bounding_box_t current_bounding_box;
	/*bool previous_bounding_box_valid;*/
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
	source_t *source;
	vector<sink_t *> sinks;
	vector<sink_t *> current_sinks;
	//point centroid;
	int nearest_rr_node;
	box current_bounding_box;
	box scheduler_bounding_box;
	box saved_scheduler_bounding_box;
} virtual_net_t;

//typedef struct simple_net_t {
	//const net_t *parent_net;
	//source_t source
//} 

typedef struct rr_node_property_t {
	t_rr_type type;
	bool inc_direction;
	int xlow;
	int ylow;
	int xhigh;
	int yhigh;
	int real_xlow;
	int real_ylow;
	int real_xhigh;
	int real_yhigh;
	float R;
	float C;
	int cost_index;
	int capacity;
	int occ;
	int recalc_occ;
	float pres_cost;
	float acc_cost;
} rr_node_property_t;

typedef struct rr_edge_property_t {
	bool buffered;
	float switch_delay;
	float R;
} rr_edge_property_t;

typedef graph_t<rr_node_property_t, rr_edge_property_t> RRGraph;
typedef vertex_t<rr_node_property_t, rr_edge_property_t> RRNode;
typedef edge_t<rr_edge_property_t> RREdge;

typedef struct route_state_t {
	int rr_node;
	const RREdge *prev_edge;
	float upstream_R;
	float delay;
	float known_cost;
	float cost;
} route_state_t;

template<typename BoundingBox>
int get_bounding_box_area(const BoundingBox &bb)
{
	assert(bb.xmax >= bb.xmin && bb.ymax >= bb.ymin);
	int area = (bb.xmax - bb.xmin + 1) * (bb.ymax - bb.ymin + 1);
	assert(area >= 0);
	return area;
}

template<typename BoundingBox>
bool box_overlap(const BoundingBox &box_a, const BoundingBox &box_b)
{
	using namespace boost::numeric;
	interval<int> a_hor(box_a.xmin, box_a.xmax);
	interval<int> b_hor(box_b.xmin, box_b.xmax);

	interval<int> a_vert(box_a.ymin, box_a.ymax);
	interval<int> b_vert(box_b.ymin, box_b.ymax);

	return overlap(a_hor, b_hor) && overlap(a_vert, b_vert);
}

bool hybrid_route(t_router_opts *opts);
bool partitioning_route(t_router_opts *opts);
bool greedy_route(t_router_opts *opts);
bool partitioning_route_bounding_box(t_router_opts *opts);

float get_timing_driven_expected_cost(const RRNode &current, const RRNode &target,
		float criticality_fac, float R_upstream);

void update_one_cost(RRGraph &g, const vector<int>::const_iterator &rr_nodes_begin, const vector<int>::const_iterator &rr_nodes_end, int delta, float pres_fac);

#include "route_tree.h"

void update_one_cost(RRGraph &g, route_tree_t &rt, const RouteTreeNode &node, int delta, float pres_fac);

void update_one_cost_internal(RRNode &rr_node, int delta, float pres_fac);

bool operator<(const route_state_t &a, const route_state_t &b);

#endif
