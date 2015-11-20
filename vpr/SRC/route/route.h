#ifndef ROUTE_H
#define ROUTE_H

#include "graph.h"

typedef struct perf_t {
	int num_heap_pushes;
} perf_t;

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
} bounding_box_t;

template<typename BoundingBox>
int get_bounding_box_area(const BoundingBox &bb)
{
	int area = (bb.xmax - bb.xmin + 1) * (bb.ymax - bb.ymin + 1);
	assert(area >= 0);
	return area;
}

typedef struct source_t {
	int rr_node;
	int x;
	int y;
} source_t;

typedef struct sink_t {
	int id;
	float criticality_fac;
	int rr_node;
	int x;
	int y;
	source_t source;
	bounding_box_t current_bounding_box;
	bounding_box_t previous_bounding_box;
} sink_t;

bool operator<(const sink_t &a, const sink_t &b);

typedef struct net_t {
	int vpr_id;
	int local_id;
	int current_local_id;
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
	//int previous_sink_index;
	/*bool previous_sink_valid;*/
	bounding_box_t current_bounding_box;
	/*bool previous_bounding_box_valid;*/
	/*bounding_box_t box;*/
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

typedef edge_t<rr_edge_property_t> RREdge;

typedef struct route_state_t {
	int rr_node;
	const RREdge *prev_edge;
	float upstream_R;
	float delay;
	float known_cost;
	float cost;
} route_state_t;

typedef struct rt_node_property_t {
	bool valid;
	bool reexpand;
	int rr_node;
	const RREdge *prev_edge;
	float upstream_R;	
	/*float upstream_R_from_route_state;*/
	float downstream_C;
	float delay;
} rt_node_property_t;

typedef struct rt_edge_property_t {
	const RREdge *rr_edge;
} rt_edge_property_t;

typedef graph_t<rr_node_property_t, rr_edge_property_t> RRGraph;
typedef vertex_t<rr_node_property_t, rr_edge_property_t> RRNode;

typedef graph_t<rt_node_property_t, rt_edge_property_t> RouteTree;
typedef vertex_t<rt_node_property_t, rt_edge_property_t> RouteTreeNode;
typedef edge_t<rt_edge_property_t> RouteTreeEdge;

typedef struct route_tree_t {
	RouteTree graph;
	std::map<int, int> rr_node_to_rt_node;
	/*map<int, vector<int>> sink_rr_node_to_path;*/
	int root;
} route_tree_t;

typedef std::vector<int> Segment;

typedef struct trace_t {
	int first_sink_rr_node;
	std::map<int, Segment> segments;
	/* for debugging */
	//int num_sources;
	std::set<int> existing_nodes;
} trace_t;

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

#endif
