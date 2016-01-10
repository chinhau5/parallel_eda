#ifndef NEW_RR_GRAPH_H
#define NEW_RR_GRAPH_H

#include "vpr_types.h"
#include "graph.h"
#include <tbb/tbb.h>

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
	tbb::spin_mutex *lock;
} rr_node_property_t;

//typedef struct rr_node_property_t {
	//t_rr_type type;
	//bool inc_direction;
	//int xlow;
	//int ylow;
	//int xhigh;
	//int yhigh;
	//int real_xlow;
	//int real_ylow;
	//int real_xhigh;
	//int real_yhigh;
	//float R;
	//float C;
	//int cost_index;
	//tbb::spin_mutex *lock;
	//int capacity;
	//int occ;
	//std::vector<int> users;
	//int recalc_occ;
	//float pres_cost;
	//float acc_cost;
//} rr_node_property_t;

typedef struct rr_edge_property_t {
	bool buffered;
	float switch_delay;
	float R;
} rr_edge_property_t;

typedef graph_t<rr_node_property_t, rr_edge_property_t> RRGraph;
typedef vertex_t<rr_node_property_t, rr_edge_property_t> RRNode;
typedef edge_t<rr_edge_property_t> RREdge;

#endif
