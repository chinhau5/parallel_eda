#ifndef TRACE_H
#define TRACE_H

#include <vector>
#include <map>
#include <set>
#include "route.h"

typedef std::vector<int> Segment;

typedef struct trace_t {
	//int first_sink_rr_node;
	std::map<int, Segment> segments;
	/* for debugging */
	//int num_sources;
	std::vector<const Segment *> paths_starting_with_source;
	std::set<int> existing_nodes;
} trace_t;

void trace_init(trace_t &trace);

bool trace_empty(const trace_t &trace);

void trace_check_if_only_path(const trace_t &trace);

bool trace_has_node(const trace_t &trace, int rr_node);

const Segment &trace_add_path(trace_t &trace, const RRGraph &g, const route_state_t *state, int sink_rr_node, int vpr_net_id);

void trace_rip_up_net(trace_t &trace, RRGraph &g, float pres_fac);

void trace_rip_up_segment(trace_t &trace, RRGraph &g, int sink_rr_node, float pres_fac);

#endif
