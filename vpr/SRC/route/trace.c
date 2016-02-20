#include "trace.h"
#include "route.h"
#include "utility.h"
#include "route_tree.h"
#include "log.h"

extern zlog_category_t *delta_log;

void trace_init(trace_t &trace)
{
	/*trace.num_sources = 0;*/
	/*trace.first_sink_rr_node = -1;*/
}

bool trace_empty(const trace_t &trace)
{
	return trace.existing_nodes.empty() && trace.paths_starting_with_source.empty() && trace.segments.empty();
}

void trace_check_if_only_path(const trace_t &trace)
{
	/*trace.segments*/
}

bool trace_has_node(const trace_t &trace, int rr_node)
{
	return trace.existing_nodes.find(rr_node) != trace.existing_nodes.end();
}

const Segment &trace_add_path(trace_t &trace, const RRGraph &g, const route_state_t *state, int sink_rr_node, int vpr_net_id)
{
	assert(trace.segments.find(sink_rr_node) == trace.segments.end());

	auto &new_segment = trace.segments.insert(make_pair(sink_rr_node, Segment())).first->second;

	char buffer[256];
	int prev_rr_node = -1;
	int current_rr_node = sink_rr_node;
	while (state[current_rr_node].prev_edge) {
		if (trace.existing_nodes.find(current_rr_node) != trace.existing_nodes.end()) {
			sprintf_rr_node(current_rr_node, buffer);
			char prev_s[256];
			sprintf_rr_node(prev_rr_node, prev_s);
			zlog_warn(delta_log, "Warning: Trying to add existing node %s to trace. Stopping traceback\n", buffer);
			/*assert(false);*/
			/*current_rr_node = prev_rr_node;*/
			break;
		}

		new_segment.push_back(current_rr_node);
		trace.existing_nodes.insert(current_rr_node);

		sprintf_rr_node(current_rr_node, buffer);
		zlog_level(delta_log, ROUTER_V2, "Net %d trace: %s\n", vpr_net_id, buffer);

		int parent_rr_node = id(get_source(g, *state[current_rr_node].prev_edge));
		prev_rr_node = current_rr_node;
		current_rr_node = parent_rr_node;
	}

	/* we're checking for various conditions: */
	const RRNode &first_node = get_vertex(g, current_rr_node);
	sprintf_rr_node(current_rr_node, buffer);
	if (first_node.properties.type == SOURCE) {
		if (trace.paths_starting_with_source.empty()) {
			/* this is a first segment */
			/*assert(trace.segments.size() == 1);*/ /* this is not true after the first iteration because we're ripping up segment by segment */
			assert(trace.existing_nodes.find(current_rr_node) == trace.existing_nodes.end());

			zlog_level(delta_log, ROUTER_V2, "Net %d first segment trace: %s\n", vpr_net_id, buffer);

			/* shouldn't reach this case actually because paths_starting_with_source is empty */
			/*if (trace.existing_nodes.find(current_rr_node) != trace.existing_nodes.end()) {*/
				/*zlog_error(delta_log, "Error: Trying to add existing node %s prev_edge %d to trace\n", buffer, state[current_rr_node].prev_edge ? 1 : 0);*/
				/*assert(false);*/
			/*}*/

			new_segment.push_back(current_rr_node);
			trace.existing_nodes.insert(current_rr_node);
			trace.paths_starting_with_source.push_back(&new_segment);
		} else { /* we have existing paths that start with SOURCE */
			if (trace.existing_nodes.find(current_rr_node) != trace.existing_nodes.end()) {
				/* this is an existing SOURCE
				 * we need to check whether it has capacity of > 1
				 * >1 : no error 
				 * <=1 : error */
				assert(any_of(trace.paths_starting_with_source.begin(), trace.paths_starting_with_source.end(),
							[&current_rr_node] (const Segment *path) -> bool {
							return path->back() == current_rr_node;
							}));

				if (first_node.properties.capacity > 1) {
					zlog_level(delta_log, ROUTER_V2, "Reconnecting to existing %s\n", buffer);
					new_segment.push_back(current_rr_node);
					trace.paths_starting_with_source.push_back(&new_segment);
				} else {
					zlog_error(delta_log, "Error: Trying to reconnect to %s with only capacity of %d\n", buffer, first_node.properties.capacity);
					assert(false);
				}
			} else {
				/* this is a new SOURCE 
				 * error because we now have mutliple sources for the net */
				zlog_error(delta_log, "Error: Trying to start the trace with another source %s when existing sources ", buffer);
				for (const auto &p : trace.paths_starting_with_source) {
					char old_source[256];
					sprintf_rr_node(p->back(), old_source);
					zlog_error(delta_log, "%s ", buffer);
				}
				zlog_error(delta_log, "exist\n");
				assert(false);
			}
		}
	} else {
		assert(first_node.properties.type == CHANX || first_node.properties.type == CHANY || first_node.properties.type == OPIN);

		if (trace.segments.size() == 1) {
			zlog_error(delta_log, "Error: First path started with %s instead of a SOURCE\n", buffer);
			assert(false);
		}
	}

	/*if (get_vertex(g, current_rr_node).properties.type == SOURCE) {*/
		/*if (trace.num_sources != 0) {*/
			/*zlog_error(delta_log, "Error: Trying to add another path that starts from SOURCE [num_sources: %d]\n", trace.num_sources);*/
			/*assert(false);*/
		/*}*/
		/*new_segment.push_back(current_rr_node);*/
		/*++trace.num_sources;*/
	/*}*/

	return new_segment;
}

void trace_rip_up_net(trace_t &trace, RRGraph &g, float pres_fac)
{
	for (const auto &segment : trace.segments) {
		update_one_cost(g, nullptr, segment.second.begin(), segment.second.end(), -1, pres_fac, false, nullptr);
	}
	trace.segments.clear();
	/*trace.first_sink_rr_node = -1;*/
	trace.paths_starting_with_source.clear();
	trace.existing_nodes.clear();
}

void trace_rip_up_segment(trace_t &trace, RRGraph &g, int sink_rr_node, float pres_fac)
{
	/*const vector<int> *path = route_tree_get_path_to_sink(rt, sink_rr_node);*/
	/*if (path) {*/
		/*update_one_cost(g, path->begin(), path->end()-1, -1, pres_fac);*/
		/*route_tree_remove_path(rt, sink_rr_node);*/
	/*}*/
	auto iter = trace.segments.find(sink_rr_node);
	if (iter != trace.segments.end()) {
		update_one_cost(g, nullptr, iter->second.begin(), iter->second.end(), -1, pres_fac, false, nullptr);
		/*if (get_vertex(g, iter->second.back()).properties.type == SOURCE) {*/
			/*--trace.num_sources;*/
			/*if (trace.num_sources != 0) {*/
				/*char buffer[256];*/
				/*sprintf_rr_node(sink_rr_node, buffer);*/
				/*zlog_error(delta_log, "Error: Ripped up a segment for %s but num_sources (%d) != 0\n", buffer, trace.num_sources);*/
				/*assert(false);*/
			/*}*/
		/*}*/

		for (const auto &node : iter->second) {
			char buffer[256];
			sprintf_rr_node(node, buffer);
			if (trace.existing_nodes.find(node) == trace.existing_nodes.end()) {
				/* this is possible because SOURCE with capacity > 1 only appears
				 * once in the set */
				const auto &prop = get_vertex(g, node).properties;
				if (prop.type == SOURCE && prop.capacity > 1 && !trace.paths_starting_with_source.empty()) {
					zlog_level(delta_log, ROUTER_V2, "Ripping up %s with capacity > 1\n", buffer);
				} else {
					zlog_error(delta_log, "Error: Ripping up non-existing node %s\n", buffer);
					assert(false);
				}
			} else {
				zlog_level(delta_log, ROUTER_V2, "Ripping up node %s from the trace\n", buffer);
				trace.existing_nodes.erase(node);
			}
		}

		if (get_vertex(g, iter->second.back()).properties.type == SOURCE) {
			auto ptr_iter = find_if(trace.paths_starting_with_source.begin(), trace.paths_starting_with_source.end(), [&sink_rr_node] (const Segment *path) {
					return path->front() == sink_rr_node;
					});
			assert(ptr_iter != trace.paths_starting_with_source.end());
			trace.paths_starting_with_source.erase(ptr_iter);
		}

		trace.segments.erase(iter);
	}
}

