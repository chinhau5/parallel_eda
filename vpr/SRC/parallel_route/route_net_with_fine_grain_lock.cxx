#include "pch.h"

#include "vpr_types.h"
#include "path_delay.h"

#include "route.h"
#include "route_tree.h"
#include "congestion.h"
#include "router.h"
#include "log.h"
#include "expand.h"
#include "router.h"

RREdge get_previous_edge(int rr_node_id, const route_state_t *state, const route_tree_t &rt, const RRGraph &g);

template<typename NodeCallback>
void get_path(int sink_rr_node_id, const route_state_t *state, const route_tree_t &rt, const RRGraph &g, int vpr_net_id, const NodeCallback &callback)
{
	int current_rr_node_id = sink_rr_node_id;
	RREdge previous_edge = get_previous_edge(current_rr_node_id, state, rt, g);

	while (valid(previous_edge)) {

		//char c[256];
		/* printing */
		//sprintf_rr_node(current_rr_node_id, c);
		//zlog_level(delta_log, ROUTER_V2, "Net %d get path: %s\n", vpr_net_id, c);
		//
		callback(previous_edge);

		current_rr_node_id = get_source(g, previous_edge);
		previous_edge = get_previous_edge(current_rr_node_id, state, rt, g);
	}

	//callback(current_rr_node_id, RRGraph::null_edge());
}

vector<const sink_t *> route_net_with_fine_grain_lock(const RRGraph &g, int vpr_id, const source_t *source, const vector<const sink_t *> &sinks, const route_parameters_t &params, route_state_t *state, congestion_locked_t *congestion, route_tree_t &rt, t_net_timing &net_timing, bool lock, perf_t *perf, lock_perf_t *lock_perf)
{
	std::priority_queue<route_state_t> heap;

	vector<int> modified;

	zlog_level(delta_log, ROUTER_V1, "Routing net %d (%lu sinks)\n", vpr_id, sinks.size());

	vector<const sink_t *> sorted_sinks = sinks;
	std::sort(begin(sorted_sinks), end(sorted_sinks), [] (const sink_t *a, const sink_t *b) -> bool {
			return a->criticality_fac > b->criticality_fac;
			});

	char buffer[256];

	sprintf_rr_node(source->rr_node, buffer);
	zlog_level(delta_log, ROUTER_V1, "Source: %s\n", buffer);

	if (route_tree_empty(rt)) {
		/* special case */
		zlog_level(delta_log, ROUTER_V2, "Empty route tree. Setting root to %s\n", buffer);

		const auto &source_rr_node = get_vertex_props(g, source->rr_node);
		RouteTreeNode root_rt_node = route_tree_add_rr_node(rt, source->rr_node, g);
		route_tree_set_node_properties(rt, root_rt_node, true, source_rr_node.R, 0.5 * source_rr_node.R * source_rr_node.C);
		route_tree_add_root(rt, source->rr_node);

		/*update_one_cost_internal(source_rr_node, 1, params.pres_fac);*/
	} else {
		assert(rt.root_rt_nodes.size() == 1);
		const auto &rt_root_p = get_vertex_props(rt.graph, rt.root_rt_nodes[0]);
		if (source && rt_root_p.rr_node != source->rr_node) {
			char root[256];
			char source_str[256];
			sprintf_rr_node(rt_root_p.rr_node, root);
			sprintf_rr_node(source->rr_node, source_str);
			zlog_error(delta_log, "Error: Non empty route tree node has root %s that is different from net source %s\n",
					root, source_str);
			assert(false);
		}
	}

	int isink = 0;
	int num_routed_sinks = 0;
	vector<const sink_t *> unrouted_sinks;
	RRNode existing_opin = RRGraph::null_vertex();
	for (const auto &sink : sorted_sinks) {
		sprintf_rr_node(sink->rr_node, buffer);
		zlog_level(delta_log, ROUTER_V1, "Net %d Sink %d: %s criticality: %g BB: %d-%d %d-%d Prev BB: %d-%d %d-%d\n", vpr_id, sink->id, buffer, sink->criticality_fac, sink->current_bounding_box.xmin, sink->current_bounding_box.xmax, sink->current_bounding_box.ymin, sink->current_bounding_box.ymax, sink->previous_bounding_box.xmin, sink->previous_bounding_box.xmax, sink->previous_bounding_box.ymin, sink->previous_bounding_box.ymax);

		const auto &sink_rr_node = get_vertex_props(g, sink->rr_node);

		route_tree_multi_root_add_to_heap(rt, g, sink->rr_node, sink->criticality_fac, params.astar_fac, sink->current_bounding_box, heap, perf);

		bool found_sink = false;
		while (!heap.empty() && !found_sink) {
			auto item = heap.top();
			heap.pop();

			if (perf) {
				++perf->num_heap_pops;
			}

			sprintf_rr_node(item.rr_node, buffer);
			const auto &v = get_vertex_props(g, item.rr_node);
			zlog_level(delta_log, ROUTER_V3, "Current: %s occ/cap: %d/%d prev: %d new_cost: %g new_known: %g new_delay: %g old_cost: %g old_known: %g old_delay: %g\n", buffer, congestion[item.rr_node].cong.occ, v.capacity, valid(item.prev_edge) ? get_source(g, item.prev_edge) : -1, item.cost, item.known_cost, item.delay, state[item.rr_node].cost, state[item.rr_node].known_cost, state[item.rr_node].delay);

			if (item.rr_node == sink->rr_node) {
				state[item.rr_node] = item;
				modified.push_back(item.rr_node);
				found_sink = true;
			} else if (item.known_cost < state[item.rr_node].known_cost && item.cost < state[item.rr_node].cost) {
				/*if (route_tree_get_rt_node(rt, item.rr_node)) {*/
					/*sprintf_rr_node(item.rr_node, buffer);*/
					/*zlog_warn(delta_log, "Warning: Found a lower cost path to existing route tree node %s\n", buffer);*/
					/*assert(false);*/
				/*}*/

				state[item.rr_node] = item;
				modified.push_back(item.rr_node);

				expand_neighbors_with_fine_grain_lock(g, item.rr_node, state, congestion, sink_rr_node, sink->criticality_fac, params.astar_fac, heap, [&g, &sink, &sink_rr_node, &existing_opin] (const RRNode &v) -> bool {

					/*if (trace_has_node(prev_trace, id(v))) {*/
						/*zlog_level(delta_log, ROUTER_V3, " existing node route tree ");*/
					/*}*/
					const auto &prop = get_vertex_props(g, v);

					if (existing_opin != RRGraph::null_vertex() && prop.type == OPIN && v != existing_opin) {
					return false;
					}

					if (prop.xhigh < sink->current_bounding_box.xmin
							|| prop.xlow > sink->current_bounding_box.xmax
							|| prop.yhigh < sink->current_bounding_box.ymin
							|| prop.ylow > sink->current_bounding_box.ymax) {
					zlog_level(delta_log, ROUTER_V3, "outside of bounding box\n");
					return false;
					}

					if (prop.type == IPIN
							&& (prop.xhigh != sink_rr_node.xhigh 
								|| prop.yhigh != sink_rr_node.yhigh)) {
					zlog_level(delta_log, ROUTER_V3, "not target IPIN\n");
					return false;
					}

					/*if (net.sinks.size() >= HIGH_FANOUT_NET_LIM) {*/
					/*if (prop.xhigh < target_x - highfanout_rlim*/
							/*|| prop.xlow > target_x + highfanout_rlim*/
							/*|| prop.yhigh < target_y - highfanout_rlim*/
							/*|| prop.ylow > target_y + highfanout_rlim) {*/
						/*return false;*/
					/*}*/
					/*}*/
					return true;
				}, lock, perf, lock_perf);
			} else {
				zlog_level(delta_log, ROUTER_V3, "\tNot expanding neighbor because known_cost %g > %g or cost %g > %g\n", item.known_cost, state[item.rr_node].known_cost, item.cost, state[item.rr_node].cost);
			}
		}

		if (!found_sink) {
			zlog_error(delta_log, "Error: Failed to find sink\n");
			unrouted_sinks.push_back(sink);
			assert(heap.empty());
		} else {
			struct s_trace *tptr, *prevptr, *temptail;

			tptr = new s_trace(); /* SINK on the end of the connection */
			tptr->index = sink->rr_node;
			tptr->iswitch = OPEN;
			tptr->next = NULL;
			temptail = tptr; /* This will become the new tail at the end */

			RouteTreeNode rt_node = route_tree_add_rr_node(rt, sink->rr_node, g);
			assert(rt_node != RouteTree::null_vertex());
			route_tree_set_node_properties(rt, rt_node, false, state[sink->rr_node].upstream_R, state[sink->rr_node].delay);
			update_one_cost_internal(sink->rr_node, g, congestion, 1, params.pres_fac, true, lock_perf);

			sprintf_rr_node(sink->rr_node, buffer);
			zlog_level(delta_log, ROUTER_V3, "Adding %s to route tree\n", buffer);

			RRNode child_rr_node = sink->rr_node;

			get_path(sink->rr_node, state, rt, g, vpr_id, [&] (const RREdge &edge) -> void
					{
					RRNode parent_rr_node = get_source(g, edge);

					const auto &parent_rr_node_p = get_vertex_props(g, parent_rr_node);

					RouteTreeNode rt_node = route_tree_add_rr_node(rt, parent_rr_node, g);
					if (rt_node != RouteTree::null_vertex()) {
						route_tree_set_node_properties(rt, rt_node, parent_rr_node_p.type != IPIN && parent_rr_node_p.type != SINK, state[parent_rr_node].upstream_R, state[parent_rr_node].delay);

						update_one_cost_internal(parent_rr_node, g, congestion, 1, params.pres_fac, true, lock_perf);
					} else if (get_vertex_props(g, parent_rr_node).type == SOURCE) {
						update_one_cost_internal(parent_rr_node, g, congestion, 1, params.pres_fac, true, lock_perf);
					}

					sprintf_rr_node(parent_rr_node, buffer);
					zlog_level(delta_log, ROUTER_V3, "Adding %s to route tree\n", buffer);

					assert(child_rr_node == get_target(g, edge));

					char buffer2[256];
					sprintf_rr_node(child_rr_node, buffer2);

					zlog_level(delta_log, ROUTER_V3, "Adding edge %s -> %s\n", buffer, buffer2);

					const RouteTreeEdge &rt_edge = route_tree_add_edge_between_rr_node(rt, parent_rr_node, child_rr_node); 
					auto &rt_edge_props = get_edge_props(rt.graph, rt_edge);
					rt_edge_props.rr_edge = edge;

					if (parent_rr_node_p.type == OPIN && existing_opin == RRGraph::null_vertex()) {
						existing_opin = parent_rr_node;
					}

					child_rr_node = parent_rr_node;

					/* linkage to vpr */
					prevptr = new s_trace();
					prevptr->index = parent_rr_node;
					prevptr->iswitch = get_edge_props(g, edge).index;
					prevptr->next = tptr;
					tptr = prevptr;

					});

			extern struct s_trace **trace_head; /* [0..(num_nets-1)] */
			extern struct s_trace **trace_tail; /* [0..(num_nets-1)] */

			if (trace_tail[vpr_id] != NULL) {
				trace_tail[vpr_id]->next = tptr; /* Traceback ends with tptr */
			} else { /* This was the first "chunk" of the net's routing */
				trace_head[vpr_id] = tptr;
			}

			trace_tail[vpr_id] = temptail;

			net_timing.delay[sink->id+1] = get_vertex_props(rt.graph, route_tree_get_rt_node(rt, sink->rr_node)).delay;

			++num_routed_sinks;
		}

		for (const auto &m : modified)  {
			state[m].known_cost = std::numeric_limits<float>::max();
			state[m].cost = std::numeric_limits<float>::max();
		}

		heap = std::priority_queue<route_state_t>();
	}

	zlog_level(delta_log, ROUTER_V1, "\n");

	return unrouted_sinks;

	/*delete [] state;*/

	/*check_route_tree(rt, net, g);*/
}
