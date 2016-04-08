
void route_net_fast(RRGraph &g, const net_t &net, const route_parameters_t &params, route_state_t *state, route_tree_t &rt, trace_t &trace, t_net_timing &net_timing, perf_t *perf)
{
	std::priority_queue<route_state_t> heap;

	vector<int> modified;

	zlog_level(delta_log, ROUTER_V1, "Routing VPR net id %d\n", net.vpr_id);

	if (route_tree_empty(rt)) {
		/* special case */
		/*update_one_cost_internal(get_vertex_props(g, net.current_source.rr_node), 1, params.pres_fac);*/

		/*route_tree_set_source(rt, get_vertex_props(g, net.source.rr_node));*/
	} else {
		RouteTreeNode rt_root = get_vertex_props(rt.graph, rt.root_rt_node_id);
		if (rt_root.rr_node != net.source.rr_node) {
			char root[256];
			char source[256];
			sprintf_rr_node(rt_root.rr_node, root);
			sprintf_rr_node(net.source.rr_node, source);
			zlog_error(delta_log, "Error: Non empty route tree node has root %s that is different from net source %s\n",
					root, source);
			assert(false);
		}
	}

	vector<sink_t> sorted_sinks = net.sinks;
	std::sort(sorted_sinks.begin(), sorted_sinks.end());

	char buffer[256];

	sprintf_rr_node(net.source.rr_node, buffer);
	zlog_level(delta_log, ROUTER_V1, "Source: %s\n", buffer);

	int isink = 0;
	for (const auto &sink : sorted_sinks) {
		sprintf_rr_node(sink.rr_node, buffer);
		zlog_level(delta_log, ROUTER_V1, "Sink %d: %s criticality: %g BB: %d-%d %d-%d\n", sink.id, buffer, sink.criticality_fac, sink.current_bounding_box.xmin, sink.current_bounding_box.xmax, sink.current_bounding_box.ymin, sink.current_bounding_box.ymax);
		route_tree_add_to_heap(rt, g, get_vertex_props(g, sink.rr_node), sink.criticality_fac, params.astar_fac, sink.current_bounding_box, heap, perf);

		bool found_sink = false;
		while (!heap.empty() && !found_sink) {
			auto item = heap.top();
			heap.pop();

			sprintf_rr_node(item.rr_node, buffer);
			zlog_level(delta_log, ROUTER_V2, "Current: %s prev: %d old_delay: %g new_delay: %g old_known: %g new_known: %g old_cost: %g new_cost: %g\n", buffer, item.prev_edge != -1 ? get_source(g, item.prev_edge) : -1, state[item.rr_node].delay, item.delay, state[item.rr_node].known_cost, item.known_cost, state[item.rr_node].cost, item.cost);

			if (item.rr_node == sink.rr_node) {
				state[item.rr_node] = item;
				modified.push_back(item.rr_node);
				found_sink = true;
			} else if (item.known_cost < state[item.rr_node].known_cost && item.cost < state[item.rr_node].cost) {
				/*if (route_tree_get_rt_node(rt, item.rr_node)) {*/
					/*sprintf_rr_node(item.rr_node, buffer);*/
					/*zlog_warn(delta_log, "Warning: Found a lower cost path to existing route tree node %s\n", buffer);*/
					/*assert(false);*/
				/*}*/
				const auto &sink_vertex = get_vertex_props(g, sink.rr_node);

				state[item.rr_node] = item;
				modified.push_back(item.rr_node);

				for (const auto &e_i : get_vertex_props(g, item.rr_node).edges) {
					auto &e = get_edge(g, e_i);
					auto &neighbor = get_vertex_props(g, get_target(g, e_i));

					char buffer[256];
					sprintf_rr_node(id(neighbor), buffer);
					zlog_level(delta_log, ROUTER_V3, "\tNeighbor: %s ", buffer);

					const auto &prop = neighbor.properties;

					if (prop.xhigh < sink.current_bounding_box.xmin
							|| prop.xlow > sink.current_bounding_box.xmax
							|| prop.yhigh < sink.current_bounding_box.ymin
							|| prop.ylow > sink.current_bounding_box.ymax) {
						zlog_level(delta_log, ROUTER_V3, "outside of bounding box\n");
						continue;
					}

					if (prop.type == IPIN
							&& (prop.xhigh != sink_vertex.xhigh 
								|| prop.yhigh != sink_vertex.yhigh)) {
						zlog_level(delta_log, ROUTER_V3, "not target IPIN\n");
						continue;
					}

					route_state_t new_item;

					new_item.rr_node = id(neighbor);
					new_item.prev_edge = id(e);

					float unbuffered_upstream_R = item.upstream_R;
					float upstream_R = e.R + neighbor.R;
					if (!e.buffered) {
						upstream_R += unbuffered_upstream_R;
					}
					new_item.upstream_R = upstream_R;

					new_item.delay = item.delay + get_delay(e, neighbor, unbuffered_upstream_R);

					float known_cost = item.known_cost + get_known_cost(g, e, sink.criticality_fac, unbuffered_upstream_R);
					new_item.known_cost = known_cost;

					float expected_cost = get_timing_driven_expected_cost(get_vertex_props(g, new_item.rr_node), sink_vertex, sink.criticality_fac, upstream_R);
					new_item.cost = known_cost + params.astar_fac * expected_cost;

					heap.push(new_item);

					if (perf) {
						++perf->num_heap_pushes;
					}

					zlog_level(delta_log, ROUTER_V3, "added with prev_edge: %X upstream_R: %g delay: %g known_cost: %g expected_cost: %g cost: %g\n", new_item.prev_edge, new_item.upstream_R, new_item.delay, new_item.known_cost, expected_cost, new_item.cost);
				}
			} else {
				zlog_level(delta_log, ROUTER_V3, "\tNot expanding neighbor because known_cost %g > %g or cost %g > %g\n", item.known_cost, state[item.rr_node].known_cost, item.cost, state[item.rr_node].cost);
			}
		}

		assert(found_sink);

		const auto &new_path = trace_add_path(trace, g, state, sink.rr_node, net.vpr_id);
		/* new_path.second-1 because we do not update the last node because it's the node in existing route tree (which was already updated previously */
		/*int new_path_first_rr_node = *(new_path.end()-1);*/
		/* if the first node is a connection to existing route, don't update its cost again */
		/*if (route_tree_get_rt_node(rt, new_path_first_rr_node)) {*/
			/*update_one_cost(g, new_path.begin(), new_path.end()-1, 1, params.pres_fac);*/
		/*} else {*/
			/*update_one_cost(g, new_path.begin(), new_path.end(), 1, params.pres_fac);*/
		/*}*/
		/*update_one_cost(g, new_path.begin(), new_path.end(), 1, params.pres_fac, false, nullptr);*/

		/* important to update route tree only after updating cost because the previous check relies on whether the
		 * first node is already in the route tree */
		route_tree_add_path(rt, g, state, new_path, net.vpr_id);

		net_timing.delay[sink.id+1] = route_tree_get_rt_node(rt, sink.rr_node)->delay;

		for (const auto &m : modified)  {
			state[m].known_cost = std::numeric_limits<float>::max();
			state[m].cost = std::numeric_limits<float>::max();
		}

		heap = std::priority_queue<route_state_t>();
	}

	/*delete [] state;*/

	/*check_route_tree(rt, net, g);*/
}

void route_net(RRGraph &g, int vpr_id, const source_t *source, const vector<sink_t *> &sinks, const route_parameters_t &params, route_state_t *state, route_tree_t &rt, trace_t &trace, const trace_t &prev_trace, t_net_timing &net_timing, perf_t *perf)
{
	std::priority_queue<route_state_t> heap;

	vector<int> modified;

	zlog_level(delta_log, ROUTER_V1, "Routing VPR net id %d\n", vpr_id);

	if (route_tree_empty(rt)) {
		/* special case */
		/*update_one_cost_internal(get_vertex_props(g, net.current_source.rr_node), 1, params.pres_fac);*/

		/*route_tree_set_source(rt, get_vertex_props(g, source->rr_node));*/
		auto &source_rr_node = get_vertex_props(g, source->rr_node);
		RouteTreeNode *root_rt_node = route_tree_add_rr_node(rt, source_rr_node);
		route_tree_set_node_properties(*root_rt_node, true, -1, source_rr_node.R, 0.5 * source_rr_node.R * source_rr_node.C);
		route_tree_set_root(rt, source->rr_node);
	} else {
		RouteTreeNode rt_root = get_vertex_props(rt.graph, rt.root_rt_node_id);
		if (rt_root.rr_node != source->rr_node) {
			char root[256];
			char source_str[256];
			sprintf_rr_node(rt_root.rr_node, root);
			sprintf_rr_node(source->rr_node, source_str);
			zlog_error(delta_log, "Error: Non empty route tree node has root %s that is different from net source %s\n",
					root, source_str);
			assert(false);
		}
	}

	vector<sink_t *> sorted_sinks = sinks;
	std::sort(sorted_sinks.begin(), sorted_sinks.end());

	char buffer[256];

	sprintf_rr_node(source->rr_node, buffer);
	zlog_level(delta_log, ROUTER_V1, "Source: %s\n", buffer);

	int isink = 0;
	for (const auto &sink : sorted_sinks) {
		sprintf_rr_node(sink->rr_node, buffer);
		zlog_level(delta_log, ROUTER_V1, "Net %d Sink %d: %s criticality: %g BB: %d-%d %d-%d Prev BB: %d-%d %d-%d\n", vpr_id, sink->id, buffer, sink->criticality_fac, sink->current_bounding_box.xmin, sink->current_bounding_box.xmax, sink->current_bounding_box.ymin, sink->current_bounding_box.ymax, sink->previous_bounding_box.xmin, sink->previous_bounding_box.xmax, sink->previous_bounding_box.ymin, sink->previous_bounding_box.ymax);
		route_tree_add_to_heap(rt, g, get_vertex_props(g, sink->rr_node), sink->criticality_fac, params.astar_fac, sink->current_bounding_box, heap, perf);

		bool found_sink = false;
		while (!heap.empty() && !found_sink) {
			auto item = heap.top();
			heap.pop();

			sprintf_rr_node(item.rr_node, buffer);
			const auto &v = get_vertex_props(g, item.rr_node);
			/*zlog_level(delta_log, ROUTER_V2, "Current: %s occ/cap: %d/%d prev: %d old_cost: %g new_cost: %g old_delay: %g new_delay: %g old_known: %g new_known: %g \n", buffer, congestion[item.rr_node].occ, v.capacity, item.prev_edge ? id(get_source(g, *item.prev_edge)) : -1, state[item.rr_node].cost, item.cost, state[item.rr_node].delay, item.delay, state[item.rr_node].known_cost, item.known_cost);*/

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
				const auto &sink_vertex = get_vertex_props(g, sink->rr_node);

				state[item.rr_node] = item;
				modified.push_back(item.rr_node);

				expand_neighbors(g, v, nullptr, nullptr, sink_vertex, sink->criticality_fac, params.astar_fac, heap, [&sink, &prev_trace, &sink_vertex] (const RRNode &v) -> bool {

					if (trace_has_node(prev_trace, id(v))) {
						zlog_level(delta_log, ROUTER_V3, " existing node route tree ");
					}
					const auto &prop = v.properties;

					if (prop.xhigh < sink->current_bounding_box.xmin
							|| prop.xlow > sink->current_bounding_box.xmax
							|| prop.yhigh < sink->current_bounding_box.ymin
							|| prop.ylow > sink->current_bounding_box.ymax) {
					zlog_level(delta_log, ROUTER_V3, "outside of bounding box\n");
					return false;
					}

					if (prop.type == IPIN
							&& (prop.xhigh != sink_vertex.xhigh 
								|| prop.yhigh != sink_vertex.yhigh)) {
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
				}, false, perf, nullptr);
			} else {
				zlog_level(delta_log, ROUTER_V3, "\tNot expanding neighbor because known_cost %g > %g or cost %g > %g\n", item.known_cost, state[item.rr_node].known_cost, item.cost, state[item.rr_node].cost);
			}
		}

		if (!found_sink) {
			zlog_error(delta_log, "Failed to find sink\n");
			assert(false);
		}

		const auto &new_path = trace_add_path(trace, g, state, sink->rr_node, vpr_id);
		/* new_path.second-1 because we do not update the last node because it's the node in existing route tree (which was already updated previously */
		/*int new_path_first_rr_node = *(new_path.end()-1);*/
		/* if the first node is a connection to existing route, don't update its cost again */
		/*if (route_tree_get_rt_node(rt, new_path_first_rr_node)) {*/
			/*update_one_cost(g, new_path.begin(), new_path.end()-1, 1, params.pres_fac);*/
		/*} else {*/
			/*update_one_cost(g, new_path.begin(), new_path.end(), 1, params.pres_fac);*/
		/*}*/
		/*update_one_cost(g, new_path.begin(), new_path.end(), 1, params.pres_fac, false, nullptr);*/

		/* important to update route tree only after updating cost because the previous check relies on whether the
		 * first node is already in the route tree */
		route_tree_add_path(rt, g, state, new_path, vpr_id);

		net_timing.delay[sink->id+1] = route_tree_get_rt_node(rt, sink->rr_node)->delay;

		for (const auto &m : modified)  {
			state[m].known_cost = std::numeric_limits<float>::max();
			state[m].cost = std::numeric_limits<float>::max();
		}

		heap = std::priority_queue<route_state_t>();
	}

	/*delete [] state;*/

	/*check_route_tree(rt, net, g);*/
}

void route_net(RRGraph &g, const net_t &net, const route_parameters_t &params, route_state_t *state, route_tree_t &rt, trace_t &trace, t_net_timing &net_timing, perf_t *perf)
{
	std::priority_queue<route_state_t> heap;

	vector<int> modified;

	zlog_level(delta_log, ROUTER_V1, "Routing VPR net id %d\n", net.vpr_id);

	if (route_tree_empty(rt)) {
		/* special case */
		/*update_one_cost_internal(get_vertex_props(g, net.current_source.rr_node), 1, params.pres_fac);*/

		/*route_tree_set_source(rt, get_vertex_props(g, net.source.rr_node));*/
	} else {
		RouteTreeNode rt_root = get_vertex_props(rt.graph, rt.root_rt_node_id);
		if (rt_root.rr_node != net.source.rr_node) {
			char root[256];
			char source[256];
			sprintf_rr_node(rt_root.rr_node, root);
			sprintf_rr_node(net.source.rr_node, source);
			zlog_error(delta_log, "Error: Non empty route tree node has root %s that is different from net source %s\n",
					root, source);
			assert(false);
		}
	}

	vector<sink_t> sorted_sinks = net.sinks;
	std::sort(sorted_sinks.begin(), sorted_sinks.end());

	char buffer[256];

	sprintf_rr_node(net.source.rr_node, buffer);
	zlog_level(delta_log, ROUTER_V1, "Source: %s\n", buffer);

	int isink = 0;
	for (const auto &sink : sorted_sinks) {
		sprintf_rr_node(sink.rr_node, buffer);
		zlog_level(delta_log, ROUTER_V1, "Sink %d: %s criticality: %g BB: %d-%d %d-%d\n", sink.id, buffer, sink.criticality_fac, sink.current_bounding_box.xmin, sink.current_bounding_box.xmax, sink.current_bounding_box.ymin, sink.current_bounding_box.ymax);
		route_tree_add_to_heap(rt, g, get_vertex_props(g, sink.rr_node), sink.criticality_fac, params.astar_fac, sink.current_bounding_box, heap, perf);

		bool found_sink = false;
		while (!heap.empty() && !found_sink) {
			auto item = heap.top();
			heap.pop();

			sprintf_rr_node(item.rr_node, buffer);
			zlog_level(delta_log, ROUTER_V2, "Current: %s prev: %d old_delay: %g new_delay: %g old_known: %g new_known: %g old_cost: %g new_cost: %g\n", buffer, item.prev_edge != -1 ? get_source(g, item.prev_edge) : -1, state[item.rr_node].delay, item.delay, state[item.rr_node].known_cost, item.known_cost, state[item.rr_node].cost, item.cost);

			if (item.rr_node == sink.rr_node) {
				state[item.rr_node] = item;
				modified.push_back(item.rr_node);
				found_sink = true;
			} else if (item.known_cost < state[item.rr_node].known_cost && item.cost < state[item.rr_node].cost) {
				/*if (route_tree_get_rt_node(rt, item.rr_node)) {*/
					/*sprintf_rr_node(item.rr_node, buffer);*/
					/*zlog_warn(delta_log, "Warning: Found a lower cost path to existing route tree node %s\n", buffer);*/
					/*assert(false);*/
				/*}*/
				const auto &sink_vertex = get_vertex_props(g, sink.rr_node);

				state[item.rr_node] = item;
				modified.push_back(item.rr_node);

				expand_neighbors(g, get_vertex_props(g, item.rr_node), nullptr, nullptr, sink_vertex, sink.criticality_fac, params.astar_fac, heap, [&net, &sink, &sink_vertex] (const RRNode &v) -> bool {
					const auto &prop = v.properties;

					if (prop.xhigh < sink.current_bounding_box.xmin
							|| prop.xlow > sink.current_bounding_box.xmax
							|| prop.yhigh < sink.current_bounding_box.ymin
							|| prop.ylow > sink.current_bounding_box.ymax) {
					zlog_level(delta_log, ROUTER_V3, "outside of bounding box\n");
					return false;
					}

					if (prop.type == IPIN
							&& (prop.xhigh != sink_vertex.xhigh 
								|| prop.yhigh != sink_vertex.yhigh)) {
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
				}, false, perf, nullptr);
			} else {
				zlog_level(delta_log, ROUTER_V3, "\tNot expanding neighbor because known_cost %g > %g or cost %g > %g\n", item.known_cost, state[item.rr_node].known_cost, item.cost, state[item.rr_node].cost);
			}
		}

		assert(found_sink);

		const auto &new_path = trace_add_path(trace, g, state, sink.rr_node, net.vpr_id);
		/* new_path.second-1 because we do not update the last node because it's the node in existing route tree (which was already updated previously */
		/*int new_path_first_rr_node = *(new_path.end()-1);*/
		/* if the first node is a connection to existing route, don't update its cost again */
		/*if (route_tree_get_rt_node(rt, new_path_first_rr_node)) {*/
			/*update_one_cost(g, new_path.begin(), new_path.end()-1, 1, params.pres_fac);*/
		/*} else {*/
			/*update_one_cost(g, new_path.begin(), new_path.end(), 1, params.pres_fac);*/
		/*}*/
		/*update_one_cost(g, new_path.begin(), new_path.end(), 1, params.pres_fac, false, nullptr);*/

		/* important to update route tree only after updating cost because the previous check relies on whether the
		 * first node is already in the route tree */
		route_tree_add_path(rt, g, state, new_path, net.vpr_id);

		net_timing.delay[sink.id+1] = route_tree_get_rt_node(rt, sink.rr_node)->delay;

		for (const auto &m : modified)  {
			state[m].known_cost = std::numeric_limits<float>::max();
			state[m].cost = std::numeric_limits<float>::max();
		}

		heap = std::priority_queue<route_state_t>();
	}

	/*delete [] state;*/

	/*check_route_tree(rt, net, g);*/
}
