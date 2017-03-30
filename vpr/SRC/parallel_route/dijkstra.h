#ifndef DIJKSTRA_H
#define DIJKSTRA_H

template<typename Edge>
bool operator<(const heap_node_t<Edge> &a, const heap_node_t<Edge> &b)
{
	return a.distance > b.distance;
}

template<typename Graph, typename Edge, typename EdgeWeightFunc, typename ExpandCheckFunc, typename Callbacks>
void dijkstra(const Graph &g, const std::vector<heap_node_t<Edge>> &sources, int sink, float *known_distance, float *distance, Edge *prev_edge, const ExpandCheckFunc &expand_node, const EdgeWeightFunc &edge_weight, Callbacks &callbacks)
{
	using Item = heap_node_t<Edge>;
	std::priority_queue<Item> heap;

	for (const auto &s : sources) {
		heap.push(s);
	}

	bool found = false;
	while (!heap.empty() && !found) {
		Item item = heap.top();
		heap.pop();

		//unsigned int xikd, xkd, xid, xd;
		//memcpy(&xikd, &item.known_distance, sizeof (xikd));
		//memcpy(&xkd, &known_distance[item.node], sizeof (xkd));
		//memcpy(&xid, &item.distance, sizeof (xid));
		//memcpy(&xd, &distance[item.node], sizeof (xd));

		//zlog_level(delta_log, ROUTER_V3, "Current: %d [kd=%g %a %X okd=%g %a %X] [d=%g %a %X od=%g %a %X] prev=%d\n",
					//item.node,
					//item.known_distance, item.known_distance, xikd,
					//known_distance[item.node], known_distance[item.node], xkd,
					//item.distance, item.distance, xid,
					//distance[item.node], distance[item.node], xd,
					//get_source(g, item.prev_edge));

		callbacks.popped_node(item);

		//if (!callbacks.expand_node(item.node)) {
			//continue;
		//}

		if (item.node == sink) {
			found = true;

			assert(item.distance < distance[item.node]);

			known_distance[item.node] = item.known_distance;
			distance[item.node] = item.distance;
			prev_edge[item.node] = item.prev_edge;

			//zlog_level(delta_log, ROUTER_V3, "Relaxing %d\n", item.node);

			callbacks.relax_node(item);
		} else if (item.known_distance < known_distance[item.node] && item.distance < distance[item.node]) {
			//zlog_level(delta_log, ROUTER_V3, "VS %g %g\n", item.distance-item.known_distance, distance[item.node]-known_distance[item.node]);
			//assert(item.known_distance < known_distance[item.node]);

			known_distance[item.node] = item.known_distance;
			distance[item.node] = item.distance;
			prev_edge[item.node] = item.prev_edge;

			//zlog_level(delta_log, ROUTER_V3, "Relaxing %d\n", item.node);

			callbacks.relax_node(item);

			for (const auto &e : get_out_edges(g, item.node)) {
				int v = get_target(g, e);

				//zlog_level(delta_log, ROUTER_V3, "\tNeighbor: %d\n", v);

				if (expand_node(v)) {
					const auto &weight = edge_weight(e);
					float kd = known_distance[item.node] + weight.first;
					float d = known_distance[item.node] + weight.second;

					//zlog_level(delta_log, ROUTER_V3, "\t[w1 %X w2 %X] [kd=%X okd=%X] [d=%X od=%X] [kd=%g okd=%g] [d=%g od=%g]\n",
					//*(unsigned int *)&weight.first, *(unsigned int *)&weight.second, *(unsigned int *)&kd, *(unsigned int *)&known_distance[v], *(unsigned int *)&d, *(unsigned int *)&distance[v], kd, known_distance[v], d, distance[v]);

					if (kd < known_distance[v] && d < distance[v]) {
						//assert(kd <= known_distance[v]);
						//zlog_level(delta_log, ROUTER_V3, "\tPushing neighbor %d [pkd %g %a pd %g %a][edge k %g %a p %g %a][kd %g %a d %g %a] to heap\n",
								//v,
								//known_distance[item.node], known_distance[item.node],
								//distance[item.node], distance[item.node],
								//weight.first, weight.first, weight.second, weight.second,
								//kd, kd, d, d);

						Item new_node { v, kd, d, e };

						callbacks.push_node(new_node);

						heap.push(new_node);
					}
				}
			}
		} 
	}

	assert(found);
}

#endif
