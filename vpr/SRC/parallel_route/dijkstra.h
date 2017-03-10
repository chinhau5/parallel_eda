#ifndef DIJKSTRA_H
#define DIJKSTRA_H

template<typename Edge>
bool operator<(const existing_source_t<Edge> &a, const existing_source_t<Edge> &b)
{
	return a.distance > b.distance;
}

template<typename Graph, typename Edge, typename EdgeWeightFunc, typename ExpandCheckFunc, typename Callbacks>
void dijkstra(const Graph &g, const std::vector<existing_source_t<Edge>> &sources, int sink, float *known_distance, float *distance, Edge *prev_edge, const ExpandCheckFunc &expand_node, const EdgeWeightFunc &edge_weight, Callbacks &callbacks)
{
	using Item = existing_source_t<Edge>;
	std::priority_queue<Item> heap;

	for (const auto &s : sources) {
		heap.push(s);
	}

	bool found = false;
	while (!heap.empty() && !found) {
		Item item = heap.top();
		heap.pop();

		zlog_level(delta_log, ROUTER_V3, "Current: %d [kd=%g okd=%g] [d=%g od=%g] prev=%d\n",
					item.node, item.known_distance, known_distance[item.node], item.distance, distance[item.node], get_source(g, item.prev_edge));
		callbacks.popped_node(item.node);

		//if (!callbacks.expand_node(item.node)) {
			//continue;
		//}

		if (item.distance < distance[item.node]) {
			assert(item.known_distance <= known_distance[item.node]);

			known_distance[item.node] = item.known_distance;
			distance[item.node] = item.distance;
			prev_edge[item.node] = item.prev_edge;

			zlog_level(delta_log, ROUTER_V3, "Relaxing %d\n", item.node);
			callbacks.relax_node(item.node, item.prev_edge);

			if (item.node != sink) {
				for (const auto &e : get_out_edges(g, item.node)) {
					int v = get_target(g, e);

					zlog_level(delta_log, ROUTER_V3, "\tNeighbor: %d\n", v);

					if (!expand_node(v)) {
						continue;
					}

					const auto &weight = edge_weight(e);
					float kd = known_distance[item.node] + weight.first;
					float d = known_distance[item.node] + weight.second;

					//zlog_level(delta_log, ROUTER_V3, "\t[w1 %X w2 %X] [kd=%X okd=%X] [d=%X od=%X] [kd=%g okd=%g] [d=%g od=%g]\n",
							//*(unsigned int *)&weight.first, *(unsigned int *)&weight.second, *(unsigned int *)&kd, *(unsigned int *)&known_distance[v], *(unsigned int *)&d, *(unsigned int *)&distance[v], kd, known_distance[v], d, distance[v]);

					if (d < distance[v]) {
						assert(kd <= known_distance[v]);

						zlog_level(delta_log, ROUTER_V3, "\tPushing neighbor %d to heap\n", v);

						callbacks.push_node(v);

						heap.push({ v, kd, d, e });
					}
				}
			} else {
				found = true;
			}
		}
	}

	assert(found);
}

#endif
