#ifndef DELTA_STEPPING_H
#define DELTA_STEPPING_H

#include "cache_graph.h"
#include "delta_stepping_common.h"

//using namespace std;

template<typename Edge, typename Callbacks>
void relax(Buckets &buckets, float delta, std::vector<int> &in_bucket, const std::vector<bool> &vertex_deleted,
		float *known_distance, float *distance, Edge *prev_edge,
		int v, float new_known_distance, float new_distance, const Edge &edge,
		Callbacks &callbacks)
{
	if (new_known_distance < known_distance[v] && new_distance < distance[v]) {
		zlog_level(delta_log, ROUTER_V3, "\t\tRelaxing %d [known %g -> %g dist %g -> %g bucket %d -> %d]\n", v, known_distance[v], new_known_distance, distance[v], new_distance, (int)floor(distance[v]/delta), (int)floor(new_distance/delta));

		bucket_remove(buckets, delta, in_bucket, vertex_deleted, v, distance);
		bucket_insert(buckets, delta, in_bucket, v, new_distance);

		known_distance[v] = new_known_distance;
		distance[v] = new_distance;
		prev_edge[v] = edge;

		callbacks.relax_node(v, edge);
	}
}

template<typename Edge, typename EdgeWeightFunc>
bool is_light_edge(const Edge &e, const EdgeWeightFunc &edge_weight, float delta)
{
	return edge_weight(e).second <= delta;
}

template<typename Edge, typename Extra>
struct heap_node_t {
	int node;
	float known_distance;
	float distance;
	Edge prev_edge;
	Extra extra;
};

template<typename Graph, typename Edge, typename EdgeWeightFunc, typename Callbacks, typename Extra>
void delta_stepping(const Graph &g, const std::vector<heap_node_t<Edge, Extra>> &sources, int sink, float delta, float *known_distance, float *distance, Edge *prev_edge, const EdgeWeightFunc &edge_weight, Callbacks &callbacks)
{
	Buckets buckets;
	std::vector<int> in_bucket(num_vertices(g), -1);
	std::vector<bool> vertex_deleted(num_vertices(g), false);

	for (const auto &s : sources) {
		relax(buckets, delta, in_bucket, vertex_deleted,
				known_distance, distance, prev_edge,
				s.node, s.known_distance, s.distance, s.prev_edge,
				callbacks);
	}

	int i;

	bool found = false;
	while ((i = bucket_get_non_empty(buckets)) >= 0 && !found) {
		std::vector<int> deleted_vertices;

		auto *bucket = buckets[i];

		zlog_level(delta_log, ROUTER_V3, "Bucket: %d (%d items)\n", i, bucket->size());

		zlog_level(delta_log, ROUTER_V3, "-- Relaxing light edges --\n");
		while (!bucket->empty()) {
			int u = bucket->back();
			bucket->pop_back();

			if (u == sink) {
				found = true;
			}

			assert(in_bucket[u] >= 0);
			in_bucket[u] = -1;

			if (!vertex_deleted[u]) {
				vertex_deleted[u] = true;
				deleted_vertices.push_back(u);
			}

			zlog_level(delta_log, ROUTER_V3, "Vertex: %d\n", u);

			for (const auto &e : get_out_edges(g, u)) {
				int neighbor = get_target(g, e);
				bool light = is_light_edge(e, edge_weight, delta);
				if (light) {
				}

				if (light && callbacks.expand_node(neighbor)) {
					//zlog_level(delta_log, ROUTER_V3, "\tneighbor: %d edge weight: %g,%g known_distance: %g\n", neighbor, edge_weight(e), distance[u]+edge_weight(e));
					const auto &weight = edge_weight(e);
					relax(buckets, delta, in_bucket, vertex_deleted,
							known_distance, distance, prev_edge,
							neighbor, known_distance[u]+weight.first, known_distance[u]+weight.second, e,
							callbacks);
				}
			}
		}

		zlog_level(delta_log, ROUTER_V3, "Deleted vertices: %d\n", deleted_vertices.size());

		if (!found) {
			zlog_level(delta_log, ROUTER_V3, "-- Relaxing heavy edges --\n");
			for (const auto &u : deleted_vertices) {
				zlog_level(delta_log, ROUTER_V3, "Vertex: %d\n", u);

				for (const auto &e : get_out_edges(g, u)) {
					int neighbor = get_target(g, e);
					if (!is_light_edge(e, edge_weight, delta) && callbacks.expand_node(neighbor)) {
						//zlog_level(delta_log, ROUTER_V3, "\tneighbor: %d edge weight: %g distance: %g\n", neighbor, edge_weight(e), distance[u]+edge_weight(e));
						const auto &weight = edge_weight(e);
						relax(buckets, delta, in_bucket, vertex_deleted,
								known_distance, distance, prev_edge,
								neighbor, known_distance[u]+weight.first, known_distance[u]+weight.second, e,
								callbacks);
					}
				}
			}
		}

		zlog_level(delta_log, ROUTER_V3, "\n");
	}

	assert(found);
}

#endif
