#ifndef DELTA_STEPPING_H
#define DELTA_STEPPING_H

#include "cache_graph.h"

using namespace std;

using Buckets = vector<vector<int> *>;

void bucket_insert(Buckets &buckets, float delta, vector<bool> &in_bucket, int v, float distance);

bool bucket_remove(Buckets &buckets, float delta, vector<bool> &in_bucket, const vector<bool> &vertex_deleted, int v, float *distance);

int bucket_get_non_empty(const Buckets &buckets);

template<typename EdgeProperties, typename Callbacks>
void relax(Buckets &buckets, float delta, vector<bool> &in_bucket, const vector<bool> &vertex_deleted,
		float *known_distance, float *distance, cache_edge_t<EdgeProperties> *prev_edge,
		int v, float new_known_distance, float new_distance, const cache_edge_t<EdgeProperties> &edge,
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

template<typename EdgeProperties, typename EdgeWeightFunc>
bool is_light_edge(const cache_edge_t<EdgeProperties> &e, const EdgeWeightFunc &edge_weight, float delta)
{
	return edge_weight(e).second <= delta;
}

template<typename EdgeProperties>
struct existing_source_t {
	int node;
	float known_distance;
	float distance;
	cache_edge_t<EdgeProperties> prev_edge;
};

template<typename VertexProperties, typename EdgeProperties, typename EdgeWeightFunc, typename Callbacks>
void delta_stepping(cache_graph_t<VertexProperties, EdgeProperties> &g, const vector<existing_source_t<EdgeProperties>> &sources, int sink, float delta, float *known_distance, float *distance, cache_edge_t<EdgeProperties> *prev_edge, const EdgeWeightFunc &edge_weight, Callbacks &callbacks)
{
	Buckets buckets;
	vector<bool> in_bucket(num_vertices(g), false);
	vector<bool> vertex_deleted(num_vertices(g), false);

	for (const auto &s : sources) {
		relax(buckets, delta, in_bucket, vertex_deleted,
				known_distance, distance, prev_edge,
				s.node, s.known_distance, s.distance, s.prev_edge,
				callbacks);
	}

	int i;

	bool found = false;
	while ((i = bucket_get_non_empty(buckets)) >= 0 && !found) {
		vector<int> deleted_vertices;

		auto *bucket = buckets[i];

		zlog_level(delta_log, ROUTER_V3, "Bucket: %d (%d items)\n", i, bucket->size());

		zlog_level(delta_log, ROUTER_V3, "-- Relaxing light edges --\n");
		while (!bucket->empty()) {
			int u = bucket->back();
			bucket->pop_back();

			if (u == sink) {
				found = true;
			}

			assert(in_bucket[u]);
			in_bucket[u] = false;

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
