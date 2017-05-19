#ifndef PARALLEL_DELTA_STEPPING_H
#define PARALLEL_DELTA_STEPPING_H

#include <tbb/tbb.h>
#include "cache_graph.h"
#include "delta_stepping_common.h"

//using namespace std;

template<typename Edge, typename Callbacks>
void relax(int v, float new_known_distance, float new_distance, const Edge &edge,
		Callbacks &callbacks,
		delta_stepping_t &ds,
		float *known_distance, float *distance, Edge *prev_edge)
{
	DEBUG_PRINTS("\t\tRelaxing %d [known %g -> %g dist %g -> %g bucket %d -> %d]\n", v, known_distance[v], new_known_distance, distance[v], new_distance, (int)floor(distance[v]/ds.delta), (int)floor(new_distance/ds.delta));

	if (new_known_distance < known_distance[v] && new_distance < distance[v]) {
		bucket_remove(ds, v, distance);
		bucket_insert(ds, v, new_distance);

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

struct convert {
	std::pair<int, int> operator()(int u)
	{
		return std::make_pair(u, 0);
	}
};

template<typename Graph, typename Edge, typename ExpandCheckFunc, typename EdgeWeightFunc, typename BucketChanged, typename Callbacks>
void relax_neighbors(const Graph &g, int u, float *known_distance, float *distance, Edge *prev_edge, const ExpandCheckFunc &expand_node, const EdgeWeightFunc &edge_weight, Callbacks &callbacks, const BucketChanged &bucket_changed, delta_stepping_t &ds, tbb::spin_mutex *locks)
{
	for (const auto &e : get_out_edges(g, u)) {
		bool light = is_light_edge(e, edge_weight, ds.delta);

		int v = get_target(g, e);

		if (light && expand_node(v)) {
			//DEBUG_PRINT("\tneighbor: %d edge weight: %g,%g known_distance: %g\n", v, edge_weight(e), distance[u]+edge_weight(e));
			const auto &weight = edge_weight(e);

			locks[v].lock();

			int old_bucket = ds.in_bucket[v];

			relax(v, known_distance[u]+weight.first, known_distance[u]+weight.second, e,
					callbacks,
					ds,
					known_distance, distance, prev_edge);

			int new_bucket = ds.in_bucket[v];

			locks[v].unlock();

			if (old_bucket != new_bucket) {
				bucket_changed(v, new_bucket);
			}

			//if (ds.in_bucket[v] == i) {
				//feeder.add(v);
			//}
		}
	}
}

template<typename Graph, typename Edge, typename ExpandCheckFunc, typename EdgeWeightFunc, typename Callbacks, typename Extra>
void delta_stepping(const Graph &g, const std::vector<heap_node_t<Edge, Extra>> &sources, int sink, float delta, float *known_distance, float *distance, Edge *prev_edge, const ExpandCheckFunc &expand_node, const EdgeWeightFunc &edge_weight, Callbacks &callbacks)
{
	delta_stepping_t ds;

	ds.in_bucket.resize(num_vertices(g), -1);
	ds.vertex_deleted.resize(num_vertices(g), false);
	ds.delta = delta;

	for (const auto &s : sources) {
		relax(s.node, s.known_distance, s.distance, s.prev_edge,
				callbacks,
				ds,
				known_distance, distance, prev_edge);
	}

	int i;

	tbb::spin_mutex *locks = new tbb::spin_mutex[num_vertices(g)];

	bool found = false;
	while ((i = bucket_get_non_empty(ds.buckets)) >= 0 && !found) {
		std::vector<int> deleted_vertices;

		auto *bucket = ds.buckets[i];

		DEBUG_PRINTS("Bucket: %d (%d items)\n", i, bucket->size());

		DEBUG_PRINT("-- Relaxing light edges --\n");
		tbb::parallel_do(std::begin(*bucket), std::end(*bucket),
				[&] (int u, tbb::parallel_do_feeder<int> &feeder) -> void {
					DEBUG_PRINTS("Vertex: %d\n", u);

					//if (u == sink) {
						//found = true;
					//}

					locks[u].lock();

					assert(ds.in_bucket[u] == i);
					/* mark as removed without actually removing it */
					ds.in_bucket[u] = -1;

					if (!ds.vertex_deleted[u]) {
						ds.vertex_deleted[u] = true;
						deleted_vertices.push_back(u);
					}

					locks[u].unlock();

					relax_neighbors(g, u, known_distance, distance, prev_edge, expand_node, edge_weight, callbacks, [&] (int v, int new_bucket) -> void { if (new_bucket == i) { feeder.add(v); } }, ds, locks);
				});

		bucket->clear();

		DEBUG_PRINTS("Deleted vertices: %d\n", deleted_vertices.size());

		//if (!found) {
			DEBUG_PRINT("-- Relaxing heavy edges --\n");
			tbb::parallel_do(std::begin(deleted_vertices), std::end(deleted_vertices),
				[&] (int u, tbb::parallel_do_feeder<int> &feeder) -> void {
					DEBUG_PRINTS("Vertex: %d\n", u);

					relax_neighbors(g, u, known_distance, distance, prev_edge, expand_node, edge_weight, callbacks, [] (int v, int new_bucket) -> void { }, ds, locks);
				});
		//}

		DEBUG_PRINT("\n");

		for (const auto &d : deleted_vertices) {
			ds.vertex_deleted[d] = false;
		}
	}

	//assert(found);
}

#endif
