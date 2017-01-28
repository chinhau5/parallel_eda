#include <tbb/tbb.h>
#include <metis.h>
#include <mlpack/methods/kmeans/kmeans.hpp>
#include <mlpack/methods/kmeans/refined_start.hpp>
#include "geometry.h"
#include "route.h"
#include "metis_partitioner.h"

using namespace std;

class ManhattanDistance {
	public:
	template<typename VecType1, typename VecType2>
		double Evaluate(const VecType1 &a, const VecType2 &b) const
		{
			return abs(a[0]-b[0])+abs(a[1]-b[1]);
		}
};

void cluster(const std::vector<virtual_net_t *> &virtual_nets, int num_clusters, arma::Col<size_t> &assignments)
{
	using namespace mlpack::kmeans;
	/*KMeans<mlpack::metric::EuclideanDistance, RefinedStart> k;*/

	/* the 3rd row is used for indexing purposes and not used to calculate
	 * the distance for clustering */
	arma::mat data(2, virtual_nets.size());
	for (int i = 0; i < data.n_cols; ++i) {
		data(0, i) = virtual_nets[i]->centroid.get<0>();
		data(1, i) = virtual_nets[i]->centroid.get<1>(); 
	}

	box bb = bg::make_inverse<box>();

	for (int i = 0; i < virtual_nets.size(); ++i) {
		bg::expand(bb, virtual_nets[i]->centroid);
	}

	const int max_tries = 100;
	const int max_iters = 1000;
	const double overclustering = 1;

	using Distance = ManhattanDistance;
	KMeans<Distance> k(max_iters, overclustering);

	tbb::spin_mutex lock;
	double min_metric = std::numeric_limits<double>::max();
	arma::Col<size_t> best_assignment;

	std::mt19937 mt(0);
	/*std::uniform_int_distribution<int> uni(0, num_clusters-1);*/
	assert(bb.max_corner().get<0>() >= bb.min_corner().get<0>());
	assert(bb.max_corner().get<1>() >= bb.min_corner().get<1>());
	std::uniform_int_distribution<int> uni_x(bb.min_corner().get<0>(), bb.max_corner().get<0>());
	std::uniform_int_distribution<int> uni_y(bb.min_corner().get<1>(), bb.max_corner().get<1>());

	assignments.set_size(virtual_nets.size());

	/*printf("sink bb %d-%d %d-%d\n", bl.x, tr.x, bl.y, tr.y);*/

	int max_sink_bb_area = std::numeric_limits<int>::min();

	vector<arma::mat> centroids(max_tries);
	for (int i = 0; i < max_tries; ++i) {
		centroids[i].set_size(2, num_clusters);

		for (int j = 0; j < num_clusters; ++j) {
			centroids[i](0, j) = uni_x(mt);
			centroids[i](1, j) = uni_y(mt);
			/*printf("centroid %d: %g,%g\n", j, centroids(0, j), centroids(1, j));*/
		}
	}

	tbb::parallel_for(tbb::blocked_range<int>(0, max_tries, 20),
			[&] (const tbb::blocked_range<int> &range) -> void {
			/*for (int j = 0; j < sinks.size(); ++j) {*/
			/*assignments[j] = uni(mt);*/
			/*assert(assignments[j] >= 0 && assignments[j] < num_clusters);*/
			/*}*/
			for (int i = range.begin(); i != range.end(); ++i) {
				k.Cluster(data, num_clusters, assignments, centroids[i], false, true);

				double total_d = 0;
				for (int j = 0; j < virtual_nets.size(); ++j) {
					int assign = assignments[j];
					double d = k.Metric().Evaluate(centroids[i].col(assign), data.col(j));
					total_d += d;

					/*printf("%g,%g <-> %g,%g is %g\n", data(0, j), data(1, j), centroids(0, assign), centroids(1, assign), d);*/
				}

				/*printf("total distance %g ", total_d);*/

				lock.lock();
				if (total_d < min_metric) {
					/*printf("is < %g", min_metric);*/

					min_metric = total_d;
					best_assignment = assignments;
				}
				lock.unlock();

				/*printf("\n");*/
			}
			});

	assignments = best_assignment;
	/*double total_d = 0;*/
	/*for (int j = 0; j < sinks.size(); ++j) {*/
	/*int assign = assignments[j];*/
	/*double d = k.Metric().Evaluate(centroids.col(assign), data.col(j));*/
	/*total_d += d;*/

	/*printf("%g,%g <-> %g,%g is %g\n", data(0, j), data(1, j), centroids(0, assign), centroids(1, assign), d);*/
	/*}*/
}

void partition_nets_by_clustering(vector<virtual_net_t *> &virtual_nets, int num_partitions, vector<vector<int>> &partitions, vector<bool> &has_interpartition_overlap)
{
	arma::Col<size_t> assignments;
	cluster(virtual_nets, num_partitions, assignments);

	partitions.resize(num_partitions);
	for (int i = 0; i < virtual_nets.size(); ++i) {
		assert(assignments[i] >= 0 && assignments[i] < partitions.size());
		partitions[assignments[i]].push_back(i);
	}

	has_interpartition_overlap.resize(virtual_nets.size());

	tbb::parallel_for(tbb::blocked_range<size_t>(0, virtual_nets.size(), 1024), [&has_interpartition_overlap, &assignments, &virtual_nets] (const tbb::blocked_range<size_t> &range) -> void {
			for (int i = range.begin(); i != range.end(); ++i) {
				bool has = false;
				for (int j = 0; j < virtual_nets.size() && !has; ++j) {
					if (assignments[i] != assignments[j] && bg::intersects(virtual_nets[i]->scheduler_bounding_box, virtual_nets[j]->scheduler_bounding_box)) {
						has = true;
					}
				}
				has_interpartition_overlap[i] = has;
			}
			});
}

void test_partition()
{
	virtual_net_t vnets[4];
	vnets[0].scheduler_bounding_box = box(point(0, 0), point(5, 5));
	vnets[1].scheduler_bounding_box = box(point(1, 1), point(5, 5));
	vnets[2].scheduler_bounding_box = box(point(5, 5), point(8, 8));
	vnets[3].scheduler_bounding_box = box(point(7, 6), point(9, 9));

	vector<pair<box, virtual_net_t *>> vnet_ptrs = { 
		{ vnets[0].scheduler_bounding_box, &vnets[0] },
		{ vnets[0].scheduler_bounding_box, &vnets[1] },
		{ vnets[0].scheduler_bounding_box, &vnets[2] },
		{ vnets[0].scheduler_bounding_box, &vnets[3] }
	};

	vector<vector<int>> partitions;
	vector<vector<int>> overlaps;
	vector<bool> has;
	partition_nets(vnet_ptrs, 2, 1, overlaps, partitions, has);
}

