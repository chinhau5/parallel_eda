#include <stdio.h>
#include <random>
#include <zlog.h>
#include <mlpack/methods/kmeans/kmeans.hpp>
#include <mlpack/methods/kmeans/refined_start.hpp>
#include "route.h"
#include "geometry.h"

using namespace mlpack::kmeans;

void test_clustering(int num_clusters, int num_points)
{
	/*std::mt19937 mt(time(NULL));*/
	std::mt19937 mt(0);
	std::uniform_int_distribution<int> uni(0, 10);
	const int dimen = 2;
	// The dataset we are clustering.
	arma::mat data(dimen, num_points);
	for (int i = 0; i < num_points; ++i) {
		for (int j = 0; j < dimen; ++j) {
			data(j, i) = uni(mt);
			printf("%g ", data(j, i));
		}
		printf("\n");
	}

	arma::Col<size_t> assignments;

	KMeans<mlpack::metric::EuclideanDistance, RefinedStart> k;
	k.Cluster(data, num_clusters, assignments);

	FILE *output = fopen("cluster.txt", "w");
	for (int col = 0; col < data.n_cols; ++col) {
		for (int row = 0; row < data.n_rows; ++row) {
			fprintf(output, "%g ", data(row, col));	
		}
		fprintf(output, "%d\n", assignments(col));
	}
	fclose(output);
}

typedef struct metric_t {
	template<typename VecType>
		double Evaluate(const VecType &a, const VecType &b) const
		{
			return sqrt((a[0]-b[0])*(a[0]-b[0])+
					(a[1]-b[1])*(a[1]-b[1]));
		}
};

class ManhattanDistance {
	public:
	template<typename VecType1, typename VecType2>
		double Evaluate(const VecType1 &a, const VecType2 &b) const
		{
			return abs(a[0]-b[0])+abs(a[1]-b[1]);
		}
};

void cluster(const std::vector<sink_t> &sinks, int &num_clusters, arma::Col<size_t> &assignments, int sink_bb_area_threshold)
{
	/*KMeans<mlpack::metric::EuclideanDistance, RefinedStart> k;*/

	/* the 3rd row is used for indexing purposes and not used to calculate
	 * the distance for clustering */
	arma::mat data(2, sinks.size());
	for (int i = 0; i < data.n_cols; ++i) {
		data(0, i) = sinks[i].x;
		data(1, i) = sinks[i].y;
	}

	box bb = bg::make_inverse<box>();

	for (int i = 0; i < sinks.size(); ++i) {
		bg::expand(bb, point(sinks[i].x, sinks[i].y));
	}

	const int max_tries = 100;
	const int max_iters = 1000;
	const double overclustering = 1;

	using Distance = ManhattanDistance;
	KMeans<Distance> k(max_iters, overclustering);

	double min_metric = std::numeric_limits<double>::max();
	arma::Col<size_t> best_assignment;

	std::mt19937 mt(0);
	/*std::uniform_int_distribution<int> uni(0, num_clusters-1);*/
	assert(bb.max_corner().get<0>() >= bb.min_corner().get<0>());
	assert(bb.max_corner().get<1>() >= bb.min_corner().get<1>());
	std::uniform_int_distribution<int> uni_x(bb.min_corner().get<0>(), bb.max_corner().get<0>());
	std::uniform_int_distribution<int> uni_y(bb.min_corner().get<1>(), bb.max_corner().get<1>());

	assignments.set_size(sinks.size());

	/*printf("sink bb %d-%d %d-%d\n", bl.x, tr.x, bl.y, tr.y);*/

	int max_sink_bb_area = std::numeric_limits<int>::min();

	do {
		arma::mat centroids(2, num_clusters);

		for (int i = 0; i < max_tries; ++i) {
			/*for (int j = 0; j < sinks.size(); ++j) {*/
			/*assignments[j] = uni(mt);*/
			/*assert(assignments[j] >= 0 && assignments[j] < num_clusters);*/
			/*}*/
			for (int j = 0; j < num_clusters; ++j) {
				centroids(0, j) = uni_x(mt);
				centroids(1, j) = uni_y(mt);
				/*printf("centroid %d: %g,%g\n", j, centroids(0, j), centroids(1, j));*/
			}
			k.Cluster(data, num_clusters, assignments, centroids, false, true);

			double total_d = 0;
			for (int j = 0; j < sinks.size(); ++j) {
				int assign = assignments[j];
				double d = k.Metric().Evaluate(centroids.col(assign), data.col(j));
				total_d += d;

				/*printf("%g,%g <-> %g,%g is %g\n", data(0, j), data(1, j), centroids(0, assign), centroids(1, assign), d);*/
			}

			/*printf("total distance %g ", total_d);*/

			if (total_d < min_metric) {
				/*printf("is < %g", min_metric);*/

				min_metric = total_d;
				best_assignment = assignments;
			}

			/*printf("\n");*/
		}
		assignments = best_assignment;
		/*double total_d = 0;*/
		/*for (int j = 0; j < sinks.size(); ++j) {*/
			/*int assign = assignments[j];*/
			/*double d = k.Metric().Evaluate(centroids.col(assign), data.col(j));*/
			/*total_d += d;*/

			/*printf("%g,%g <-> %g,%g is %g\n", data(0, j), data(1, j), centroids(0, assign), centroids(1, assign), d);*/
		/*}*/

		max_sink_bb_area = std::numeric_limits<int>::min();
		for (int i = 0; i < num_clusters; ++i) {
			box sink_bb = bg::make_inverse<box>();
			for (int j = 0; j < sinks.size(); ++j) {
				if (assignments[j] == i) {
					bg::expand(sink_bb, point(sinks[j].x, sinks[j].y));
				}
			}
			bg::subtract_value(sink_bb.min_corner(), 1);
			max_sink_bb_area = std::max((int)bg::area(sink_bb), max_sink_bb_area);
			/*printf("cur max_sink_bb_area: %d\n", max_sink_bb_area);*/
		}

		if (max_sink_bb_area > sink_bb_area_threshold) {
			num_clusters = std::min((int)sinks.size(), num_clusters + 1);
		}
	} while (max_sink_bb_area > sink_bb_area_threshold);
}

void print_cluster(const char *filename, const vector<sink_t> &sinks, const arma::Col<size_t> &assignments)
{
	assert(sinks.size() == assignments.n_elem);
	FILE *file = fopen(filename, "w");
	for (int i = 0; i < sinks.size(); ++i) {
		fprintf(file, "%d %d %d\n", sinks[i].x, sinks[i].y, assignments(i));	
	}
	fclose(file);
}

void print_virtual_nets(const vector<net_t> &nets)
{
	FILE *file = fopen("virtual_nets.txt", "w");
	for (const auto &net : nets) {
		fprintf(file, "Net %d \n", net.local_id);
		for (const auto &virtual_net : net.virtual_nets) {
			fprintf(file, "Virtual net %d Sinks: \n", virtual_net->id);
			for (const auto &sink : virtual_net->sinks) {
				fprintf(file, "%d \n", sink->rr_node);
			}
		}
		fprintf(file, "\n");
	}
	fclose(file);
}

void create_clustered_virtual_nets(vector<net_t> &nets, int num_nodes_per_cluster, int sink_bb_area_threshold, vector<vector<virtual_net_t>> &virtual_nets)
{
	virtual_nets.resize(nets.size());

	extern zlog_category_t *delta_log;

	for (auto &net : nets) {
		vector<virtual_net_t> current_virtual_nets;
		int num_clusters = ceil(sqrt((float)net.sinks.size()));
		/*if (num_clusters > 1000000000) {*/
			/*assert(false);*/
		/*if (num_clusters > 2) {*/
			arma::Col<size_t> assignments;
			cluster(net.sinks, num_clusters, assignments, sink_bb_area_threshold);

			int num_virtual_nets = 0;
			for (int clus = 0; clus < num_clusters; ++clus) {
				virtual_net_t virtual_net;

				virtual_net.id = clus;
				virtual_net.routed = false;
				virtual_net.source = &net.source;

				/*virtual_net.centroid.x = 0;*/
				/*virtual_net.centroid.y = 0;*/
				multi_point mp;
				bg::append(mp, point(net.source.x, net.source.y));

				/*int num_sinks_in_cluster = 0;*/
				for (int sink = 0; sink < net.sinks.size(); ++sink) {
					if (assignments[sink] == clus) {
						/*virtual_net.centroid.x += net.sinks[sink].x;*/
						/*virtual_net.centroid.y += net.sinks[sink].y;*/
						virtual_net.sinks.push_back(&net.sinks[sink]);
						/*++num_sinks_in_cluster;*/
						++num_virtual_nets;
						zlog_level(delta_log, ROUTER_V3, "Net %d cluster %d sink %d %d,%d\n", net.vpr_id, clus, sink, net.sinks[sink].x, net.sinks[sink].y);

						bg::append(mp, point(net.sinks[sink].x, net.sinks[sink].y));
					}
				}
				/*virtual_net.centroid.x /= num_sinks_in_cluster;*/
				/*virtual_net.centroid.y /= num_sinks_in_cluster;*/
				/*bg::centroid(mp, virtual_net.centroid);*/
				bg::envelope(mp, virtual_net.current_bounding_box);
				/*zlog_level(delta_log, ROUTER_V3, "Net %d cluster %d centroid %d,%d\n", net.vpr_id, clus, virtual_net.centroid.get<0>(), virtual_net.centroid.get<1>());*/
				virtual_net.nearest_rr_node = -1;

				current_virtual_nets.push_back(virtual_net);
			}
			assert(num_virtual_nets == net.sinks.size());
			char filename[256];
			sprintf(filename, "/Volumes/DATA/clusters/net_%d_cluster.txt", net.vpr_id);
			print_cluster(filename, net.sinks, assignments);
		/*} else {*/
			/*for (int sink = 0; sink < net.sinks.size(); ++sink) {*/
				/*virtual_net_t virtual_net;*/

				/*virtual_net.id = sink;*/
				/*virtual_net.routed = false;*/
				/*virtual_net.source = &net.source;*/
				/*virtual_net.sinks.push_back(&net.sinks[sink]);*/
				/*[>virtual_net.centroid.set<0>(net.sinks[sink].x);<]*/
				/*[>virtual_net.centroid.set<1>(net.sinks[sink].y);<]*/
				/*virtual_net.nearest_rr_node = -1;*/

				/*current_virtual_nets.push_back(virtual_net);*/
			/*}*/
		/*}*/

		virtual_nets[net.local_id] = std::move(current_virtual_nets);
	}	

	for (int i = 0; i < nets.size(); ++i) {
		for (auto &virtual_net : virtual_nets[i]) {
			nets[i].virtual_nets.push_back(&virtual_net);
		}
	}
}

/*void */
