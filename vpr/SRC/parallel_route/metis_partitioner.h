#ifndef PARTITION_H
#define PARTITION_H

#include "route.h"
#include "metis.h"

template<typename Graph>
void partition_graph(const Graph &g, int num_partitions, float ubvec, std::vector<int> &_part)
{
	idx_t *adjncy, *xadj;
	xadj = new idx_t[num_vertices(g)+1];
	adjncy = new idx_t[2*num_edges(g)];

	std::vector<std::vector<int>> redundant_edges(num_vertices(g));
	for (const auto &e : get_edges(g)) {
		int from = get_source(g, e);
		int to = get_target(g, e);
	
		redundant_edges[to].push_back(from);
	}

	int edge = 0;
	for (int i = 0; i < num_vertices(g); ++i) {
		xadj[i] = edge;
		const auto &v = get_vertex_props(g, i);
		for (const auto &e : get_out_edges(g, i)) {
			adjncy[edge] = get_target(g, e);
			++edge;
		}
		for (const auto &v : redundant_edges[i]) {
			adjncy[edge] = v;
			++edge;
		}
		xadj[i+1] = edge;
	}
	assert(edge == num_edges(g)*2);

	idx_t nvtxs = num_vertices(g);
	idx_t ncon = 1;
	idx_t nparts = num_partitions;
	//real_t ubvec = 2;
	idx_t objval;
	idx_t *part = _part.data();
	idx_t options[METIS_NOPTIONS];
	idx_t *vwgt = new idx_t[nvtxs];

	for (int i = 0; i < nvtxs; ++i) {
		vwgt[i] = get_vertex_props(g, i).weight;
	}
	METIS_SetDefaultOptions(options);
	//options[METIS_OPTION_DBGLVL] = METIS_DBG_INFO;
	options[METIS_OPTION_NUMBERING] = 0;
	assert(METIS_PartGraphRecursive(&nvtxs, &ncon, xadj, adjncy, vwgt, NULL, NULL, &nparts, NULL, &ubvec, options, &objval, part) == METIS_OK);
	//assert(METIS_PartGraphRecursive(&nvtxs, &ncon, xadj, adjncy, vwgt, NULL, NULL, &nparts, NULL, NULL, options, &objval, part) == METIS_OK);
//idx t *nvtxs, idx t *ncon, idx t *xadj, idx t *adjncy, idx t *vwgt, idx t *vsize, idx t *adjwgt, idx t *nparts, real t *tpwgts, real t ubvec, idx t *options, idx t *objval, idx t *part
	printf("ubvec: %g edgecut: %d num_edges: %lu percentage: %g\n", ubvec, objval, num_edges(g), (float)objval/num_edges(g)*100);
	//partitions.resize(num_partitions);
	//for (int i = 0; i < num_partitions; ++i) {
		//assert(partitions[i].empty());
	//}
	//for (int i = 0; i < virtual_nets.size(); ++i) {
		//assert(part[i] >= 0 && part[i] < partitions.size());
		//partitions[part[i]].push_back(i);
	//}

	delete [] xadj;
	delete [] adjncy;
	delete [] vwgt;
}

template<typename Net>
void partition_nets(std::vector<std::pair<box, Net>> &virtual_nets, int num_partitions, float ubvec, std::vector<std::vector<int>> &overlaps, std::vector<std::vector<int>> &partitions, std::vector<bool> &has_interpartition_overlap)
{
	overlaps.resize(virtual_nets.size());
	tbb::atomic<unsigned long> num_edges = 0;
	tbb::parallel_for(tbb::blocked_range<size_t>(0, virtual_nets.size(), 1024), [&virtual_nets, &overlaps, &num_edges] (const tbb::blocked_range<size_t> &range) -> void {
			for (int i = range.begin(); i != range.end(); ++i) {
			for (int j = 0; j < virtual_nets.size(); ++j) {
			if (i != j && bg::intersects(virtual_nets[i].first, virtual_nets[j].first)) {
			overlaps[i].push_back(j);
			++num_edges;
			}
			}
			}
			});

	idx_t *adjncy, *xadj;
	xadj = new idx_t[virtual_nets.size()+1];
	adjncy = new idx_t[2*num_edges];

	int edge = 0;
	for (int i = 0; i < virtual_nets.size(); ++i) {
		xadj[i] = edge;
		xadj[i+1] = edge + overlaps[i].size();
		for (int j = 0; j < overlaps[i].size(); ++j) {
			adjncy[edge+j] = overlaps[i][j];
		}
		edge = xadj[i+1];
	}

	idx_t nvtxs = virtual_nets.size();
	idx_t ncon = 1;
	idx_t nparts = num_partitions;
	//real_t ubvec = 2;
	idx_t objval;
	idx_t *part = new idx_t[virtual_nets.size()];
	idx_t options[METIS_NOPTIONS];
	idx_t *vwgt = new idx_t[nvtxs];

	for (int i = 0; i < nvtxs; ++i) {
		int area = bg::area(virtual_nets[i].first);
		assert(area > 0);
		vwgt[i] = area;
	}
	METIS_SetDefaultOptions(options);
	//options[METIS_OPTION_DBGLVL] = METIS_DBG_INFO;
	options[METIS_OPTION_NUMBERING] = 0;
	assert(METIS_PartGraphRecursive(&nvtxs, &ncon, xadj, adjncy, vwgt, NULL, NULL, &nparts, NULL, &ubvec, options, &objval, part) == METIS_OK);
	//assert(METIS_PartGraphRecursive(&nvtxs, &ncon, xadj, adjncy, vwgt, NULL, NULL, &nparts, NULL, NULL, options, &objval, part) == METIS_OK);
//idx t *nvtxs, idx t *ncon, idx t *xadj, idx t *adjncy, idx t *vwgt, idx t *vsize, idx t *adjwgt, idx t *nparts, real t *tpwgts, real t ubvec, idx t *options, idx t *objval, idx t *part
	printf("ubvec: %g edgecut: %d num_edges: %lu percentage: %g\n", ubvec, objval, num_edges, (float)objval/num_edges*100);
	partitions.resize(num_partitions);
	for (int i = 0; i < num_partitions; ++i) {
		assert(partitions[i].empty());
	}
	for (int i = 0; i < virtual_nets.size(); ++i) {
		assert(part[i] >= 0 && part[i] < partitions.size());
		partitions[part[i]].push_back(i);
	}

	has_interpartition_overlap.resize(virtual_nets.size());

	tbb::parallel_for(tbb::blocked_range<size_t>(0, virtual_nets.size(), 1024), [&has_interpartition_overlap, &overlaps, &part] (const tbb::blocked_range<size_t> &range) -> void {
			for (int i = range.begin(); i != range.end(); ++i) {
				bool has = std::any_of(begin(overlaps[i]), end(overlaps[i]), [&part, &i] (int other) -> bool {
							return part[i] != part[other];
						});
				has_interpartition_overlap[i] = has;
			}
			});

	delete [] xadj;
	delete [] adjncy;
	delete [] part;
	delete [] vwgt;
}

template<typename Net>
void partition_nets_overlap_area_metric(std::vector<std::pair<box, Net>> &virtual_nets, int num_partitions, float ubvec, std::vector<std::vector<std::pair<int, int>>> &overlaps, std::vector<std::vector<int>> &partitions, std::vector<bool> &has_interpartition_overlap)
{
	overlaps.resize(virtual_nets.size());
	tbb::atomic<int> num_edges = 0;
	tbb::parallel_for(tbb::blocked_range<size_t>(0, virtual_nets.size(), 1024), [&virtual_nets, &overlaps, &num_edges] (const tbb::blocked_range<size_t> &range) -> void {
			for (int i = range.begin(); i != range.end(); ++i) {
			for (int j = 0; j < virtual_nets.size(); ++j) {
			box intersection;
			bg::intersection(virtual_nets[i].first, virtual_nets[j].first, intersection);
			bg::add_value(intersection.max_corner(), 1);
			int overlap_area = bg::area(intersection);
			if (i != j && bg::intersects(virtual_nets[i].first, virtual_nets[j].first)) {
			overlaps[i].push_back(std::make_pair(j, overlap_area));
			++num_edges;
			}
			}
			}
			});

	idx_t *adjncy, *xadj, *adjwgt;
	xadj = new idx_t[virtual_nets.size()+1];
	adjncy = new idx_t[2*num_edges];
	adjwgt = new idx_t[2*num_edges];

	int edge = 0;
	for (int i = 0; i < virtual_nets.size(); ++i) {
		xadj[i] = edge;
		xadj[i+1] = edge + overlaps[i].size();
		for (int j = 0; j < overlaps[i].size(); ++j) {
			adjncy[edge+j] = overlaps[i][j].first;
			adjwgt[edge+j] = overlaps[i][j].second;
		}
		edge = xadj[i+1];
	}

	idx_t nvtxs = virtual_nets.size();
	idx_t ncon = 1;
	idx_t nparts = num_partitions;
	//real_t ubvec = 2;
	idx_t objval;
	idx_t *part = new idx_t[virtual_nets.size()];
	idx_t options[METIS_NOPTIONS];
	idx_t *vwgt = new idx_t[nvtxs];

	for (int i = 0; i < nvtxs; ++i) {
		int area = bg::area(virtual_nets[i].first);
		assert(area > 0);
		vwgt[i] = area;
	}
	METIS_SetDefaultOptions(options);
	//options[METIS_OPTION_DBGLVL] = METIS_DBG_INFO;
	options[METIS_OPTION_NUMBERING] = 0;
	assert(METIS_PartGraphRecursive(&nvtxs, &ncon, xadj, adjncy, vwgt, NULL, adjwgt, &nparts, NULL, &ubvec, options, &objval, part) == METIS_OK);
	//assert(METIS_PartGraphRecursive(&nvtxs, &ncon, xadj, adjncy, vwgt, NULL, NULL, &nparts, NULL, NULL, options, &objval, part) == METIS_OK);
//idx t *nvtxs, idx t *ncon, idx t *xadj, idx t *adjncy, idx t *vwgt, idx t *vsize, idx t *adjwgt, idx t *nparts, real t *tpwgts, real t ubvec, idx t *options, idx t *objval, idx t *part
	partitions.resize(num_partitions);
	for (int i = 0; i < num_partitions; ++i) {
		assert(partitions[i].empty());
	}
	for (int i = 0; i < virtual_nets.size(); ++i) {
		assert(part[i] >= 0 && part[i] < partitions.size());
		partitions[part[i]].push_back(i);
	}

	has_interpartition_overlap.resize(virtual_nets.size());

	tbb::parallel_for(tbb::blocked_range<size_t>(0, virtual_nets.size(), 1024), [&has_interpartition_overlap, &overlaps, &part] (const tbb::blocked_range<size_t> &range) -> void {
			for (int i = range.begin(); i != range.end(); ++i) {
				bool has = std::any_of(begin(overlaps[i]), end(overlaps[i]), [&part, &i] (const std::pair<int, int> &other) -> bool {
							return part[i] != part[other.first];
						});
				has_interpartition_overlap[i] = has;
			}
			});

	delete [] xadj;
	delete [] adjncy;
	delete [] adjwgt;
	delete [] part;
	delete [] vwgt;
}

void partition_nets_by_clustering(std::vector<virtual_net_t *> &virtual_nets, int num_partitions, std::vector<std::vector<int>> &partitions, std::vector<bool> &has_interpartition_overlap);

#endif
