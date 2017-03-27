#include "pch.h"

#include <mpi.h>

#include "vpr_types.h"
#include "path_delay.h"
#include "rr_graph.h"

#include "log.h"
#include "barrier.h"
#include "graph.h"
#include "fast_graph.h"
#include "filtered_graph.h"
#include "route.h"
#include "route_tree.h"
#include "trace.h"
#include "old_misr.h"
#include "geometry.h"
#include "quadtree.h"
#include "utility.h"
#include "net_cluster.h"
#include "args.h"
#include "init.h"
#include "router.h"
#include "congestion.h"
#include "metis_partitioner.h"
#include "fm.h"
#include "rr_graph_partitioner.h"

using namespace std;

void get_stats(const RRGraph &g)
{
	for (const auto &v : get_vertices(g)) {
		const auto &ver = get_vertex_props(g, v);
		if (ver.type == CHANX || ver.type == CHANY) {
		}
	}
}

template<typename VertexProperties, typename EdgeProperties, typename EdgeFilter>
void custom_minimum_spanning_tree(const graph_t<VertexProperties, EdgeProperties> &g, const EdgeFilter &edge_valid, vector<int> &edges)
{
	vector<bool> in_tree(num_vertices(g), false);

	for (const auto &e : get_edges(g)) {
		if (!edge_valid(e)) {
			continue;
		}

		int from = get_source(g, e);
		int to = get_target(g, e);
		if (!in_tree[from] || !in_tree[to]) {
			edges.push_back(e);
			in_tree[from] = true;
			in_tree[to] = true;
		}
	}

	//assert(all_of(begin(in_tree), end(in_tree), [] (bool val) -> bool { return val; }));
}

void dump_rr_graph(const RRGraph &g, const char *filename)
{
	FILE *dump = fopen(filename, "w");
	char buffer[256];
	for (const auto &v : get_vertices(g)) {
		sprintf_rr_node(v, buffer);

		const auto &rr_node_p = get_vertex_props(g, v);

		fprintf(dump, "%s Capacity: %d Num out edges: %d\n", buffer, rr_node_p.capacity, num_out_edges(g, v));

		for (const auto &e : get_out_edges(g, v)) {
			int to = get_target(g, e);

			sprintf_rr_node(v, buffer);
			fprintf(dump, "\t%s -> ", buffer);

			sprintf_rr_node(to, buffer);
			fprintf(dump, "%s\n", buffer);
		}
	}
	fclose(dump);
}

static void dump_edges(const RRGraph &g, const char *filename)
{
	FILE *dump = fopen(filename, "w");
	char buffer[256];
	for (const auto &e : get_edges(g)) {
		int from = get_source(g, e);
		int to = get_target(g, e);

		sprintf_rr_node(from, buffer);
		fprintf(dump, "%s -> ", buffer);

		sprintf_rr_node(to, buffer);
		fprintf(dump, "%s\n", buffer);
	}
	fclose(dump);
}

int get_direction(t_rr_type type, bool inc_direction)
{
	int dir;
	if (type == CHANY) {
		if (inc_direction) {
			dir = 0;
		} else {
			dir = 2;
		}
	} else {
		assert(type == CHANX);

		if (inc_direction) {
			dir = 1;
		} else {
			dir = 3;
		}
	}
	return dir;
}

bool starts_at(const rr_node_property_t &rr_node, int x, int y)
{
	bool starts;
	assert(rr_node.type == CHANX || rr_node.type == CHANY);

	if (rr_node.inc_direction) {
		starts = rr_node.xlow == x && rr_node.ylow == y;
	} else {
		starts = rr_node.xhigh == x && rr_node.yhigh == y;
	}
	return starts;
}

void partition_graph(const RRGraph &g, int num_partitions, const vector<int> &vertex_weights, float ubvec, vector<int> &partition)
{
	idx_t *adjncy, *xadj;
	xadj = new idx_t[num_vertices(g)+1];
	adjncy = new idx_t[num_edges(g)*2];

	vector<vector<int>> redundant_edges(num_vertices(g));
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
	idx_t *local_partition = new idx_t[num_vertices(g)];
	idx_t options[METIS_NOPTIONS];
	idx_t *vwgt = new idx_t[nvtxs];

	assert(vertex_weights.size() == nvtxs);
	for (int i = 0; i < nvtxs; ++i) {
		//vwgt[i] = num_out_edges(g, get_vertex_props(g, i));
		vwgt[i] = vertex_weights[i];
	}
	METIS_SetDefaultOptions(options);
	options[METIS_OPTION_DBGLVL] = 511;
	options[METIS_OPTION_NUMBERING] = 0;
	assert(METIS_PartGraphRecursive(&nvtxs, &ncon, xadj, adjncy, vwgt, NULL, NULL, &nparts, NULL, &ubvec, options, &objval, local_partition) == METIS_OK);
	//assert(METIS_PartGraphRecursive(&nvtxs, &ncon, xadj, adjncy, vwgt, NULL, NULL, &nparts, NULL, NULL, options, &objval, part) == METIS_OK);
//idx t *nvtxs, idx t *ncon, idx t *xadj, idx t *adjncy, idx t *vwgt, idx t *vsize, idx t *adjwgt, idx t *nparts, real t *tpwgts, real t ubvec, idx t *options, idx t *objval, idx t *part
	unsigned long num_e = num_edges(g);
	printf("ubvec: %g edgecut: %d num_edges: %lu percentage: %g\n", ubvec, objval, num_e, (float)objval/num_e*100);

	partition.resize(num_vertices(g));

	for (int i = 0; i < num_vertices(g); ++i) {
		assert(local_partition[i] >= 0 && local_partition[i] < partition.size());
		partition[i] = local_partition[i];
	}

	//has_interpartition_overlap.resize(virtual_nets.size());

	//tbb::parallel_for(tbb::blocked_range<size_t>(0, virtual_nets.size(), 1024), [&has_interpartition_overlap, &overlaps, &part] (const tbb::blocked_range<size_t> &range) -> void {
			//for (int i = range.begin(); i != range.end(); ++i) {
				//bool has = std::any_of(begin(overlaps[i]), end(overlaps[i]), [&part, &i] (int other) -> bool {
							//return part[i] != part[other];
						//});
				//has_interpartition_overlap[i] = has;
			//}
			//});

	delete [] xadj;
	delete [] adjncy;
	delete [] local_partition;
	//delete [] vwgt;
}

template<typename Graph>
void topological_sort(const Graph &g, vector<int> &sorted)
{
	vector<int> num_incoming_edges(num_vertices(g), 0);
	for (const auto &e : get_edges(g)) {
		int to = get_target(g, e);
		++num_incoming_edges[to];
	}
	queue<int> s;
	for (int i = 0; i < num_incoming_edges.size(); ++i) {
		if (num_incoming_edges[i] == 0) {
			s.push(i);
		}
	}
	while (!s.empty()) {
		int item = s.front(); s.pop();
		sorted.push_back(item);
		for (const auto &e : get_out_edges(g, item)) {
			int to = get_target(g, e);
			if (--num_incoming_edges[to] == 0) {
				s.push(to);
			}
		}
	}
	assert(all_of(begin(num_incoming_edges), end(num_incoming_edges), [] (int val) -> bool { return val == 0; }));
}

void test_partition_graph()
{
	vector<int> partitions;
	RRGraph g;
	add_vertex(g, 3);
	add_edge(g, 0, 1);
	add_edge(g, 0, 2);
	vector<int> weights = { 4, 1, 1 };
	partition_graph(g, 2, weights, 1.001, partitions);
	for (int i = 0; i < num_vertices(g); ++i) {
		printf("v %d part %d\n", i, partitions[i]);
	}
}

void test_topo()
{
	graph_t<int, int> g;
	add_vertex(g, 5);
	add_edge(g, 0, 2);
	add_edge(g, 1, 2);
	add_edge(g, 2, 3);
	add_edge(g, 2, 4);
	add_edge(g, 1, 4);
	vector<int> sorted;
	topological_sort(g, sorted);
	assert(sorted[0] == 0);
	assert(sorted[1] == 1);
	assert(sorted[2] == 2);
	assert(sorted[3] == 3);
	assert(sorted[4] == 4);
}

void test_fast_graph()
{
	fast_graph_t<int, int> fg;
	add_vertex(fg, 3);
	int num_edges = 0;
	auto b = begin(get_edges(fg));
	auto e = end(get_edges(fg));
	for (auto iter = b; iter != e; ++iter) {
	//for (const auto &e : get_edges(fg)) {
		++num_edges;
	}
	assert(num_edges == 0);
	add_edge(fg, 0, 1);
	add_edge(fg, 0, 2);
	vector<fast_edge_t<int>> edges;
	for (const auto &e : get_edges(fg)) {
		edges.push_back(e);
	}
	assert(edges.size() == 2);
	assert(get_source(fg, edges[0]) == 0);
	assert(get_target(fg, edges[0]) == 1);
	assert(get_source(fg, edges[1]) == 0);
	assert(get_target(fg, edges[1]) == 2);
	auto e1 = get_edge(fg, 0, 1);
	assert(get_source(fg, e1) == 0);
	assert(get_target(fg, e1) == 1);
	assert(has_edge(fg, 0, 1));
	assert(has_edge(fg, 0, 2));
	assert(!has_edge(fg, 1, 2));

	auto out = get_out_edges(fg, 0);
	auto outb = begin(out);
	auto oute = end(out);
	assert(get_source(fg, *outb) == 0);
	assert(get_target(fg, *outb) == 1);
	assert(get_source(fg, *(outb+1)) == 0);
	assert(get_target(fg, *(outb+1)) == 2);
	assert(outb+2 == oute);

	get_edge_props(fg, *outb) = 5;
	get_vertex_props(fg, 0) = 4;
}

void test_fm()
{
	struct imbalance_t {
		bool is_imbalanced(int v, int to, int &imba) const
		{
			imba = 0;
			return false;
		}
		
		void move(int v, int to)
		{
		}

		void add(int v, int to)
		{
		}
	};
	RRGraph g;
	add_vertex(g, 4);

	add_edge(g, 0, 1);
	add_edge(g, 1, 0);

	add_edge(g, 0, 2);
	add_edge(g, 2, 0);

	fm<RRGraph, imbalance_t> f;
	imbalance_t imba;
	vector<int> initial_partition = { 0, 0, 0, 0 };
	f.init(g, initial_partition, imba);
	f.run();
}

//template<typename ValidVertex2>
//struct connected_components_visitor_t {
	//vector<int> visited_edges;
	////vector<int> visited_nodes;
	//const ValidVertex2 &valid_vertex;

	//connected_components_visitor_t(const ValidVertex2 &valid_vertex) :
		//valid_vertex(valid_vertex)
	//{
	//}

	//void tree_edge(int e, const RRGraph &g)
	//{
		//if (valid_vertex(get_vertex_props(g, get_source(g, e))) && valid_vertex(get_vertex_props(g, get_target(g, e)))) {
			//visited_edges.push_back(e);
		//}
	//}

	//void examine_edge(int e, const RRGraph &g)
	//{
	//}

	//void discover_vertex(int v, const RRGraph &g)
	//{
		////if (valid_vertex(get_vertex_props(g, v))) {
			////visited_nodes.push_back(v);
		////}
	//}

	//void examine_vertex(int v, const RRGraph &g)
	//{
	//}
//};

struct connected_components_visitor_t {
	vector<RREdge> visited_edges;
	vector<int> visited_nodes;

	template<typename Graph>
	bool tree_edge(const RREdge &e, const Graph &g)
	{
		assert(find(begin(visited_edges), end(visited_edges), e) == end(visited_edges));
		visited_edges.push_back(e);
		return true;
	}

	template<typename Graph>
	void examine_edge(const RREdge &e, const Graph &g)
	{
	}

	template<typename Graph>
	bool discover_vertex(int v, const Graph &g)
	{
		assert(find(begin(visited_nodes), end(visited_nodes), v) == end(visited_nodes));
		visited_nodes.push_back(v);
		return true;
	}

	template<typename Graph>
	void examine_vertex(int v, const Graph &g)
	{
	}
};

template<template<typename, typename> class Graph, typename VertexProperties, typename EdgeProperties>
void connected_components_convert_to_undirected(const Graph<VertexProperties, EdgeProperties> &in, graph_t<VertexProperties, EdgeProperties> &out)
{
	add_vertex(out, num_vertices(in));
	for (const auto &v : get_vertices(in)) {
		get_vertex_props(out, v)= get_vertex_props(in, v);
	}

	for (const auto &e : get_edges(in)) {
		int from = get_source(in, e);
		int to = get_target(in, e);
		/* add reverse edge */
		get_edge_props(out, add_edge(out, from, to)) = get_edge_props(in, e);
		get_edge_props(out, add_edge(out, to, from)) = get_edge_props(in, e);
	}
}

template<typename Graph>
void connected_components(const Graph &g, vector<vector<int>> &components)
{
	vector<VertexColor> color(num_vertices(g), VertexColor::WHITE);

	for (const auto &v : get_vertices(g)) {
		if (color[v] == VertexColor::WHITE) {
			connected_components_visitor_t visitor;
			bfs(g, { v }, color, visitor);
			assert(visitor.visited_nodes.size() == visitor.visited_edges.size()+1);

			components.emplace_back(visitor.visited_nodes);
		}
	}
}

//template<typename VertexProperties, typename EdgeProperties, typename ValidVertex>
//void connected_components(const graph_t<VertexProperties, EdgeProperties> &g, const ValidVertex &valid_vertex, vector<vector<int>> &components)
//{
	//graph_t<VertexProperties, EdgeProperties> temp_g = g;
	////add_vertex(temp_g, num_vertices(g));
	//for (const auto &e : get_edges(g)) {
		//int from = get_source(g, e);
		//int to = get_target(g, e);
		//[> add reverse edge <]
		//add_edge(temp_g, to, from);
	//}

	//vector<VertexColor> color(num_vertices(temp_g), VertexColor::WHITE);

	//for (const auto &v : get_vertices(temp_g)) {
		//if (valid_vertex(get_vertex_props(g, v)) && color[v] == VertexColor::WHITE) {

			//connected_components_visitor_t<ValidVertex> visitor(valid_vertex);
			//bfs(temp_g, { v }, color, visitor);
			//set<int> visited_nodes;
			//for (const auto &e : visitor.visited_edges) {
				//visited_nodes.insert(get_source(temp_g, e));
				//visited_nodes.insert(get_target(temp_g, e));
			//}
			//assert(visited_nodes.size() == visitor.visited_edges.size()+1);

			//components.emplace_back(vector<int>(begin(visited_nodes), end(visited_nodes)));
		//}
	//}
//}

void test_filter_graph()
{
	//graph_t<int, int> g;

	//add_vertex(g, 4);

	//add_edge(g, 0, 1);
	//add_edge(g, 0, 2);
	//add_edge(g, 1, 3);

	//[> no filter <]
	//auto fg1 = make_filtered_graph(g, [] (unsigned long v) -> bool { return true; }, [] (unsigned long e) -> bool { return true; });

	//const auto &viter1 = get_vertices(fg1);
	//vector<unsigned long> vs1(begin(viter1), end(viter1));

	//assert(vs1.size() == 4);
	//assert(vs1[0] == 0);
	//assert(vs1[1] == 1);
	//assert(vs1[2] == 2);
	//assert(vs1[3] == 3);

	//const auto &eiter1 = get_edges(fg1);
	//vector<int> es1(begin(eiter1), end(eiter1));

	//assert(es1.size() == 3);
	//assert(get_source(fg1, es1[0]) == 0 && get_target(fg1, es1[0]) == 1);
	//assert(get_source(fg1, es1[1]) == 0 && get_target(fg1, es1[1]) == 2);
	//assert(get_source(fg1, es1[2]) == 1 && get_target(fg1, es1[2]) == 3);

	//[> filtered <]
	//vector<bool> vf1 = { true, false, false, true };
	//vector<bool> ef1 = { true, false, true };

	//auto fg2 = make_filtered_graph(g, [&vf1] (unsigned long v) -> bool { return vf1[v]; }, [&ef1] (unsigned long e) -> bool { return ef1[e]; });

	//const auto &viter2 = get_vertices(fg2);
	//vector<unsigned long> vs2(begin(viter2), end(viter2));

	//assert(vs2.size() == 2);
	//assert(vs2[0] == 0);
	//assert(vs2[1] == 3);

	//const auto &eiter2 = get_edges(fg2);
	//vector<int> es2(begin(eiter2), end(eiter2));

	//assert(es2.size() == 2);
	//assert(get_source(fg2, es2[0]) == 0 && get_target(fg2, es2[0]) == 1);
	//assert(get_source(fg2, es2[1]) == 1 && get_target(fg2, es2[1]) == 3);
}

void test_connected_components()
{
	RRGraph g;
	add_vertex(g, 6);

	add_edge(g, 0, 1);
	add_edge(g, 1, 0);

	add_edge(g, 1, 2);
	add_edge(g, 2, 1);

	add_edge(g, 3, 4);
	add_edge(g, 4, 3);
	
	add_edge(g, 4, 5);
	add_edge(g, 5, 4);

	vector<vector<int>> components;
	connected_components(g, components);

	assert(components.size() == 2);
	assert(components[0] == vector<int>({ 0, 1, 2 }));
	assert(components[1] == vector<int>({ 3, 4, 5 }));

	RRGraph g2;
	add_vertex(g2, 6);

	add_edge(g2, 1, 0);
	add_edge(g2, 0, 1);

	add_edge(g2, 1, 2);
	add_edge(g2, 2, 1);

	add_edge(g2, 5, 4);
	add_edge(g2, 4, 5);

	add_edge(g2, 3, 4);
	add_edge(g2, 4, 3);

	vector<vector<int>> components2;
	connected_components(g2, components2);

	assert(components2.size() == 2);
	assert(components2[0] == vector<int>({ 0, 1, 2 }));
	assert(components2[1] == vector<int>({ 3, 4, 5 }));
}

int get_rr_node_index(int x, int y, t_rr_type rr_type, int ptc, t_ivec *** L_rr_node_indices);

void init_channel_only_graph(RRGraph &channel_g)
{
	extern s_rr_node *rr_node;
	extern int num_rr_nodes;
	extern t_rr_indexed_data *rr_indexed_data;
	extern struct s_switch_inf *switch_inf;

	add_vertex(channel_g, num_rr_nodes);

	for (int i = 0; i < num_vertices(channel_g); ++i) {
		auto &v = get_vertex_props(channel_g, i);
		v.type = rr_node[i].type;
		v.inc_direction = rr_node[i].direction == INC_DIRECTION;
		v.xlow = rr_node[i].xlow;
		v.ylow = rr_node[i].ylow;
		v.xhigh = rr_node[i].xhigh;
		v.yhigh = rr_node[i].yhigh;

		extern struct s_grid_tile **grid;
		auto type = &grid[v.xlow][v.ylow];

		//v.real_xlow = rr_node[i].xlow;
		//v.real_ylow = rr_node[i].ylow;
		//v.real_xhigh = rr_node[i].xhigh;
		//v.real_yhigh = rr_node[i].ylow + type->offset;
		v.R = rr_node[i].R;
		v.C = rr_node[i].C;
		v.cost_index = rr_node[i].cost_index;
		v.capacity = rr_node[i].capacity;

		char buffer[256];
		sprintf_rr_node(i, buffer);
		//zlog_debug(rr_log, "%s: real_xlow: %d real_xhigh: %d real_ylow: %d real_yhigh: %d\n", buffer, v.real_xlow, v.real_xhigh, v.real_ylow, v.real_yhigh);

		for (int j = 0; j < rr_node[i].num_edges; ++j) {
			int neighbor_id = rr_node[i].edges[j];

			//if ((rr_node[i].type == CHANX || rr_node[i].type == CHANY) &&
					//(rr_node[neighbor_id].type == CHANX || rr_node[neighbor_id].type == CHANY) &&
					//get_track_domain(rr_node[neighbor_id].ptc_num, num_partitions) != get_track_domain(rr_node[i].ptc_num, num_partitions)) {
				//continue;
			//}
			if ((rr_node[i].type != CHANX && rr_node[i].type != CHANY) ||
					(rr_node[neighbor_id].type != CHANX && rr_node[neighbor_id].type != CHANY)) {
				continue;
			}

			auto &e = add_edge(channel_g, i, neighbor_id);

			int si = rr_node[i].switches[j];

			auto &e_p = get_edge_props(channel_g, e);
			
			//e_p.buffered = switch_inf[si].buffered; 
			//e_p.switch_delay = switch_inf[si].Tdel; 
			//e_p.R = switch_inf[si].R; 
		}
	}
	//printf("RR graph num vertices: %d\n", num_vertices(channel_g));
	//printf("RR graph num edges: %d\n", num_edges(channel_g));

	//dump_rr_graph(channel_g, "/Volumes/DATA/channel_rr_graph.txt");
}

std::pair<int, int> get_node_start(int inode);

struct visitor_2_t {
	vector<VertexColor> *color;
	set<int> *all_visited_nodes;
	vector<RREdge> visited_edges;
	vector<int> visited_nodes;
	map<tuple<int, int, t_rr_type, bool>, int> visited_channels;
	map<tuple<int, int, t_rr_type, bool>, int> channel_distance;
	int num_nodes_discovered;
	vector<int> distance;
	vector<int> pred;

	visitor_2_t(int num_vertices) :
		num_nodes_discovered(0), distance(num_vertices, std::numeric_limits<int>::max()), pred(num_vertices, -1)
	{
	}

	bool tree_edge(const RREdge &e, const RRGraph &g)
	{
		int source_id = get_source(g, e);
		assert((*color)[source_id] == VertexColor::GRAY);
		if (find(begin(visited_nodes), end(visited_nodes), source_id) != end(visited_nodes)) {
			if (record_source(g, get_target(g, e), distance[source_id]+1, source_id)) {
				visited_edges.push_back(e);
			}
		}
		return true;
	}

	void examine_edge(const RREdge &e, const RRGraph &g)
	{
	}

	bool discover_vertex(int v, const RRGraph &g)
	{
		int d = distance[v];
		if (d < std::numeric_limits<int>::max()) {
			char buffer[256];
			sprintf_rr_node(v, buffer);
			for (int i = 0; i < d; ++i) {
				zlog_level(delta_log, ROUTER_V3, "\t");
			}
			zlog_level(delta_log, ROUTER_V3, "Current: %s\n", buffer);
		}
		++num_nodes_discovered;
		return true;
	}

	void examine_vertex(int v, const RRGraph &g)
	{
	}

	bool record_source(const RRGraph &g, int v, int d, int p)
	{
		assert((*color)[v] == VertexColor::WHITE);

		const auto &ver = get_vertex_props(g, v);

		const auto &start = get_node_start(v);
		const auto &key = make_tuple(start.first, start.second, ver.type, ver.inc_direction);

		auto d_iter = channel_distance.find(key);
		auto iter = visited_channels.find(key);

		bool record;
		if (iter == visited_channels.end()) {
			visited_channels.insert(make_pair(key, 1));
			assert(d_iter == channel_distance.end());
			channel_distance.insert(make_pair(key, d));
			record = true;
		} else {
			record = d_iter->second == d;
			if (record) {
				++iter->second;
			}
		}

		if (record) {
			assert(all_visited_nodes->find(v) == all_visited_nodes->end());
			assert(find(begin(visited_nodes), end(visited_nodes), v) == end(visited_nodes));
			visited_nodes.push_back(v);
		}

		distance[v] = d;
		pred[v] = p;

		return record;
	}

	void print_path(int node)
	{
		int current = node;
		char buffer[256];
		sprintf_rr_node(current, buffer);
		printf("Path to %s\n", buffer);

		while (current != -1) {
			sprintf_rr_node(current, buffer);
			printf("\t%s\n", buffer);
			current = pred[current];
		}
	}
};

struct visitor_t {
	vector<VertexColor> *color;
	set<int> *all_visited_nodes;
	vector<RREdge> visited_edges;
	vector<int> visited_nodes;
	set<tuple<int, int, t_rr_type, bool>> visited_channels;
	int num_nodes_discovered;

	visitor_t() :
		num_nodes_discovered(0)
	{
	}

	bool tree_edge(const RREdge &e, const RRGraph &g)
	{
		assert((*color)[get_source(g, e)] == VertexColor::GRAY);
		if (find(begin(visited_nodes), end(visited_nodes), get_source(g, e)) != end(visited_nodes)) {
			bool recorded_to = record(g, get_target(g, e));
			if (recorded_to) {
				visited_edges.push_back(e);
			}
		}
		return true;
	}

	void examine_edge(const RREdge &e, const RRGraph &g)
	{
	}

	bool discover_vertex(int v, const RRGraph &g)
	{
		char buffer[256];
		sprintf_rr_node(v, buffer);
		zlog_level(delta_log, ROUTER_V3, "Current: %s\n", buffer);
		++num_nodes_discovered;
		return true;
	}

	void examine_vertex(int v, const RRGraph &g)
	{
	}

	bool record(const RRGraph &g, int v)
	{
		const auto &ver = get_vertex_props(g, v);
		const auto &start = get_node_start(v);
		const auto &key = make_tuple(start.first, start.second, ver.type, ver.inc_direction);
		bool recorded;
		if (visited_channels.find(key) == visited_channels.end()) {
			visited_channels.insert(key);
			assert((*color)[v] != VertexColor::BLACK);
			assert(find(begin(visited_nodes), end(visited_nodes), v) == end(visited_nodes));
			assert(all_visited_nodes->find(v) == all_visited_nodes->end());
			visited_nodes.push_back(v);
			recorded = true;
		} else {
			recorded = false;
		}
		return recorded;
	}

	bool record_source(const RRGraph &g, int v, int d, int p)
	{	
		const auto &ver = get_vertex_props(g, v);
		const auto &start = get_node_start(v);
		const auto &key = make_tuple(start.first, start.second, ver.type, ver.inc_direction);

		if (visited_channels.find(key) == visited_channels.end()) {
			visited_channels.insert(key);
		} 

		assert((*color)[v] != VertexColor::BLACK);
		assert(all_visited_nodes->find(v) == all_visited_nodes->end());
		assert(find(begin(visited_nodes), end(visited_nodes), v) == end(visited_nodes));
		visited_nodes.push_back(v);

		return true;
	}
};

template<typename Visitor>
void bfs_checked(const RRGraph &g, const vector<unsigned long> &nodes, vector<VertexColor> &color, Visitor &visitor, set<int> &all_visited_nodes, set<int> &all_visited_edges)
{
	visitor.color = &color;
	visitor.all_visited_nodes = &all_visited_nodes;

	for (const auto &node : nodes) {
		assert(visitor.record_source(g, node, 0, -1));
	}

	bfs(g, nodes, color, visitor);

	for (const auto &v : visitor.visited_nodes) {
		assert(all_visited_nodes.find(v) == all_visited_nodes.end());
		all_visited_nodes.insert(v);
	}

	for (const auto &e : visitor.visited_edges) {
		assert(all_visited_edges.find(e) == all_visited_edges.end());
		all_visited_edges.insert(e);
	}

	std::fill(begin(color), end(color), VertexColor::WHITE);
	for (const auto &n : all_visited_nodes) {
		color[n] = VertexColor::BLACK;
	}

	/* check for duplicates */
	set<int> visited_nodes_set(begin(visitor.visited_nodes), end(visitor.visited_nodes));
	assert(visited_nodes_set.size() == visitor.visited_nodes.size());
	/* make sure our starting points are inside the visited nodes */
	for (const auto &node : nodes) {
		assert(visited_nodes_set.find(node) != visited_nodes_set.end());
	}
	/* generate visited node set from edges just to double check */
	set<int> visited_nodes_from_edges;
	for (const auto &node : nodes) {
		visited_nodes_from_edges.insert(node);
	}
	for (const auto &e : visitor.visited_edges) {
		visited_nodes_from_edges.insert(get_source(g, e));
		visited_nodes_from_edges.insert(get_target(g, e));
	}
	assert(visited_nodes_set == visited_nodes_from_edges);
	/* check that sum of visited channels is equal to num of visited nodes */
	int total_channel_visit_count = 0;
	for (const auto &item : visitor.visited_channels) {
		total_channel_visit_count += item.second;
	}
	assert(total_channel_visit_count == visitor.visited_nodes.size());
}

template<typename Graph>
void routability(const Graph &g)
{
	struct routability_visitor_t {
		vector<set<int>> visited_sinks;
		vector<set<int>> link_to_previous_visited_sinks;
		vector<int> first_visited_source;
		int current_source;
		bool pending;
		int pending_to;

		routability_visitor_t(int num_vertices)
			: visited_sinks(num_vertices), link_to_previous_visited_sinks(num_vertices), first_visited_source(num_vertices, -1), current_source(-1)
		{
		}

		bool tree_edge(const RREdge &e, const RRGraph &g)
		{
			int to = get_target(g, e);
			if (pending) {
				assert(pending_to == to);
				assert(link_to_previous_visited_sinks[current_source].find(first_visited_source[to]) != 
						link_to_previous_visited_sinks[current_source].end());
				link_to_previous_visited_sinks[current_source].erase(first_visited_source[to]);
				pending = false;
			}
			return true;
		}

		void examine_edge(const RREdge &e, const RRGraph &g)
		{
			int to = get_target(g, e);
			if (first_visited_source[to] != -1 && first_visited_source[to] != current_source) {
				if (link_to_previous_visited_sinks[current_source].find(first_visited_source[to]) == 
link_to_previous_visited_sinks[current_source].end()) {
					link_to_previous_visited_sinks[current_source].insert(first_visited_source[to]);
					pending = true;
					pending_to = to;
				} else {
					pending = false;
				}
			} else {
				pending = false;
			}
		}

		void discover_vertex(int v, const RRGraph &g)
		{
			assert(first_visited_source[v] == -1);
			first_visited_source[v] = current_source;

			const auto &ver = get_vertex_props(g, v);
			if (ver.type == SINK) {
				assert(visited_sinks[current_source].find(v) == visited_sinks[current_source].end());
				visited_sinks[current_source].insert(v);
				//for (const auto &links : link_to_previous_visited_sinks[current_source]) {

				//}
			}
		}

		void examine_vertex(int v, const RRGraph &g)
		{
		}
	};

	int num_sinks = 0;
	for (const auto &v : get_vertices(g)) {
		const auto &ver = get_vertex_props(g, v);

		if (ver.type == SINK) {
			++num_sinks;
		}
	}

	vector<VertexColor> color(num_vertices(g), VertexColor::WHITE);
	routability_visitor_t visitor(num_vertices(g));

	for (const auto &v : get_vertices(g)) {
		const auto &ver = get_vertex_props(g, v);

		if (ver.type == SOURCE) {
			visitor.current_source = v;
			bfs(g, { static_cast<unsigned long>(v) }, color, visitor);

			set<int> all_visited_sinks;
			for (const auto &s : visitor.visited_sinks[v]) {
				all_visited_sinks.insert(s);
			}
			for (const auto &link : visitor.link_to_previous_visited_sinks[v]) {
				for (const auto &s : visitor.visited_sinks[link]) {
					assert(all_visited_sinks.find(s) == all_visited_sinks.end());
					all_visited_sinks.insert(s);
				}
			}

			struct simple_sink_visitor_t {
				set<int> visited_sinks;

				bool tree_edge(const RREdge &e, const RRGraph &g)
				{
					//int to = get_target(g, e);
					//if (pending) {
					//assert(pending_to == to);
					//link_to_previous_visited_sinks[current_source].insert(first_visited_source[to]);
					//pending = false;
					//}
					return true;
				}

				void examine_edge(const RREdge &e, const RRGraph &g)
				{
				}

				void discover_vertex(int v, const RRGraph &g)
				{
					const auto &ver = get_vertex_props(g, v);
					if (ver.type == SINK) {
						assert(visited_sinks.find(v) == visited_sinks.end());
						visited_sinks.insert(v);
					}
				}

				void examine_vertex(int v, const RRGraph &g)
				{
				}
			};

			simple_sink_visitor_t sv;
			vector<VertexColor> temp_color(num_vertices(g), VertexColor::WHITE);
			bfs(g, { static_cast<unsigned long>(v) }, temp_color, sv);

			assert(sv.visited_sinks == all_visited_sinks);

			//printf("Source %d visited %d/%d (%g) sinks\n", v, all_visited_sinks.size(), num_sinks, 100.0*all_visited_sinks.size()/num_sinks);
		}
	}
}

void init_partitioned_graph_6(int num_partitions, const RRGraph &channel_with_interior_g, const RRGraph &channel_without_interior_g, const RRGraph &orig_g, vector<RRGraph *> &graphs, int source_x, int source_y)
{
	//extern s_rr_node *rr_node;
	//extern int num_rr_nodes;
	//extern t_rr_indexed_data *rr_indexed_data;
	//extern struct s_switch_inf *switch_inf;

	//extern int nx, ny;
	//extern t_ivec ***rr_node_indices;
	//extern t_grid_tile **grid;
	//extern int *chan_width_x, *chan_width_y;

	//assert(num_vertices(orig_g) == num_vertices(channel_with_interior_g));
	//assert(num_vertices(orig_g) == num_vertices(channel_without_interior_g));
	
	//int num_visited_edges = 0;
	//int ptc = 0;
	////int source_x = (nx+2)/2;
	////int source_y = (ny+2)/2;
	////int source_x = 0;
	////int source_y = 1;
	//int i = 0;
	//t_type_ptr type = grid[source_x][source_y].type;

	////while (ptc < type->num_class && type->class_inf[ptc].type != DRIVER) {
		////++ptc;
	////}

	////while (num_visited_edges < num_edges(orig_g) && ptc < type->num_class) {
		////int source_rr_node_id = get_rr_node_index(source_x, source_y, SOURCE, ptc, rr_node_indices);
		////
	//set<int> all_visited_nodes;
	//set<int> all_visited_edges;
	//vector<VertexColor> color(num_vertices(channel_without_interior_g), VertexColor::WHITE);
	//vector<int> partition(num_vertices(orig_g), -1);
	//vector<vector<int>> partition_nodes;
	//vector<map<tuple<int, int, t_rr_type, bool>, int>> partition_visited_channels;
	//int part = 0;
	//for (int pass = 0; pass < 2; ++pass) {
		//int width = pass == 0 ? chan_width_x[source_x] : chan_width_y[source_y];
		//int real_num_tracks;
		//if (pass == 0) {
			//real_num_tracks = rr_node_indices[CHANX][source_y][source_x].nelem;
		//} else {
			//real_num_tracks = rr_node_indices[CHANY][source_x][source_y].nelem;
		//}
		////int num_tracks = rr_node_indices[pass == 0 ? CHANX : CHANY][source_x][source_y].nelem;
		////assert(real_num_tracks == num_tracks);
		//if (real_num_tracks == 0) {
			//continue;
		//}
		//vector<int> starting_wires;
		//for (int i = 0; i < width; ++i) {
			//int source_rr_node_id = get_rr_node_index(source_x, source_y, pass == 0 ? CHANX : CHANY, i, rr_node_indices);
			//const auto &source = get_vertex_props(channel_without_interior_g, source_rr_node_id);

			//if ([>starts_at(source, source_x, source_y) && <]color[source_rr_node_id] == VertexColor::WHITE) {
				//starting_wires.push_back(source_rr_node_id);
			//}
		//}
		//int inc = ceil((float)starting_wires.size()/num_partitions);
		////int inc = ceil((float)real_num_tracks/num_partitions);
		//for (int i = 0, batch = 0; i < starting_wires.size(); ++batch, ++part) {
		////for (int i = 0, batch = 0; i < real_num_tracks; ++batch, ++part) {
			//vector<unsigned long> sources;
			//for (int j = 0; j < inc && i < starting_wires.size(); ++j, ++i) {
			////for (int j = 0; j < inc && i < real_num_tracks; ++j, ++i) {
				//int source_rr_node_id = starting_wires[i];
				////int source_rr_node_id = get_rr_node_index(source_x, source_y, pass == 0 ? CHANX : CHANY, i, rr_node_indices); 
				//assert(color[source_rr_node_id] == VertexColor::WHITE);
				//sources.push_back(source_rr_node_id);
				//char buffer[256];
				//sprintf_rr_node(source_rr_node_id, buffer);
				//printf("first stage batch %d source %d: %s\n", batch, j, buffer);
			//}

			//visitor_2_t visitor(num_vertices(channel_without_interior_g));
			//bfs_checked(channel_without_interior_g, sources, color, visitor, all_visited_nodes, all_visited_edges);
			//visitor.print_path(89280);

			////set<tuple<int, int, t_rr_type, bool>> visited_channels_new;
			////for (const auto &n : visitor.visited_nodes) {
				////const auto &start = get_node_start(n);
				////const auto &node = get_vertex_props(channel_without_interior_g, n);
				////const auto &key = make_tuple(start.first, start.second, node.type, node.inc_direction);
				////assert(visited_channels_new.find(key) == visited_channels_new.end());
				////visited_channels_new.insert(key);
			////}

			//char filename[256];
			//sprintf(filename, "/Volumes/DATA/visited_pass_%d_batch_%d.txt", pass, batch);
			//FILE *file = fopen(filename, "w");
			//for (const auto &v : visitor.visited_nodes) {
				//const auto &start = get_node_start(v);
				//fprintf(file, "%d %d\n", start.first, start.second);
			//}
			//fclose(file);

			//partition_nodes.push_back(visitor.visited_nodes);
			//for (const auto &v : visitor.visited_nodes) {
				//partition[v] = part;
			//}

			//int num_unvisited_nodes = (nx+2)*(ny+2)*4 - visitor.visited_nodes.size();
			////for (int x = 0; x < nx+2; ++x) {
				////for (int y = 0; y < ny+2; ++y) {
					////if (visited_channels_new.find(make_tuple(x, y, CHANX, true)) == visited_channels_new.end()) {
						//////printf("Did not visit (%d,%d) CHANX increasing\n", x, y);
						////++num_unvisited_nodes;
					////}
					////if (visited_channels_new.find(make_tuple(x, y, CHANX, false)) == visited_channels_new.end()) {
						//////printf("Did not visit (%d,%d) CHANX decreasing\n", x, y);
						////++num_unvisited_nodes;
					////}
					////if (visited_channels_new.find(make_tuple(x, y, CHANY, true)) == visited_channels_new.end()) {
						//////printf("Did not visit (%d,%d) CHANY increasing\n", x, y);
						////++num_unvisited_nodes;
					////}
					////if (visited_channels_new.find(make_tuple(x, y, CHANY, false)) == visited_channels_new.end()) {
						//////printf("Did not visit (%d,%d) CHANY decreasing\n", x, y);
						////++num_unvisited_nodes;
					////}
				////}
			////}

			//printf("num_nodes_discovered = %d\n", visitor.num_nodes_discovered);
			//printf("num_visited_nodes = %lu\n", visitor.visited_nodes.size());
			//printf("num_unvisited_nodes = %d\n", num_unvisited_nodes);
			//printf("num_visited_edge = %lu\n", visitor.visited_edges.size());

			//num_visited_edges += visitor.visited_edges.size();	
			////while (++ptc < type->num_class && type->class_inf[ptc].type != DRIVER) {
			////}
		//}
	//}

	//int num_channels = 0;
	//for (const auto &v : get_vertices(orig_g)) {
		//const auto &ver = get_vertex_props(orig_g, v);
		//if (ver.type == CHANX || ver.type == CHANY) {
			//++num_channels;
		//}
	//}
	//printf("%d,%d num_all_visited_nodes = %lu num_channels = %d\n", source_x, source_y, all_visited_nodes.size(), num_channels);

	//char filename[256];
	//sprintf(filename, "/Volumes/DATA/unvisited.txt");
	//FILE *file = fopen(filename, "w");
	//int num_unvisited_nodes = 0;
	//for (const auto &v : get_vertices(orig_g)) {
		//const auto &ver = get_vertex_props(orig_g, v);
		//if ((ver.type == CHANX || ver.type == CHANY) && all_visited_nodes.find(v) == all_visited_nodes.end()) {
			//const auto &start = get_node_start(v);
			//fprintf(file, "%d %d\n", start.first, start.second);
			//++num_unvisited_nodes;
		//}
	//}
	//fclose(file);

	//sprintf(filename, "/Volumes/DATA/visited.txt");
	//file = fopen(filename, "w");
	//for (const auto &v : all_visited_nodes) {
		//const auto &start = get_node_start(v);
		//fprintf(file, "%d %d\n", start.first, start.second);
	//}
	//fclose(file);

	//for (int i = 0; i < partition_nodes.size(); ++i) {
		//for (const auto &v : partition_nodes[i]) {
			//assert(partition[v] == i);
		//}
	//}

	//for (const auto &e : get_edges(channel_with_interior_g)) {
		//int from = get_source(channel_with_interior_g, e);
		//int to = get_target(channel_with_interior_g, e);
		//if (!has_edge(channel_without_interior_g, from, to) || all_visited_edges.find(get_edge(channel_without_interior_g, from, to)) ==
					//all_visited_edges.end()) {
			//if (!(partition[from] == -1 && partition[to] == -1)) {
				//char buffer[256];
				////sprintf_rr_node
				//assert(false);
			//}
			//if (partition[from] == -1 || partition[to] == -1) {
			//} else if (partition[from] == partition[to]) {
				////partition_nodes.
			//}
		//}
	//}

	////for (const auto &v : get_vertices(channel_without_interior_g)) {
		////if ((v.type == CHANX || v.type == CHANY)) {
			////if (all_visited_nodes.find(id(v)) == all_visited_nodes.end()) {
				////visitor_t visitor;
				////bfs_checked(channel_without_interior_g, { v }, color, visitor, all_visited_nodes);

				////int num_unvisited_nodes = (nx+2)*(ny+2)*4 - visitor.visited_nodes.size();

				////printf("source: %d\n", id(v));
				////printf("num_nodes_discovered = %d\n", visitor.num_nodes_discovered);
				////printf("num_visited_nodes = %lu\n", visitor.visited_nodes.size());
				////printf("num_unvisited_nodes = %d\n", num_unvisited_nodes);
				////printf("num_visited_edge = %lu\n", visitor.visited_edges.size());
			////} else {
				////assert(color[id(v)] == VertexColor::BLACK);
			////}
		////}
	////}

	//return;

	//vector<int> partitions;
	//int num_chan_edges = 0;
	//int num_misc_edges = -1;

	//for (int i = 0; i < num_partitions; ++i) {
		//RRGraph &new_g = *new RRGraph;
		//add_vertex(new_g, num_vertices(channel_with_interior_g));

		//for (int i = 0; i < num_vertices(new_g); ++i) {
			//get_vertex_props(new_g, i) = get_vertex_props(orig_g, i);
		//}

		//int num_intra_partiton_interior_edges = 0;
		//int num_interior_edges = 0;
		//for (const auto &e : get_edges(channel_with_interior_g)) {
			//int from = get_source(channel_with_interior_g, e);
			//int to = get_target(channel_with_interior_g, e);

			//if (!has_edge(channel_without_interior_g, from, to)) {
				//++num_interior_edges;
			//}

			//if (partitions[from] == i && partitions[from] == partitions[to]) {
				//if (!has_edge(channel_without_interior_g, from, to)) {
					//++num_intra_partiton_interior_edges;
				//}

				//const RREdge &new_e = add_edge(new_g, from, to);
				//get_edge_props(new_g, new_e) = get_edge_props(channel_with_interior_g, e);
			//}
		//}

		//num_chan_edges += num_edges(new_g);

		//printf("graph %d num intra partition interior edges/num interior edges = %d/%d (%g)\n", i, num_intra_partiton_interior_edges, num_interior_edges, (float)num_intra_partiton_interior_edges/num_interior_edges*100);

		//vector<vector<int>> components;
		//vector<bool> vf(num_vertices(new_g));
		//vector<bool> ef(num_edges(new_g));
		////connected_components(make_filtered_graph(new_g, &vf, &ef), components);
		//int total_num_nodes_in_components = 0;
		//for (const auto &comp : components) {
			//total_num_nodes_in_components += comp.size();
		//}
		////assert(total_num_nodes_in_components == num_vertices(new_g));
		//printf("graph %d num connected components = %lu average num nodes = %g\n", i, components.size(), (float)total_num_nodes_in_components/components.size());

		//char filename[256];
		//char buffer[256];
		//sprintf(filename, "/Volumes/DATA/components/graph_%d_components.txt", i);
		//FILE *file = fopen(filename, "w");
		//for (int j = 0; j < components.size(); ++j) {
			//fprintf(file, "Component %d Size %lu\n", j, components[j].size());
			//for (int k = 0; k < components[j].size(); ++k) {
				//sprintf_rr_node(components[j][k], buffer);
				////if (get_vertex_props(orig_g, components[j][k]).type != CHANX && 
////get_vertex_props(orig_g, components[j][k]).type != CHANY) {
					////assert(false);
				////}
				//fprintf(file, "\t%s Part: %d\n", buffer, partitions[components[j][k]]);
			//}
		//}
		//fclose(file);

		//int temp_num_misc_edges = 0;
		//for (const auto &e : get_edges(orig_g)) {
			//int from = get_source(orig_g, e);
			//int to = get_target(orig_g, e);
			//const auto &from_p = get_vertex_props(orig_g, from);
			//const auto &to_p = get_vertex_props(orig_g, to);

			//if ((from_p.type == CHANX || from_p.type == CHANY) && 
					//(to_p.type == CHANX || to_p.type == CHANY)) {
				//continue;
			//}

			//const RREdge &new_e = add_edge(new_g, from, to);
			//get_edge_props(new_g, new_e) = get_edge_props(orig_g, e);

			//++temp_num_misc_edges;
		//}

		//if (num_misc_edges == -1) {
			//num_misc_edges = temp_num_misc_edges;
		//}

		//assert(num_edges(new_g) <= num_edges(orig_g));

		//graphs.push_back(&new_g);
	//}

	//assert(num_misc_edges + num_chan_edges <= num_edges(orig_g));
	//printf("num_chan_edges: %d total_chan_edges: %d percentage: %g\n", num_chan_edges, num_edges(orig_g)-num_misc_edges, num_chan_edges*100.0/(num_edges(orig_g)-num_misc_edges));
	////assert(all_of(begin(edge_valid), end(edge_valid), [] (bool val) -> bool { return !val; }));
	////assert(num_edges_spanned == num_edges(g));
}

void init_partitioned_graph_5(int num_partitions, const RRGraph &channel_with_interior_g, const RRGraph &channel_without_interior_g, const RRGraph &orig_g, vector<RRGraph *> &graphs, int source_x, int source_y)
{
	extern s_rr_node *rr_node;
	extern int num_rr_nodes;
	extern t_rr_indexed_data *rr_indexed_data;
	extern struct s_switch_inf *switch_inf;

	extern int nx, ny;
	extern t_ivec ***rr_node_indices;
	extern t_grid_tile **grid;
	extern int *chan_width_x, *chan_width_y;

	assert(num_vertices(orig_g) == num_vertices(channel_with_interior_g));
	assert(num_vertices(orig_g) == num_vertices(channel_without_interior_g));
	
	int num_visited_edges = 0;
	int ptc = 0;
	//int source_x = (nx+2)/2;
	//int source_y = (ny+2)/2;
	//int source_x = 0;
	//int source_y = 1;
	int i = 0;
	t_type_ptr type = grid[source_x][source_y].type;

	struct visitor_t {
		vector<VertexColor> *color;
		set<int> *all_visited_nodes;
		vector<RREdge> visited_edges;
		vector<int> visited_nodes;
		set<tuple<int, int, t_rr_type, bool>> visited_channels;
		int num_nodes_discovered;

		visitor_t() :
			num_nodes_discovered(0)
		{
		}

		bool tree_edge(const RREdge &e, const RRGraph &g)
		{
			assert((*color)[get_source(g, e)] == VertexColor::GRAY);
			if (find(begin(visited_nodes), end(visited_nodes), get_source(g, e)) != end(visited_nodes)) {
				bool recorded_to = record(g, get_target(g, e));
				if (recorded_to) {
					visited_edges.push_back(e);
				}
			}
			return true;
		}

		void examine_edge(const RREdge &e, const RRGraph &g)
		{
		}

		void discover_vertex(int v, const RRGraph &g)
		{
			char buffer[256];
			sprintf_rr_node(v, buffer);
			zlog_level(delta_log, ROUTER_V3, "Current: %s\n", buffer);
			++num_nodes_discovered;
		}

		void examine_vertex(int v, const RRGraph &g)
		{
		}

		bool record(const RRGraph &g, int v)
		{
			const auto &ver = get_vertex_props(g, v);
			const auto &start = get_node_start(v);
			const auto &key = make_tuple(start.first, start.second, ver.type, ver.inc_direction);
			bool recorded;
			if (visited_channels.find(key) == visited_channels.end()) {
				visited_channels.insert(key);
				assert((*color)[v] != VertexColor::BLACK);
				assert(find(begin(visited_nodes), end(visited_nodes), v) == end(visited_nodes));
				assert(all_visited_nodes->find(v) == all_visited_nodes->end());
				visited_nodes.push_back(v);
				recorded = true;
			} else {
				recorded = false;
			}
			return recorded;
		}
	};


	//while (ptc < type->num_class && type->class_inf[ptc].type != DRIVER) {
		//++ptc;
	//}

	//while (num_visited_edges < num_edges(orig_g) && ptc < type->num_class) {
		//int source_rr_node_id = get_rr_node_index(source_x, source_y, SOURCE, ptc, rr_node_indices);
		//
	set<int> all_visited_nodes;
	vector<VertexColor> color(num_vertices(channel_without_interior_g), VertexColor::WHITE);
	for (int pass = 0; pass < 2; ++pass) {
		int width = pass == 0 ? chan_width_x[source_x] : chan_width_y[source_y];
		int real_num_tracks;
		if (pass == 0) {
			real_num_tracks = rr_node_indices[CHANX][source_y][source_x].nelem;
		} else {
			real_num_tracks = rr_node_indices[CHANY][source_x][source_y].nelem;
		}
		//int num_tracks = rr_node_indices[pass == 0 ? CHANX : CHANY][source_x][source_y].nelem;
		//assert(real_num_tracks == num_tracks);
		if (real_num_tracks == 0) {
			continue;
		}
		for (int i = 0; i < width; ++i) {
			int source_rr_node_id = get_rr_node_index(source_x, source_y, pass == 0 ? CHANX : CHANY, i, rr_node_indices);
			const auto &source = get_vertex_props(channel_without_interior_g, source_rr_node_id);

			if (!starts_at(source, source_x, source_y) || color[source_rr_node_id] != VertexColor::WHITE) {
				continue;
			}

			char buffer[256];
			sprintf_rr_node(source_rr_node_id, buffer);
			printf("first stage source %d: %s\n", i, buffer);

			visitor_t visitor;
			visitor.color = &color;
			visitor.all_visited_nodes = &all_visited_nodes;
			assert(visitor.record(channel_without_interior_g, source_rr_node_id));
			bfs(channel_without_interior_g, { static_cast<unsigned long>(source_rr_node_id) }, color, visitor);

			for (const auto &v : visitor.visited_nodes) {
				assert(all_visited_nodes.find(v) == all_visited_nodes.end());
				all_visited_nodes.insert(v);
			}

			std::fill(begin(color), end(color), VertexColor::WHITE);
			for (const auto &n : all_visited_nodes) {
				color[n] = VertexColor::BLACK;
			}

			set<int> debug_visited_nodes(begin(visitor.visited_nodes), end(visitor.visited_nodes));
			assert(debug_visited_nodes.size() == visitor.visited_nodes.size());

			assert(debug_visited_nodes.find(source_rr_node_id) != debug_visited_nodes.end());

			set<int> debug_visited_nodes_2;
			debug_visited_nodes_2.insert(source_rr_node_id);
			for (const auto &e : visitor.visited_edges) {
				debug_visited_nodes_2.insert(get_source(channel_without_interior_g, e));
				debug_visited_nodes_2.insert(get_target(channel_without_interior_g, e));
			}
			assert(debug_visited_nodes == debug_visited_nodes_2);

			set<tuple<int, int, t_rr_type, bool>> visited_channels_new;
			for (const auto &n : visitor.visited_nodes) {
				const auto &start = get_node_start(n);
				const auto &node = get_vertex_props(channel_without_interior_g, n);
				const auto &key = make_tuple(start.first, start.second, node.type, node.inc_direction);
				assert(visited_channels_new.find(key) == visited_channels_new.end());
				visited_channels_new.insert(key);
			}

			char filename[256];
			sprintf(filename, "/Volumes/DATA/visited_pass_%d_track_%d.txt", pass, i);
			FILE *file = fopen(filename, "w");
			for (const auto &v : visitor.visited_nodes) {
				const auto &start = get_node_start(v);
				fprintf(file, "%d %d\n", start.first, start.second);
			}
			fclose(file);

			int num_unvisited_nodes = 0;
			for (int x = 0; x < nx+2; ++x) {
				for (int y = 0; y < ny+2; ++y) {
					if (visited_channels_new.find(make_tuple(x, y, CHANX, true)) == visited_channels_new.end()) {
						//printf("Did not visit (%d,%d) CHANX increasing\n", x, y);
						++num_unvisited_nodes;
					}
					if (visited_channels_new.find(make_tuple(x, y, CHANX, false)) == visited_channels_new.end()) {
						//printf("Did not visit (%d,%d) CHANX decreasing\n", x, y);
						++num_unvisited_nodes;
					}
					if (visited_channels_new.find(make_tuple(x, y, CHANY, true)) == visited_channels_new.end()) {
						//printf("Did not visit (%d,%d) CHANY increasing\n", x, y);
						++num_unvisited_nodes;
					}
					if (visited_channels_new.find(make_tuple(x, y, CHANY, false)) == visited_channels_new.end()) {
						//printf("Did not visit (%d,%d) CHANY decreasing\n", x, y);
						++num_unvisited_nodes;
					}
				}
			}

			printf("num_nodes_discovered = %d\n", visitor.num_nodes_discovered);
			printf("num_visited_nodes = %lu\n", visitor.visited_nodes.size());
			printf("num_unvisited_nodes = %d\n", num_unvisited_nodes);
			printf("num_visited_edge = %lu\n", visitor.visited_edges.size());

			num_visited_edges += visitor.visited_edges.size();	
			//while (++ptc < type->num_class && type->class_inf[ptc].type != DRIVER) {
			//}
		}
	}

	int num_channels = 0;
	for (const auto &v : get_vertices(orig_g)) {
		const auto &ver = get_vertex_props(orig_g, v);
		if (ver.type == CHANX || ver.type == CHANY) {
			++num_channels;
		}
	}
	printf("%d,%d num_all_visited_nodes = %lu num_channels = %d\n", source_x, source_y, all_visited_nodes.size(), num_channels);

	char filename[256];
	sprintf(filename, "/Volumes/DATA/unvisited.txt");
	FILE *file = fopen(filename, "w");
	int num_unvisited_nodes = 0;
	for (const auto &v : get_vertices(orig_g)) {
		const auto &ver = get_vertex_props(orig_g, v);
		if ((ver.type == CHANX || ver.type == CHANY) && all_visited_nodes.find(v) == all_visited_nodes.end()) {
			const auto &start = get_node_start(v);
			fprintf(file, "%d %d\n", start.first, start.second);
			++num_unvisited_nodes;
		}
	}
	fclose(file);

	sprintf(filename, "/Volumes/DATA/visited.txt");
	file = fopen(filename, "w");
	for (const auto &v : all_visited_nodes) {
		const auto &start = get_node_start(v);
		fprintf(file, "%d %d\n", start.first, start.second);
	}
	fclose(file);

	for (const auto &v : get_vertices(channel_without_interior_g)) {
		const auto &ver = get_vertex_props(channel_without_interior_g, v);
		if ((ver.type == CHANX || ver.type == CHANY)) {
			if (all_visited_nodes.find(v) == all_visited_nodes.end()) {
				visitor_t visitor;
				visitor.color = &color;
				visitor.all_visited_nodes = &all_visited_nodes;
				assert(visitor.record(channel_without_interior_g, v));
				bfs(channel_without_interior_g, { v }, color, visitor);

				for (const auto &v : visitor.visited_nodes) {
					assert(all_visited_nodes.find(v) == all_visited_nodes.end());
					all_visited_nodes.insert(v);
				}

				std::fill(begin(color), end(color), VertexColor::WHITE);
				for (const auto &n : all_visited_nodes) {
					color[n] = VertexColor::BLACK;
				}

				set<int> debug_visited_nodes(begin(visitor.visited_nodes), end(visitor.visited_nodes));
				assert(debug_visited_nodes.size() == visitor.visited_nodes.size());

				assert(debug_visited_nodes.find(v) != debug_visited_nodes.end());

				set<int> debug_visited_nodes_2;
				debug_visited_nodes_2.insert(v);
				for (const auto &e : visitor.visited_edges) {
					debug_visited_nodes_2.insert(get_source(channel_without_interior_g, e));
					debug_visited_nodes_2.insert(get_target(channel_without_interior_g, e));
				}
				assert(debug_visited_nodes == debug_visited_nodes_2);


				int num_unvisited_nodes = (nx+2)*(ny+2)*4 - visitor.visited_nodes.size();

				printf("source: %d\n", v);
				printf("num_nodes_discovered = %d\n", visitor.num_nodes_discovered);
				printf("num_visited_nodes = %lu\n", visitor.visited_nodes.size());
				printf("num_unvisited_nodes = %d\n", num_unvisited_nodes);
				printf("num_visited_edge = %lu\n", visitor.visited_edges.size());
			} else {
				assert(color[v] == VertexColor::BLACK);
			}
		}
	}

	return;

	vector<int> partitions;
	//partition_graph(channel_without_interior_g, num_partitions, partitions);

	int num_chan_edges = 0;
	int num_misc_edges = -1;

	for (int i = 0; i < num_partitions; ++i) {
		RRGraph &new_g = *new RRGraph;
		add_vertex(new_g, num_vertices(channel_with_interior_g));

		for (int i = 0; i < num_vertices(new_g); ++i) {
			get_vertex_props(new_g, i) = get_vertex_props(orig_g, i);
		}

		int num_intra_partiton_interior_edges = 0;
		int num_interior_edges = 0;
		for (const auto &e : get_edges(channel_with_interior_g)) {
			int from = get_source(channel_with_interior_g, e);
			int to = get_target(channel_with_interior_g, e);

			if (!has_edge(channel_without_interior_g, from, to)) {
				++num_interior_edges;
			}

			if (partitions[from] == i && partitions[from] == partitions[to]) {
				if (!has_edge(channel_without_interior_g, from, to)) {
					++num_intra_partiton_interior_edges;
				}

				const RREdge &new_e = add_edge(new_g, from, to);
				get_edge_props(new_g, new_e) = get_edge_props(channel_with_interior_g, e);
			}
		}

		num_chan_edges += num_edges(new_g);

		printf("graph %d num intra partition interior edges/num interior edges = %d/%d (%g)\n", i, num_intra_partiton_interior_edges, num_interior_edges, (float)num_intra_partiton_interior_edges/num_interior_edges*100);

		vector<vector<int>> components;
		vector<bool> vf(num_vertices(new_g));
		vector<bool> ef(num_edges(new_g));
		//connected_components(make_filtered_graph(new_g, &vf, &ef), components);
		int total_num_nodes_in_components = 0;
		for (const auto &comp : components) {
			total_num_nodes_in_components += comp.size();
		}
		//assert(total_num_nodes_in_components == num_vertices(new_g));
		printf("graph %d num connected components = %lu average num nodes = %g\n", i, components.size(), (float)total_num_nodes_in_components/components.size());

		char filename[256];
		char buffer[256];
		sprintf(filename, "/Volumes/DATA/components/graph_%d_components.txt", i);
		FILE *file = fopen(filename, "w");
		for (int j = 0; j < components.size(); ++j) {
			fprintf(file, "Component %d Size %lu\n", j, components[j].size());
			for (int k = 0; k < components[j].size(); ++k) {
				sprintf_rr_node(components[j][k], buffer);
				//if (get_vertex_props(orig_g, components[j][k]).type != CHANX && 
//get_vertex_props(orig_g, components[j][k]).type != CHANY) {
					//assert(false);
				//}
				fprintf(file, "\t%s Part: %d\n", buffer, partitions[components[j][k]]);
			}
		}
		fclose(file);

		int temp_num_misc_edges = 0;
		for (const auto &e : get_edges(orig_g)) {
			int from = get_source(orig_g, e);
			int to = get_target(orig_g, e);
			const auto &from_p = get_vertex_props(orig_g, from);
			const auto &to_p = get_vertex_props(orig_g, to);

			if ((from_p.type == CHANX || from_p.type == CHANY) && 
					(to_p.type == CHANX || to_p.type == CHANY)) {
				continue;
			}

			const RREdge &new_e = add_edge(new_g, from, to);
			get_edge_props(new_g, new_e) = get_edge_props(orig_g, e);

			++temp_num_misc_edges;
		}

		if (num_misc_edges == -1) {
			num_misc_edges = temp_num_misc_edges;
		}

		assert(num_edges(new_g) <= num_edges(orig_g));

		graphs.push_back(&new_g);
	}

	assert(num_misc_edges + num_chan_edges <= num_edges(orig_g));
	printf("num_chan_edges: %d total_chan_edges: %d percentage: %g\n", num_chan_edges, num_edges(orig_g)-num_misc_edges, num_chan_edges*100.0/(num_edges(orig_g)-num_misc_edges));
	//assert(all_of(begin(edge_valid), end(edge_valid), [] (bool val) -> bool { return !val; }));
	//assert(num_edges_spanned == num_edges(g));
}

void init_partitioned_graph_4(int num_partitions, RRGraph &channel_with_interior_g, RRGraph &channel_without_interior_g, const RRGraph &orig_g, vector<RRGraph *> &graphs, vector<vector<vector<int>>> &all_partition_components)
{
	//extern s_rr_node *rr_node;
	//extern int num_rr_nodes;
	//extern t_rr_indexed_data *rr_indexed_data;
	//extern struct s_switch_inf *switch_inf;

	//extern int nx, ny;
	//extern t_ivec ***rr_node_indices;
	//extern t_grid_tile **grid;
	//extern int *chan_width_x, *chan_width_y;

	//assert(num_vertices(orig_g) == num_vertices(channel_with_interior_g));
	//assert(num_vertices(orig_g) == num_vertices(channel_without_interior_g));

	//vector<int> partition;
	//vector<int> weights(num_vertices(channel_without_interior_g));
	//for (const auto &v : get_vertices(channel_without_interior_g)) {
		//const auto &ver = get_vertex_props(channel_without_interior_g, v);
		//if (ver.type == CHANX || ver.type == CHANY) {
			//weights[v] = 1;
		//} else {
			//weights[v] = 0;
		//}
	//}
	//partition_graph(channel_without_interior_g, num_partitions, weights, 1.001, partition);
	//vector<int> num_nodes(num_partitions, 0);
	//vector<int> num_non_chan_nodes(num_partitions, 0);
	//for (int i = 0; i < num_vertices(channel_without_interior_g); ++i) {
		//const auto &v = get_vertex_props(channel_without_interior_g, i);
		//if (v.type != CHANX && v.type != CHANY) {
			//++num_non_chan_nodes[partition[i]];
		//}
		//++num_nodes[partition[i]];
	//}
	//int total_num_nodes = 0;
	//vector<int> num_chan_nodes(num_partitions);
	//for (int i = 0; i < num_partitions; ++i) {
		//printf("part %d num nodes %d num_chan_nodes: %d num_non_chan_nodes: %d\n", i, num_nodes[i], num_nodes[i] - num_non_chan_nodes[i], num_non_chan_nodes[i]);
		//total_num_nodes += num_nodes[i];
		//num_chan_nodes[i] = num_nodes[i] - num_non_chan_nodes[i];
	//}
	//assert(total_num_nodes == num_vertices(channel_without_interior_g));

	//for (int i = 0; i < num_vertices(orig_g); ++i) {
		//const auto &ver = get_vertex_props(orig_g, i);
		//if (ver.type != CHANX && ver.type != CHANY) {
			//partition[i] = -1;
		//}
	//}

	//int num_chan_edges = 0;
	//int num_non_chan_edges = -1;

	

	//for (int i = 0; i < num_partitions; ++i) {
		//RRGraph &new_g = *new RRGraph;
		//add_vertex(new_g, num_vertices(channel_with_interior_g));

		//for (int i = 0; i < num_vertices(new_g); ++i) {
			//get_vertex_props(new_g, i) = get_vertex_props(orig_g, i);
		//}

		//int num_intra_partiton_interior_edges = 0;
		//int num_interior_edges = 0;
		//for (const auto &e : get_edges(channel_with_interior_g)) {
			//int from = get_source(channel_with_interior_g, e);
			//int to = get_target(channel_with_interior_g, e);

			//if (!has_edge(channel_without_interior_g, from, to)) {
				//++num_interior_edges;
			//}

			//if (partition[from] == i && partition[from] == partition[to]) {
				//if (!has_edge(channel_without_interior_g, from, to)) {
					//++num_intra_partiton_interior_edges;
				//}

				//const RREdge &new_e = add_edge(new_g, from, to);
				//get_edge_props(new_g, new_e) = get_edge_props(channel_with_interior_g, e);
			//}
		//}

		//num_chan_edges += num_edges(new_g);

		//printf("graph %d num intra partition interior edges/num interior edges = %d/%d (%g)\n", i, num_intra_partiton_interior_edges, num_interior_edges, (float)num_intra_partiton_interior_edges/num_interior_edges*100);

		//int local_num_non_chan_edges = 0;
		//for (const auto &e : get_edges(orig_g)) {
			//int from = get_source(orig_g, e);
			//int to = get_target(orig_g, e);
			//const auto &from_p = get_vertex_props(orig_g, from);
			//const auto &to_p = get_vertex_props(orig_g, to);

			//if ((from_p.type == CHANX || from_p.type == CHANY) && 
					//(to_p.type == CHANX || to_p.type == CHANY)) {
				//continue;
			//}

			//const RREdge &new_e = add_edge(new_g, from, to);
			//get_edge_props(new_g, new_e) = get_edge_props(orig_g, e);

			//++local_num_non_chan_edges;
		//}

		//if (num_non_chan_edges == -1) {
			//num_non_chan_edges = local_num_non_chan_edges;
		//}

		//assert(num_edges(new_g) <= num_edges(orig_g));

		//graphs.push_back(&new_g);

		//vector<vector<int>> &components = all_partition_components[i];
		////connected_components(new_g, [&partition, &i] (const RRNode &v) -> bool { return partition[id(v)] == i && (v.type == CHANX || v.type == CHANY); }, components);
		////connected_components(new_g, [&partition, &i] (const RRNode &v) -> bool { return partition[id(v)] == i || (v.type != CHANX && v.type != CHANY); }, components);
		//RRGraph temp_g;
		//connected_components_convert_to_undirected(new_g, temp_g);

		//int expected_total_num_nodes_in_components = 0;	
		//vector<bool> vf(num_vertices(temp_g));
		//for (const auto &v : get_vertices(temp_g)) {
			//const auto &ver = get_vertex_props(temp_g, v);
			//if (ver.type != CHANX && ver.type != CHANY) {
				//assert(partition[v] == -1);
			//}
			//vf[v] = partition[v] == i || (ver.type != CHANX && ver.type != CHANY);
			//if (vf[v]) {
				//++expected_total_num_nodes_in_components;
			//}
		//}
		//vector<bool> ef(num_edges(temp_g));
		//for (const auto &e : get_edges(temp_g)) {
			//int from = get_source(temp_g, e);
			//int to = get_target(temp_g, e);
			//ef[e] = vf[from] && vf[to];
		//}

		////connected_components(make_filtered_graph(temp_g, &vf, &ef), components);
		//int total_num_nodes_in_components = 0;
		//for (const auto &comp : components) {
			//total_num_nodes_in_components += comp.size();
		//}
		////assert(total_num_nodes_in_components == num_chan_nodes[i]);
		//assert(total_num_nodes_in_components == expected_total_num_nodes_in_components);
		//printf("graph %d num connected components = %lu average num nodes = %g\n", i, components.size(), (float)total_num_nodes_in_components/components.size());

		////assert(total_num_nodes_in_components == num_vertices(orig_g));

		//char filename[256];
		//char buffer[256];
		//sprintf(filename, "/Volumes/DATA/components/graph_%d_components.txt", i);
		//FILE *file = fopen(filename, "w");
		//for (int j = 0; j < components.size(); ++j) {
			//fprintf(file, "Component %d Size %lu\n", j, components[j].size());
			//for (int k = 0; k < components[j].size(); ++k) {
				//sprintf_rr_node(components[j][k], buffer);
				////if (get_vertex_props(orig_g, components[j][k]).type != CHANX && 
////get_vertex_props(orig_g, components[j][k]).type != CHANY) {
					////assert(false);
				////}
				//fprintf(file, "\t%s Part: %d\n", buffer, partition[components[j][k]]);
			//}
		//}
		//fclose(file);

	//}

	//assert(num_non_chan_edges + num_chan_edges <= num_edges(orig_g));
	//printf("num_chan_edges_after_partition: %d total_chan_edges: %d percentage: %g\n", num_chan_edges, num_edges(orig_g)-num_non_chan_edges, num_chan_edges*100.0/(num_edges(orig_g)-num_non_chan_edges));
	////assert(all_of(begin(edge_valid), end(edge_valid), [] (bool val) -> bool { return !val; }));
	////assert(num_edges_spanned == num_edges(g));
}

void init_partitioned_graph_3(int num_partitions, RRGraph &g, vector<RRGraph *> &graphs)
{
	extern s_rr_node *rr_node;
	extern int num_rr_nodes;
	extern t_rr_indexed_data *rr_indexed_data;
	extern struct s_switch_inf *switch_inf;

	set<RREdge> duplicate_edges;
	RRGraph orig_g;

	add_vertex(g, num_rr_nodes);
	add_vertex(orig_g, num_rr_nodes);

	for (int i = 0; i < num_vertices(g); ++i) {
		auto &v = get_vertex_props(g, i);
		v.type = rr_node[i].type;
		v.inc_direction = rr_node[i].direction == INC_DIRECTION;
		v.xlow = rr_node[i].xlow;
		v.ylow = rr_node[i].ylow;
		v.xhigh = rr_node[i].xhigh;
		v.yhigh = rr_node[i].yhigh;

		extern struct s_grid_tile **grid;
		auto type = &grid[v.xlow][v.ylow];

		//v.real_xlow = rr_node[i].xlow;
		//v.real_ylow = rr_node[i].ylow;
		//v.real_xhigh = rr_node[i].xhigh;
		//v.real_yhigh = rr_node[i].ylow + type->offset;
		v.R = rr_node[i].R;
		v.C = rr_node[i].C;
		v.cost_index = rr_node[i].cost_index;
		v.capacity = rr_node[i].capacity;

		char buffer[256];
		sprintf_rr_node(i, buffer);
		//zlog_debug(rr_log, "%s: real_xlow: %d real_xhigh: %d real_ylow: %d real_yhigh: %d\n", buffer, v.real_xlow, v.real_xhigh, v.real_ylow, v.real_yhigh);

		for (int j = 0; j < rr_node[i].num_edges; ++j) {
			int neighbor_id = rr_node[i].edges[j];

			//if ((rr_node[i].type == CHANX || rr_node[i].type == CHANY) &&
					//(rr_node[neighbor_id].type == CHANX || rr_node[neighbor_id].type == CHANY) &&
					//get_track_domain(rr_node[neighbor_id].ptc_num, num_partitions) != get_track_domain(rr_node[i].ptc_num, num_partitions)) {
				//continue;
			//}
			if ((rr_node[i].type != CHANX && rr_node[i].type != CHANY) ||
					(rr_node[neighbor_id].type != CHANX && rr_node[neighbor_id].type != CHANY)) {
				continue;
			}

			auto &e = add_edge(g, i, neighbor_id);

			int si = rr_node[i].switches[j];
			
			auto &e_p = get_edge_props(g, e);
			
			//e_p.buffered = switch_inf[si].buffered; 
			//e_p.switch_delay = switch_inf[si].Tdel; 
			//e_p.R = switch_inf[si].R; 

			auto &orig_e = add_edge(orig_g, i, neighbor_id);
			get_edge_props(orig_g, orig_e) = get_edge_props(g, e);

			auto &e2 = add_edge(g, neighbor_id, i);

			get_edge_props(orig_g, e2) = get_edge_props(g, e);

			duplicate_edges.insert(e2);
		}
	}
	printf("RR graph num vertices: %d\n", num_vertices(orig_g));
	printf("RR graph num edges: %d\n", num_edges(orig_g));

	dump_rr_graph(orig_g, "/Volumes/DATA/rr_graph.txt");

	extern int nx, ny;
	extern t_ivec ***rr_node_indices;
	extern t_grid_tile **grid;
	extern int *chan_width_x, *chan_width_y;

	vector<int> partitions;
	//partition_graph(g, num_partitions, partitions);

	for (int i = 0; i < num_partitions; ++i) {
		RRGraph new_g;
		add_vertex(new_g, num_vertices(g));

		for (const auto &e : get_edges(g)) {
			if (duplicate_edges.find(e) != end(duplicate_edges)) {
				continue;
			}

			int from = get_source(g, e);
			int to = get_target(g, e);

			if (partitions[from] == i && partitions[from] == partitions[to]) {
				const RREdge &new_e = add_edge(new_g, from, to);
				get_edge_props(new_g, new_e) = get_edge_props(g, e);
			}
		}

		char filename[256];

		sprintf(filename, "/Volumes/DATA/rr_graph_%d.txt", i);
		dump_rr_graph(new_g, filename);

		sprintf(filename, "/Volumes/DATA/rr_graph_%d_edges.txt", i);
		dump_edges(new_g, filename);
	}
	//assert(all_of(begin(edge_valid), end(edge_valid), [] (bool val) -> bool { return !val; }));
	//assert(num_edges_spanned == num_edges(g));
}

void init_partitioned_graph_2(int num_partitions, RRGraph &g, vector<RRGraph *> &graphs)
{
	extern s_rr_node *rr_node;
	extern int num_rr_nodes;
	extern t_rr_indexed_data *rr_indexed_data;
	extern struct s_switch_inf *switch_inf;

	add_vertex(g, num_rr_nodes);
	for (int i = 0; i < num_vertices(g); ++i) {
		auto &v = get_vertex_props(g, i);
		v.type = rr_node[i].type;
		v.inc_direction = rr_node[i].direction == INC_DIRECTION;
		v.xlow = rr_node[i].xlow;
		v.ylow = rr_node[i].ylow;
		v.xhigh = rr_node[i].xhigh;
		v.yhigh = rr_node[i].yhigh;

		extern struct s_grid_tile **grid;
		auto type = &grid[v.xlow][v.ylow];

		//v.real_xlow = rr_node[i].xlow;
		//v.real_ylow = rr_node[i].ylow;
		//v.real_xhigh = rr_node[i].xhigh;
		//v.real_yhigh = rr_node[i].ylow + type->offset;
		v.R = rr_node[i].R;
		v.C = rr_node[i].C;
		v.cost_index = rr_node[i].cost_index;
		v.capacity = rr_node[i].capacity;

		char buffer[256];
		sprintf_rr_node(i, buffer);
		//zlog_debug(rr_log, "%s: real_xlow: %d real_xhigh: %d real_ylow: %d real_yhigh: %d\n", buffer, v.real_xlow, v.real_xhigh, v.real_ylow, v.real_yhigh);

		for (int j = 0; j < rr_node[i].num_edges; ++j) {
			int neighbor = rr_node[i].edges[j];

			//if ((rr_node[i].type == CHANX || rr_node[i].type == CHANY) &&
					//(rr_node[neighbor].type == CHANX || rr_node[neighbor].type == CHANY) &&
					//get_track_domain(rr_node[neighbor].ptc_num, num_partitions) != get_track_domain(rr_node[i].ptc_num, num_partitions)) {
				//continue;
			//}

			auto &e = add_edge(g, i, neighbor);

			int si = rr_node[i].switches[j];
			
			auto &e_p = get_edge_props(g, e);
			
			//e_p.buffered = switch_inf[si].buffered; 
			//e_p.switch_delay = switch_inf[si].Tdel; 
			//e_p.R = switch_inf[si].R; 
		}
	}
	printf("RR graph num vertices: %d\n", num_vertices(g));
	printf("RR graph num edges: %d\n", num_edges(g));

	dump_rr_graph(g, "/Volumes/DATA/rr_graph.txt");

	extern int nx, ny;
	extern t_ivec ***rr_node_indices;
	extern t_grid_tile **grid;
	extern int *chan_width_x, *chan_width_y;
	
	vector<bool> visited(num_vertices(g), false);
	int num_visited_edges = 0;
	int ptc = 0;
	int source_x = (nx+2)/2;
	int source_y = (ny+2)/2;
	int num_srcs = rr_node_indices[SOURCE][source_x][source_y].nelem;
	int i = 0;
	t_type_ptr type = grid[source_x][source_y].type;

	//while (ptc < type->num_class && type->class_inf[ptc].type != DRIVER) {
		//++ptc;
	//}

	//while (num_visited_edges < num_edges(g) && ptc < type->num_class) {
		//int source_rr_node_id = get_rr_node_index(source_x, source_y, SOURCE, ptc, rr_node_indices);
	for (int pass = 0; pass < 2; ++pass) {
		int width = pass == 0 ? chan_width_x[source_x] : chan_width_y[source_y];
		for (int i = 0; i < width; ++i) {
			int source_rr_node_id = get_rr_node_index(source_x, source_y, pass == 0 ? CHANX : CHANY, i, rr_node_indices);

			printf("source: %d\n", source_rr_node_id);

			const auto &source = get_vertex_props(g, source_rr_node_id);
			if (!starts_at(source, source_x, source_y)) {
				continue;
			}

			vector<RREdge> visited_edges;
			//bfs(g, source, pass == 0, source.inc_direction, visited, visited_edges);

			num_visited_edges += visited_edges.size();	
			//while (++ptc < type->num_class && type->class_inf[ptc].type != DRIVER) {
			//}

			printf("visited nodes size: %lu\n", visited_edges.size());

			RRGraph new_g;
			add_vertex(new_g, num_vertices(g));

			for (const auto &e_id : visited_edges) {
				const auto &e_p = get_edge_props(g, e_id);
				const RREdge &new_e = add_edge(new_g, get_source(g, e_id), get_target(g, e_id));
				auto &new_e_p = get_edge_props(new_g, new_e);
				new_e_p = e_p;
			}
			char filename[256];
			sprintf(filename, "/Volumes/DATA/rr_graph_type_%d_track_%d.txt", pass, i);
			dump_rr_graph(new_g, filename);
		}
	}
	//assert(all_of(begin(edge_valid), end(edge_valid), [] (bool val) -> bool { return !val; }));
	//assert(num_edges_spanned == num_edges(g));
}

void init_partitioned_graph(int num_partitions, RRGraph &g, vector<RRGraph *> &graphs)
{
	//extern s_rr_node *rr_node;
	//extern int num_rr_nodes;
	//extern t_rr_indexed_data *rr_indexed_data;
	//extern struct s_switch_inf *switch_inf;

	//add_vertex(g, num_rr_nodes);
	//for (int i = 0; i < num_vertices(g); ++i) {
		//auto &v = get_vertex_props(g, i);
		//v.type = rr_node[i].type;
		//v.inc_direction = rr_node[i].direction == INC_DIRECTION;
		//v.xlow = rr_node[i].xlow;
		//v.ylow = rr_node[i].ylow;
		//v.xhigh = rr_node[i].xhigh;
		//v.yhigh = rr_node[i].yhigh;

		//extern struct s_grid_tile **grid;
		//auto type = &grid[v.xlow][v.ylow];

		////v.real_xlow = rr_node[i].xlow;
		////v.real_ylow = rr_node[i].ylow;
		////v.real_xhigh = rr_node[i].xhigh;
		////v.real_yhigh = rr_node[i].ylow + type->offset;
		//v.R = rr_node[i].R;
		//v.C = rr_node[i].C;
		//v.cost_index = rr_node[i].cost_index;
		//v.capacity = rr_node[i].capacity;

		//char buffer[256];
		//sprintf_rr_node(i, buffer);
		////zlog_debug(rr_log, "%s: real_xlow: %d real_xhigh: %d real_ylow: %d real_yhigh: %d\n", buffer, v.real_xlow, v.real_xhigh, v.real_ylow, v.real_yhigh);

		//for (int j = 0; j < rr_node[i].num_edges; ++j) {
			//int neighbor = rr_node[i].edges[j];

			////if ((rr_node[i].type == CHANX || rr_node[i].type == CHANY) &&
					////(rr_node[neighbor].type == CHANX || rr_node[neighbor].type == CHANY) &&
					////get_track_domain(rr_node[neighbor].ptc_num, num_partitions) != get_track_domain(rr_node[i].ptc_num, num_partitions)) {
				////continue;
			////}

			//auto &e = add_edge(g, i, neighbor);

			//int si = rr_node[i].switches[j];
			
			//auto &e_p = get_edge_props(g, e);
			
			//e_p.buffered = switch_inf[si].buffered; 
			//e_p.switch_delay = switch_inf[si].Tdel; 
			//e_p.R = switch_inf[si].R; 
		//}
	//}
	//printf("RR graph num vertices: %d\n", num_vertices(g));
	//printf("RR graph num edges: %d\n", num_edges(g));

	//dump_rr_graph(g, "rr_graph.txt");
	
	//vector<bool> edge_valid(num_edges(g), true);
	//vector<vector<int>> trees;
	//int num_edges_spanned = 0;
	//while (num_edges_spanned < num_edges(g)) {
		//trees.push_back(vector<int>());

		//auto &edges = trees.back();
		//custom_minimum_spanning_tree(g, [&] (int e) -> bool {
				//const auto &from = get_vertex_props(g, get_source(g, e));
				//const auto &to = get_vertex_props(g, get_target(g, e));
				//return edge_valid[e] && 
				//(from.type == CHANX || from.type == CHANY) &&
				 //(to.type == CHANX || to.type == CHANY);

				//}, edges);
		//num_edges_spanned += edges.size();
		//if (edges.size() == 0) {
			//break;
		//}
		//printf("mst edges size: %lu\n", edges.size());
		//for (const auto &e : edges) {
			//edge_valid[e] = false;
		//}
	//}
	////assert(all_of(begin(edge_valid), end(edge_valid), [] (bool val) -> bool { return !val; }));
	////assert(num_edges_spanned == num_edges(g));

	//for (int i = 0; i < trees.size(); ++i) {
		//RRGraph new_g;
		//add_vertex(new_g, num_vertices(g));

		//for (auto e : trees[i]) {
			//RREdge &new_e = add_edge(new_g, get_source(g, e), get_target(g, e));
			//const auto &edge = get_edge(g, e);
			//new_e.properties = edge.properties;
		//}
		//char filename[256];
		//sprintf(filename, "rr_graph_%d.txt", i);
		//dump_rr_graph(new_g, filename);
	//}

	//int low = 0;
	//int high = trees.size();
	//int inc = ceil((float)trees.size()/(num_partitions*2));
	//vector<int> used(trees.size(), 0);
	//while (low < high) {
		//RRGraph &new_g = *new RRGraph;
		//add_vertex(new_g, num_vertices(g));
		//for (int i = 0; i < num_vertices(g); ++i) {
			//get_vertex_props(new_g, i) = get_vertex_props(g, i);
		//}
		//for (int i = low; i < std::min(low+inc, high); ++i) {
			//assert(used[i] == 0);
			//++used[i];
			//for (auto e : trees[i]) {
				//RREdge &new_e = add_edge(new_g, get_source(g, e), get_target(g, e));
				//const auto &edge = get_edge(g, e);
				//new_e.properties = edge.properties;
			//}
		//}
		//low += inc;
		//for (int i = std::max(high-inc, low+1); i < high; ++i) {
			//assert(used[i] == 0);
			//++used[i];
			//for (auto e : trees[i]) {
				//RREdge &new_e = add_edge(new_g, get_source(g, e), get_target(g, e));
				//const auto &edge = get_edge(g, e);
				//new_e.properties = edge.properties;
			//}
		//}
		//high -= inc;
		//graphs.push_back(&new_g);
		//assert(num_edges(new_g) > 0);
	//}

	//assert(all_of(begin(used), end(used), [] (int val) -> bool { return val == 1; }));
	//assert(graphs.size() == num_partitions);

	//for (const auto e : get_edges(g)) {
		//if (get_vertex_props(g, get_source(g, e)).type == SOURCE || get_vertex_props(g, get_target(g, e)).type == SINK) {
			//for (auto &new_g : graphs) {
				//int from = get_source(g, e);
				//int to = get_target(g, e);

				//if (!has_edge(*new_g, from, to)) {
					//RREdge &new_e = add_edge(*new_g, from, to);
					//const auto &edge = get_edge(g, e);
					//new_e.properties = edge.properties;
				//}
			//}
		//}
	//}

	////for (int i = 0; i < graphs.size(); ++i) {
		////char filename[256];
		////sprintf(filename, "rr_graph_%d.txt", i);
		////dump_rr_graph(*graphs[i], filename);
	////}
}

vector<vector<FILE *>> delta_log_files;
vector<vector<FILE *>> missing_edge_log_files;
FILE *ss_log_file = nullptr;

static void log_impl(zlog_msg_t *msg, FILE *&file)
{
	if (!file) {
		char filename[256];
		assert(false);
		//sprintf(filename, "%s%s", LOG_PATH_PREFIX, msg->path);

		file = fopen(filename, "w");
		if (!file) {
			perror(nullptr);
			assert(false);
		}
	}
	fprintf(file, "%s", msg->buf);
	fflush(file);
}

int ss_log_output(zlog_msg_t *msg)
{
	int tid;
	int iter;

	log_impl(msg, ss_log_file);

	return 0;
}

int delta_log_output(zlog_msg_t *msg)
{
	int tid;
	int iter;

	if (sscanf(msg->path, "iter_%d_tid_%d.log", &iter, &tid) != 2) {
		return 0;
	}

	//concurrent_log_impl(msg, delta_log_files, iter, tid);

	return 0;
}

int missing_edge_log_output(zlog_msg_t *msg)
{
	int tid;
	int iter;

	if (sscanf(msg->path, "missing_edge_iter_%d_tid_%d.log", &iter, &tid) != 2) {
		return 0;
	}

	//concurrent_log_impl(msg, missing_edge_log_files, iter, tid);

	return 0;
}

bool is_channel(const rr_node_property_t &node)
{
	return node.type == CHANX || node.type == CHANY;
}

void update_pseudo_sources(pseudo_net_t &pnet, unrouted_t &unrouted, const RRGraph &orig_g, const vector<int> &pid, int num_partitions, int &net_next_pid)
{
	pnet.pseudo_sources.clear();

	do {
		zlog_level(delta_log, ROUTER_V3, "Getting pseudo sources in partition %d\n", net_next_pid);

		set<int> added_pseudo_sources;
		for (auto &b : unrouted.boundary_nodes) {
			RRNode bnode = b.rr_node;
			zlog_level(delta_log, ROUTER_V3, "Looking at neighbors of boundary node %d\n", bnode);

			for (const auto &e : get_out_edges(orig_g, bnode)) {
				int to = get_target(orig_g, e);
				const auto &to_ver = get_vertex_props(orig_g, to);
				//if () {
				//assert(to_ver.type != SOURCE);
				//if (to_ver.type == OPIN) {
				//pnet->pseudo_sources.push_back(to);
				//pnet->pseudo_sources_prev_edge.push_back(e);
				//}
				//} else if (pid[to] == to_pid) {
				if (pid[to] == -1) {
					assert(to_ver.type == IPIN);
				} else if (pid[to] == net_next_pid) {
					assert(to_ver.type == CHANX || to_ver.type == CHANY); 
					if (added_pseudo_sources.find(to) == added_pseudo_sources.end()) {
						pseudo_source_t psource;
						psource.node = to;
						psource.prev_edge = e;
						psource.bnode = &b;
						pnet.pseudo_sources.emplace_back(psource);
						added_pseudo_sources.insert(to);
						zlog_level(delta_log, ROUTER_V3, "\tAdding node %d\n", to);
					} else {
						zlog_level(delta_log, ROUTER_V3, "\tNode %d already added\n", to);
					}
				} else {
					zlog_level(delta_log, ROUTER_V3, "\tNode %d is in partition %d\n", to, pid[to]);
				}
			}
		} 

		if (pnet.pseudo_sources.empty()) {
			net_next_pid = (net_next_pid+1) % num_partitions;
		}
	} while (pnet.pseudo_sources.empty());
}

pseudo_net_t *get_pseudo_net(unrouted_t &unrouted, net_t *net, const RRGraph &orig_g, const vector<int> &pid, int &net_next_pid, int num_partitions)
{
	pseudo_net_t *pnet = new pseudo_net_t;

	pnet->net = net;

	for (auto &u : unrouted.unrouted_sinks) {
		pseudo_sink_t psink;

		psink.sink = u;
		psink.pseudo_source = nullptr;

		pnet->pseudo_sinks.emplace_back(psink);
	}

	update_pseudo_sources(*pnet, unrouted, orig_g, pid, num_partitions, net_next_pid);

	return pnet;	
}

void get_sinks_to_route(net_t *net, const route_tree_t &rt, const vector<sink_t *> &unroutable_sinks, vector<sink_t *> &sinks_to_route)
{
	sinks_to_route.reserve(net->sinks.size());
	for (auto &sink : net->sinks) {
		RouteTreeNode sink_rt_node = route_tree_get_rt_node(rt, sink.rr_node);
		bool routable = find(begin(unroutable_sinks), end(unroutable_sinks), &sink) == end(unroutable_sinks);

		if (sink_rt_node == RouteTree::null_vertex()) {
			if (routable) {
				sinks_to_route.push_back(&sink);
			}
		} else {
			assert(!get_vertex_props(rt.graph, sink_rt_node).pending_rip_up);
		}
	}
}

void get_sinks_to_route(net_t *net, const route_tree_t &rt, vector<sink_t *> &sinks)
{
	sinks.reserve(net->sinks.size());
	for (auto &sink : net->sinks) {
		RouteTreeNode sink_rt_node = route_tree_get_rt_node(rt, sink.rr_node);
		if (sink_rt_node == RouteTree::null_vertex()) {
			sinks.push_back(&sink);
		} else {
			assert(!get_vertex_props(rt.graph, sink_rt_node).pending_rip_up);
		}
	}
}

void init_route_structs(const RRGraph &g, const vector<net_t> &nets, const vector<net_t> &global_nets, route_state_t **states, congestion_t **congestion, vector<route_tree_t> &route_trees, t_net_timing **net_timing)
{
	*states = new route_state_t[num_vertices(g)];
	for (int i = 0; i < num_vertices(g); ++i) {
		(*states)[i].rr_node = -1;
		(*states)[i].known_cost = std::numeric_limits<float>::max();
		(*states)[i].cost = std::numeric_limits<float>::max();
		(*states)[i].prev_edge = RRGraph::null_edge();
		(*states)[i].upstream_R = -1;
		(*states)[i].delay = std::numeric_limits<float>::max();
	}

	*congestion = new congestion_t[num_vertices(g)];
    for (int i = 0; i < num_vertices(g); ++i) {
        (*congestion)[i].acc_cost = 1;
        (*congestion)[i].pres_cost = 1;
        (*congestion)[i].occ = 0;
    }

	route_trees.resize(nets.size());
    for (int i = 0; i < nets.size(); ++i) {
        route_tree_init(route_trees[i]);
    }

    *net_timing = new t_net_timing[nets.size()+global_nets.size()];
    init_net_timing(nets, global_nets, *net_timing);
}

void init_route_structs_mpi(const RRGraph &g, const vector<net_t> &nets, const vector<net_t> &global_nets, route_state_t **states, congestion_t **congestion, MPI_Win *win, vector<route_tree_t> &route_trees, t_net_timing **net_timing)
{
	*states = new route_state_t[num_vertices(g)];
	for (int i = 0; i < num_vertices(g); ++i) {
		(*states)[i].rr_node = -1;
		(*states)[i].known_cost = std::numeric_limits<float>::max();
		(*states)[i].cost = std::numeric_limits<float>::max();
		(*states)[i].prev_edge = RRGraph::null_edge();
		(*states)[i].upstream_R = -1;
		(*states)[i].delay = std::numeric_limits<float>::max();
	}

	assert(MPI_Win_allocate(sizeof(congestion_t)*num_vertices(g), 1, MPI_INFO_NULL, MPI_COMM_WORLD, congestion, win) == MPI_SUCCESS);
    for (int i = 0; i < num_vertices(g); ++i) {
        (*congestion)[i].acc_cost = 1;
        (*congestion)[i].pres_cost = 1;
        (*congestion)[i].occ = 0;
    }

	route_trees.resize(nets.size());
    for (int i = 0; i < nets.size(); ++i) {
        route_tree_init(route_trees[i]);
    }

    *net_timing = new t_net_timing[nets.size()+global_nets.size()];
    init_net_timing(nets, global_nets, *net_timing);
}

void init_route_structs_locked(const RRGraph &g, const vector<net_t> &nets, const vector<net_t> &global_nets, int num_threads, vector<route_state_t *> &states, congestion_locked_t **congestion, vector<route_tree_t> &route_trees, t_net_timing **net_timing)
{
    states.resize(num_threads);
    for (int i = 0; i < num_threads; ++i) {
        states[i] = new route_state_t[num_vertices(g)];
        for (int j = 0; j < num_vertices(g); ++j) {
            states[i][j].rr_node = -1;
            states[i][j].known_cost = std::numeric_limits<float>::max();
            states[i][j].cost = std::numeric_limits<float>::max();
            states[i][j].prev_edge = RRGraph::null_edge();
            states[i][j].upstream_R = -1;
            states[i][j].delay = std::numeric_limits<float>::max();
        }
    }

    *congestion = new congestion_locked_t[num_vertices(g)];
    for (int i = 0; i < num_vertices(g); ++i) {
        (*congestion)[i].cong.acc_cost = 1;
        (*congestion)[i].cong.pres_cost = 1;
        (*congestion)[i].cong.occ = 0;
    }

	route_trees.resize(nets.size());
    for (int i = 0; i < nets.size(); ++i) {
        route_tree_init(route_trees[i]);
    }

    *net_timing = new t_net_timing[nets.size()+global_nets.size()];
    init_net_timing(nets, global_nets, *net_timing);
}

static MPI_Datatype net_timing_dt;
static MPI_Datatype recalc_occ_dt;
static MPI_Datatype pn_send_dt;
static MPI_Datatype net_sink_dt;
static MPI_Datatype net_dt;

MPI_Datatype get_net_timing_dt()
{
	return net_timing_dt;
}

void init_datatypes()
{
	int bl[] = { 1, 1 };
	MPI_Aint disp[] = { offsetof(t_net_timing, delay), sizeof(congestion_t) };
	MPI_Datatype dts[] = { MPI_INT, MPI_UB };

	assert(MPI_Type_create_struct(1, bl, disp, dts, &net_timing_dt) == MPI_SUCCESS);
	assert(MPI_Type_commit(&net_timing_dt) == MPI_SUCCESS);

	int recalc_bl[] = { 1, 1, 1 };
	MPI_Aint recalc_disp[] = { offsetof(congestion_t, recalc_occ), 0, sizeof(congestion_t) };
	MPI_Datatype recalc_dts[] = { MPI_INT, MPI_LB, MPI_UB };

	assert(MPI_Type_create_struct(3, recalc_bl, recalc_disp, recalc_dts, &recalc_occ_dt) == MPI_SUCCESS);
	assert(MPI_Type_commit(&recalc_occ_dt) == MPI_SUCCESS);

	//assert(MPI_Type_contiguous(4, MPI_INT, &pn_send_dt) == MPI_SUCCESS);
	//assert(MPI_Type_commit(&pn_send_dt) == MPI_SUCCESS);
	//MPI_Aint lb, extent;
	//MPI_Type_get_extent(pn_send_dt, &lb, &extent);
	//assert(lb == 0 && extent == sizeof(path_node_send_t));
	//
	int pn_send_bl[] = { 2, 2 };
	MPI_Aint pn_send_disp[] = { offsetof(path_node_send_t, rr_node_id), offsetof(path_node_send_t, delay) };
	MPI_Datatype pn_send_dts[] = { MPI_INT, MPI_FLOAT };

	assert(MPI_Type_create_struct(2, pn_send_bl, pn_send_disp, pn_send_dts, &pn_send_dt) == MPI_SUCCESS);
	assert(MPI_Type_commit(&pn_send_dt) == MPI_SUCCESS);
	MPI_Aint lb, extent;
	MPI_Type_get_extent(pn_send_dt, &lb, &extent);
	assert(lb == 0 && extent == sizeof(path_node_send_t));

	int net_bl[] = {
		1,
		1, 1, 1,
		1, 1, 1, 1 };
	MPI_Aint net_disp[] = {
		offsetof(net_t, vpr_id),

		/* source */
		offsetof(net_t, source)+offsetof(source_t, rr_node),
		offsetof(net_t, source)+offsetof(source_t, x),
		offsetof(net_t, source)+offsetof(source_t, y),
		
		/* bounding box */
		offsetof(net_t, bounding_box)+offsetof(bounding_box_t, xmin),
		offsetof(net_t, bounding_box)+offsetof(bounding_box_t, xmax),
		offsetof(net_t, bounding_box)+offsetof(bounding_box_t, ymin),
		offsetof(net_t, bounding_box)+offsetof(bounding_box_t, ymax)
	};
	MPI_Aint net_disp_2[8];
	MPI_Aint base;
	net_t tmp_net;
	MPI_Get_address(&tmp_net, &base);
	MPI_Get_address(&tmp_net.vpr_id, &net_disp_2[0]);
	MPI_Get_address(&tmp_net.source.rr_node, &net_disp_2[1]);
	MPI_Get_address(&tmp_net.source.x, &net_disp_2[2]);
	MPI_Get_address(&tmp_net.source.y, &net_disp_2[3]);
	//MPI_Get_address(&tmp_net.bounding_box.xmin, &net_disp_2[4]);
	//MPI_Get_address(&tmp_net.bounding_box.xmax, &net_disp_2[5]);
	//MPI_Get_address(&tmp_net.bounding_box.ymin, &net_disp_2[6]);
	//MPI_Get_address(&tmp_net.bounding_box.ymax, &net_disp_2[7]);
	/* TODO: fix the bounding box offsets */
	assert(false);
	for (int i = 0; i < 8; ++i) {
		net_disp_2[i] -= base;
		assert(net_disp_2[i] == net_disp[i]);
	}
	MPI_Datatype net_dts[] = {
		MPI_INT,
		MPI_INT, MPI_INT, MPI_INT,
		MPI_INT, MPI_INT, MPI_INT, MPI_INT };

	MPI_Datatype net_tmp_dt;
	assert(MPI_Type_create_struct(8, net_bl, net_disp, net_dts, &net_tmp_dt) == MPI_SUCCESS);
	assert(MPI_Type_commit(&net_tmp_dt) == MPI_SUCCESS);
	assert(MPI_Type_create_resized(net_tmp_dt, net_disp[0], sizeof(net_t), &net_dt) == MPI_SUCCESS);
	assert(MPI_Type_commit(&net_dt) == MPI_SUCCESS);

	int net_sink_bl[] = {
		1, 1, 1,
		1, 1, 1, 1
	};
	MPI_Aint net_sink_disp[] = {
		offsetof(sink_t, rr_node),
		offsetof(sink_t, x),
		offsetof(sink_t, y),

		offsetof(sink_t, current_bounding_box)+offsetof(bounding_box_t, xmin),
		offsetof(sink_t, current_bounding_box)+offsetof(bounding_box_t, xmax),
		offsetof(sink_t, current_bounding_box)+offsetof(bounding_box_t, ymin),
		offsetof(sink_t, current_bounding_box)+offsetof(bounding_box_t, ymax)
	};
	MPI_Aint net_sink_disp_2[7];
	sink_t tmp_sink;
	MPI_Get_address(&tmp_sink, &base);
	MPI_Get_address(&tmp_sink.rr_node, &net_sink_disp_2[0]);
	MPI_Get_address(&tmp_sink.x, &net_sink_disp_2[1]);
	MPI_Get_address(&tmp_sink.y, &net_sink_disp_2[2]);
	MPI_Get_address(&tmp_sink.current_bounding_box.xmin, &net_sink_disp_2[3]);
	MPI_Get_address(&tmp_sink.current_bounding_box.xmax, &net_sink_disp_2[4]);
	MPI_Get_address(&tmp_sink.current_bounding_box.ymin, &net_sink_disp_2[5]);
	MPI_Get_address(&tmp_sink.current_bounding_box.ymax, &net_sink_disp_2[6]);
	for (int i = 0; i < 7; ++i) {
		net_sink_disp_2[i] -= base;
		assert(net_sink_disp_2[i] == net_sink_disp[i]);
	}
	MPI_Datatype net_sink_dts[] = {
		MPI_INT, MPI_INT, MPI_INT,
		MPI_INT, MPI_INT, MPI_INT, MPI_INT };

	MPI_Datatype net_sink_tmp_dt;
	assert(MPI_Type_create_struct(7, net_sink_bl, net_sink_disp, net_sink_dts, &net_sink_tmp_dt) == MPI_SUCCESS);
	assert(MPI_Type_commit(&net_sink_tmp_dt) == MPI_SUCCESS);
	assert(MPI_Type_create_resized(net_sink_tmp_dt, net_sink_disp[0], sizeof(sink_t), &net_sink_dt) == MPI_SUCCESS);
	assert(MPI_Type_commit(&net_sink_dt) == MPI_SUCCESS);
}

void recalc_sum(congestion_t *in, congestion_t *inout, int *len, MPI_Datatype *dptr)
{
	for (int i = 0; i < *len; ++i) {
		inout[i].recalc_occ = in[i].recalc_occ + inout[i].recalc_occ;
	}
}

void init_displ_nets(const vector<vector<net_t*>> &partitions, int **recvcounts, int **displs)
{
	*recvcounts = new int[partitions.size()];
	for (int pi = 0; pi < partitions.size(); ++pi) {
		(*recvcounts)[pi] = partitions[pi].size();
	}

	*displs = new int[partitions.size()];
	for (int i = 0; i < partitions.size(); ++i) {
		(*displs)[i] = 0;
		for (int j = 0; j < i; ++j) {
			(*displs)[i] += (*recvcounts)[j];
		}
	}
}

void init_displ(const vector<vector<net_t*>> &partitions, int **recvcounts, int **displs)
{
	*recvcounts = new int[partitions.size()];
	for (int i = 0; i < partitions.size(); ++i) {
		(*recvcounts)[i] = 0;
	}

	for (int pi = 0; pi < partitions.size(); ++pi) {
		for (auto &net : partitions[pi]) {
			(*recvcounts)[pi] += net->sinks.size();
		}
	}

	*displs = new int[partitions.size()];
	for (int i = 0; i < partitions.size(); ++i) {
		(*displs)[i] = 0;
		for (int j = 0; j < i; ++j) {
			(*displs)[i] += (*recvcounts)[j];
		}
	}
}

void init_displ(int num_procs, int current_level, const vector<pair<box, net_t *>> &nets_to_route, int initial_num_procs, int **recvcounts, int **displs)
{
	*recvcounts = new int[num_procs];
	for (int i = 0; i < num_procs; ++i) {
		(*recvcounts)[i] = 0;
	}
	for (int pi = 0; pi < num_procs; ++pi) {
		for (int i = pi*pow(2, current_level); i < nets_to_route.size(); i += initial_num_procs) {
			for (int j = 0; j < pow(2, current_level) && i+j < nets_to_route.size(); ++j) {
				net_t *net = nets_to_route[i+j].second;
				(*recvcounts)[pi] += net->sinks.size();
			}
		}
	}

	*displs = new int[num_procs];
	for (int i = 0; i < num_procs; ++i) {
		(*displs)[i] = 0;
		for (int j = 0; j < i; ++j) {
			(*displs)[i] += (*recvcounts)[j];
		}
	}
}

void sync_net_delay(const vector<vector<net_t *>> &partitions, int procid, int num_procs, int *recvcounts, int *displs, MPI_Comm comm, t_net_timing *net_timing)
{
	assert(procid >= 0 && procid < partitions.size());
	assert(partitions.size() == num_procs);

	int num_delays = 0;
	for (auto &net : partitions[procid]) {
		num_delays += net->sinks.size();
	}

	float *delays = new float[num_delays];
	int idx = 0;
	for (auto &net : partitions[procid]) {
		zlog_level(delta_log, ROUTER_V3, "Net vpd id %d send delays:\n", net->vpr_id);
		for (int k = 1; k <= net->sinks.size(); ++k) {
			delays[idx] = net_timing[net->vpr_id].delay[k];
			zlog_level(delta_log, ROUTER_V3, "\t%g\n", delays[idx]);
			++idx;
		}
	}

	assert(idx == num_delays);
	assert(idx == recvcounts[procid]);

	int total_num_sinks = 0;
	for (int i = 0; i < num_procs; ++i) {
		total_num_sinks += recvcounts[i];
	}

	float *all_delays = new float[total_num_sinks];

	MPI_Allgatherv(delays, recvcounts[procid], MPI_FLOAT, all_delays, recvcounts, displs, MPI_FLOAT, comm);

	for (int pi = 0; pi < num_procs; ++pi) {
		idx = 0;
		for (auto &net : partitions[pi]) {
			zlog_level(delta_log, ROUTER_V3, "Net vpr id %d recv delays:\n", net->vpr_id);
			for (int k = 1; k <= net->sinks.size(); ++k) {
				if (pi == procid) {
					assert(net_timing[net->vpr_id].delay[k] == all_delays[displs[pi] + idx]);
				} else {
					net_timing[net->vpr_id].delay[k] = all_delays[displs[pi] + idx];
				}
				zlog_level(delta_log, ROUTER_V3, "\t%g\n", net_timing[net->vpr_id].delay[k]);
				++idx;
			}
		}
		assert(idx == recvcounts[pi]);
	}

	delete [] delays;
	delete [] all_delays;
}

void sync_net_delay(const vector<pair<box, net_t *>> &nets_to_route, int procid, int num_procs, int initial_num_procs, int *recvcounts, int *displs, int current_level, MPI_Comm comm, t_net_timing *net_timing)
{
	int num_delays = 0;
	for (int i = procid*pow(2, current_level); i < nets_to_route.size(); i += initial_num_procs) {
		for (int j = 0; j < pow(2, current_level) && i+j < nets_to_route.size(); ++j) {
			net_t *net = nets_to_route[i+j].second;
			num_delays += net->sinks.size();
		}
	}

	float *delays = new float[num_delays];
	int idx = 0;
	for (int i = procid*pow(2, current_level); i < nets_to_route.size(); i += initial_num_procs) {
		for (int j = 0; j < pow(2, current_level) && i+j < nets_to_route.size(); ++j) {
			net_t *net = nets_to_route[i+j].second;

			zlog_level(delta_log, ROUTER_V3, "Net index %d vpd id %d send delays:\n", i+j, net->vpr_id);
			for (int k = 1; k <= net->sinks.size(); ++k) {
				delays[idx] = net_timing[net->vpr_id].delay[k];
				zlog_level(delta_log, ROUTER_V3, "\t%g\n", delays[idx]);
				++idx;
			}
		}
	}

	assert(idx == num_delays);
	assert(idx == recvcounts[procid]);

	int total_num_sinks = 0;
	for (int i = 0; i < num_procs; ++i) {
		total_num_sinks += recvcounts[i];
	}

	float *all_delays = new float[total_num_sinks];

	MPI_Allgatherv(delays, recvcounts[procid], MPI_FLOAT, all_delays, recvcounts, displs, MPI_FLOAT, comm);

	for (int pi = 0; pi < num_procs; ++pi) {
		idx = 0;
		for (int i = pi*pow(2, current_level); i < nets_to_route.size(); i += initial_num_procs) {
			for (int j = 0; j < pow(2, current_level) && i+j < nets_to_route.size(); ++j) {
				net_t *net = nets_to_route[i+j].second;
				zlog_level(delta_log, ROUTER_V3, "Net index %d vpr id %d recv delays:\n", i+j, net->vpr_id);
				for (int k = 1; k <= net->sinks.size(); ++k) {
					if (pi == procid) {
						assert(net_timing[net->vpr_id].delay[k] == all_delays[displs[pi] + idx]);
					} else {
						net_timing[net->vpr_id].delay[k] = all_delays[displs[pi] + idx];
					}
					zlog_level(delta_log, ROUTER_V3, "\t%g\n", net_timing[net->vpr_id].delay[k]);
					++idx;
				}
			}
		}
		assert(idx == recvcounts[pi]);
	}

	delete [] delays;
	delete [] all_delays;
}

void recv_route_tree(const net_t *net, const RRGraph &g, vector<route_tree_t> &route_trees, t_net_timing *net_timing, int from_procid, MPI_Comm comm)
{
	assert(route_tree_empty(route_trees[net->local_id]));

	MPI_Status status;

	MPI_Probe(from_procid, net->local_id, comm, &status);

	int count;
	MPI_Get_count(&status, pn_send_dt, &count);

	assert(count > 0);
	assert(status.MPI_TAG == net->local_id);

	vector<path_node_send_t> recv(count);
	MPI_Recv(recv.data(), count, pn_send_dt, status.MPI_SOURCE, status.MPI_TAG, comm, MPI_STATUS_IGNORE);

	char buffer[256];
	zlog_level(delta_log, ROUTER_V3, "Recv route tree of net %d\n", net->vpr_id);

	int num_sources = 0;
	int num_sinks = 0;
	RRNode source_rr_node = RRGraph::null_vertex();
	for (int i = 0; i < count; ++i) {
		int rr_node = recv[i].rr_node_id;

		//RREdge rr_edge_to_parent;
		//if (recv[i].parent_rr_node == RRGraph::null_vertex()) {
			//rr_edge_to_parent = RRGraph::null_edge();
		//} else {
			//rr_edge_to_parent = get_edge(g, recv[i].parent_rr_node, rr_node);
		//}

		const auto &rr_node_p = get_vertex_props(g, rr_node);

		RouteTreeNode rt_node = route_tree_add_rr_node(route_trees[net->local_id], rr_node, g);
		assert(rt_node != RouteTree::null_vertex());

		route_tree_set_node_properties(route_trees[net->local_id], rt_node, rr_node_p.type != IPIN && rr_node_p.type != SINK, recv[i].upstream_R, recv[i].delay);

		if (rr_node_p.type == SOURCE) {
			route_tree_add_root(route_trees[net->local_id], rr_node);
			source_rr_node = rr_node;
			++num_sources;
		} else if (rr_node_p.type == SINK) {
			++num_sinks;

			int num_matches = 0;
			int match = -1;
			for (int j = 0; j < net->sinks.size(); ++j) {
				if (net->sinks[j].rr_node == rr_node) {
					match = j;
					++num_matches;
				}
			}
			assert(num_matches == 1);

			net_timing[net->vpr_id].delay[net->sinks[match].id+1] = recv[i].delay;
		}

		//update_one_cost_internal(rr_node, g, congestion, 1, pres_fac);

		sprintf_rr_node(rr_node, buffer);
		zlog_level(delta_log, ROUTER_V3, "\t%s parent: %d delay: %g upstream_R: %g\n", buffer, recv[i].parent_rr_node, recv[i].delay, recv[i].upstream_R);
	}

	assert(num_sources == 1);

	/* set up edges */
	for (int i = 0; i < count; ++i) {
		int rr_node = recv[i].rr_node_id;

		if (recv[i].parent_rr_node != RRGraph::null_vertex()) {
			route_tree_add_edge_between_rr_node(route_trees[net->local_id], recv[i].parent_rr_node, rr_node);
		} 
	}

	/* hacky fix for titan benchmark */
	//RouteTreeNode source_rt_node = route_tree_get_rt_node(route_trees[net->local_id], source_rr_node);
	//assert(source_rt_node != RouteTree::null_vertex());
	//int delta = num_out_edges(route_trees[net->local_id].graph, source_rt_node)-1;
	//update_one_cost_internal(source_rr_node, g, congestion, delta, pres_fac);
}

void recv_route_tree(net_t *net, const RRGraph &g, vector<vector<sink_t *>> &routed_sinks, vector<route_tree_t> &route_trees, t_net_timing *net_timing, int from_procid, MPI_Comm comm)
{
	MPI_Status status;
	MPI_Probe(from_procid, net->local_id, comm, &status);

	int num_routed_sinks;
	MPI_Get_count(&status, MPI_INT, &num_routed_sinks);

	assert(num_routed_sinks > 0);

	int *routed_sink_ids = new int[num_routed_sinks];
	MPI_Recv(routed_sink_ids, num_routed_sinks, MPI_INT, status.MPI_SOURCE, status.MPI_TAG, comm, MPI_STATUS_IGNORE);

	assert(routed_sinks[net->local_id].empty());

	for (int i = 0; i < num_routed_sinks; ++i) {
		int num_matches = 0;
		int match = -1;
		for (int j = 0; j < net->sinks.size(); ++j) {
			if (net->sinks[j].id == routed_sink_ids[i]) {
				match = j;
				++num_matches;
			}
		}
		assert(num_matches == 1);
		routed_sinks[net->local_id].push_back(&net->sinks[match]);
	}

	delete [] routed_sink_ids;

	if (routed_sinks[net->local_id].empty()) {
		zlog_level(delta_log, ROUTER_V3, "Recving, net %d has no routed sinks\n", net->vpr_id);
	}

	assert(route_tree_empty(route_trees[net->local_id]));

	MPI_Probe(from_procid, net->local_id, comm, &status);

	int count;
	MPI_Get_count(&status, pn_send_dt, &count);

	assert(count > 0);
	assert(status.MPI_TAG == net->local_id);

	path_node_send_t *recv = new path_node_send_t[count];
	MPI_Recv(recv, count, pn_send_dt, status.MPI_SOURCE, status.MPI_TAG, comm, MPI_STATUS_IGNORE);

	char buffer[256];
	zlog_level(delta_log, ROUTER_V3, "Recv route tree of net %d\n", net->vpr_id);

	int num_sources = 0;
	int num_sinks = 0;
	RRNode source_rr_node = RRGraph::null_vertex();
	for (int i = 0; i < count; ++i) {
		int rr_node = recv[i].rr_node_id;

		RREdge rr_edge_to_parent;
		if (recv[i].parent_rr_node == RRGraph::null_vertex()) {
			rr_edge_to_parent = RRGraph::null_edge();
		} else {
			rr_edge_to_parent = get_edge(g, recv[i].parent_rr_node, rr_node);
		}

		const auto &rr_node_p = get_vertex_props(g, rr_node);

		RouteTreeNode rt_node = route_tree_add_rr_node(route_trees[net->local_id], rr_node, g);
		assert(rt_node != RouteTree::null_vertex());

		route_tree_set_node_properties(route_trees[net->local_id], rt_node, rr_node_p.type != IPIN && rr_node_p.type != SINK, recv[i].upstream_R, recv[i].delay);

		if (rr_node_p.type == SOURCE) {
			route_tree_add_root(route_trees[net->local_id], rr_node);
			source_rr_node = rr_node;
			++num_sources;
		} else if (rr_node_p.type == SINK) {
			++num_sinks;

			auto iter = find_if(begin(routed_sinks[net->local_id]), end(routed_sinks[net->local_id]), [&rr_node] (const sink_t *sink) -> bool { return sink->rr_node == rr_node; });
			assert(iter != end(routed_sinks[net->local_id]));

			net_timing[net->vpr_id].delay[(*iter)->id+1] = recv[i].delay;
		}

		//update_one_cost_internal(rr_node, g, congestion, 1, pres_fac);

		sprintf_rr_node(rr_node, buffer);
		zlog_level(delta_log, ROUTER_V3, "\t%s parent: %d delay: %g upstream_R: %g\n", buffer, recv[i].parent_rr_node, recv[i].delay, recv[i].upstream_R);
	}

	assert(num_sources == 1);
	assert(num_sinks == routed_sinks[net->local_id].size());

	/* set up edges */
	for (int i = 0; i < count; ++i) {
		int rr_node = recv[i].rr_node_id;

		if (recv[i].parent_rr_node != RRGraph::null_vertex()) {
			route_tree_add_edge_between_rr_node(route_trees[net->local_id], recv[i].parent_rr_node, rr_node);
		} 
	}

	/* hacky fix for titan benchmark */
	//RouteTreeNode source_rt_node = route_tree_get_rt_node(route_trees[net->local_id], source_rr_node);
	//assert(source_rt_node != RouteTree::null_vertex());
	//int delta = num_out_edges(route_trees[net->local_id].graph, source_rt_node)-1;
	//update_one_cost_internal(source_rr_node, g, congestion, delta, pres_fac);

	delete [] recv;
}

void send_route_tree(const net_t *net, const vector<route_tree_t> &route_trees, int to_procid, MPI_Comm comm)
{
	zlog_level(delta_log, ROUTER_V3, "Sending route tree of net %d\n", net->vpr_id);
	
	int num_rt_nodes = 0;
	for (const auto &rt_node : route_tree_get_nodes(route_trees[net->local_id])) {
		++num_rt_nodes;
	}

	vector<path_node_send_t> send(num_rt_nodes);

	int i = 0;
	for (const auto &rt_node : route_tree_get_nodes(route_trees[net->local_id])) {
		const auto &rt_node_p = get_vertex_props(route_trees[net->local_id].graph, rt_node);

		send[i].rr_node_id = rt_node_p.rr_node;
		if (valid(rt_node_p.rt_edge_to_parent)) {
			send[i].parent_rr_node = get_vertex_props(route_trees[net->local_id].graph, get_source(route_trees[net->local_id].graph, rt_node_p.rt_edge_to_parent)).rr_node;
		} else {
			send[i].parent_rr_node = RRGraph::null_vertex();
		}
		send[i].delay = rt_node_p.delay;
		send[i].upstream_R = rt_node_p.upstream_R;

		char buffer[256];
		sprintf_rr_node(send[i].rr_node_id, buffer);
		
		zlog_level(delta_log, ROUTER_V3, "\t%s parent: %d delay: %g upstream_R: %g\n", buffer, send[i].parent_rr_node, send[i].delay, send[i].upstream_R);
		
		++i;
	}

	MPI_Send(send.data(), num_rt_nodes, pn_send_dt, to_procid, net->local_id, comm);
}

void send_route_tree(const net_t *net, const vector<vector<sink_t *>> &routed_sinks, const vector<route_tree_t> &route_trees, int to_procid, MPI_Comm comm)
{
	if (routed_sinks[net->local_id].empty()) {
		zlog_level(delta_log, ROUTER_V3, "Sending, net %d has no routed sinks\n", net->vpr_id);
	}

	int *routed_sink_ids = new int[routed_sinks[net->local_id].size()];
	for (int i = 0; i < routed_sinks[net->local_id].size(); ++i) {
		routed_sink_ids[i] = routed_sinks[net->local_id][i]->id;
	}
	MPI_Send(routed_sink_ids, routed_sinks[net->local_id].size(), MPI_INT, to_procid, net->local_id, comm);
	delete [] routed_sink_ids;

	zlog_level(delta_log, ROUTER_V3, "Sending route tree of net %d\n", net->vpr_id);
	
	int num_rt_nodes = 0;
	for (const auto &rt_node : route_tree_get_nodes(route_trees[net->local_id])) {
		++num_rt_nodes;
	}

	path_node_send_t *send = new path_node_send_t[num_rt_nodes];

	int i = 0;
	for (const auto &rt_node : route_tree_get_nodes(route_trees[net->local_id])) {
		const auto &rt_node_p = get_vertex_props(route_trees[net->local_id].graph, rt_node);

		send[i].rr_node_id = rt_node_p.rr_node;
		if (valid(rt_node_p.rt_edge_to_parent)) {
			send[i].parent_rr_node = get_vertex_props(route_trees[net->local_id].graph, get_source(route_trees[net->local_id].graph, rt_node_p.rt_edge_to_parent)).rr_node;
		} else {
			send[i].parent_rr_node = RRGraph::null_vertex();
		}
		send[i].delay = rt_node_p.delay;
		send[i].upstream_R = rt_node_p.upstream_R;

		char buffer[256];
		sprintf_rr_node(send[i].rr_node_id, buffer);
		
		zlog_level(delta_log, ROUTER_V3, "\t%s parent: %d delay: %g upstream_R: %g\n", buffer, send[i].parent_rr_node, send[i].delay, send[i].upstream_R);
		
		++i;
	}

	MPI_Send(send, num_rt_nodes, pn_send_dt, to_procid, net->local_id, comm);

	delete [] send;
}


//void recv_route_tree(net_t *net, const RRGraph &g, vector<vector<sink_t *>> &routed_sinks, route_state_t *states, congestion_t *congestion, float pres_fac, vector<route_tree_t> &route_trees, t_net_timing *net_timing, int from_procid, MPI_Comm comm)
//{
	//MPI_Status status;
	//MPI_Probe(from_procid, net->local_id, comm, &status);

	//int num_routed_sinks;
	//MPI_Get_count(&status, MPI_INT, &num_routed_sinks);

	//int *routed_sink_ids = new int[num_routed_sinks];
	//MPI_Recv(routed_sink_ids, num_routed_sinks, MPI_INT, status.MPI_SOURCE, status.MPI_TAG, comm, MPI_STATUS_IGNORE);

	//assert(routed_sinks[net->local_id].empty());
	//for (int i = 0; i < num_routed_sinks; ++i) {
		//int num_matches = 0;
		//int match = -1;
		//for (int j = 0; j < net->sinks.size(); ++j) {
			//if (net->sinks[j].id == routed_sink_ids[i]) {
				//match = j;
				//++num_matches;
			//}
		//}
		//assert(num_matches == 1);
		//routed_sinks[net->local_id].push_back(&net->sinks[match]);
	//}

	//delete [] routed_sink_ids;

	//if (routed_sinks[net->local_id].empty()) {
		//zlog_level(delta_log, ROUTER_V3, "Recving, net %d has no routed sinks\n", net->vpr_id);
	//}

	//assert(route_tree_empty(route_trees[net->local_id]));

	//for (int i = 0; i < routed_sinks[net->local_id].size(); ++i) {
		//const auto &rs = routed_sinks[net->local_id][i];

		//MPI_Status status;
		//MPI_Probe(from_procid, net->local_id, comm, &status);

		//int count;
		//MPI_Get_count(&status, pn_send_dt, &count);

		//assert(count > 0);
		//assert(status.MPI_TAG == net->local_id);

		//path_node_send_t *recv = new path_node_send_t[count];
		//MPI_Recv(recv, count, pn_send_dt, status.MPI_SOURCE, status.MPI_TAG, comm, MPI_STATUS_IGNORE);

		//char buffer[256];
		//sprintf_rr_node(rs->rr_node, buffer);
		//zlog_level(delta_log, ROUTER_V3, "Recv path of net %d sink %s\n", net->vpr_id, buffer);

		//const auto &path = make_shared<vector<path_node_t>>(count);
		//for (int j = 0; j < count; ++j) {
			//int rr_node = recv[j].rr_node_id;

			//(*path)[j].rr_node_id = rr_node;
			//(*path)[j].update_cost = recv[j].update_cost ? true : false;
			//if (recv[j].prev_edge_a == -1 && recv[j].prev_edge_b == -1) {
				//(*path)[j].prev_edge = RRGraph::null_edge();
			//} else {
				//(*path)[j].prev_edge = get_edge(g, recv[j].prev_edge_a, recv[j].prev_edge_b);
			//}

			//states[rr_node].delay = recv[j].delay;
			//states[rr_node].upstream_R = recv[j].upstream_R;

			//sprintf_rr_node(rr_node, buffer);
			//zlog_level(delta_log, ROUTER_V3, "\t%s update: %d prev_a: %d prev_b: %d delay: %g upstream_R: %g\n", buffer, recv[j].update_cost, recv[j].prev_edge_a, recv[j].prev_edge_b, recv[j].delay, recv[j].upstream_R);
		//}

		//if (i == 0) {
			//const auto &rr_node_p = get_vertex_props(g, recv[count-1].rr_node_id);
			//assert(rr_node_p.type == SOURCE);

			//RouteTreeNode rt_node = route_tree_add_rr_node(route_trees[net->local_id], recv[count-1].rr_node_id, g);
			//assert(rt_node != RouteTree::null_vertex());

			//route_tree_set_node_properties(get_vertex_props(route_trees[net->local_id].graph, rt_node), rr_node_p.type != IPIN && rr_node_p.type != SINK, (*path)[count-1].prev_edge, recv[count-1].upstream_R, recv[count-1].delay);

			//route_tree_add_root(route_trees[net->local_id], recv[count-1].rr_node_id);
		//}

		//route_tree_add_path(route_trees[net->local_id], path, g, states, false);

		//vector<RRNode> added_nodes;
		//for (const auto &n : *path) {
			//if (n.update_cost) {
				//added_nodes.push_back(n.rr_node_id);
			//}
		//}

		//update_one_cost(g, congestion, added_nodes.begin(), added_nodes.end(), 1, pres_fac);

		//assert(get_vertex_props(g, (*path)[0].rr_node_id).type == SINK);
		//assert(rs->rr_node == (*path)[0].rr_node_id);

		//net_timing[net->vpr_id].delay[rs->id+1] = states[rs->rr_node].delay;

		//delete [] recv;
	//}
//}

//void send_route_tree(net_t *net, const RRGraph &g, const vector<vector<sink_t *>> &routed_sinks, const vector<route_tree_t> &route_trees, int to_procid, MPI_Comm comm)
//{
	//if (routed_sinks[net->local_id].empty()) {
		//zlog_level(delta_log, ROUTER_V3, "Sending, net %d has no routed sinks\n", net->vpr_id);
	//}

	//int *routed_sink_ids = new int[routed_sinks[net->local_id].size()];
	//for (int i = 0; i < routed_sinks[net->local_id].size(); ++i) {
		//routed_sink_ids[i] = routed_sinks[net->local_id][i]->id;
	//}
	//MPI_Send(routed_sink_ids, routed_sinks[net->local_id].size(), MPI_INT, to_procid, net->local_id, comm);
	//delete [] routed_sink_ids;

	//for (const auto &rs : routed_sinks[net->local_id]) {
		//const auto &path = route_tree_get_path(route_trees[net->local_id], rs->rr_node);

		//char buffer[256];
		//sprintf_rr_node(rs->rr_node, buffer);
		//zlog_level(delta_log, ROUTER_V3, "Sending path of net %d sink %s\n", net->vpr_id, buffer);

		//path_node_send_t *send = new path_node_send_t[path->size()];
		//for (int i = 0; i < path->size(); ++i) {
			//send[i].rr_node_id = (*path)[i].rr_node_id;
			//send[i].update_cost = (*path)[i].update_cost ? 1 : 0;
			//send[i].prev_edge_a = get_source(g, (*path)[i].prev_edge);
			//send[i].prev_edge_b = get_target(g, (*path)[i].prev_edge);

			//RouteTreeNode rt_node = route_tree_get_rt_node(route_trees[net->local_id], send[i].rr_node_id);
			//const auto &rt_node_p = get_vertex_props(route_trees[net->local_id].graph, rt_node);
			//send[i].delay = rt_node_p.delay;
			//send[i].upstream_R = rt_node_p.upstream_R;

			//sprintf_rr_node(send[i].rr_node_id, buffer);
			//zlog_level(delta_log, ROUTER_V3, "\t%s update: %d prev_a: %d prev_b: %d delay: %g upstream_R: %g\n", buffer, send[i].update_cost, send[i].prev_edge_a, send[i].prev_edge_b, send[i].delay, send[i].upstream_R);
		//}

		//MPI_Send(send, path->size(), pn_send_dt, to_procid, net->local_id, comm);

		//delete [] send;
	//}
//}
//
void sync(congestion_t *congestion, const RRGraph &g, float pres_fac, int this_pid, int num_procs, MPI_Comm comm);

void free_circuit();
void free_cb(t_pb *pb);

static void free_blocks()
{
	extern struct s_block *block;
	extern int num_blocks;

	if (block != NULL) {
		for (int i = 0; i < num_blocks; i++) {
			if (block[i].pb != NULL) {
				free_cb(block[i].pb);
				free(block[i].pb);
			}
			free(block[i].nets);
			free(block[i].name);
		}
	}
	free(block);
	block = NULL;
}

void recv_sinks(vector<net_t> &nets, MPI_Comm comm)
{
	for (int i = 0; i < nets.size(); ++i) {
		zlog_level(delta_log, ROUTER_V3, "Recving net vpr id %d source rr %d x %d y %d bb %d %d %d %d\n", nets[i].vpr_id, nets[i].source.rr_node, nets[i].source.x, nets[i].source.y, bg::get<bg::min_corner, 0>(nets[i].bounding_box), bg::get<bg::max_corner, 0>(nets[i].bounding_box), bg::get<bg::min_corner, 1>(nets[i].bounding_box), bg::get<bg::max_corner, 1>(nets[i].bounding_box));

		int num_sinks;
		MPI_Bcast(&num_sinks, 1, MPI_INT, 0, comm);

		nets[i].sinks.resize(num_sinks);
		MPI_Bcast(nets[i].sinks.data(), num_sinks, net_sink_dt, 0, comm);

		zlog_level(delta_log, ROUTER_V3, "Net %d sink capacity %d\n", nets[i].vpr_id, nets[i].sinks.capacity()*sizeof(sink_t));

		for (int j = 0; j < num_sinks; ++j) {
			zlog_level(delta_log, ROUTER_V3, "\tRecving net sink rr %d x %d y %d bb %d %d %d %d\n", nets[i].sinks[j].rr_node, nets[i].sinks[j].x, nets[i].sinks[j].y, nets[i].sinks[j].current_bounding_box.xmin, nets[i].sinks[j].current_bounding_box.xmax, nets[i].sinks[j].current_bounding_box.ymin, nets[i].sinks[j].current_bounding_box.ymax);

			nets[i].sinks[j].id = j;
			nets[i].sinks[j].criticality_fac = std::numeric_limits<float>::max();
		}
	}
}

void send_sinks(vector<net_t> &nets, MPI_Comm comm)
{
	for (int i = 0; i < nets.size(); ++i) {
		int num_sinks = nets[i].sinks.size();

		assert(num_sinks > 0);

		zlog_level(delta_log, ROUTER_V3, "Sending net vpr id %d source rr %d x %d y %d bb %d %d %d %d\n", nets[i].vpr_id, nets[i].source.rr_node, nets[i].source.x, nets[i].source.y, bg::get<bg::min_corner, 0>(nets[i].bounding_box), bg::get<bg::max_corner, 0>(nets[i].bounding_box), bg::get<bg::min_corner, 1>(nets[i].bounding_box), bg::get<bg::max_corner, 1>(nets[i].bounding_box));

		MPI_Bcast(&num_sinks, 1, MPI_INT, 0, comm);
		MPI_Bcast(nets[i].sinks.data(), num_sinks, net_sink_dt, 0, comm);

		for (int j = 0; j < num_sinks; ++j) {
			zlog_level(delta_log, ROUTER_V3, "\tSending net sink rr %d x %d y %d bb %d %d %d %d\n", nets[i].sinks[j].rr_node, nets[i].sinks[j].x, nets[i].sinks[j].y, nets[i].sinks[j].current_bounding_box.xmin, nets[i].sinks[j].current_bounding_box.xmax, nets[i].sinks[j].current_bounding_box.ymin, nets[i].sinks[j].current_bounding_box.ymax);
		}
	}
}

void sync_nets(vector<net_t> &nets, vector<net_t> &global_nets, int procid, MPI_Comm comm)
{
	if (procid == 0) {
		int num_global_nets = global_nets.size();
		MPI_Bcast(&num_global_nets, 1, MPI_INT, 0, comm);
		MPI_Bcast(global_nets.data(), global_nets.size(), net_dt, 0, comm);

		printf("sent num_global_nets: %d\n", num_global_nets);
		
		int num_nets = nets.size();
		MPI_Bcast(&num_nets, 1, MPI_INT, 0, comm);
		MPI_Bcast(nets.data(), nets.size(), net_dt, 0, comm);

		printf("sent num_nets: %d\n", num_nets);

		printf("sent nets\n");

		send_sinks(nets, comm);
		send_sinks(global_nets, comm);
	} else {
		int num_global_nets;
		MPI_Bcast(&num_global_nets, 1, MPI_INT, 0, comm);
		global_nets.resize(num_global_nets);

		printf("recv num_global_nets: %d\n", num_global_nets);

		MPI_Bcast(global_nets.data(), global_nets.size(), net_dt, 0, comm);

		int num_nets;
		MPI_Bcast(&num_nets, 1, MPI_INT, 0, comm);
		nets.resize(num_nets);

		printf("recv num_nets: %d\n", num_nets);

		MPI_Bcast(nets.data(), nets.size(), net_dt, 0, comm);

		printf("recv nets\n");

		int local_id = 0;
		for (int i = 0; i < nets.size(); ++i) {
			nets[i].local_id = local_id++;
		}
		for (int i = 0; i < global_nets.size(); ++i) {
			global_nets[i].local_id = -1;
		}

		recv_sinks(nets, comm);
		recv_sinks(global_nets, comm);
	}
}

void sync_recalc_occ(congestion_t *congestion, int num_vertices, int procid, int num_procs, MPI_Comm comm)
{
	vector<int> all_recalc_occ(num_vertices);
	for (int i = 0; i < num_vertices; ++i) {
		all_recalc_occ[i] = congestion[i].recalc_occ;
	}

	assert(MPI_Allreduce(MPI_IN_PLACE, all_recalc_occ.data(), num_vertices, MPI_INT, MPI_SUM, comm) == MPI_SUCCESS);

	for (int i = 0; i < num_vertices; ++i) {
		congestion[i].recalc_occ = all_recalc_occ[i];
	}
}

int get_num_interpartition_nets(const vector<net_t> &nets, int num_partitions);

