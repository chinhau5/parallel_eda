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
#include "scheduler.h"
#include "geometry.h"
#include "quadtree.h"
#include "utility.h"
#include "cluster.h"
#include "args.h"
#include "init.h"
#include "router.h"
#include "congestion.h"
#include "partition.h"
#include "fm.h"

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

static void dump_rr_graph(const RRGraph &g, const char *filename)
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

enum class VertexColor {
	WHITE,
	GRAY,
	BLACK
};

template<typename Graph, typename Visitor>
//void bfs(const RRGraph &g, const RRNode &rr_node, bool horizontal, bool inc_direction, vector<bool> &visited, vector<int> &visited_edges) 
void bfs(const Graph &g, const vector<unsigned long> &rr_nodes, vector<VertexColor> &color, Visitor &visitor)
{
	//struct queue_item {
		//const RRNode *node;
		//const RREdge *edge;
	//};

	std::queue<int> s;
	for (auto rr_node : rr_nodes) {
		assert(color[rr_node] == VertexColor::WHITE);
		visitor.discover_vertex(rr_node, g);
		s.push(rr_node);
		color[rr_node] = VertexColor::GRAY;
	}

	//extern int nx, ny;
	//vector<vector<map<pair<int, int>, int>>> visit_state(nx+2, vector<map<pair<int, int>, int>>(ny+2));
	//for (int x = 0; x < nx+2; ++x) {
		//for (int y = 0; y < ny+2; ++y) {
			//visit_state[x][y]
		//}
	//}
	//
	//char buffer[256];
	//extern zlog_category_t *delta_log;

	while (!s.empty()) {
		int current = s.front();
		s.pop();
		visitor.examine_vertex(current, g);
		//int rr_node_id = id(*current);
		//sprintf_rr_node(rr_node_id, buffer);
		//zlog_level(delta_log, ROUTER_V3, "Current: %s Add: %d\n", buffer, item.add ? 1 : 0);

		//if (get_vertex_props(g, rr_node_id).type != CHANX && get_vertex_props(g, rr_node_id).type != CHANY) {assert(false);}
		//
		//bool visit_current = visit(*current);
		//if (item.add && visit_current) {
		//visited_nodes.push_back(rr_node_id);

		//if (item.edge) {
		//visited_edges.push_back(id(*item.edge));
		//}
		//}

		for (auto e : get_out_edges(g, current)) {
			visitor.examine_edge(e, g);
			int neighbor = get_target(g, e);
			//if (!((current.type == CHANX || current.type == CHANY) &&
			//(neighbor.type == CHANX || neighbor.type == CHANY))
			//|| current.inc_direction == neighbor.inc_direction) {
			//s.push({ &neighbor, &e });
			//}

			VertexColor vc = color[neighbor];
			if (vc == VertexColor::WHITE) {
				if (visitor.tree_edge(e, g)) {
					visitor.discover_vertex(neighbor, g);
					s.push(neighbor);
					color[neighbor] = VertexColor::GRAY;
				}
				//sprintf_rr_node(id(neighbor), buffer);
				//zlog_level(delta_log, ROUTER_V3, "\tNeighbor: %s Parent add: %d Parent visit: %d\n", buffer, item.add ? 1 : 0, visit_current ? 1 : 0);
			} else if (vc == VertexColor::GRAY) {
			} else {
				assert(vc == VertexColor::BLACK);
				//visitor.black_target(
			}

			//if (current.type == CHANX || current.type == CHANY) {
			//if ((horizontal && neighbor.type == CHANX) || (!horizontal && neighbor.type == CHANY)) {
			//if (neighbor.inc_direction == inc_direction) {
			//s.push({ &neighbor, &e });
			//} else {
			//printf("Not pushing");
			//}
			//} else {
			//s.push({ &neighbor, &e });
			//}
			//} else {
			//s.push({ &neighbor, &e });
			//}
		} 
		color[current] = VertexColor::BLACK;
	}
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

int get_rr_node_index(int x, int y, t_rr_type rr_type, int ptc,
		t_ivec *** L_rr_node_indices);

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
			
			e_p.buffered = switch_inf[si].buffered; 
			e_p.switch_delay = switch_inf[si].Tdel; 
			e_p.R = switch_inf[si].R; 
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
struct partition_vertex_predicate_t {
	const Graph &g;
	const vector<int> &pid;
	int p;
	
	partition_vertex_predicate_t(const Graph &g, const vector<int> &pid, int p)
		: g(g), pid(pid), p(p)
	{
	}

	bool operator()(unsigned long v) const {
		const auto &ver = get_vertex_props(g, v);
		return (ver.type == CHANX || ver.type == CHANY) && pid[v] == p;
	}
};

template<typename Predicate>
struct subgraph_vertex_predicate_t {
	const Predicate &vertex_valid;

	subgraph_vertex_predicate_t(const Predicate &vertex_valid) :
		vertex_valid(vertex_valid)
	{
	}

	bool operator()(unsigned long v) const
	{
		return vertex_valid(v);
	}
};

template<typename Graph, typename Predicate>
struct subgraph_edge_predicate_t {
	const Graph &g;
	const Predicate &vertex_valid;

	subgraph_edge_predicate_t(const Graph &g, const Predicate &vertex_valid) : g(g), vertex_valid(vertex_valid) {}

	bool operator()(const RREdge &e) const
	{
		int from = get_source(g, e);
		int to = get_target(g, e);

		return vertex_valid(from) && vertex_valid(to);
	}
};

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

class rr_graph_partitioner {
	private:
		int num_tracks;

	public:
		RRGraph orig_g;
		RRGraph undirected_orig_g;
		RRGraph channel_with_interior_g;
		RRGraph channel_without_interior_g;
		int num_levels;
		vector<int> result_pid;
        vector<vector<int>> result_pid_by_level;
		vector<vector<RRNode>> sink_in_nodes;
		vector<vector<RRNode>> ipin_in_nodes;

	private:
		vector<vector<vector<int>>> num_starting_tracks; // [type][x][y]
		vector<int> real_track_num; // [rr_node]

		int num_partitions;
		//vector<int> initial_partition;
		//int partition_index_offset;
		vector<vector<vector<vector<int>>>> num_tracks_in_partition; //[type][part][x][y]

		int total_num_chans;

	public:
		void init_graphs(t_router_opts *opts, struct s_det_routing_arch *det_routing_arch, t_direct_inf *directs, int num_directs, t_segment_inf *segment_inf, t_timing_inf *timing_inf)
		{
			extern int num_types;
			extern struct s_type_descriptor *type_descriptors;
			extern int nx, ny;
			extern struct s_grid_tile **grid;

			//free_rr_graph();

			//int warnings;

			//build_rr_graph(GRAPH_UNIDIR, num_types, type_descriptors, nx, ny, grid,
					//opts->fixed_channel_width, NULL, det_routing_arch->switch_block_type,
					//det_routing_arch->Fs, det_routing_arch->num_segment,
					//det_routing_arch->num_switch, segment_inf,
					//det_routing_arch->global_route_switch,
					//det_routing_arch->delayless_switch, *timing_inf,
					//det_routing_arch->wire_to_ipin_switch, opts->base_cost_type,
					//directs, num_directs, FALSE,
					//&warnings);

			//init_channel_only_graph(channel_with_interior_g);

			//dump_rr_graph(channel_with_interior_g, "/Volumes/DATA/rr_graph_channel_with_interior.txt");

			init_graph(orig_g, sink_in_nodes, ipin_in_nodes);

			init_sprintf_rr_node(&orig_g);

			dump_rr_graph(orig_g, "/Volumes/DATA/rr_graph.txt");

			undirected_orig_g = orig_g;

			for (const auto &e : get_edges(orig_g)) {
				int from = get_source(orig_g, e);
				int to = get_target(orig_g, e);
				add_edge(undirected_orig_g, to, from);
			}

			//free_rr_graph();
			//for (int i = 0; i < det_routing_arch->num_segment; ++i) {
				//for (int j = 0; j < segment_inf[i].sb_len; ++j) {
					//if (j != 0 && j != segment_inf[i].sb_len-1) {
						//segment_inf[i].sb[j] = FALSE;
					//}
				//}
			//}
			//build_rr_graph(GRAPH_UNIDIR, num_types, type_descriptors, nx, ny, grid,
					//opts->fixed_channel_width, NULL, det_routing_arch->switch_block_type,
					//det_routing_arch->Fs, det_routing_arch->num_segment,
					//det_routing_arch->num_switch, segment_inf,
					//det_routing_arch->global_route_switch,
					//det_routing_arch->delayless_switch, *timing_inf,
					//det_routing_arch->wire_to_ipin_switch, opts->base_cost_type,
					//directs, num_directs, FALSE,
					//&warnings);

			//init_channel_only_graph(channel_without_interior_g);

			//dump_rr_graph(channel_without_interior_g, "/Volumes/DATA/rr_graph_channel_without_interior.txt");
		}

		void load_real_track_numbers()
		{
			extern int nx, ny;
			extern t_ivec ***rr_node_indices;

			real_track_num.resize(num_vertices(orig_g), -1);

			num_starting_tracks.resize(2);
			for (int type = 0; type < 2; ++type) {
				num_starting_tracks[type].resize(nx+2);
				for (int x = 0; x < nx+1; ++x) {
					num_starting_tracks[type][x].resize(ny+2, 0);
					for (int y = 0; y < ny+1; ++y) {
						int num_tracks;
						if (type == 0) {
							num_tracks = rr_node_indices[CHANX][y][x].nelem;
						} else {
							num_tracks = rr_node_indices[CHANY][x][y].nelem;
						}
						for (int i = 0; i < num_tracks; ++i) {
							int rr_node = get_rr_node_index(x, y, type == 0 ? CHANX : CHANY, i, rr_node_indices);
							const auto &source = get_vertex_props(channel_without_interior_g, rr_node);

							if (starts_at(source, x, y)) {
								assert(real_track_num[rr_node] == -1);
								real_track_num[rr_node] = num_starting_tracks[type][x][y]++;
							} 
						}
					}
				}
			}
			int total_num_tracks = 0;
			for (int type = 0; type < 2; ++type) {
				for (int x = 0; x < nx+1; ++x) {
					for (int y = 0; y < ny+1; ++y) {
						total_num_tracks += num_starting_tracks[type][x][y];
					}
				}
			}
			assert(total_num_tracks == total_num_chans);
		}

		void init(t_router_opts *opts, struct s_det_routing_arch *det_routing_arch, t_direct_inf *directs, int num_directs, t_segment_inf *segment_inf, t_timing_inf *timing_inf)
		{
			extern int *chan_width_x, *chan_width_y;
			assert(chan_width_x[0] == chan_width_y[0]);
			num_tracks = chan_width_x[0];
			init_graphs(opts, det_routing_arch, directs, num_directs, segment_inf, timing_inf);
			total_num_chans = 0;
			for (const auto &v : get_vertices(orig_g)) {
				const auto &ver = get_vertex_props(orig_g, v);
				if (ver.type == CHANX || ver.type == CHANY) {
					++total_num_chans;
				}
			}
			//load_real_track_numbers();
			alloc_num_tracks_in_partition();
			result_pid.resize(num_vertices(orig_g));
		}

		//void load_initial_partitions()
		//{
			//initial_partition.resize(num_vertices(orig_g));
			//for (const auto &v : get_vertices(orig_g)) {
				//const auto &ver = get_vertex_props(orig_g, v);
				//if (ver.type == CHANX || ver.type == CHANY) {
					//assert(real_track_num[v] >= 0 && real_track_num[v] < num_tracks);
					//initial_partition[v] = real_track_num[v] % 2;
					//int type = ver.type == CHANX ? 0 : 1;
					//int x, y;
					//std::tie(x, y) = get_node_start(v);

					//++num_tracks_in_partition[type][initial_partition[v]][x][y];
				//}
			//}
		//}

		void alloc_num_tracks_in_partition()
		{
			extern int nx, ny;

			num_tracks_in_partition.resize(2);
			for (int type = 0; type < 2; ++type) {
				num_tracks_in_partition[type].resize(2);
				for (int part = 0; part < 2; ++part) {
					num_tracks_in_partition[type][part].resize(nx+2);
					for (int x = 0; x < nx+2; ++x) {
						num_tracks_in_partition[type][part][x].resize(ny+2);
					}
				}
			}
		}

		void init_num_tracks_in_partition()
		{
			extern int nx, ny;

			for (int type = 0; type < 2; ++type) {
				for (int part = 0; part < 2; ++part) {
					for (int x = 0; x < nx+2; ++x) {
						for (int y = 0; y < ny+2; ++y) {
							num_tracks_in_partition[type][part][x][y] = 0;
						}
					}
				}
			}
		}

		//void load_initial_partitions(int num_partitions)
		//{
		//initial_partition.resize(num_vertices(orig_g));
		//for (const auto &v : get_vertices(orig_g)) {
		//const auto &ver = get_vertex_props(orig_g, v);
		//if (ver.type == CHANX || ver.type == CHANY) {
		//assert(real_track_num[v] >= 0 && real_track_num[v] < num_tracks);
		//initial_partition[v] = real_track_num[v] % num_partitions;
		//int type = ver.type == CHANX ? 0 : 1;
		//int x, y;
		//std::tie(x, y) = get_node_start(v);

		//++num_tracks_in_partition[type][initial_partition[v]][x][y];
		//}
		//}
		//}

		int type_index(int v) const
		{
			const auto &ver = get_vertex_props(orig_g, v);
			assert(ver.type == CHANX || ver.type == CHANY);
			return ver.type == CHANX ? 0 : 1;
		}

		bool is_imbalanced(int v, int to, int &imba) const
		{
			int from = (to+1) % 2;
			int x, y;
			std::tie(x, y) = get_node_start(v);

			assert(from == 0 || from == 1);
			assert(to == 0 || to == 1);

			int average = num_starting_tracks[type_index(v)][x][y] / num_partitions;
			int after = num_tracks_in_partition[type_index(v)][from][x][y] - 1; 
			imba = abs(after-average);

			assert(after >= 0);

			//return after == 0;
			return imba > 3;
		}

		void move(int v, int to)
		{
			int from = (to+1) % 2;
			int x, y;
			std::tie(x, y) = get_node_start(v);

			assert(from == 0 || from == 1);
			assert(to == 0 || to == 1);

			num_tracks_in_partition[type_index(v)][to][x][y]++;
			num_tracks_in_partition[type_index(v)][from][x][y]--;
		}

		void add(int v, int to)
		{
			int x, y;
			std::tie(x, y) = get_node_start(v);

			assert(to == 0 || to == 1);
			num_tracks_in_partition[type_index(v)][to][x][y]++;
		}

		//void recursive_bipartition(int start, int end)
		//{
		//int current_num_partitions = end-start+1;
		//assert(current_num_partitions % 2 == 0);
		//int temp = current_num_partitions / 2;
		//if (temp > 1) {
		//recursive_bipartition(start, start+temp-1);
		//recursive_bipartition(start+temp, end);
		//} else {
		//partition_index_offset = start;
		//assert(current_num_partitions == 2);
		//auto current_g = make_filtered_graph(orig_g, 
		//[&] (unsigned long v) { 
		//const auto &ver = get_vertex_props(orig_g, v);
		//return (ver.type == CHANX || ver.type == CHANY) && (initial_partition[v] == start || initial_partition[v] == end);
		//},
		//[&] (unsigned long e) { 
		//const auto &from = get_vertex_props(orig_g, get_source(orig_g, e));
		//const auto &to = get_vertex_props(orig_g, get_target(orig_g, e));
		//return (from.type == CHANX || from.type == CHANY) && (initial_partition[id(from)] == start || initial_partition[id(from)] == end) &&
		//(to.type == CHANX || to.type == CHANY) && (initial_partition[id(to)] == start || initial_partition[id(to)] == end);
		//});

		//fm<decltype(current_g), rr_graph_partitioner> fm;
		//vector<int> adjusted_initial_partition(num_vertices(current_g), -1);
		//for (const auto &v : get_vertices(current_g)) {
		//adjusted_initial_partition[v] = initial_partition[v]-partition_index_offset;
		//}

		//fm.init(current_g, adjusted_initial_partition, *this);
		//fm.run();
		//}
		//}

		template<typename Graph, typename Valid>
			decltype(auto) make_subgraph(Graph &g, const Valid &vertex_valid)
			{
				return make_filtered_graph(g, subgraph_vertex_predicate_t<Valid>(vertex_valid), subgraph_edge_predicate_t<typename Graph::base, Valid>(g, vertex_valid));
			}

		template<typename Graph>
			void grow_initial_partition_2(const Graph &g, vector<int> &initial_partition)
			{
				struct visitor_t {
					set<int> visited_edges;
					set<int> visited_nodes;
					set<tuple<int, int, t_rr_type, bool>> visited_channels;
					int num_nodes_discovered;

					visitor_t() :
						num_nodes_discovered(0)
					{
					}

					bool tree_edge(const RREdge &e, const Graph &g)
					{
						bool expand = false;
						assert(false);
						//if (visited_nodes.find(get_source(g, e)) != end(visited_nodes)) {
							//bool recorded_to = record(g, get_target(g, e));
							//if (recorded_to) {
								//assert(visited_edges.find(e) == end(visited_edges));
								//visited_edges.insert(e);
								//expand = true;
							//}
						//}
						return expand;
					}

					void examine_edge(const RREdge &e, const Graph &g)
					{
					}

					bool discover_vertex(int v, const Graph &g)
					{
						char buffer[256];
						sprintf_rr_node(v, buffer);
						zlog_level(delta_log, ROUTER_V3, "Current: %s\n", buffer);
						++num_nodes_discovered;
					}

					void examine_vertex(int v, const Graph &g)
					{
					}

					bool record_impl(const Graph &g, int v)
					{
						const auto &ver = get_vertex_props(g, v);
						const auto &start = get_node_start(v);
						const auto &key = make_tuple(start.first, start.second, ver.type, ver.inc_direction);
						bool recorded;
						if (visited_channels.find(key) == visited_channels.end()) {
							visited_channels.insert(key);
							recorded = true;
						} else {
							recorded = false;
						}
						return recorded;
					}

					bool record(const Graph &g, int v)
					{
						bool recorded = record_impl(g, v);
						if (recorded) {
							assert(visited_nodes.find(v) == end(visited_nodes));
							visited_nodes.insert(v);
						}
						return recorded;
					}

					bool record_source(const Graph &g, int v)
					{	
						record_impl(g, v);

						assert(visited_nodes.find(v) == end(visited_nodes));
						visited_nodes.insert(v);

						return true;
					}
				};

				using clock = std::chrono::high_resolution_clock;

				vector<VertexColor> color(num_vertices(g), VertexColor::WHITE);
				int num_nodes = 0;
				set<int> all_visited_nodes;
				vector<vector<int>> visited_nodes;
				for (const auto &v : get_vertices(g)) {
					if (color[v] == VertexColor::WHITE) {
						visitor_t visitor;
						visitor.record_source(g, v);

						auto start = clock::now();
						bfs(g, { static_cast<unsigned long>(v) }, color, visitor);
						auto time = clock::now()-start;

						start = clock::now();
						assert(visitor.visited_nodes.size() == visitor.visited_edges.size()+1);
						assert(visitor.visited_nodes.find(v) != end(visitor.visited_nodes));
						for (const auto &v : visitor.visited_nodes) {
							//initial_partition[v] = current_pid;
							assert(color[v] == VertexColor::BLACK);
							assert(all_visited_nodes.find(v) == all_visited_nodes.end());
							all_visited_nodes.insert(v);
						}
						visited_nodes.emplace_back(begin(visitor.visited_nodes), end(visitor.visited_nodes));
						//for (const auto &v : get_vertices(g)) {
							//if (all_visited_nodes.find(v) == end(all_visited_nodes)) {
								//assert(color[v] == VertexColor::WHITE);
							//} else {
								//assert(color[v] == VertexColor::BLACK);
							//}
						//}
						auto update_time = clock::now()-start;

						//num_nodes_in_partition[current_pid] += visitor.visited_nodes.size();
						//current_pid = (current_pid+1) % 2;

						//printf("ver: %d time = %lld ms update time = %lld ms visited nodes = %lu\n", v, std::chrono::duration_cast<std::chrono::milliseconds>(time).count(), std::chrono::duration_cast<std::chrono::milliseconds>(update_time).count(), visitor.visited_nodes.size());
						//if (visitor.visited_nodes.size() == 1) {
							//printf("
						//}
					}
					++num_nodes;
				}
				std::sort(begin(visited_nodes), end(visited_nodes), [] (const auto &a, const auto &b) -> bool { return a.size()>b.size(); });
				vector<int> num_nodes_in_partition(2, 0);
				int current_pid = 0;
				for (const auto &nodes : visited_nodes) {
					num_nodes_in_partition[current_pid] += nodes.size();
					for (const auto &v : nodes) {
						initial_partition[v] = current_pid;
					}
					current_pid = (current_pid+1) % 2;
				}
				assert(num_nodes_in_partition[0] + num_nodes_in_partition[1] == num_nodes);
				printf("P0 Size = %d P1 Size = %d Imbalance: %d\n", num_nodes_in_partition[0], num_nodes_in_partition[1], num_nodes_in_partition[0]-num_nodes_in_partition[1]);
			}

		template<typename Graph>
			void grow_initial_partition_3(const Graph &g, vector<int> &initial_partition)
			{
				struct visitor_t {
					set<RREdge> visited_edges;
					set<int> visited_nodes;
					set<tuple<int, int, t_rr_type, bool>> visited_channels;
					int num_nodes_discovered;
					const vector<int> &pid;
					bool reachable_partition[2];

					visitor_t(const vector<int> &pid) :
						num_nodes_discovered(0) , pid(pid)
					{
						std::fill(begin(reachable_partition), end(reachable_partition), false);
						assert(!reachable_partition[0]);
						assert(!reachable_partition[1]);
					}

					bool tree_edge(const RREdge &e, const Graph &g)
					{
						bool expand = false;
						if (visited_nodes.find(get_source(g, e)) != end(visited_nodes)) {
							bool recorded_to = record(g, get_target(g, e));
							if (recorded_to) {
								assert(visited_edges.find(e) == end(visited_edges));
								visited_edges.insert(e);
								expand = true;
							}
						}
						return expand;
					}

					void examine_edge(const RREdge &e, const Graph &g)
					{
						int to = get_target(g, e);
						int to_pid = pid[to];
						if (to_pid != -1) {
							assert(to_pid == 0 || to_pid == 1);
							if (!reachable_partition[to_pid]) {
								reachable_partition[to_pid] = true;
							}
						}
					}

					void discover_vertex(int v, const Graph &g)
					{
						char buffer[256];
						sprintf_rr_node(v, buffer);
						zlog_level(delta_log, ROUTER_V3, "Current: %s\n", buffer);
						++num_nodes_discovered;
					}

					void examine_vertex(int v, const Graph &g)
					{
					}

					bool record_impl(const Graph &g, int v)
					{
						const auto &ver = get_vertex_props(g, v);
						const auto &start = get_node_start(v);
						const auto &key = make_tuple(start.first, start.second, ver.type, ver.inc_direction);
						bool recorded;
						if (visited_channels.find(key) == visited_channels.end()) {
							visited_channels.insert(key);
							recorded = true;
						} else {
							recorded = false;
						}
						return recorded;
					}

					bool record(const Graph &g, int v)
					{
						bool recorded = record_impl(g, v);
						if (recorded) {
							assert(visited_nodes.find(v) == end(visited_nodes));
							visited_nodes.insert(v);
						}
						return recorded;
					}

					bool record_source(const Graph &g, int v)
					{	
						assert(record_impl(g, v));

						assert(visited_nodes.find(v) == end(visited_nodes));
						visited_nodes.insert(v);

						return true;
					}
				};

				using clock = std::chrono::high_resolution_clock;

				vector<VertexColor> color(num_vertices(g), VertexColor::WHITE);
				int num_nodes = 0;
				set<int> all_visited_nodes;
				vector<int> num_nodes_in_partition(2, 0);
				//vector<int> sorted;
				//topological_sort(g, sorted);
				for (const auto &v : get_vertices(g)) {
					if (color[v] == VertexColor::WHITE) {
						visitor_t visitor(initial_partition);
						visitor.record_source(g, v);

						auto start = clock::now();
						bfs(g, { static_cast<unsigned long>(v) }, color, visitor);
						auto time = clock::now()-start;

						assert(visitor.visited_nodes.size() == visitor.visited_edges.size()+1);
						assert(visitor.visited_nodes.find(v) != end(visitor.visited_nodes));
						int target_partition;
						if (num_nodes_in_partition[0] == 0) {
							target_partition = 0;
						} else if (num_nodes_in_partition[1] == 0) {
							target_partition = 1;
						} else {
							if ((visitor.reachable_partition[0] && visitor.reachable_partition[1]) ||
									(!visitor.reachable_partition[0] && !visitor.reachable_partition[1])) {
								if (num_nodes_in_partition[0] < num_nodes_in_partition[1]) {
									target_partition = 0;
								} else {
									target_partition = 1;
								}
							} else {
								if (visitor.reachable_partition[0]) {
									target_partition = 0;
								} else {
									assert(visitor.reachable_partition[1]); 
									target_partition = 1;
								} 
							}
						}
						for (const auto &v : visitor.visited_nodes) {
							assert(color[v] == VertexColor::BLACK);
							assert(all_visited_nodes.find(v) == all_visited_nodes.end());
							initial_partition[v] = target_partition;
							all_visited_nodes.insert(v);
						}
						num_nodes_in_partition[target_partition] += visitor.visited_nodes.size();
					}
					++num_nodes;
				}
				assert(num_nodes_in_partition[0] + num_nodes_in_partition[1] == num_nodes);
				printf("P0 Size = %d P1 Size = %d Imbalance: %d\n", num_nodes_in_partition[0], num_nodes_in_partition[1], num_nodes_in_partition[0]-num_nodes_in_partition[1]);
			}

		template<typename Graph>
			void grow_initial_partition(const Graph &g, vector<int> &initial_partition)
			{
				struct visitor_t {
					set<RREdge> visited_edges;
					set<int> visited_nodes;
					set<tuple<int, int, t_rr_type, bool>> visited_channels;
					int num_nodes_discovered;
					const vector<int> &pid;
					bool reachable_partition[2];

					visitor_t(const vector<int> &pid) :
						num_nodes_discovered(0) , pid(pid)
					{
						std::fill(begin(reachable_partition), end(reachable_partition), false);
						assert(!reachable_partition[0]);
						assert(!reachable_partition[1]);
					}

					bool tree_edge(const RREdge &e, const Graph &g)
					{
						bool expand = false;
						if (visited_nodes.find(get_source(g, e)) != end(visited_nodes)) {
							bool recorded_to = record(g, get_target(g, e));
							if (recorded_to) {
								assert(visited_edges.find(e) == end(visited_edges));
								visited_edges.insert(e);
								expand = true;
							}
						}
						return expand;
					}

					void examine_edge(const RREdge &e, const Graph &g)
					{
						int to = get_target(g, e);
						int to_pid = pid[to];
						if (to_pid != -1) {
							assert(to_pid == 0 || to_pid == 1);
							if (!reachable_partition[to_pid]) {
								reachable_partition[to_pid] = true;
							}
						}
					}

					void discover_vertex(int v, const Graph &g)
					{
						char buffer[256];
						sprintf_rr_node(v, buffer);
						zlog_level(delta_log, ROUTER_V3, "Current: %s\n", buffer);
						++num_nodes_discovered;
					}

					void examine_vertex(int v, const Graph &g)
					{
					}

					bool record_impl(const Graph &g, int v)
					{
						const auto &ver = get_vertex_props(g, v);
						const auto &start = get_node_start(v);
						const auto &key = make_tuple(start.first, start.second, ver.type, ver.inc_direction);
						bool recorded;
						if (visited_channels.find(key) == visited_channels.end()) {
							visited_channels.insert(key);
							recorded = true;
						} else {
							recorded = false;
						}
						return recorded;
					}

					bool record(const Graph &g, int v)
					{
						bool recorded = record_impl(g, v);
						if (recorded) {
							assert(visited_nodes.find(v) == end(visited_nodes));
							visited_nodes.insert(v);
						}
						return recorded;
					}

					bool record_source(const Graph &g, int v)
					{	
						assert(record_impl(g, v));

						assert(visited_nodes.find(v) == end(visited_nodes));
						visited_nodes.insert(v);

						return true;
					}
				};

				using clock = std::chrono::high_resolution_clock;

				vector<VertexColor> color(num_vertices(g), VertexColor::WHITE);
				int num_nodes = 0;
				set<int> all_visited_nodes;
				vector<int> num_nodes_in_partition(2, 0);
				//vector<int> sorted;
				//topological_sort(g, sorted);
				for (const auto &v : get_vertices(g)) {
					if (color[v] == VertexColor::WHITE) {
						visitor_t visitor(initial_partition);
						visitor.record_source(g, v);

						auto start = clock::now();
						bfs(g, { static_cast<unsigned long>(v) }, color, visitor);
						auto time = clock::now()-start;

						assert(visitor.visited_nodes.size() == visitor.visited_edges.size()+1);
						assert(visitor.visited_nodes.find(v) != end(visitor.visited_nodes));
						int target_partition;
						if (num_nodes_in_partition[0] == 0) {
							target_partition = 0;
						} else if (num_nodes_in_partition[1] == 0) {
							target_partition = 1;
						} else {
							if ((visitor.reachable_partition[0] && visitor.reachable_partition[1]) ||
									(!visitor.reachable_partition[0] && !visitor.reachable_partition[1])) {
								if (num_nodes_in_partition[0] < num_nodes_in_partition[1]) {
									target_partition = 0;
								} else {
									target_partition = 1;
								}
							} else {
								if (visitor.reachable_partition[0]) {
									target_partition = 0;
								} else {
									assert(visitor.reachable_partition[1]); 
									target_partition = 1;
								} 
							}
						}
						for (const auto &v : visitor.visited_nodes) {
							assert(color[v] == VertexColor::BLACK);
							assert(all_visited_nodes.find(v) == all_visited_nodes.end());
							initial_partition[v] = target_partition;
							all_visited_nodes.insert(v);
						}
						num_nodes_in_partition[target_partition] += visitor.visited_nodes.size();
					}
					++num_nodes;
				}
				assert(num_nodes_in_partition[0] + num_nodes_in_partition[1] == num_nodes);
				printf("P0 Size = %d P1 Size = %d Imbalance: %d\n", num_nodes_in_partition[0], num_nodes_in_partition[1], num_nodes_in_partition[0]-num_nodes_in_partition[1]);
			}

		template<typename Graph>
			void recursive_bipartition(const Graph &g, const Graph &undirected_g, int level, int node)
			{
				init_num_tracks_in_partition();

				vector<int> initial_partition(num_vertices(g), -1);

				/* get initial partition by BFS */
				grow_initial_partition(g, initial_partition);

				/* refine partition */
				//fm<Graph, rr_graph_partitioner> fm;
				//fm.init(undirected_g, initial_partition, *this);
				//fm.run();

				//const vector<int> &pid = fm.get_pid();
				const vector<int> &pid = initial_partition;

                for (const auto &v : get_vertices(g)) {
                    assert(result_pid_by_level[level][v] == -1);
                    result_pid_by_level[level][v] = 2*node+pid[v];
                }

				if (level > 0) {
					//auto valid_vertex_0 = [] (unsigned long v) -> bool { return true; };
					//auto valid_vertex_1 = [] (unsigned long v) -> bool { return true; };

					for (int i = 0; i < 2; ++i) {
						recursive_bipartition(
								make_subgraph(g.base_graph(), partition_vertex_predicate_t<typename Graph::base>(g.base_graph(), pid, i)),
								make_subgraph(undirected_g.base_graph(), partition_vertex_predicate_t<typename Graph::base>(undirected_g.base_graph(), pid, i)),
								level-1, 2*node+i);
					}
				} else {
					for (const auto &v : get_vertices(g)) {
						assert(result_pid[v] == -1);
						assert(pid[v] == 0 || pid[v] == 1);
						result_pid[v] = 2*node+pid[v];
					}
					//printf("last level index: %d %d\n", 2*node, 2*node+1);
				}
			}

		void partition_without_ipin(int _num_partitions, vector<RRGraph *> &graphs)
		{
			assert(std::pow(2, std::log2(_num_partitions)) == _num_partitions);

			num_partitions = _num_partitions;
            num_levels = std::log2(num_partitions)+1;

            result_pid_by_level.resize(num_levels, vector<int>(num_vertices(orig_g), -1));

			for (const auto &v : get_vertices(orig_g)) {
				const auto &props = get_vertex_props(orig_g, v);
				if (props.type == CHANX || props.type == CHANY) {
					result_pid_by_level[num_levels-1][v] = 0;
				} else {
					assert(result_pid_by_level[num_levels-1][v] == -1);
				}
			}

			std::fill(begin(result_pid), end(result_pid), -1);

			if (num_partitions == 1) {
				assert(num_levels == 1);
				for (const auto &v : get_vertices(orig_g)) {
					const auto &props = get_vertex_props(orig_g, v);
					if (props.type == CHANX || props.type == CHANY) {
						result_pid[v] = 0;
					} else {
						assert(result_pid[v] == -1);
					}
				}
			} else {
				vector<int> pid(num_vertices(orig_g), 0);

				recursive_bipartition(
						make_subgraph(orig_g.base_graph(), partition_vertex_predicate_t<typename decltype(orig_g)::base>(orig_g.base_graph(), pid, 0)),
						make_subgraph(undirected_orig_g.base_graph(), partition_vertex_predicate_t<typename decltype(undirected_orig_g)::base>(undirected_orig_g.base_graph(), pid, 0)), num_levels-2, 0);
			}

			for (int level = 0; level < result_pid_by_level.size(); ++level) {
				int max_pid = 0;
				for (int i = 0; i < result_pid_by_level[level].size(); ++i) {
					max_pid = std::max(max_pid, result_pid_by_level[level][i]);
				}
				assert(max_pid == pow(2, result_pid_by_level.size()-level-1)-1);
				vector<int> num_nodes(max_pid+1, 0);
				for (const auto &v : get_vertices(orig_g)) {
					const auto &props = get_vertex_props(orig_g, v);
					if (props.type == CHANX || props.type == CHANY) {
						++num_nodes[result_pid_by_level[level][v]];
					} else {
						assert(result_pid_by_level[level][v] == -1);
					}
				}
				for (int i = 0; i < num_nodes.size(); ++i) {
					printf("level %d p%d num nodes: %d\n", level, i, num_nodes[i]);
				}
			}

			//int num_nodes_partitioned = 0;
			//vector<int> num_nodes_in_partition(num_partitions, 0);
			//for (const auto &p : result_pid) {
				//if (p != -1) {
					//assert(p >= 0 && p < num_partitions);
					//++num_nodes_in_partition[p];
					//++num_nodes_partitioned;
				//}
			//}
			//int num_ipins = 0;
			//for (const auto &v : get_vertices(orig_g)) {
				//const auto &props = get_vertex_props(orig_g, v);
				//if (props.type == IPIN) {
					//++num_ipins;
				//}
			//}
			//assert(num_nodes_partitioned == total_num_chans + num_ipins_partitioned);

			//for (int i = 0; i < num_partitions; ++i) {
				//printf("Num nodes in partition %d = %d\n", i, num_nodes_in_partition[i]);
				//char filename[256];
				//sprintf(filename, "/Volumes/DATA/graph_part_%d.txt", i);
				//FILE *file = fopen(filename, "w");
				//for (const auto &v : get_vertices(orig_g)) {
					//if (result_pid[v] != i) {
						//continue;
					//}
					//const auto &ver = get_vertex_props(orig_g, v);
					//if ((ver.type == CHANX || ver.type == CHANY)) {
						//const auto &start = get_node_start(v);
						//fprintf(file, "%d %d\n", start.first, start.second);
					//}
				//}
				//fclose(file);

				//RRGraph &new_g = *new RRGraph;

				//add_vertex(new_g, num_vertices(orig_g));
				//for (const auto &v : get_vertices(orig_g)) {
					//get_vertex_props(new_g, v) = get_vertex_props(orig_g, v);
				//}

				//for (const auto &e : get_edges(orig_g)) {
					//int from = get_source(orig_g, e);
					//int to = get_target(orig_g, e);

					//const auto &from_ver = get_vertex_props(orig_g, from);
					//const auto &to_ver = get_vertex_props(orig_g, to);

					//if ((from_ver.type == CHANX || from_ver.type == CHANY)
							//&& (to_ver.type == CHANX || to_ver.type == CHANY)) {
						//assert(result_pid[from] != -1 && result_pid[to] != -1);
						//if (result_pid[from] == i && result_pid[from] == result_pid[to]) {
							//const RREdge &new_e = add_edge(new_g, from, to);
							//get_edge_props(new_g, new_e) = get_edge_props(orig_g, e);
						//}
					//} else {
						//assert(result_pid[from] == -1 || result_pid[to] == -1);
						//const RREdge &new_e = add_edge(new_g, from, to);
						//get_edge_props(new_g, new_e) = get_edge_props(orig_g, e);
					//}
				//}

				//graphs.push_back(&new_g);
			//}
		}

		void partition(int _num_partitions, vector<RRGraph *> &graphs)
		{
			assert(std::pow(2, std::log2(_num_partitions)) == _num_partitions);

			num_partitions = _num_partitions;
            num_levels = std::log2(num_partitions)+1;

            result_pid_by_level.resize(num_levels, vector<int>(num_vertices(orig_g), -1));

			for (const auto &v : get_vertices(orig_g)) {
				const auto &props = get_vertex_props(orig_g, v);
				if (props.type == CHANX || props.type == CHANY) {
					result_pid_by_level[num_levels-1][v] = 0;
				} else {
					assert(result_pid_by_level[num_levels-1][v] == -1);
				}
			}

			std::fill(begin(result_pid), end(result_pid), -1);

			if (num_partitions == 1) {
				assert(num_levels == 1);
				for (const auto &v : get_vertices(orig_g)) {
					const auto &props = get_vertex_props(orig_g, v);
					if (props.type == CHANX || props.type == CHANY || props.type == IPIN) {
						result_pid[v] = 0;
					} else {
						assert(result_pid[v] == -1);
					}
				}
			} else {
				vector<int> pid(num_vertices(orig_g), 0);

				recursive_bipartition(
						make_subgraph(orig_g.base_graph(), partition_vertex_predicate_t<typename decltype(orig_g)::base>(orig_g.base_graph(), pid, 0)),
						make_subgraph(undirected_orig_g.base_graph(), partition_vertex_predicate_t<typename decltype(undirected_orig_g)::base>(undirected_orig_g.base_graph(), pid, 0)), num_levels-2, 0);
			}

			int num_ipins_partitioned = 0;
			for (int level = 0; level < result_pid_by_level.size(); ++level) {
				int cur_num_partitions = pow(2, result_pid_by_level.size()-level-1);

				vector<vector<int>> sink_partition_connect_count(num_vertices(orig_g), vector<int>(cur_num_partitions, 0));

				for (const auto &v : get_vertices(orig_g)) {
					const auto &rr_node_p = get_vertex_props(orig_g, v);

					if (rr_node_p.type == IPIN) {
						assert(num_out_edges(orig_g, v) == 1);

						const auto &e = *begin(get_out_edges(orig_g, v));

						const auto &sink_rr_node = get_target(orig_g, e);

						assert(get_vertex_props(orig_g, sink_rr_node).type == SINK);

						if (!ipin_in_nodes[v].empty()) {
							int min_partition = -1;
							int min_count = std::numeric_limits<int>::max();
							vector<int> num_in_nodes_per_partition(cur_num_partitions, 0);
							for (const auto &in_node : ipin_in_nodes[v]) {
								const auto &in_node_p = get_vertex_props(orig_g, in_node);

								assert(in_node_p.type == CHANX || in_node_p.type == CHANY);

								int part = result_pid_by_level[level][in_node];

								assert(part >= 0 && part < cur_num_partitions);

								++num_in_nodes_per_partition[part];

								if (sink_partition_connect_count[sink_rr_node][part] < min_count) {
									min_partition = part;
									min_count = sink_partition_connect_count[sink_rr_node][part];
								}
							}
							assert(min_partition != -1);

							sink_partition_connect_count[sink_rr_node][min_partition] += num_in_nodes_per_partition[min_partition];
							assert(result_pid_by_level[level][v] == -1);
							result_pid_by_level[level][v] = min_partition;

							++num_ipins_partitioned;
						}
					}
				}
			}

			for (int level = 0; level < result_pid_by_level.size(); ++level) {
				int max_pid = 0;
				for (int i = 0; i < result_pid_by_level[level].size(); ++i) {
					max_pid = std::max(max_pid, result_pid_by_level[level][i]);
				}
				assert(max_pid == pow(2, result_pid_by_level.size()-level-1)-1);
				vector<int> num_nodes(max_pid+1, 0);
				for (const auto &v : get_vertices(orig_g)) {
					const auto &props = get_vertex_props(orig_g, v);
					if (props.type == CHANX || props.type == CHANY || props.type == IPIN) {
						int part = result_pid_by_level[level][v];
						if (part == -1) {
							assert(props.type == IPIN);
						} else {
							++num_nodes[part];
						}
					} else {
						assert(result_pid_by_level[level][v] == -1);
					}
				}
				for (int i = 0; i < num_nodes.size(); ++i) {
					printf("level %d p%d num nodes: %d\n", level, i, num_nodes[i]);
				}
			}

			//int num_nodes_partitioned = 0;
			//vector<int> num_nodes_in_partition(num_partitions, 0);
			//for (const auto &p : result_pid) {
				//if (p != -1) {
					//assert(p >= 0 && p < num_partitions);
					//++num_nodes_in_partition[p];
					//++num_nodes_partitioned;
				//}
			//}
			//int num_ipins = 0;
			//for (const auto &v : get_vertices(orig_g)) {
				//const auto &props = get_vertex_props(orig_g, v);
				//if (props.type == IPIN) {
					//++num_ipins;
				//}
			//}
			//assert(num_nodes_partitioned == total_num_chans + num_ipins_partitioned);

			//for (int i = 0; i < num_partitions; ++i) {
				//printf("Num nodes in partition %d = %d\n", i, num_nodes_in_partition[i]);
				//char filename[256];
				//sprintf(filename, "/Volumes/DATA/graph_part_%d.txt", i);
				//FILE *file = fopen(filename, "w");
				//for (const auto &v : get_vertices(orig_g)) {
					//if (result_pid[v] != i) {
						//continue;
					//}
					//const auto &ver = get_vertex_props(orig_g, v);
					//if ((ver.type == CHANX || ver.type == CHANY)) {
						//const auto &start = get_node_start(v);
						//fprintf(file, "%d %d\n", start.first, start.second);
					//}
				//}
				//fclose(file);

				//RRGraph &new_g = *new RRGraph;

				//add_vertex(new_g, num_vertices(orig_g));
				//for (const auto &v : get_vertices(orig_g)) {
					//get_vertex_props(new_g, v) = get_vertex_props(orig_g, v);
				//}

				//for (const auto &e : get_edges(orig_g)) {
					//int from = get_source(orig_g, e);
					//int to = get_target(orig_g, e);

					//const auto &from_ver = get_vertex_props(orig_g, from);
					//const auto &to_ver = get_vertex_props(orig_g, to);

					//if ((from_ver.type == CHANX || from_ver.type == CHANY)
							//&& (to_ver.type == CHANX || to_ver.type == CHANY)) {
						//assert(result_pid[from] != -1 && result_pid[to] != -1);
						//if (result_pid[from] == i && result_pid[from] == result_pid[to]) {
							//const RREdge &new_e = add_edge(new_g, from, to);
							//get_edge_props(new_g, new_e) = get_edge_props(orig_g, e);
						//}
					//} else {
						//assert(result_pid[from] == -1 || result_pid[to] == -1);
						//const RREdge &new_e = add_edge(new_g, from, to);
						//get_edge_props(new_g, new_e) = get_edge_props(orig_g, e);
					//}
				//}

				//graphs.push_back(&new_g);
			//}
		}
};

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
			
			e_p.buffered = switch_inf[si].buffered; 
			e_p.switch_delay = switch_inf[si].Tdel; 
			e_p.R = switch_inf[si].R; 

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
			
			e_p.buffered = switch_inf[si].buffered; 
			e_p.switch_delay = switch_inf[si].Tdel; 
			e_p.R = switch_inf[si].R; 
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
		sprintf(filename, "%s%s", LOG_PATH_PREFIX, msg->path);

		file = fopen(filename, "w");
		if (!file) {
			perror(nullptr);
			assert(false);
		}
	}
	fprintf(file, "%s", msg->buf);
	fflush(file);
}

static int ss_log_output(zlog_msg_t *msg)
{
	int tid;
	int iter;

	log_impl(msg, ss_log_file);

	return 0;
}

static int delta_log_output(zlog_msg_t *msg)
{
	int tid;
	int iter;

	if (sscanf(msg->path, "iter_%d_tid_%d.log", &iter, &tid) != 2) {
		return 0;
	}

	concurrent_log_impl(msg, delta_log_files, iter, tid);

	return 0;
}

static int missing_edge_log_output(zlog_msg_t *msg)
{
	int tid;
	int iter;

	if (sscanf(msg->path, "missing_edge_iter_%d_tid_%d.log", &iter, &tid) != 2) {
		return 0;
	}

	concurrent_log_impl(msg, missing_edge_log_files, iter, tid);

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
	MPI_Datatype net_dts[] = {
		MPI_INT,
		MPI_INT, MPI_INT, MPI_INT,
		MPI_INT, MPI_INT, MPI_INT, MPI_INT };

	MPI_Datatype net_tmp_dt;
	assert(MPI_Type_create_struct(8, net_bl, net_disp, net_dts, &net_tmp_dt) == MPI_SUCCESS);
	assert(MPI_Type_commit(&net_tmp_dt) == MPI_SUCCESS);
	assert(MPI_Type_create_resized(net_tmp_dt, 0, sizeof(net_t), &net_dt) == MPI_SUCCESS);
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
	MPI_Datatype net_sink_dts[] = {
		MPI_INT, MPI_INT, MPI_INT,
		MPI_INT, MPI_INT, MPI_INT, MPI_INT };

	MPI_Datatype net_sink_tmp_dt;
	assert(MPI_Type_create_struct(7, net_sink_bl, net_sink_disp, net_sink_dts, &net_sink_tmp_dt) == MPI_SUCCESS);
	assert(MPI_Type_commit(&net_sink_tmp_dt) == MPI_SUCCESS);
	assert(MPI_Type_create_resized(net_sink_tmp_dt, 0, sizeof(sink_t), &net_sink_dt) == MPI_SUCCESS);
	assert(MPI_Type_commit(&net_sink_dt) == MPI_SUCCESS);
}

void recalc_sum(congestion_t *in, congestion_t *inout, int *len, MPI_Datatype *dptr)
{
	for (int i = 0; i < *len; ++i) {
		inout[i].recalc_occ = in[i].recalc_occ + inout[i].recalc_occ;
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

void recv_route_tree(net_t *net, const RRGraph &g, vector<vector<sink_t *>> &routed_sinks, route_state_t *states, congestion_t *congestion, float pres_fac, vector<route_tree_t> &route_trees, t_net_timing *net_timing, int from_procid, MPI_Comm comm)
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

		route_tree_set_node_properties(route_trees[net->local_id], rt_node, rr_node_p.type != IPIN && rr_node_p.type != SINK, rr_edge_to_parent, recv[i].upstream_R, recv[i].delay);

		if (rr_node_p.type == SOURCE) {
			route_tree_add_root(route_trees[net->local_id], rr_node);
			source_rr_node = rr_node;
			++num_sources;
		} else if (rr_node_p.type == SINK) {
			++num_sinks;

			auto iter = find_if(begin(routed_sinks[net->local_id]), end(routed_sinks[net->local_id]), [&rr_node] (const auto &sink) -> bool { return sink->rr_node == rr_node; });
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
	RouteTreeNode source_rt_node = route_tree_get_rt_node(route_trees[net->local_id], source_rr_node);
	assert(source_rt_node != RouteTree::null_vertex());
	int delta = num_out_edges(route_trees[net->local_id].graph, source_rt_node)-1;
	update_one_cost_internal(source_rr_node, g, congestion, delta, pres_fac);

	delete [] recv;
}

void send_route_tree(net_t *net, const RRGraph &g, const vector<vector<sink_t *>> &routed_sinks, const vector<route_tree_t> &route_trees, int to_procid, MPI_Comm comm)
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
		if (valid(rt_node_p.rr_edge_to_parent)) {
			send[i].parent_rr_node = get_source(g, rt_node_p.rr_edge_to_parent);
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
void sync(congestion_t *congestion, int this_pid, int num_procs, MPI_Comm comm);

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

void sync_nets(vector<net_t> &nets, vector<net_t> &global_nets, int procid, MPI_Comm comm)
{
	if (procid == 0) {
		int num_global_nets = global_nets.size();
		MPI_Bcast(&num_global_nets, 1, MPI_INT, 0, comm);

		printf("sent num_global_nets: %d\n", num_global_nets);
		
		int num_nets = nets.size();
		MPI_Bcast(&num_nets, 1, MPI_INT, 0, comm);

		printf("sent num_nets: %d\n", num_nets);

		MPI_Bcast(nets.data(), nets.size(), net_dt, 0, comm);

		printf("sent nets\n");

		for (int i = 0; i < nets.size(); ++i) {
			int num_sinks = nets[i].sinks.size();

			zlog_level(delta_log, ROUTER_V3, "Sending net vpr id %d source rr %d x %d y %d bb %d %d %d %d\n", nets[i].vpr_id, nets[i].source.rr_node,nets[i].source.x, nets[i].source.y, nets[i].bounding_box.xmin, nets[i].bounding_box.xmax, nets[i].bounding_box.ymin, nets[i].bounding_box.ymax);
			MPI_Bcast(&num_sinks, 1, MPI_INT, 0, comm);
			MPI_Bcast(nets[i].sinks.data(), num_sinks, net_sink_dt, 0, comm);

			for (int j = 0; j < num_sinks; ++j) {
				zlog_level(delta_log, ROUTER_V3, "\tSending net sink rr %d x %d y %d bb %d %d %d %d\n", nets[i].sinks[j].rr_node, nets[i].sinks[j].x, nets[i].sinks[j].y, nets[i].sinks[j].current_bounding_box.xmin, nets[i].sinks[j].current_bounding_box.xmax, nets[i].sinks[j].current_bounding_box.ymin, nets[i].sinks[j].current_bounding_box.ymax);
			}
		}
	} else {
		int num_global_nets;
		MPI_Bcast(&num_global_nets, 1, MPI_INT, 0, comm);
		global_nets.resize(num_global_nets);

		printf("recv num_global_nets: %d\n", num_global_nets);

		int num_nets;
		MPI_Bcast(&num_nets, 1, MPI_INT, 0, comm);
		nets.resize(num_nets);

		printf("recv num_nets: %d\n", num_nets);

		MPI_Bcast(nets.data(), nets.size(), net_dt, 0, comm);

		printf("recv nets\n");

		int local_id = 0;

		for (int i = 0; i < num_nets; ++i) {
			zlog_level(delta_log, ROUTER_V3, "Recving net vpr id %d source rr %d x %d y %d bb %d %d %d %d\n", nets[i].vpr_id, nets[i].source.rr_node,nets[i].source.x, nets[i].source.y, nets[i].bounding_box.xmin, nets[i].bounding_box.xmax, nets[i].bounding_box.ymin, nets[i].bounding_box.ymax);

			nets[i].local_id = local_id;

			int num_sinks;
			MPI_Bcast(&num_sinks, 1, MPI_INT, 0, comm);
			nets[i].sinks.resize(num_sinks);
			MPI_Bcast(nets[i].sinks.data(), num_sinks, net_sink_dt, 0, comm);

			for (int j = 0; j < num_sinks; ++j) {
				zlog_level(delta_log, ROUTER_V3, "\tRecving net sink rr %d x %d y %d bb %d %d %d %d\n", nets[i].sinks[j].rr_node, nets[i].sinks[j].x, nets[i].sinks[j].y, nets[i].sinks[j].current_bounding_box.xmin, nets[i].sinks[j].current_bounding_box.xmax, nets[i].sinks[j].current_bounding_box.ymin, nets[i].sinks[j].current_bounding_box.ymax);

				nets[i].sinks[j].id = j;
				nets[i].sinks[j].criticality_fac = std::numeric_limits<float>::max();
			}

			++local_id;
		}
	}
}

void sync_recalc_occ(congestion_t *congestion, int num_vertices, int procid, int num_procs, MPI_Comm comm)
{
	if (procid == 0) {
		int num_recvs = 0;
		//int rr_node_pid = partitioner.result_pid_by_level[current_level][i];
		//int from_pid = rr_node_pid == -1 ? 0 : rr_node_pid;

		//if (from_pid != procid) {
		int *all_recalc_occ = new int[num_vertices];

		for (int from_pid = 1; from_pid < num_procs; ++from_pid) {
			int count;

			//do { 
			//MPI_Status status;

			//MPI_Recv(all_recalc_occ, num_vertices, MPI_INT, from_pid, from_pid, cur_comm, &status);
			//MPI_Get_count(&status, MPI_INT, &count);
			//} while (count <= 0);

			MPI_Status status;
			assert(MPI_Recv(all_recalc_occ, num_vertices, MPI_INT, from_pid, from_pid, comm, &status) == MPI_SUCCESS);
			MPI_Get_count(&status, MPI_INT, &count);

			assert(count == num_vertices);

			for (int i = 0; i < num_vertices; ++i) {
				//zlog_level(delta_log, ROUTER_V3, "%d Recvd %d recalc_occ %d\n", from_pid, i, recalc_occ);
				congestion[i].recalc_occ += all_recalc_occ[i];
			}
		}

		delete [] all_recalc_occ;
	} else {
		int num_sends = 0;

		int *all_recalc_occ = new int[num_vertices];

		for (int i = 0; i < num_vertices; ++i) {
			//int pid = partitioner.result_pid_by_level[current_level][i];
			//if (pid == procid) {
			zlog_level(delta_log, ROUTER_V3, "Sending %d recalc_occ %d\n", i, congestion[i].recalc_occ);
			all_recalc_occ[i] = congestion[i].recalc_occ;
			++num_sends;
			//}
		}

		assert(MPI_Send(all_recalc_occ, num_vertices, MPI_INT, 0, procid, comm) == MPI_SUCCESS);

		delete [] all_recalc_occ;

		printf("[%d] Num sends: %d\n", procid, num_sends);
	}
}

bool mpi_spatial_route_flat(t_router_opts *opts, struct s_det_routing_arch det_routing_arch, t_direct_inf *directs, int num_directs, t_segment_inf *segment_inf, t_timing_inf timing_inf)
{
    using std::chrono::duration_cast;
    using std::chrono::nanoseconds;
    using clock = std::chrono::high_resolution_clock;

	int initial_num_procs, num_procs, initial_procid, procid;
	MPI_Comm cur_comm = MPI_COMM_WORLD;

    MPI_Comm_size(cur_comm, &initial_num_procs);
    MPI_Comm_rank(cur_comm, &initial_procid);

	num_procs = initial_num_procs;
	procid = initial_procid;

	init_congestion_mpi_datatype();
	init_datatypes();

    init_logging();
    zlog_set_record("custom_output", delta_log_output);
    zlog_set_record("missing_edge", missing_edge_log_output);
    zlog_set_record("ss", ss_log_output);
    delta_log_files.resize(opts->max_router_iterations);
    for (int i = 0; i < opts->max_router_iterations; ++i) {
        delta_log_files[i].resize(num_procs, nullptr);
    }
    missing_edge_log_files.resize(opts->max_router_iterations);
    for (int i = 0; i < opts->max_router_iterations; ++i) {
        missing_edge_log_files[i].resize(num_procs, nullptr);
    }

    //fclose(fopen(LOG_PATH_PREFIX"iter__tid_.log", "w"));

    test_fast_graph();
    test_topo();
    test_fm();
    test_filter_graph();
    test_partition_graph();
    test_connected_components();

    vector<RRGraph *> graphs;
    rr_graph_partitioner partitioner;
    partitioner.init(opts, &det_routing_arch, directs, num_directs, segment_inf, &timing_inf);
    partitioner.partition_without_ipin(1, graphs);

    //for (const auto &g : graphs) {
        //routability(*g);
    //}
    //
    //RRGraph combined_g;
    //add_vertex(combined_g, num_vertices(partitioner.orig_g));

    //for (const auto &g : graphs) {
        //for (const auto &e : get_edges(*g)) {
            //int from = get_source(*g, e);
            //int to = get_target(*g, e);
            //const auto &from_ver = get_vertex_props(*g, from);
            //const auto &to_ver = get_vertex_props(*g, to);

            //if (is_channel(from_ver) && is_channel(to_ver)) {
                //assert(!has_edge(combined_g, from, to));
                //add_edge(combined_g, from, to);
            //}
        //}
    //}

    //for (const auto &e : get_edges(partitioner.orig_g)) {
        //int from = get_source(partitioner.orig_g, e);
        //int to = get_target(partitioner.orig_g, e);
        //const auto &from_ver = get_vertex_props(partitioner.orig_g, from);
        //const auto &to_ver = get_vertex_props(partitioner.orig_g, to);
        //if (!is_channel(from_ver) || !is_channel(to_ver)) {
            //assert(!has_edge(combined_g, from, to));
            //add_edge(combined_g, from, to);
        //}
    //}

    //printf("Combined/Orig graph has %d/%d (%g) edges.\n", num_edges(combined_g), num_edges(partitioner.orig_g), 100.0*num_edges(combined_g)/num_edges(partitioner.orig_g));

    //goto lol;

    //extern int num_types;
    //extern struct s_type_descriptor *type_descriptors;
    //extern int nx, ny;
    //extern struct s_grid_tile **grid;

    //free_rr_graph();

    //int warnings;

    //build_rr_graph(GRAPH_UNIDIR, num_types, type_descriptors, nx, ny, grid,
            //opts->fixed_channel_width, NULL, det_routing_arch.switch_block_type,
            //det_routing_arch.Fs, det_routing_arch.num_segment,
            //det_routing_arch.num_switch, segment_inf,
            //det_routing_arch.global_route_switch,
            //det_routing_arch.delayless_switch, timing_inf,
            //det_routing_arch.wire_to_ipin_switch, opts->base_cost_type,
            //directs, num_directs, FALSE,
            //&warnings);

    //RRGraph channel_with_interior_g;
    //init_channel_only_graph(channel_with_interior_g);

    //dump_rr_graph(channel_with_interior_g, "/Volumes/DATA/rr_graph_channel_with_interior.txt");

    //RRGraph orig_g;
    //init_graph(orig_g);

    //dump_rr_graph(orig_g, "/Volumes/DATA/rr_graph.txt");

    //free_rr_graph();
    //for (int i = 0; i < det_routing_arch.num_segment; ++i) {
        //for (int j = 0; j < segment_inf[i].sb_len; ++j) {
            //if (j != 0 && j != segment_inf[i].sb_len-1) {
                //segment_inf[i].sb[j] = FALSE;
            //}
        //}
    //}
    //build_rr_graph(GRAPH_UNIDIR, num_types, type_descriptors, nx, ny, grid,
            //opts->fixed_channel_width, NULL, det_routing_arch.switch_block_type,
            //det_routing_arch.Fs, det_routing_arch.num_segment,
            //det_routing_arch.num_switch, segment_inf,
            //det_routing_arch.global_route_switch,
            //det_routing_arch.delayless_switch, timing_inf,
            //det_routing_arch.wire_to_ipin_switch, opts->base_cost_type,
            //directs, num_directs, FALSE,
            //&warnings);

    //RRGraph channel_without_interior_g;
    //init_channel_only_graph(channel_without_interior_g);

    //dump_rr_graph(channel_without_interior_g, "/Volumes/DATA/rr_graph_channel_without_interior.txt");

    //vector<vector<vector<int>>> all_partition_components(opts->num_threads);
    ////for (int x = 0; x < nx+1; ++x) {
        ////for (int y = 0; y < ny+1; ++y) {
            ////init_partitioned_graph_5(opts->num_threads, channel_with_interior_g, channel_without_interior_g, orig_g, graphs, x, y);
        ////}
    ////}
    ////init_partitioned_graph_4(opts->num_threads, channel_with_interior_g, channel_without_interior_g, orig_g, graphs, 14, 14);
    //init_partitioned_graph_4(opts->num_threads, channel_with_interior_g, channel_without_interior_g, orig_g, graphs, all_partition_components);

    //for (int i = 0; i < graphs.size(); ++i) {
        //char filename[256];

        //sprintf(filename, "/Volumes/DATA/rr_graph_%d.txt", i);
        //dump_rr_graph(*graphs[i], filename);

        //sprintf(filename, "/Volumes/DATA/rr_graph_%d_edges.txt", i);
        //dump_edges(*graphs[i], filename);
    //}
	//
	zlog_put_mdc("iter", "0");

	char buffer[256];
	sprintf(buffer, "%d", procid);
	zlog_put_mdc("tid", buffer);

    vector<net_t> nets;
    vector<net_t> global_nets;
    //init_nets(nets, global_nets, opts->bb_factor, partitioner.sink_in_nodes);	
	if (procid == 0) {
		printf("[%d] initializing nets\n", procid);
		init_nets(nets, global_nets, opts->bb_factor, opts->large_bb);	
	}
	printf("[%d] syncing nets\n", procid);
	sync_nets(nets, global_nets, procid, cur_comm);

	printf("[%d] done initializing nets\n", procid);

    vector<net_t *> nets_ptr(nets.size());
    for (int i = 0; i < nets.size(); ++i) {
        nets_ptr[i] = &nets[i];
    }
    extern s_bb *route_bb;
    sort(nets_ptr.begin(), nets_ptr.end(), [] (const net_t *a, const net_t *b) -> bool {
            return get_bounding_box_area(a->bounding_box) > get_bounding_box_area(b->bounding_box);
            });
    int rank = 0;
    for (auto &net : nets_ptr) {
        net->bb_area_rank = rank++;
    }

	route_state_t *states;
	congestion_t *congestion;
	vector<route_tree_t> route_trees;
	t_net_timing *net_timing;
	init_route_structs(partitioner.orig_g, nets, global_nets, &states, &congestion, route_trees, &net_timing);

    route_parameters_t params;
    params.criticality_exp = opts->criticality_exp;
    params.astar_fac = opts->astar_fac;
    params.max_criticality = opts->max_criticality;
    params.bend_cost = opts->bend_cost;
	params.pres_fac = opts->first_iter_pres_fac; /* Typically 0 -> ignore cong. */

	if (procid != 0) {
		free_circuit();
	}
	free_rr_graph();

    vector<pair<box, net_t *>> nets_to_route;
    
    for (auto &net : nets) {
        box b = bg::make_inverse<box>();

        bg::expand(b, point(net.source.x, net.source.y));
        for (const auto &sink : net.sinks) {
            bg::expand(b, point(sink.x, sink.y));
        }
        bg::subtract_value(b.min_corner(), 1+opts->bb_factor);
        bg::add_value(b.max_corner(), opts->bb_factor);

        nets_to_route.push_back(make_pair(b, &net));
    }
    std::sort(begin(nets_to_route), end(nets_to_route), [] (const pair<box, net_t *> &a, const pair<box, net_t *> &b) -> bool {
            return a.second->sinks.size() > b.second->sinks.size();
            });

    bool routed = false;

    int iter;
    float crit_path_delay;
	int current_level = 0;

	vector<vector<sink_t *>> unroutable_sinks(nets.size());
	//vector<vector<sink_t *>> fixed_sinks(nets.size());
	vector<vector<sink_t *>> routed_sinks(nets.size());
	bool has_unroutable_sinks = false;
	unsigned long prev_num_overused_nodes = std::numeric_limits<unsigned long>::max();

    clock::duration total_greedy_route_time = clock::duration::zero();
    clock::duration total_update_cost_time = clock::duration::zero();
    clock::duration total_analyze_timing_time = clock::duration::zero();
    clock::duration total_iter_time = clock::duration::zero();

    for (iter = 0; iter < opts->max_router_iterations && !routed; ++iter) {
        clock::duration greedy_route_time = clock::duration::zero();
        clock::duration update_cost_time = clock::duration::zero();
        clock::duration partitioning_time = clock::duration::zero();
        clock::duration analyze_timing_time = clock::duration::zero();
        clock::duration iter_time = clock::duration::zero();

        //for (int i = 0; i < opts->num_threads; ++i) {
            //sprintf(buffer, LOG_PATH_PREFIX"iter_%d_tid_%d.log", iter, i);
            //fclose(fopen(buffer, "w"));
        //}

		MPI_Barrier(cur_comm);
        
        auto iter_start = clock::now();

        //auto route_start = clock::now();

        //tbb::enumerable_thread_specific<state_t *> state_tls;
        //
		vector<perf_t> perfs(num_procs);
        vector<int> thread_num_nets_routed(num_procs, 0);
        vector<int> thread_num_nets_to_route(num_procs, 0);
        vector<int> thread_num_sinks_routed(num_procs, 0);
        vector<int> thread_num_sinks_to_route(num_procs, 0);
		vector<int> thread_num_interpartition_sinks(num_procs, 0);
        bool has_interpartition_sinks = false;
        vector<vector<interpartition_sink_t>> interpartition_sinks(nets.size()); 
		//vector<vector<RRNode>> net_sinks(nets.size());

		sprintf(buffer, "%d", iter);
		zlog_put_mdc("iter", buffer);

		sprintf(buffer, "%d", procid);
		zlog_put_mdc("tid", buffer);

        zlog_info(delta_log, "Routing iteration: %d\n", iter);
        printf("Routing iteration: %d\n", iter);

		vector<ongoing_transaction_t> transactions;

        auto greedy_route_start = clock::now();

		perfs[procid].num_heap_pushes = 0;
		perfs[procid].num_heap_pops = 0;
		perfs[procid].num_neighbor_visits = 0;

		for (int i = procid*pow(2, current_level); i < nets_to_route.size(); i += initial_num_procs) {
			for (int j = 0; j < pow(2, current_level) && i+j < nets_to_route.size(); ++j) {
				zlog_level(delta_log, ROUTER_V1, "Net index: %d\n", i+j);

				net_t *net = nets_to_route[i+j].second;

				update_sink_criticalities(*net, net_timing[net->vpr_id], params);

				//printf("Routing net %d\n", net->vpr_id);
				//vector<sink_t *> temp_routed_sinks = routed_sinks[net->local_id];

				//vector<RRNode> sinks_to_mark;
				//for (const auto &sink : temp_routed_sinks) {
					//bool fixed = find(begin(fixed_sinks[net->local_id]), end(fixed_sinks[net->local_id]), sink) != end(fixed_sinks[net->local_id]);
					//if (!fixed) {
						//sinks_to_mark.push_back(sink->rr_node);
					//} else {
						//routed_sinks[net->local_id].push_back(sink);
					//}
				//}

				//route_tree_mark_paths_to_be_ripped(route_trees[net->local_id], partitioner.orig_g, partitioner.result_pid_by_level[current_level], procid, sinks_to_mark);
				if (opts->rip_up_always) {
					route_tree_mark_all_nodes_to_be_ripped(route_trees[net->local_id], partitioner.orig_g);
				} else {
					route_tree_mark_congested_nodes_to_be_ripped(route_trees[net->local_id], partitioner.orig_g, congestion);
				}

				route_tree_rip_up_marked_mpi_send_recv(route_trees[net->local_id], partitioner.orig_g, congestion, params.pres_fac, procid, num_procs, cur_comm, transactions);

				routed_sinks[net->local_id].clear();
				for (auto &sink : net->sinks) {
					if (route_tree_get_rt_node(route_trees[net->local_id], sink.rr_node) != RouteTree::null_vertex()) {
						routed_sinks[net->local_id].push_back(&sink);
					}
				}

				vector<sink_t *> sinks;	
				get_sinks_to_route(net, route_trees[net->local_id], unroutable_sinks[net->local_id], sinks);

				//local_perf.total_rip_up_time += clock::now()-rip_up_start;

				//auto route_start = clock::now();

				if (!sinks.empty()) {
					int previous_num_unroutable_sinks = unroutable_sinks[net->local_id].size();

					int previous_num_routed_sinks = routed_sinks[net->local_id].size();

					route_net_mpi_send_recv(partitioner.orig_g, partitioner.result_pid_by_level[0], procid, net->vpr_id, &net->source, sinks, params, states, congestion, num_procs, cur_comm, transactions, route_trees[net->local_id], net_timing[net->vpr_id], routed_sinks[net->local_id], unroutable_sinks[net->local_id], &perfs[procid]);

					assert(routed_sinks[net->local_id].size() + unroutable_sinks[net->local_id].size() == net->sinks.size());

					if (!unroutable_sinks[net->local_id].empty() && previous_num_unroutable_sinks > 0) {
						assert(previous_num_unroutable_sinks == unroutable_sinks[net->local_id].size());
					}

					if (!has_unroutable_sinks) {
						has_unroutable_sinks = !unroutable_sinks[net->local_id].empty();
					}

					++thread_num_nets_routed[procid];
					++thread_num_nets_to_route[procid];

					thread_num_sinks_to_route[procid] += sinks.size();
					thread_num_sinks_routed[procid] += routed_sinks[net->local_id].size() - previous_num_routed_sinks;
				} else {
					assert(routed_sinks[net->local_id].size() + unroutable_sinks[net->local_id].size() == net->sinks.size());
					zlog_level(delta_log, ROUTER_V2, "Net %d has no sinks in the current iteration because there are %lu/%lu non-routable/all sinks\n", net->vpr_id, unroutable_sinks[net->local_id].size(), net->sinks.size());
				}

				//local_perf.total_route_time += clock::now()-rip_up_start;
			}
		}

		MPI_Barrier(cur_comm);

		sync(congestion, procid, num_procs, cur_comm);

		MPI_Barrier(cur_comm);

		//greedy_end_time = clock::now();

        greedy_route_time = clock::now()-greedy_route_start;

		//if (has_interpartition_sinks) {
			//for (const auto &net_interpartition_sinks : interpartition_sinks) {
				//for (const auto &is : net_interpartition_sinks) {
					//vector<RRNode> added_nodes;
					//for (const auto &node : is.path) {
						//if (node.update_cost) {
							//added_nodes.push_back(node.rr_node_id);
						//}
					//}
					//update_one_cost(partitioner.orig_g, congestion, added_nodes.begin(), added_nodes.end(), 1, params.pres_fac, false, nullptr);
				//}
			//}
		//}

        //if (greedy_rip_up_all) {
            //next_greedy_rip_up_iter += greedy_rip_up_all_period;
            //++greedy_rip_up_all_period;
            //prev_num_overused_nodes = std::numeric_limits<unsigned long>::max();
        //}

        //route_time = clock::now()-route_start;

        iter_time = clock::now()-iter_start;


		int total_num_sinks_to_route;
		if (procid == 0) {
			int total_num_nets_to_route = thread_num_nets_to_route[0];
			int total_num_nets_routed = thread_num_nets_routed[0];
			int total_num_sinks_routed = thread_num_sinks_routed[0];
			total_num_sinks_to_route = thread_num_sinks_to_route[0];

			for (int i = 1; i < num_procs; ++i) {
				MPI_Recv(&thread_num_nets_to_route[i], 1, MPI_INT, i, i, cur_comm, MPI_STATUS_IGNORE);
				MPI_Recv(&thread_num_nets_routed[i], 1, MPI_INT, i, i, cur_comm, MPI_STATUS_IGNORE);
				MPI_Recv(&thread_num_sinks_to_route[i], 1, MPI_INT, i, i, cur_comm, MPI_STATUS_IGNORE);
				MPI_Recv(&thread_num_sinks_routed[i], 1, MPI_INT, i, i, cur_comm, MPI_STATUS_IGNORE);

				total_num_nets_to_route += thread_num_nets_to_route[i];
				total_num_nets_routed += thread_num_nets_routed[i];
				total_num_sinks_routed += thread_num_sinks_routed[i];
				total_num_sinks_to_route += thread_num_sinks_to_route[i];
			}
			assert(total_num_nets_to_route == total_num_nets_routed);

			printf("num nets routed: ");
			for (int i = 0; i < num_procs; ++i) {
				printf("%d/%d (%g), ", thread_num_nets_routed[i], thread_num_nets_to_route[i], thread_num_nets_routed[i]*100.0/thread_num_nets_to_route[i]);
			}
			printf("\n");

			printf("num sinks routed: ");
			for (int i = 0; i < num_procs; ++i) {
				printf("%d/%d (%g), ", thread_num_sinks_routed[i], thread_num_sinks_to_route[i], thread_num_sinks_routed[i]*100.0/thread_num_sinks_to_route[i]);
			}
			printf("\n");

			printf("Total num nets routed: %d/%d (%g)\n", total_num_nets_routed, total_num_nets_to_route, total_num_nets_routed*100.0/total_num_nets_to_route);
			printf("Total num sinks routed: %d/%d (%g)\n", total_num_sinks_routed, total_num_sinks_to_route, (total_num_sinks_routed)*100.0/total_num_sinks_to_route);

			for (int i = 1; i < num_procs; ++i) {
				MPI_Recv(&perfs[i].num_heap_pushes, 1, MPI_INT, i, i, cur_comm, MPI_STATUS_IGNORE);
				MPI_Recv(&perfs[i].num_heap_pops, 1, MPI_INT, i, i, cur_comm, MPI_STATUS_IGNORE);
				MPI_Recv(&perfs[i].num_neighbor_visits, 1, MPI_INT, i, i, cur_comm, MPI_STATUS_IGNORE);
			}

			printf("num_heap_pushes: ");
			for (int i = 0; i < num_procs; ++i) {
				printf("%lu ", perfs[i].num_heap_pushes);
			}
			printf("\n");

			printf("num_heap_pops: ");
			for (int i = 0; i < num_procs; ++i) {
				printf("%lu ", perfs[i].num_heap_pops);
			}
			printf("\n");

			printf("num_neighbor_visits: ");
			for (int i = 0; i < num_procs; ++i) {
				printf("%lu ", perfs[i].num_neighbor_visits);
			}
			printf("\n");

			unsigned long total_num_heap_pushes = 0;
			unsigned long total_num_heap_pops = 0;
			unsigned long total_num_neighbor_visits = 0;

			for (int i = 0; i < num_procs; ++i) {
				total_num_heap_pushes += perfs[i].num_heap_pushes;
				total_num_heap_pops += perfs[i].num_heap_pops;
				total_num_neighbor_visits += perfs[i].num_neighbor_visits;
			}

			printf("total_num_heap_pushes: %lu\n", total_num_heap_pushes);
			printf("total_num_heap_pops: %lu\n", total_num_heap_pops); 
			printf("total_num_neighbor_visits: %lu\n", total_num_neighbor_visits);
		} else {
			MPI_Send(&thread_num_nets_to_route[procid], 1, MPI_INT, 0, procid, cur_comm);
			MPI_Send(&thread_num_nets_routed[procid], 1, MPI_INT, 0, procid, cur_comm);
			MPI_Send(&thread_num_sinks_to_route[procid], 1, MPI_INT, 0, procid, cur_comm);
			MPI_Send(&thread_num_sinks_routed[procid], 1, MPI_INT, 0, procid, cur_comm);

			MPI_Send(&perfs[procid].num_heap_pushes, 1, MPI_INT, 0, procid, cur_comm);
			MPI_Send(&perfs[procid].num_heap_pops, 1, MPI_INT, 0, procid, cur_comm);
			MPI_Send(&perfs[procid].num_neighbor_visits, 1, MPI_INT, 0, procid, cur_comm);
		}

		/* checking */
        for (int i = 0; i < num_vertices(partitioner.orig_g); ++i) {
            congestion[i].recalc_occ = 0; 
        }

		for (int i = procid*pow(2, current_level); i < nets_to_route.size(); i += initial_num_procs) {
			for (int j = 0; j < pow(2, current_level) && i+j < nets_to_route.size(); ++j) {
				net_t *net = nets_to_route[i+j].second;

				check_route_tree(route_trees[net->local_id], *net, routed_sinks[net->local_id], partitioner.orig_g);
				recalculate_occ(route_trees[net->local_id], partitioner.orig_g, congestion);
			}
        }

		//int *recalc_occ = new int[num_vertices(partitioner.orig_g)];
		//for (int i = 0; i < num_vertices(partitioner.orig_g); ++i) {
			//recalc_occ[i] = congestion[i].recalc_occ;	
		//}
		//if (procid == 0) {
			//int *recv_recalc_occ = new int[num_vertices(partitioner.orig_g)];
			//MPI_Reduce(recalc_occ, recv_recalc_occ, num_vertices(partitioner.orig_g), MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);
			//for (int i = 0; i < num_vertices(partitioner.orig_g); ++i) {
				//congestion[i].recalc_occ = recv_recalc_occ[i];
			//}
			//delete [] recv_recalc_occ;
		//} else {
			//MPI_Reduce(recalc_occ, nullptr, num_vertices(partitioner.orig_g), MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);
		//}
		//delete [] recalc_occ;
		
		//MPI_Op recalc_sum_op;
		//MPI_Op_create((MPI_User_function *)recalc_sum, 1, &recalc_sum_op);	
		//if (procid == 0) {
			//MPI_Reduce(MPI_IN_PLACE, congestion, num_vertices(partitioner.orig_g), recalc_occ_dt, recalc_sum_op, 0, MPI_COMM_WORLD);
		//} else {
			//MPI_Reduce(congestion, nullptr, num_vertices(partitioner.orig_g), recalc_occ_dt, recalc_sum_op, 0, MPI_COMM_WORLD);
		//}
		//MPI_Op_free(&recalc_sum_op);

		sync_recalc_occ(congestion, num_vertices(partitioner.orig_g), procid, num_procs, cur_comm);

		if (procid == 0) {
			bool valid = true;
			for (int i = 0; i < num_vertices(partitioner.orig_g); ++i) {
				sprintf_rr_node(i, buffer);
				if (congestion[i].recalc_occ != congestion[i].occ) {
					zlog_error(delta_log, "Node %s occ mismatch, recalc: %d original: %d\n", buffer, congestion[i].recalc_occ, congestion[i].occ);
					valid = false;
				}
			}
			assert(valid);
		}

        int overused_total_bb_rank = 0;
        int num_congested_nets = 0;

		for (int i = procid*pow(2, current_level); i < nets_to_route.size(); i += initial_num_procs) {
			for (int j = 0; j < pow(2, current_level) && i+j < nets_to_route.size(); ++j) {
				net_t *net = nets_to_route[i+j].second;

				vector<int> overused_rr_node;
				assert(route_trees[net->local_id].root_rt_nodes.size() == 1);
				get_overused_nodes(route_trees[net->local_id], route_trees[net->local_id].root_rt_nodes[0], partitioner.orig_g, congestion, overused_rr_node);
				if (!overused_rr_node.empty()) {
					zlog_level(delta_log, ROUTER_V1, "Net %d bb_rank %d has %lu overused nodes:\n", net->vpr_id, net->bb_area_rank, overused_rr_node.size());
					for (const auto &item : overused_rr_node) {
						zlog_level(delta_log, ROUTER_V1, "%d ", item);
					}
					zlog_level(delta_log, ROUTER_V1, "\n");
					overused_total_bb_rank += net->bb_area_rank;
					++num_congested_nets;
				}
			}
		}

        zlog_level(delta_log, ROUTER_V1, "Average overused net bb rank: %g\n", (float)overused_total_bb_rank/num_congested_nets);
        zlog_level(delta_log, ROUTER_V1, "Num congested nets: %d\n", num_congested_nets);

        iter_start = clock::now();

        auto analyze_timing_start = clock::now();

		int *recvcounts = nullptr;
		int *displs = nullptr;

		init_displ(num_procs, current_level, nets_to_route, initial_num_procs, &recvcounts, &displs);

		int num_crits = 0;
		for (int i = 0; i < num_procs; ++i) {
			num_crits += recvcounts[i];
		}

		sync_net_delay(nets_to_route, procid, num_procs, initial_num_procs, recvcounts, displs, current_level, cur_comm, net_timing);

		float *all_crits;
		int idx;

		if (procid == 0) {
			crit_path_delay = analyze_timing(net_timing);

			all_crits = new float[num_crits];

			idx = 0;
			for (int pi = 0; pi < num_procs; ++pi) {
				idx = 0;
				for (int i = pi*pow(2, current_level); i < nets_to_route.size(); i += initial_num_procs) {
					for (int j = 0; j < pow(2, current_level) && i+j < nets_to_route.size(); ++j) {
						net_t *net = nets_to_route[i+j].second;
						zlog_level(delta_log, ROUTER_V3, "Net %d send crits:\n", net->vpr_id);
						for (int k = 1; k <= net->sinks.size(); ++k) {
							all_crits[displs[pi] + idx] = net_timing[net->vpr_id].timing_criticality[k];
							zlog_level(delta_log, ROUTER_V3, "\t%g\n", all_crits[displs[pi] + idx]);
							++idx;
						}
					}
				}
				assert(idx == recvcounts[pi]);
			}
		} else {
			all_crits = nullptr;
		}

		float *crits = new float[recvcounts[procid]];

		MPI_Scatterv(all_crits, recvcounts, displs, MPI_FLOAT, crits, recvcounts[procid], MPI_FLOAT, 0, cur_comm);

		idx = 0;
		for (int i = procid*pow(2, current_level); i < nets_to_route.size(); i += initial_num_procs) {
			for (int j = 0; j < pow(2, current_level) && i+j < nets_to_route.size(); ++j) {
				net_t *net = nets_to_route[i+j].second;
				zlog_level(delta_log, ROUTER_V3, "Net %d recv crits:\n", net->vpr_id);

				for (int s = 1; s <= net->sinks.size(); ++s) {
					net_timing[net->vpr_id].timing_criticality[s] = crits[idx];
					zlog_level(delta_log, ROUTER_V3, "\t%g\n", net_timing[net->vpr_id].timing_criticality[s]);
					++idx;
				}
			}
		}

		if (procid == 0) {
			delete [] all_crits;
		}
		delete [] crits;

        analyze_timing_time = clock::now()-analyze_timing_start;


		int m_routed = (feasible_routing(partitioner.orig_g, congestion) && !has_unroutable_sinks) ? 1 : 0;
		int reduced_routed;
		MPI_Allreduce(&m_routed, &reduced_routed, 1, MPI_INT, MPI_LAND, cur_comm);

        zlog_level(delta_log, ROUTER_V1, "m_routed: %d reduced_routed: %d\n", m_routed, reduced_routed);

        if (reduced_routed) {
            //dump_route(*current_traces_ptr, "route.txt");
			routed = true;	
        } else {
            unsigned long num_overused_nodes = 0;
			vector<int> overused_nodes_by_type(NUM_RR_TYPES, 0);
            for (int i = 0; i < num_vertices(partitioner.orig_g); ++i) {
                if (congestion[i].occ > get_vertex_props(partitioner.orig_g, i).capacity) {
					const auto &v_p = get_vertex_props(partitioner.orig_g, i);
					++overused_nodes_by_type[v_p.type];

                    ++num_overused_nodes;
                }
            }

			static const char *name_type[] = { "SOURCE", "SINK", "IPIN", "OPIN",
				"CHANX", "CHANY", "INTRA_CLUSTER_EDGE" };
            zlog_info(delta_log, "Num overused nodes: %lu/%d (%.2f)\n", num_overused_nodes, num_vertices(partitioner.orig_g), num_overused_nodes*100.0/num_vertices(partitioner.orig_g));
			for (int i = 0; i < overused_nodes_by_type.size(); ++i) {
				zlog_info(delta_log, "\t%s: %d (%g)\n", name_type[i], overused_nodes_by_type[i], overused_nodes_by_type[i]*100.0/num_overused_nodes);
			}

			int *all_overused_nodes_by_type = new int[num_procs*NUM_RR_TYPES];
			int *overused_nodes_by_type_send = new int[NUM_RR_TYPES];
			for (int i = 0; i < overused_nodes_by_type.size(); ++i) {
				overused_nodes_by_type_send[i] = overused_nodes_by_type[i];
			}

			MPI_Gather(overused_nodes_by_type_send, NUM_RR_TYPES, MPI_INT, all_overused_nodes_by_type, NUM_RR_TYPES, MPI_INT, 0, cur_comm);

			unsigned long *all_num_overused_nodes = new unsigned long[num_procs];
			MPI_Gather(&num_overused_nodes, 1, MPI_UNSIGNED_LONG, all_num_overused_nodes, 1, MPI_UNSIGNED_LONG, 0, cur_comm);

			if (procid == 0) {
				printf("Num overused nodes: ");
				for (int i = 0; i < num_procs; ++i) {
					printf("%lu/%d (%.2f) ", all_num_overused_nodes[i], num_vertices(partitioner.orig_g), all_num_overused_nodes[i]*100.0/num_vertices(partitioner.orig_g));
				}
				printf("\n");

				for (int i = 0; i < overused_nodes_by_type.size(); ++i) {
					printf("\t%s: ", name_type[i]);
					for (int j = 0; j < num_procs; ++j) {
						printf("%d (%g) ", all_overused_nodes_by_type[j*NUM_RR_TYPES+i], all_overused_nodes_by_type[j*NUM_RR_TYPES+i]*100.0/all_num_overused_nodes[j]);
					}
					printf("\n");
				}
			} 

			delete [] all_overused_nodes_by_type;
			delete [] overused_nodes_by_type_send;
			delete [] all_num_overused_nodes;

			int not_decreasing = (num_overused_nodes > prev_num_overused_nodes && iter > 10) ? 1 : 0;
			//int not_decreasing = current_level+1 < partitioner.result_pid_by_level.size(); [> testing <]
			//int not_decreasing = current_level+1 <= std::log2(initial_num_procs); [> testing <]
			int reduced_not_decreasing;
			MPI_Allreduce(&not_decreasing, &reduced_not_decreasing, 1, MPI_INT, MPI_LOR, cur_comm);

			prev_num_overused_nodes = num_overused_nodes;

			zlog_level(delta_log, ROUTER_V1, "not_decreasing: %d reduced_not_decreasing: %d\n", not_decreasing, reduced_not_decreasing);

			if (reduced_not_decreasing) {
				/* need to send route tree over */
				if (procid % 2 == 0) {
					/*receiver*/
					for (int i = (procid+1)*pow(2, current_level); i < nets_to_route.size(); i += initial_num_procs) {
						for (int j = 0; j < pow(2, current_level) && i+j < nets_to_route.size(); ++j) {
							net_t *net = nets_to_route[i+j].second;

							zlog_level(delta_log, ROUTER_V3, "Recving net index %d from %d\n", i+j, procid+1);

							recv_route_tree(net, partitioner.orig_g, routed_sinks, states, congestion, params.pres_fac, route_trees, net_timing, procid+1, cur_comm);
						}
					}
				} else {
					/*sender*/
					for (int i = procid*pow(2, current_level); i < nets_to_route.size(); i += initial_num_procs) {
						for (int j = 0; j < pow(2, current_level) && i+j < nets_to_route.size(); ++j) {
							net_t *net = nets_to_route[i+j].second;

							zlog_level(delta_log, ROUTER_V3, "Sending net index %d from %d\n", i+j, procid);

							send_route_tree(net, partitioner.orig_g, routed_sinks, route_trees, procid-1, cur_comm);
						}
					}
					
				}

				assert(num_procs % 2 == 0);
				num_procs /= 2;
				MPI_Comm new_comm;
				MPI_Comm_split(cur_comm, procid%2, procid, &new_comm);

				if (procid % 2 == 1) {
					/* early exit code */
					break;
				}

				cur_comm = new_comm;
				MPI_Comm_rank(cur_comm, &procid);

				++current_level;

				//assert(current_level < partitioner.result_pid_by_level.size());

				printf("[%d] Transitioned to level %d at iteration %d\n", procid, current_level, iter);
				zlog_level(delta_log, ROUTER_V1, "Transitioned to level %d at iteration %d\n", current_level, iter);
				zlog_level(delta_log, ROUTER_V1, "New pid %d for initial pid %d\n", procid, initial_procid);

				for (int i = procid*pow(2, current_level); i < nets_to_route.size(); i += initial_num_procs) {
					for (int j = 0; j < pow(2, current_level) && i+j < nets_to_route.size(); ++j) {
						net_t *net = nets_to_route[i+j].second;

						unroutable_sinks[net->local_id].clear();

						//for (const auto &rs : routed_sinks[net->local_id]) {
							//assert(find(begin(fixed_sinks[net->local_id]), end(fixed_sinks[net->local_id]), rs) == end(fixed_sinks[net->local_id]));

							//zlog_level(delta_log, ROUTER_V3, "Fixing net %d sink %d\n", net->vpr_id, rs->id);

							//fixed_sinks[net->local_id].push_back(rs);
						//}
					}
				}

				prev_num_overused_nodes = std::numeric_limits<unsigned long>::max();
				has_unroutable_sinks = false;

				for (int i = 0; i < num_vertices(partitioner.orig_g); ++i) {
					congestion[i].recalc_occ = 0; 
				}

				for (int i = procid*pow(2, current_level); i < nets_to_route.size(); i += initial_num_procs) {
					for (int j = 0; j < pow(2, current_level) && i+j < nets_to_route.size(); ++j) {
						net_t *net = nets_to_route[i+j].second;

						if (!routed_sinks[net->local_id].empty()) {
							zlog_level(delta_log, ROUTER_V3, "Checking net index %d\n", i+j);

							check_route_tree(route_trees[net->local_id], *net, routed_sinks[net->local_id], partitioner.orig_g);
							recalculate_occ(route_trees[net->local_id], partitioner.orig_g, congestion);
						} else {
							zlog_level(delta_log, ROUTER_V3, "Not checking net index %d because of empty route tree\n", i+j);
						}
					}
				}

				sync_recalc_occ(congestion, num_vertices(partitioner.orig_g),  procid, num_procs, cur_comm);

				if (procid == 0) {
					bool valid = true;
					for (int i = 0; i < num_vertices(partitioner.orig_g); ++i) {
						sprintf_rr_node(i, buffer);
						if (congestion[i].recalc_occ != congestion[i].occ) {
							zlog_error(delta_log, "Node %s occ mismatch, recalc: %d original: %d\n", buffer, congestion[i].recalc_occ, congestion[i].occ);
							valid = false;
						}
					}
					assert(valid);
				}
			} 

			//if ((float)total_num_interpartition_sinks/total_num_sinks_to_route > 0.1) {
				//++current_level;
				//assert(current_level < partitioner.result_pid_by_level.size());
			//}

            auto update_cost_start = clock::now();

            if (iter == 0) {
                params.pres_fac = opts->initial_pres_fac;
                //update_costs_mpi(partitioner.orig_g[0], partitioner.result_pid_by_level[current_level], procid, congestion, win, params.pres_fac, 0);
				update_costs(partitioner.orig_g, congestion, params.pres_fac, 0);
            } else {
                params.pres_fac *= opts->pres_fac_mult;

				/* Avoid overflow for high iteration counts, even if acc_cost is big */
                params.pres_fac = std::min(params.pres_fac, static_cast<float>(HUGE_POSITIVE_FLOAT / 1e5));

                //update_costs_mpi(partitioner.orig_g[0], partitioner.result_pid_by_level[current_level], procid, congestion, win, params.pres_fac, opts->acc_fac);
				update_costs(partitioner.orig_g, congestion, params.pres_fac, opts->acc_fac);
            }

            update_cost_time = clock::now()-update_cost_start;
        }

        iter_time += clock::now()-iter_start;

        //total_route_time += route_time;
        total_greedy_route_time += greedy_route_time;
        total_update_cost_time += update_cost_time;
        total_analyze_timing_time += analyze_timing_time;
        total_iter_time += iter_time;

		if (procid == 0) {
			printf("Iteration time: %g s.\n", duration_cast<nanoseconds>(iter_time).count() / 1e9);
			printf("\tRoute time: %g s.\n", duration_cast<nanoseconds>(greedy_route_time).count() / 1e9);
			printf("\tUpdate cost time: %g s.\n", duration_cast<nanoseconds>(update_cost_time).count() / 1e9);
			printf("\tAnalyze timing time: %g s.\n", duration_cast<nanoseconds>(analyze_timing_time).count() / 1e9);
			printf("Critical path: %g ns\n", crit_path_delay);
		}

        //clock::time_point greedy_earliest_end_time = *std::min_element(begin(greedy_end_time), begin(greedy_end_time)+num_procs);

        //printf("greedy wait time: ");
        //for (int i = 0; i < num_procs; ++i) {
			//printf("%g (%g), ", duration_cast<nanoseconds>(greedy_end_time[i]-greedy_earliest_end_time).count() / 1e9, 100.0*(greedy_end_time[i]-greedy_earliest_end_time)/greedy_route_time);
        //}
		//printf("\n");
    }

    //sprintf(buffer, "%s_run_%d_rr_graph_occ.txt", s_circuit_name, run);
    //dump_rr_graph_occ(congestion, num_vertices(g), buffer);
	if (initial_procid == 0) {
		if (routed) {
			printf("Routed in %d iterations. Total iteration time: %g\n", iter, duration_cast<nanoseconds>(total_iter_time).count() / 1e9);
			printf("\tTotal route time: %g s.\n", duration_cast<nanoseconds>(total_greedy_route_time).count() / 1e9);
			printf("\t\tTotal greedy route time: %g s.\n", duration_cast<nanoseconds>(total_greedy_route_time).count() / 1e9);
			printf("\tTotal update cost time: %g s.\n", duration_cast<nanoseconds>(total_update_cost_time).count() / 1e9);
			printf("\tTotal analyze timing time: %g s.\n", duration_cast<nanoseconds>(total_analyze_timing_time).count() / 1e9);

			printf("Final critical path: %g ns\n", crit_path_delay);
		} else {
			printf("Failed to route in %d iterations. Total iteration time: %g\n", opts->max_router_iterations, duration_cast<nanoseconds>(total_iter_time).count() / 1e9);
			printf("\tTotal route time: %g s.\n", duration_cast<nanoseconds>(total_greedy_route_time).count() / 1e9);
			printf("\t\tTotal greedy route time: %g s.\n", duration_cast<nanoseconds>(total_greedy_route_time).count() / 1e9);
			printf("\tTotal update cost time: %g s.\n", duration_cast<nanoseconds>(total_update_cost_time).count() / 1e9);
			printf("\tTotal analyze timing time: %g s.\n", duration_cast<nanoseconds>(total_analyze_timing_time).count() / 1e9);
		}
	} 

	MPI_Barrier(MPI_COMM_WORLD);

	MPI_Finalize();
	exit(0);

    //delete_graph(g);
    delete_net_timing(nets, global_nets, net_timing);	
    delete [] congestion;
    delete [] net_timing;
    delete [] states;

    return routed;
}
//
bool mpi_spatial_route_simulated(t_router_opts *opts, struct s_det_routing_arch det_routing_arch, t_direct_inf *directs, int num_directs, t_segment_inf *segment_inf, t_timing_inf timing_inf)
{
    using std::chrono::duration_cast;
    using std::chrono::nanoseconds;
    using clock = std::chrono::high_resolution_clock;

	int num_procs = opts->num_threads;

    init_logging();
    zlog_set_record("custom_output", delta_log_output);
    zlog_set_record("missing_edge", missing_edge_log_output);
    zlog_set_record("ss", ss_log_output);
    delta_log_files.resize(opts->max_router_iterations);
    for (int i = 0; i < opts->max_router_iterations; ++i) {
        delta_log_files[i].resize(num_procs, nullptr);
    }
    missing_edge_log_files.resize(opts->max_router_iterations);
    for (int i = 0; i < opts->max_router_iterations; ++i) {
        missing_edge_log_files[i].resize(num_procs, nullptr);
    }

    //fclose(fopen(LOG_PATH_PREFIX"iter__tid_.log", "w"));

    test_fast_graph();
    test_topo();
    test_fm();
    test_filter_graph();
    test_partition_graph();
    test_connected_components();

    vector<RRGraph *> graphs;
    rr_graph_partitioner partitioner;
    partitioner.init(opts, &det_routing_arch, directs, num_directs, segment_inf, &timing_inf);
    partitioner.partition_without_ipin(num_procs, graphs);
    //for (const auto &g : graphs) {
        //routability(*g);
    //}
    //
    //RRGraph combined_g;
    //add_vertex(combined_g, num_vertices(partitioner.orig_g));

    //for (const auto &g : graphs) {
        //for (const auto &e : get_edges(*g)) {
            //int from = get_source(*g, e);
            //int to = get_target(*g, e);
            //const auto &from_ver = get_vertex_props(*g, from);
            //const auto &to_ver = get_vertex_props(*g, to);

            //if (is_channel(from_ver) && is_channel(to_ver)) {
                //assert(!has_edge(combined_g, from, to));
                //add_edge(combined_g, from, to);
            //}
        //}
    //}

    //for (const auto &e : get_edges(partitioner.orig_g)) {
        //int from = get_source(partitioner.orig_g, e);
        //int to = get_target(partitioner.orig_g, e);
        //const auto &from_ver = get_vertex_props(partitioner.orig_g, from);
        //const auto &to_ver = get_vertex_props(partitioner.orig_g, to);
        //if (!is_channel(from_ver) || !is_channel(to_ver)) {
            //assert(!has_edge(combined_g, from, to));
            //add_edge(combined_g, from, to);
        //}
    //}

    //printf("Combined/Orig graph has %d/%d (%g) edges.\n", num_edges(combined_g), num_edges(partitioner.orig_g), 100.0*num_edges(combined_g)/num_edges(partitioner.orig_g));

    //goto lol;

    //extern int num_types;
    //extern struct s_type_descriptor *type_descriptors;
    //extern int nx, ny;
    //extern struct s_grid_tile **grid;

    //free_rr_graph();

    //int warnings;

    //build_rr_graph(GRAPH_UNIDIR, num_types, type_descriptors, nx, ny, grid,
            //opts->fixed_channel_width, NULL, det_routing_arch.switch_block_type,
            //det_routing_arch.Fs, det_routing_arch.num_segment,
            //det_routing_arch.num_switch, segment_inf,
            //det_routing_arch.global_route_switch,
            //det_routing_arch.delayless_switch, timing_inf,
            //det_routing_arch.wire_to_ipin_switch, opts->base_cost_type,
            //directs, num_directs, FALSE,
            //&warnings);

    //RRGraph channel_with_interior_g;
    //init_channel_only_graph(channel_with_interior_g);

    //dump_rr_graph(channel_with_interior_g, "/Volumes/DATA/rr_graph_channel_with_interior.txt");

    //RRGraph orig_g;
    //init_graph(orig_g);

    //dump_rr_graph(orig_g, "/Volumes/DATA/rr_graph.txt");

    //free_rr_graph();
    //for (int i = 0; i < det_routing_arch.num_segment; ++i) {
        //for (int j = 0; j < segment_inf[i].sb_len; ++j) {
            //if (j != 0 && j != segment_inf[i].sb_len-1) {
                //segment_inf[i].sb[j] = FALSE;
            //}
        //}
    //}
    //build_rr_graph(GRAPH_UNIDIR, num_types, type_descriptors, nx, ny, grid,
            //opts->fixed_channel_width, NULL, det_routing_arch.switch_block_type,
            //det_routing_arch.Fs, det_routing_arch.num_segment,
            //det_routing_arch.num_switch, segment_inf,
            //det_routing_arch.global_route_switch,
            //det_routing_arch.delayless_switch, timing_inf,
            //det_routing_arch.wire_to_ipin_switch, opts->base_cost_type,
            //directs, num_directs, FALSE,
            //&warnings);

    //RRGraph channel_without_interior_g;
    //init_channel_only_graph(channel_without_interior_g);

    //dump_rr_graph(channel_without_interior_g, "/Volumes/DATA/rr_graph_channel_without_interior.txt");

    //vector<vector<vector<int>>> all_partition_components(opts->num_threads);
    ////for (int x = 0; x < nx+1; ++x) {
        ////for (int y = 0; y < ny+1; ++y) {
            ////init_partitioned_graph_5(opts->num_threads, channel_with_interior_g, channel_without_interior_g, orig_g, graphs, x, y);
        ////}
    ////}
    ////init_partitioned_graph_4(opts->num_threads, channel_with_interior_g, channel_without_interior_g, orig_g, graphs, 14, 14);
    //init_partitioned_graph_4(opts->num_threads, channel_with_interior_g, channel_without_interior_g, orig_g, graphs, all_partition_components);

    //for (int i = 0; i < graphs.size(); ++i) {
        //char filename[256];

        //sprintf(filename, "/Volumes/DATA/rr_graph_%d.txt", i);
        //dump_rr_graph(*graphs[i], filename);

        //sprintf(filename, "/Volumes/DATA/rr_graph_%d_edges.txt", i);
        //dump_edges(*graphs[i], filename);
    //}

    vector<net_t> nets;
    vector<net_t> global_nets;
    //init_nets(nets, global_nets, opts->bb_factor, partitioner.sink_in_nodes);	
    init_nets(nets, global_nets, opts->bb_factor, opts->large_bb);	

    vector<net_t *> nets_ptr(nets.size());
    for (int i = 0; i < nets.size(); ++i) {
        nets_ptr[i] = &nets[i];
    }
    extern s_bb *route_bb;
    sort(nets_ptr.begin(), nets_ptr.end(), [] (const net_t *a, const net_t *b) -> bool {
            return get_bounding_box_area(route_bb[a->vpr_id]) > get_bounding_box_area(route_bb[b->vpr_id]);
            });
    int rank = 0;
    for (auto &net : nets_ptr) {
        net->bb_area_rank = rank++;
    }

	vector<route_state_t *> states;
	congestion_locked_t *congestion;
	vector<route_tree_t> route_trees;
	t_net_timing *net_timing;
	init_route_structs_locked(partitioner.orig_g, nets, global_nets, opts->num_threads, states, &congestion, route_trees, &net_timing);

    route_parameters_t params;
    params.criticality_exp = opts->criticality_exp;
    params.astar_fac = opts->astar_fac;
    params.max_criticality = opts->max_criticality;
    params.bend_cost = opts->bend_cost;
	params.pres_fac = opts->first_iter_pres_fac; /* Typically 0 -> ignore cong. */

    char buffer[256];

    vector<pair<box, net_t *>> nets_to_route;
    
    for (auto &net : nets) {
        box b = bg::make_inverse<box>();

        bg::expand(b, point(net.source.x, net.source.y));
        for (const auto &sink : net.sinks) {
            bg::expand(b, point(sink.x, sink.y));
        }
        bg::subtract_value(b.min_corner(), 1+opts->bb_factor);
        bg::add_value(b.max_corner(), opts->bb_factor);

        nets_to_route.push_back(make_pair(b, &net));
    }
    std::sort(begin(nets_to_route), end(nets_to_route), [] (const pair<box, net_t *> &a, const pair<box, net_t *> &b) -> bool {
            return a.second->sinks.size() > b.second->sinks.size();
            });

    bool routed = false;

    int iter;
    float crit_path_delay;
	int current_level = 0;

	vector<vector<sink_t *>> unroutable_sinks(nets.size());
	//vector<vector<sink_t *>> fixed_sinks(nets.size());
	vector<vector<sink_t *>> routed_sinks(nets.size());
	bool has_unroutable_sinks = false;
	unsigned long prev_num_overused_nodes = std::numeric_limits<unsigned long>::max();

    clock::duration total_greedy_route_time = clock::duration::zero();
    clock::duration total_update_cost_time = clock::duration::zero();
    clock::duration total_analyze_timing_time = clock::duration::zero();
    clock::duration total_iter_time = clock::duration::zero();

    for (iter = 0; iter < opts->max_router_iterations && !routed; ++iter) {
        clock::duration greedy_route_time = clock::duration::zero();
        clock::duration update_cost_time = clock::duration::zero();
        clock::duration partitioning_time = clock::duration::zero();
        clock::duration analyze_timing_time = clock::duration::zero();
        clock::duration iter_time = clock::duration::zero();

        //for (int i = 0; i < opts->num_threads; ++i) {
            //sprintf(buffer, LOG_PATH_PREFIX"iter_%d_tid_%d.log", iter, i);
            //fclose(fopen(buffer, "w"));
        //}

        auto iter_start = clock::now();

        //auto route_start = clock::now();

        //tbb::enumerable_thread_specific<state_t *> state_tls;
        //
		vector<perf_t> perfs(num_procs);
        vector<int> thread_num_nets_routed(num_procs, 0);
        vector<int> thread_num_nets_to_route(num_procs, 0);
        vector<int> thread_num_sinks_routed(num_procs, 0);
        vector<int> thread_num_sinks_to_route(num_procs, 0);
		//vector<vector<RRNode>> net_sinks(nets.size());
		auto greedy_route_start = clock::now();

		tbb::parallel_for(tbb::blocked_range<int>(0, opts->num_threads, 1),
				[&] (const tbb::blocked_range<int> &range) {

				assert(range.end() - range.begin() == 1);

				int procid = range.begin();

				char local_buffer[256];

				sprintf(local_buffer, "%d", iter);
				zlog_put_mdc("iter", local_buffer);

				sprintf(local_buffer, "%d", procid);
				zlog_put_mdc("tid", local_buffer);

				zlog_info(delta_log, "Routing iteration: %d\n", iter);
				printf("Routing iteration: %d\n", iter);


				perf_t local_perf;

				local_perf.num_heap_pushes = 0;
				local_perf.num_heap_pops = 0;
				local_perf.num_neighbor_visits = 0;

				int i = procid;
				while (i < nets_to_route.size()) {
					net_t *net = nets_to_route[i].second;

					update_sink_criticalities(*net, net_timing[net->vpr_id], params);

					//printf("Routing net %d\n", net->vpr_id);
					//vector<sink_t *> temp_routed_sinks = routed_sinks[net->local_id];

					//vector<RRNode> sinks_to_mark;
					//for (const auto &sink : temp_routed_sinks) {
					//bool fixed = find(begin(fixed_sinks[net->local_id]), end(fixed_sinks[net->local_id]), sink) != end(fixed_sinks[net->local_id]);
					//if (!fixed) {
					//sinks_to_mark.push_back(sink->rr_node);
					//} else {
					//routed_sinks[net->local_id].push_back(sink);
					//}
					//}

					//route_tree_mark_paths_to_be_ripped(route_trees[net->local_id], partitioner.orig_g, partitioner.result_pid_by_level[current_level], procid, sinks_to_mark);
					route_tree_mark_congested_nodes_to_be_ripped(route_trees[net->local_id], partitioner.orig_g, congestion);

					route_tree_rip_up_marked(route_trees[net->local_id], partitioner.orig_g, congestion, params.pres_fac, true, nullptr);

					routed_sinks[net->local_id].clear();
					for (auto &sink : net->sinks) {
						if (route_tree_get_rt_node(route_trees[net->local_id], sink.rr_node) != RouteTree::null_vertex()) {
							routed_sinks[net->local_id].push_back(&sink);
						}
					}

					vector<sink_t *> sinks;	
					get_sinks_to_route(net, route_trees[net->local_id], unroutable_sinks[net->local_id], sinks);

					//local_perf.total_rip_up_time += clock::now()-rip_up_start;

					//auto route_start = clock::now();

					if (!sinks.empty()) {
						int previous_num_unroutable_sinks = unroutable_sinks[net->local_id].size();

						int previous_num_routed_sinks = routed_sinks[net->local_id].size();

						route_net_with_partitioned_fine_grain_lock(partitioner.orig_g, partitioner.result_pid_by_level[current_level], procid, net->vpr_id, &net->source, sinks, params, states[procid], congestion, route_trees[net->local_id], net_timing[net->vpr_id], routed_sinks[net->local_id], unroutable_sinks[net->local_id], true, &local_perf, nullptr);

						assert(routed_sinks[net->local_id].size() + unroutable_sinks[net->local_id].size() == net->sinks.size());

						if (!unroutable_sinks[net->local_id].empty() && previous_num_unroutable_sinks > 0) {
							assert(previous_num_unroutable_sinks == unroutable_sinks[net->local_id].size());
						}

						if (!has_unroutable_sinks) {
							has_unroutable_sinks = !unroutable_sinks[net->local_id].empty();
						}

						++thread_num_nets_routed[procid];
						++thread_num_nets_to_route[procid];

						thread_num_sinks_to_route[procid] += sinks.size();
						thread_num_sinks_routed[procid] += routed_sinks[net->local_id].size() - previous_num_routed_sinks;
					} else {
						assert(routed_sinks[net->local_id].size() + unroutable_sinks[net->local_id].size() == net->sinks.size());
						zlog_level(delta_log, ROUTER_V2, "Net %d has no sinks in the current iteration because there are %lu/%lu non-routable/all sinks\n", net->vpr_id, unroutable_sinks[net->local_id].size(), net->sinks.size());
					}

				//local_perf.total_route_time += clock::now()-rip_up_start;
					i += num_procs;
				}

				perfs[procid] = local_perf;

				});

		//greedy_end_time = clock::now();

        greedy_route_time = clock::now()-greedy_route_start;

		//if (has_interpartition_sinks) {
			//for (const auto &net_interpartition_sinks : interpartition_sinks) {
				//for (const auto &is : net_interpartition_sinks) {
					//vector<RRNode> added_nodes;
					//for (const auto &node : is.path) {
						//if (node.update_cost) {
							//added_nodes.push_back(node.rr_node_id);
						//}
					//}
					//update_one_cost(partitioner.orig_g, congestion, added_nodes.begin(), added_nodes.end(), 1, params.pres_fac, false, nullptr);
				//}
			//}
		//}

        //if (greedy_rip_up_all) {
            //next_greedy_rip_up_iter += greedy_rip_up_all_period;
            //++greedy_rip_up_all_period;
            //prev_num_overused_nodes = std::numeric_limits<unsigned long>::max();
        //}

        //route_time = clock::now()-route_start;

        iter_time = clock::now()-iter_start;


		int total_num_nets_to_route = thread_num_nets_to_route[0];
		int total_num_nets_routed = thread_num_nets_routed[0];
		int total_num_sinks_routed = thread_num_sinks_routed[0];
		int total_num_sinks_to_route = thread_num_sinks_to_route[0];

		for (int i = 1; i < num_procs; ++i) {
			total_num_nets_to_route += thread_num_nets_to_route[i];
			total_num_nets_routed += thread_num_nets_routed[i];
			total_num_sinks_routed += thread_num_sinks_routed[i];
			total_num_sinks_to_route += thread_num_sinks_to_route[i];
		}
		assert(total_num_nets_to_route == total_num_nets_routed);

		printf("num nets routed: ");
		for (int i = 0; i < num_procs; ++i) {
			printf("%d/%d (%g), ", thread_num_nets_routed[i], thread_num_nets_to_route[i], thread_num_nets_routed[i]*100.0/thread_num_nets_to_route[i]);
		}
		printf("\n");

		printf("num sinks routed: ");
		for (int i = 0; i < num_procs; ++i) {
			printf("%d/%d (%g), ", thread_num_sinks_routed[i], thread_num_sinks_to_route[i], thread_num_sinks_routed[i]*100.0/thread_num_sinks_to_route[i]);
		}
		printf("\n");

		printf("Total num nets routed: %d/%d (%g)\n", total_num_nets_routed, total_num_nets_to_route, total_num_nets_routed*100.0/total_num_nets_to_route);
		printf("Total num sinks routed: %d/%d (%g)\n", total_num_sinks_routed, total_num_sinks_to_route, (total_num_sinks_routed)*100.0/total_num_sinks_to_route);

		printf("num_heap_pushes: ");
		for (int i = 0; i < num_procs; ++i) {
			printf("%lu ", perfs[i].num_heap_pushes);
		}
		printf("\n");

		printf("num_heap_pops: ");
		for (int i = 0; i < num_procs; ++i) {
			printf("%lu ", perfs[i].num_heap_pops);
		}
		printf("\n");

		printf("num_neighbor_visits: ");
		for (int i = 0; i < num_procs; ++i) {
			printf("%lu ", perfs[i].num_neighbor_visits);
		}
		printf("\n");

		unsigned long total_num_heap_pushes = 0;
		unsigned long total_num_heap_pops = 0;
		unsigned long total_num_neighbor_visits = 0;

		for (int i = 0; i < num_procs; ++i) {
			total_num_heap_pushes += perfs[i].num_heap_pushes;
			total_num_heap_pops += perfs[i].num_heap_pops;
			total_num_neighbor_visits += perfs[i].num_neighbor_visits;
		}

		printf("total_num_heap_pushes: %lu\n", total_num_heap_pushes);
		printf("total_num_heap_pops: %lu\n", total_num_heap_pops); 
		printf("total_num_neighbor_visits: %lu\n", total_num_neighbor_visits);

		/* checking */
        for (int i = 0; i < num_vertices(partitioner.orig_g); ++i) {
            congestion[i].cong.recalc_occ = 0; 
        }

		for (int i = 0; i < nets_to_route.size(); ++i) {
			net_t *net = nets_to_route[i].second;

			check_route_tree(route_trees[net->local_id], *net, routed_sinks[net->local_id], partitioner.orig_g);
			recalculate_occ(route_trees[net->local_id], partitioner.orig_g, congestion);
        }

		//int *recalc_occ = new int[num_vertices(partitioner.orig_g)];
		//for (int i = 0; i < num_vertices(partitioner.orig_g); ++i) {
			//recalc_occ[i] = congestion[i].recalc_occ;	
		//}
		//if (procid == 0) {
			//int *recv_recalc_occ = new int[num_vertices(partitioner.orig_g)];
			//MPI_Reduce(recalc_occ, recv_recalc_occ, num_vertices(partitioner.orig_g), MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);
			//for (int i = 0; i < num_vertices(partitioner.orig_g); ++i) {
				//congestion[i].recalc_occ = recv_recalc_occ[i];
			//}
			//delete [] recv_recalc_occ;
		//} else {
			//MPI_Reduce(recalc_occ, nullptr, num_vertices(partitioner.orig_g), MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);
		//}
		//delete [] recalc_occ;
		
		//MPI_Op recalc_sum_op;
		//MPI_Op_create((MPI_User_function *)recalc_sum, 1, &recalc_sum_op);	
		//if (procid == 0) {
			//MPI_Reduce(MPI_IN_PLACE, congestion, num_vertices(partitioner.orig_g), recalc_occ_dt, recalc_sum_op, 0, MPI_COMM_WORLD);
		//} else {
			//MPI_Reduce(congestion, nullptr, num_vertices(partitioner.orig_g), recalc_occ_dt, recalc_sum_op, 0, MPI_COMM_WORLD);
		//}
		//MPI_Op_free(&recalc_sum_op);

		//if (procid == 0) {
			//int num_recvs = 0;
			//for (int i = 0; i < num_vertices(*graphs[0]); ++i) {
				////int rr_node_pid = partitioner.result_pid_by_level[current_level][i];
				////int from_pid = rr_node_pid == -1 ? 0 : rr_node_pid;

				////if (from_pid != procid) {
				//for (int from_pid = 1; from_pid < num_procs; ++from_pid) {
					//int recalc_occ;
					//MPI_Recv(&recalc_occ, 1, MPI_INT, from_pid, i, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
					//zlog_level(delta_log, ROUTER_V3, "%d Recvd %d recalc_occ %d\n", from_pid, i, recalc_occ);
					//congestion[i].recalc_occ += recalc_occ;

					//++num_recvs;
				//}
			//}
			//printf("Num recvs: %d\n", num_recvs);
		//} else {
			//int num_sends = 0;

			//for (int i = 0; i < num_vertices(*graphs[0]); ++i) {
				////int pid = partitioner.result_pid_by_level[current_level][i];
				////if (pid == procid) {
					//zlog_level(delta_log, ROUTER_V3, "Sending %d recalc_occ %d\n", i, congestion[i].recalc_occ);
					//MPI_Send(&congestion[i].recalc_occ, 1, MPI_INT, 0, i, MPI_COMM_WORLD);

					//++num_sends;
				////}
			//}

			//printf("[%d] Num sends: %d\n", procid, num_sends);
		//}

		bool valid = true;
		for (int i = 0; i < num_vertices(partitioner.orig_g); ++i) {
			sprintf_rr_node(i, buffer);
			if (congestion[i].cong.recalc_occ != congestion[i].cong.occ) {
				zlog_error(delta_log, "Node %s occ mismatch, recalc: %d original: %d\n", buffer, congestion[i].cong.recalc_occ, congestion[i].cong.occ);
				valid = false;
			}
		}
		assert(valid);

        int overused_total_bb_rank = 0;
        int num_congested_nets = 0;

		for (int i = 0; i < nets_to_route.size(); ++i) {
			net_t *net = nets_to_route[i].second;

			vector<int> overused_rr_node;
			assert(route_trees[net->local_id].root_rt_nodes.size() == 1);
			get_overused_nodes(route_trees[net->local_id], route_trees[net->local_id].root_rt_nodes[0], partitioner.orig_g, congestion, overused_rr_node);
			if (!overused_rr_node.empty()) {
				zlog_level(delta_log, ROUTER_V1, "Net %d bb_rank %d has %lu overused nodes:\n", net->vpr_id, net->bb_area_rank, overused_rr_node.size());
				for (const auto &item : overused_rr_node) {
					zlog_level(delta_log, ROUTER_V1, "%d ", item);
				}
				zlog_level(delta_log, ROUTER_V1, "\n");
				overused_total_bb_rank += net->bb_area_rank;
				++num_congested_nets;
			}
		}

        zlog_level(delta_log, ROUTER_V1, "Average overused net bb rank: %g\n", (float)overused_total_bb_rank/num_congested_nets);
        zlog_level(delta_log, ROUTER_V1, "Num congested nets: %d\n", num_congested_nets);

        iter_start = clock::now();

        auto analyze_timing_start = clock::now();

		crit_path_delay = analyze_timing(net_timing);

		//float *all_crits = nullptr;

		//if (procid == 0) {


			//all_crits = all_delays;

			//idx = 0;
			//for (int i = 0; i < num_procs; ++i) {
				//idx = 0;
				//for (int j = i; j < nets_to_route.size(); j += num_procs) {
					//net_t *net = nets_to_route[j].second;
					//zlog_level(delta_log, ROUTER_V3, "Net %d send crits:\n", net->vpr_id);
					//for (int k = 1; k <= net->sinks.size(); ++k) {
						//all_crits[displs[i] + idx] = net_timing[net->vpr_id].timing_criticality[k];
						//zlog_level(delta_log, ROUTER_V3, "\t%g\n", all_crits[displs[i] + idx]);
						//++idx;
					//}
				//}
				//assert(idx == recvcounts[i]);
			//}
		//} else {
			//all_crits = nullptr;
		//}

		//float *crits = delays;

		//MPI_Scatterv(all_crits, recvcounts, displs, MPI_FLOAT, crits, recvcounts[procid], MPI_FLOAT, 0, cur_comm);

		//idx = 0;
		//for (int i = procid; i < nets_to_route.size(); i += num_procs) {
			//net_t *net = nets_to_route[i].second;
			//zlog_level(delta_log, ROUTER_V3, "Net %d recv crits:\n", net->vpr_id);
			//for (int j = 1; j <= net->sinks.size(); ++j) {
				//net_timing[net->vpr_id].timing_criticality[j] = crits[idx];
				//zlog_level(delta_log, ROUTER_V3, "\t%g\n", net_timing[net->vpr_id].timing_criticality[j]);
				//++idx;
			//}
		//}

        analyze_timing_time = clock::now()-analyze_timing_start;

        if (feasible_routing(partitioner.orig_g, congestion) && !has_unroutable_sinks) {
            //dump_route(*current_traces_ptr, "route.txt");
			routed = true;	
        } else {
            unsigned long num_overused_nodes = 0;
			vector<int> overused_nodes_by_type(NUM_RR_TYPES, 0);
            for (int i = 0; i < num_vertices(partitioner.orig_g); ++i) {
                if (congestion[i].cong.occ > get_vertex_props(partitioner.orig_g, i).capacity) {
					const auto &v_p = get_vertex_props(partitioner.orig_g, i);
					++overused_nodes_by_type[v_p.type];

                    ++num_overused_nodes;
                }
            }

			static const char *name_type[] = { "SOURCE", "SINK", "IPIN", "OPIN",
				"CHANX", "CHANY", "INTRA_CLUSTER_EDGE" };
            zlog_info(delta_log, "Num overused nodes: %lu/%d (%.2f)\n", num_overused_nodes, num_vertices(partitioner.orig_g), num_overused_nodes*100.0/num_vertices(partitioner.orig_g));
			for (int i = 0; i < overused_nodes_by_type.size(); ++i) {
				zlog_info(delta_log, "\t%s: %d (%g)\n", name_type[i], overused_nodes_by_type[i], overused_nodes_by_type[i]*100.0/num_overused_nodes);
			}

            printf("Num overused nodes: %lu/%d (%.2f)\n", num_overused_nodes, num_vertices(partitioner.orig_g), num_overused_nodes*100.0/num_vertices(partitioner.orig_g));
			for (int i = 0; i < overused_nodes_by_type.size(); ++i) {
				printf("\t%s: %d (%g)\n", name_type[i], overused_nodes_by_type[i], overused_nodes_by_type[i]*100.0/num_overused_nodes);
			}

			if (num_overused_nodes > prev_num_overused_nodes && iter > 10) {
				++current_level;

				assert(current_level < partitioner.result_pid_by_level.size());

				printf("Transitioned to level %d at iteration %d\n", current_level, iter);
				zlog_level(delta_log, ROUTER_V1, "Transitioned to level %d at iteration %d\n", current_level, iter);

				for (int i = 0; i < nets_to_route.size(); ++i) {
					net_t *net = nets_to_route[i].second;

					unroutable_sinks[net->local_id].clear();

					//for (const auto &rs : routed_sinks[net->local_id]) {
					//assert(find(begin(fixed_sinks[net->local_id]), end(fixed_sinks[net->local_id]), rs) == end(fixed_sinks[net->local_id]));

					//zlog_level(delta_log, ROUTER_V3, "Fixing net %d sink %d\n", net->vpr_id, rs->id);

					//fixed_sinks[net->local_id].push_back(rs);
					//}
				}

				prev_num_overused_nodes = std::numeric_limits<unsigned long>::max();
				has_unroutable_sinks = false;
				num_procs /= 2;
			} else {
				prev_num_overused_nodes = num_overused_nodes;
			}

			//if ((float)total_num_interpartition_sinks/total_num_sinks_to_route > 0.1) {
				//++current_level;
				//assert(current_level < partitioner.result_pid_by_level.size());
			//}

            auto update_cost_start = clock::now();

            if (iter == 0) {
                params.pres_fac = opts->initial_pres_fac;
                //update_costs_mpi(partitioner.orig_g[0], partitioner.result_pid_by_level[current_level], procid, congestion, win, params.pres_fac, 0);
				update_costs(partitioner.orig_g, congestion, params.pres_fac, 0);
            } else {
                params.pres_fac *= opts->pres_fac_mult;

				/* Avoid overflow for high iteration counts, even if acc_cost is big */
                params.pres_fac = std::min(params.pres_fac, static_cast<float>(HUGE_POSITIVE_FLOAT / 1e5));

                //update_costs_mpi(partitioner.orig_g[0], partitioner.result_pid_by_level[current_level], procid, congestion, win, params.pres_fac, opts->acc_fac);
				update_costs(partitioner.orig_g, congestion, params.pres_fac, opts->acc_fac);
            }

            update_cost_time = clock::now()-update_cost_start;
        }

        iter_time += clock::now()-iter_start;

        //total_route_time += route_time;
        total_greedy_route_time += greedy_route_time;
        total_update_cost_time += update_cost_time;
        total_analyze_timing_time += analyze_timing_time;
        total_iter_time += iter_time;

		printf("Iteration time: %g s.\n", duration_cast<nanoseconds>(iter_time).count() / 1e9);
		printf("\tRoute time: %g s.\n", duration_cast<nanoseconds>(greedy_route_time).count() / 1e9);
		printf("\tUpdate cost time: %g s.\n", duration_cast<nanoseconds>(update_cost_time).count() / 1e9);
		printf("\tAnalyze timing time: %g s.\n", duration_cast<nanoseconds>(analyze_timing_time).count() / 1e9);
		printf("Critical path: %g ns\n", crit_path_delay);

        //clock::time_point greedy_earliest_end_time = *std::min_element(begin(greedy_end_time), begin(greedy_end_time)+num_procs);

        //printf("greedy wait time: ");
        //for (int i = 0; i < num_procs; ++i) {
			//printf("%g (%g), ", duration_cast<nanoseconds>(greedy_end_time[i]-greedy_earliest_end_time).count() / 1e9, 100.0*(greedy_end_time[i]-greedy_earliest_end_time)/greedy_route_time);
        //}
		//printf("\n");
    }

    //sprintf(buffer, "%s_run_%d_rr_graph_occ.txt", s_circuit_name, run);
    //dump_rr_graph_occ(congestion, num_vertices(g), buffer);
	if (routed) {
		printf("Routed in %d iterations. Total iteration time: %g\n", iter, duration_cast<nanoseconds>(total_iter_time).count() / 1e9);
		printf("\tTotal route time: %g s.\n", duration_cast<nanoseconds>(total_greedy_route_time).count() / 1e9);
		printf("\t\tTotal greedy route time: %g s.\n", duration_cast<nanoseconds>(total_greedy_route_time).count() / 1e9);
		printf("\tTotal update cost time: %g s.\n", duration_cast<nanoseconds>(total_update_cost_time).count() / 1e9);
		printf("\tTotal analyze timing time: %g s.\n", duration_cast<nanoseconds>(total_analyze_timing_time).count() / 1e9);

		printf("Final critical path: %g ns\n", crit_path_delay);
	} else {
		printf("Failed to route in %d iterations. Total iteration time: %g\n", opts->max_router_iterations, duration_cast<nanoseconds>(total_iter_time).count() / 1e9);
		printf("\tTotal route time: %g s.\n", duration_cast<nanoseconds>(total_greedy_route_time).count() / 1e9);
		printf("\t\tTotal greedy route time: %g s.\n", duration_cast<nanoseconds>(total_greedy_route_time).count() / 1e9);
		printf("\tTotal update cost time: %g s.\n", duration_cast<nanoseconds>(total_update_cost_time).count() / 1e9);
		printf("\tTotal analyze timing time: %g s.\n", duration_cast<nanoseconds>(total_analyze_timing_time).count() / 1e9);
	}

    //delete_graph(g);
    delete_net_timing(nets, global_nets, net_timing);	
    delete [] congestion;
    delete [] net_timing;

    return routed;
}

bool mpi_spatial_route_partitioned(t_router_opts *opts, struct s_det_routing_arch det_routing_arch, t_direct_inf *directs, int num_directs, t_segment_inf *segment_inf, t_timing_inf timing_inf)
{
    using std::chrono::duration_cast;
    using std::chrono::nanoseconds;
    using clock = std::chrono::high_resolution_clock;

	int initial_num_procs, num_procs, initial_procid, procid;
	MPI_Comm cur_comm = MPI_COMM_WORLD;

    MPI_Comm_size(cur_comm, &initial_num_procs);
    MPI_Comm_rank(cur_comm, &initial_procid);

	num_procs = initial_num_procs;
	procid = initial_procid;

	init_congestion_mpi_datatype();
	init_datatypes();

    init_logging();
    zlog_set_record("custom_output", delta_log_output);
    zlog_set_record("missing_edge", missing_edge_log_output);
    zlog_set_record("ss", ss_log_output);
    delta_log_files.resize(opts->max_router_iterations);
    for (int i = 0; i < opts->max_router_iterations; ++i) {
        delta_log_files[i].resize(num_procs, nullptr);
    }
    missing_edge_log_files.resize(opts->max_router_iterations);
    for (int i = 0; i < opts->max_router_iterations; ++i) {
        missing_edge_log_files[i].resize(num_procs, nullptr);
    }

    //fclose(fopen(LOG_PATH_PREFIX"iter__tid_.log", "w"));

    test_fast_graph();
    test_topo();
    test_fm();
    test_filter_graph();
    test_partition_graph();
    test_connected_components();

    vector<RRGraph *> graphs;
    rr_graph_partitioner partitioner;
    partitioner.init(opts, &det_routing_arch, directs, num_directs, segment_inf, &timing_inf);
    partitioner.partition(num_procs, graphs);
    //for (const auto &g : graphs) {
        //routability(*g);
    //}
    //
    //RRGraph combined_g;
    //add_vertex(combined_g, num_vertices(partitioner.orig_g));

    //for (const auto &g : graphs) {
        //for (const auto &e : get_edges(*g)) {
            //int from = get_source(*g, e);
            //int to = get_target(*g, e);
            //const auto &from_ver = get_vertex_props(*g, from);
            //const auto &to_ver = get_vertex_props(*g, to);

            //if (is_channel(from_ver) && is_channel(to_ver)) {
                //assert(!has_edge(combined_g, from, to));
                //add_edge(combined_g, from, to);
            //}
        //}
    //}

    //for (const auto &e : get_edges(partitioner.orig_g)) {
        //int from = get_source(partitioner.orig_g, e);
        //int to = get_target(partitioner.orig_g, e);
        //const auto &from_ver = get_vertex_props(partitioner.orig_g, from);
        //const auto &to_ver = get_vertex_props(partitioner.orig_g, to);
        //if (!is_channel(from_ver) || !is_channel(to_ver)) {
            //assert(!has_edge(combined_g, from, to));
            //add_edge(combined_g, from, to);
        //}
    //}

    //printf("Combined/Orig graph has %d/%d (%g) edges.\n", num_edges(combined_g), num_edges(partitioner.orig_g), 100.0*num_edges(combined_g)/num_edges(partitioner.orig_g));

    //goto lol;

    //extern int num_types;
    //extern struct s_type_descriptor *type_descriptors;
    //extern int nx, ny;
    //extern struct s_grid_tile **grid;

    //free_rr_graph();

    //int warnings;

    //build_rr_graph(GRAPH_UNIDIR, num_types, type_descriptors, nx, ny, grid,
            //opts->fixed_channel_width, NULL, det_routing_arch.switch_block_type,
            //det_routing_arch.Fs, det_routing_arch.num_segment,
            //det_routing_arch.num_switch, segment_inf,
            //det_routing_arch.global_route_switch,
            //det_routing_arch.delayless_switch, timing_inf,
            //det_routing_arch.wire_to_ipin_switch, opts->base_cost_type,
            //directs, num_directs, FALSE,
            //&warnings);

    //RRGraph channel_with_interior_g;
    //init_channel_only_graph(channel_with_interior_g);

    //dump_rr_graph(channel_with_interior_g, "/Volumes/DATA/rr_graph_channel_with_interior.txt");

    //RRGraph orig_g;
    //init_graph(orig_g);

    //dump_rr_graph(orig_g, "/Volumes/DATA/rr_graph.txt");

    //free_rr_graph();
    //for (int i = 0; i < det_routing_arch.num_segment; ++i) {
        //for (int j = 0; j < segment_inf[i].sb_len; ++j) {
            //if (j != 0 && j != segment_inf[i].sb_len-1) {
                //segment_inf[i].sb[j] = FALSE;
            //}
        //}
    //}
    //build_rr_graph(GRAPH_UNIDIR, num_types, type_descriptors, nx, ny, grid,
            //opts->fixed_channel_width, NULL, det_routing_arch.switch_block_type,
            //det_routing_arch.Fs, det_routing_arch.num_segment,
            //det_routing_arch.num_switch, segment_inf,
            //det_routing_arch.global_route_switch,
            //det_routing_arch.delayless_switch, timing_inf,
            //det_routing_arch.wire_to_ipin_switch, opts->base_cost_type,
            //directs, num_directs, FALSE,
            //&warnings);

    //RRGraph channel_without_interior_g;
    //init_channel_only_graph(channel_without_interior_g);

    //dump_rr_graph(channel_without_interior_g, "/Volumes/DATA/rr_graph_channel_without_interior.txt");

    //vector<vector<vector<int>>> all_partition_components(opts->num_threads);
    ////for (int x = 0; x < nx+1; ++x) {
        ////for (int y = 0; y < ny+1; ++y) {
            ////init_partitioned_graph_5(opts->num_threads, channel_with_interior_g, channel_without_interior_g, orig_g, graphs, x, y);
        ////}
    ////}
    ////init_partitioned_graph_4(opts->num_threads, channel_with_interior_g, channel_without_interior_g, orig_g, graphs, 14, 14);
    //init_partitioned_graph_4(opts->num_threads, channel_with_interior_g, channel_without_interior_g, orig_g, graphs, all_partition_components);

    //for (int i = 0; i < graphs.size(); ++i) {
        //char filename[256];

        //sprintf(filename, "/Volumes/DATA/rr_graph_%d.txt", i);
        //dump_rr_graph(*graphs[i], filename);

        //sprintf(filename, "/Volumes/DATA/rr_graph_%d_edges.txt", i);
        //dump_edges(*graphs[i], filename);
    //}

    vector<net_t> nets;
    vector<net_t> global_nets;
    //init_nets(nets, global_nets, opts->bb_factor, partitioner.sink_in_nodes);	
    init_nets(nets, global_nets, opts->bb_factor, opts->large_bb);	

    vector<net_t *> nets_ptr(nets.size());
    for (int i = 0; i < nets.size(); ++i) {
        nets_ptr[i] = &nets[i];
    }
    extern s_bb *route_bb;
    sort(nets_ptr.begin(), nets_ptr.end(), [] (const net_t *a, const net_t *b) -> bool {
            return get_bounding_box_area(route_bb[a->vpr_id]) > get_bounding_box_area(route_bb[b->vpr_id]);
            });
    int rank = 0;
    for (auto &net : nets_ptr) {
        net->bb_area_rank = rank++;
    }

	route_state_t *states;
	congestion_t *congestion;
	vector<route_tree_t> route_trees;
	t_net_timing *net_timing;
	init_route_structs(partitioner.orig_g, nets, global_nets, &states, &congestion, route_trees, &net_timing);

    route_parameters_t params;
    params.criticality_exp = opts->criticality_exp;
    params.astar_fac = opts->astar_fac;
    params.max_criticality = opts->max_criticality;
    params.bend_cost = opts->bend_cost;
	params.pres_fac = opts->first_iter_pres_fac; /* Typically 0 -> ignore cong. */

    char buffer[256];

    vector<pair<box, net_t *>> nets_to_route;
    
    for (auto &net : nets) {
        box b = bg::make_inverse<box>();

        bg::expand(b, point(net.source.x, net.source.y));
        for (const auto &sink : net.sinks) {
            bg::expand(b, point(sink.x, sink.y));
        }
        bg::subtract_value(b.min_corner(), 1+opts->bb_factor);
        bg::add_value(b.max_corner(), opts->bb_factor);

        nets_to_route.push_back(make_pair(b, &net));
    }
    std::sort(begin(nets_to_route), end(nets_to_route), [] (const pair<box, net_t *> &a, const pair<box, net_t *> &b) -> bool {
            return a.second->sinks.size() > b.second->sinks.size();
            });

    bool routed = false;

    int iter;
    float crit_path_delay;
	int current_level = 0;

	vector<vector<sink_t *>> unroutable_sinks(nets.size());
	//vector<vector<sink_t *>> fixed_sinks(nets.size());
	vector<vector<sink_t *>> routed_sinks(nets.size());
	bool has_unroutable_sinks = false;
	unsigned long prev_num_overused_nodes = std::numeric_limits<unsigned long>::max();

    clock::duration total_greedy_route_time = clock::duration::zero();
    clock::duration total_update_cost_time = clock::duration::zero();
    clock::duration total_analyze_timing_time = clock::duration::zero();
    clock::duration total_iter_time = clock::duration::zero();

    for (iter = 0; iter < opts->max_router_iterations && !routed; ++iter) {
        clock::duration greedy_route_time = clock::duration::zero();
        clock::duration update_cost_time = clock::duration::zero();
        clock::duration partitioning_time = clock::duration::zero();
        clock::duration analyze_timing_time = clock::duration::zero();
        clock::duration iter_time = clock::duration::zero();

        //for (int i = 0; i < opts->num_threads; ++i) {
            //sprintf(buffer, LOG_PATH_PREFIX"iter_%d_tid_%d.log", iter, i);
            //fclose(fopen(buffer, "w"));
        //}

		MPI_Barrier(cur_comm);
        
        auto iter_start = clock::now();

        //auto route_start = clock::now();

        //tbb::enumerable_thread_specific<state_t *> state_tls;
        //
		vector<perf_t> perfs(num_procs);
        vector<int> thread_num_nets_routed(num_procs, 0);
        vector<int> thread_num_nets_to_route(num_procs, 0);
        vector<int> thread_num_sinks_routed(num_procs, 0);
        vector<int> thread_num_sinks_to_route(num_procs, 0);
		vector<int> thread_num_interpartition_sinks(num_procs, 0);
        bool has_interpartition_sinks = false;
        vector<vector<interpartition_sink_t>> interpartition_sinks(nets.size()); 
		//vector<vector<RRNode>> net_sinks(nets.size());

		sprintf(buffer, "%d", iter);
		zlog_put_mdc("iter", buffer);

		sprintf(buffer, "%d", procid);
		zlog_put_mdc("tid", buffer);

        zlog_info(delta_log, "Routing iteration: %d\n", iter);
        printf("Routing iteration: %d\n", iter);

        auto greedy_route_start = clock::now();

		perfs[procid].num_heap_pushes = 0;
		perfs[procid].num_heap_pops = 0;
		perfs[procid].num_neighbor_visits = 0;

		for (int i = procid*pow(2, current_level); i < nets_to_route.size(); i += initial_num_procs) {
			for (int j = 0; j < pow(2, current_level) && i+j < nets_to_route.size(); ++j) {
				zlog_level(delta_log, ROUTER_V1, "Net index: %d\n", i+j);

				net_t *net = nets_to_route[i+j].second;

				update_sink_criticalities(*net, net_timing[net->vpr_id], params);

				//printf("Routing net %d\n", net->vpr_id);
				//vector<sink_t *> temp_routed_sinks = routed_sinks[net->local_id];

				//vector<RRNode> sinks_to_mark;
				//for (const auto &sink : temp_routed_sinks) {
					//bool fixed = find(begin(fixed_sinks[net->local_id]), end(fixed_sinks[net->local_id]), sink) != end(fixed_sinks[net->local_id]);
					//if (!fixed) {
						//sinks_to_mark.push_back(sink->rr_node);
					//} else {
						//routed_sinks[net->local_id].push_back(sink);
					//}
				//}

				//route_tree_mark_paths_to_be_ripped(route_trees[net->local_id], partitioner.orig_g, partitioner.result_pid_by_level[current_level], procid, sinks_to_mark);
				route_tree_mark_congested_nodes_to_be_ripped(route_trees[net->local_id], partitioner.orig_g, congestion);

				route_tree_rip_up_marked(route_trees[net->local_id], partitioner.orig_g, congestion, params.pres_fac);

				routed_sinks[net->local_id].clear();
				for (auto &sink : net->sinks) {
					if (route_tree_get_rt_node(route_trees[net->local_id], sink.rr_node) != RouteTree::null_vertex()) {
						routed_sinks[net->local_id].push_back(&sink);
					}
				}

				vector<sink_t *> sinks;	
				get_sinks_to_route(net, route_trees[net->local_id], unroutable_sinks[net->local_id], sinks);

				//local_perf.total_rip_up_time += clock::now()-rip_up_start;

				//auto route_start = clock::now();

				if (!sinks.empty()) {
					int previous_num_unroutable_sinks = unroutable_sinks[net->local_id].size();

					int previous_num_routed_sinks = routed_sinks[net->local_id].size();

					//route_net_mpi_send_recv(partitioner.orig_g, partitioner.result_pid_by_level[current_level], procid, net->vpr_id, &net->source, sinks, params, states, congestion, num_procs, cur_comm, transactions, route_trees[net->local_id], net_timing[net->vpr_id], routed_sinks[net->local_id], unroutable_sinks[net->local_id], &perfs[procid]);

					assert(routed_sinks[net->local_id].size() + unroutable_sinks[net->local_id].size() == net->sinks.size());

					if (!unroutable_sinks[net->local_id].empty() && previous_num_unroutable_sinks > 0) {
						assert(previous_num_unroutable_sinks == unroutable_sinks[net->local_id].size());
					}

					if (!has_unroutable_sinks) {
						has_unroutable_sinks = !unroutable_sinks[net->local_id].empty();
					}

					++thread_num_nets_routed[procid];
					++thread_num_nets_to_route[procid];

					thread_num_sinks_to_route[procid] += sinks.size();
					thread_num_sinks_routed[procid] += routed_sinks[net->local_id].size() - previous_num_routed_sinks;
				} else {
					assert(routed_sinks[net->local_id].size() + unroutable_sinks[net->local_id].size() == net->sinks.size());
					zlog_level(delta_log, ROUTER_V2, "Net %d has no sinks in the current iteration because there are %lu/%lu non-routable/all sinks\n", net->vpr_id, unroutable_sinks[net->local_id].size(), net->sinks.size());
				}

				//local_perf.total_route_time += clock::now()-rip_up_start;
			}
		}

		//greedy_end_time = clock::now();
		MPI_Barrier(cur_comm);

        greedy_route_time = clock::now()-greedy_route_start;

		//if (has_interpartition_sinks) {
			//for (const auto &net_interpartition_sinks : interpartition_sinks) {
				//for (const auto &is : net_interpartition_sinks) {
					//vector<RRNode> added_nodes;
					//for (const auto &node : is.path) {
						//if (node.update_cost) {
							//added_nodes.push_back(node.rr_node_id);
						//}
					//}
					//update_one_cost(partitioner.orig_g, congestion, added_nodes.begin(), added_nodes.end(), 1, params.pres_fac, false, nullptr);
				//}
			//}
		//}

        //if (greedy_rip_up_all) {
            //next_greedy_rip_up_iter += greedy_rip_up_all_period;
            //++greedy_rip_up_all_period;
            //prev_num_overused_nodes = std::numeric_limits<unsigned long>::max();
        //}

        //route_time = clock::now()-route_start;

        iter_time = clock::now()-iter_start;


		int total_num_sinks_to_route;
		if (procid == 0) {
			int total_num_nets_to_route = thread_num_nets_to_route[0];
			int total_num_nets_routed = thread_num_nets_routed[0];
			int total_num_sinks_routed = thread_num_sinks_routed[0];
			total_num_sinks_to_route = thread_num_sinks_to_route[0];

			for (int i = 1; i < num_procs; ++i) {
				MPI_Recv(&thread_num_nets_to_route[i], 1, MPI_INT, i, i, cur_comm, MPI_STATUS_IGNORE);
				MPI_Recv(&thread_num_nets_routed[i], 1, MPI_INT, i, i, cur_comm, MPI_STATUS_IGNORE);
				MPI_Recv(&thread_num_sinks_to_route[i], 1, MPI_INT, i, i, cur_comm, MPI_STATUS_IGNORE);
				MPI_Recv(&thread_num_sinks_routed[i], 1, MPI_INT, i, i, cur_comm, MPI_STATUS_IGNORE);

				total_num_nets_to_route += thread_num_nets_to_route[i];
				total_num_nets_routed += thread_num_nets_routed[i];
				total_num_sinks_routed += thread_num_sinks_routed[i];
				total_num_sinks_to_route += thread_num_sinks_to_route[i];
			}
			assert(total_num_nets_to_route == total_num_nets_routed);

			printf("num nets routed: ");
			for (int i = 0; i < num_procs; ++i) {
				printf("%d/%d (%g), ", thread_num_nets_routed[i], thread_num_nets_to_route[i], thread_num_nets_routed[i]*100.0/thread_num_nets_to_route[i]);
			}
			printf("\n");

			printf("num sinks routed: ");
			for (int i = 0; i < num_procs; ++i) {
				printf("%d/%d (%g), ", thread_num_sinks_routed[i], thread_num_sinks_to_route[i], thread_num_sinks_routed[i]*100.0/thread_num_sinks_to_route[i]);
			}
			printf("\n");

			printf("Total num nets routed: %d/%d (%g)\n", total_num_nets_routed, total_num_nets_to_route, total_num_nets_routed*100.0/total_num_nets_to_route);
			printf("Total num sinks routed: %d/%d (%g)\n", total_num_sinks_routed, total_num_sinks_to_route, (total_num_sinks_routed)*100.0/total_num_sinks_to_route);

			for (int i = 1; i < num_procs; ++i) {
				MPI_Recv(&perfs[i].num_heap_pushes, 1, MPI_INT, i, i, cur_comm, MPI_STATUS_IGNORE);
				MPI_Recv(&perfs[i].num_heap_pops, 1, MPI_INT, i, i, cur_comm, MPI_STATUS_IGNORE);
				MPI_Recv(&perfs[i].num_neighbor_visits, 1, MPI_INT, i, i, cur_comm, MPI_STATUS_IGNORE);
			}

			printf("num_heap_pushes: ");
			for (int i = 0; i < num_procs; ++i) {
				printf("%lu ", perfs[i].num_heap_pushes);
			}
			printf("\n");

			printf("num_heap_pops: ");
			for (int i = 0; i < num_procs; ++i) {
				printf("%lu ", perfs[i].num_heap_pops);
			}
			printf("\n");

			printf("num_neighbor_visits: ");
			for (int i = 0; i < num_procs; ++i) {
				printf("%lu ", perfs[i].num_neighbor_visits);
			}
			printf("\n");

			unsigned long total_num_heap_pushes = 0;
			unsigned long total_num_heap_pops = 0;
			unsigned long total_num_neighbor_visits = 0;

			for (int i = 0; i < num_procs; ++i) {
				total_num_heap_pushes += perfs[i].num_heap_pushes;
				total_num_heap_pops += perfs[i].num_heap_pops;
				total_num_neighbor_visits += perfs[i].num_neighbor_visits;
			}

			printf("total_num_heap_pushes: %lu\n", total_num_heap_pushes);
			printf("total_num_heap_pops: %lu\n", total_num_heap_pops); 
			printf("total_num_neighbor_visits: %lu\n", total_num_neighbor_visits);
		} else {
			MPI_Send(&thread_num_nets_to_route[procid], 1, MPI_INT, 0, procid, cur_comm);
			MPI_Send(&thread_num_nets_routed[procid], 1, MPI_INT, 0, procid, cur_comm);
			MPI_Send(&thread_num_sinks_to_route[procid], 1, MPI_INT, 0, procid, cur_comm);
			MPI_Send(&thread_num_sinks_routed[procid], 1, MPI_INT, 0, procid, cur_comm);

			MPI_Send(&perfs[procid].num_heap_pushes, 1, MPI_INT, 0, procid, cur_comm);
			MPI_Send(&perfs[procid].num_heap_pops, 1, MPI_INT, 0, procid, cur_comm);
			MPI_Send(&perfs[procid].num_neighbor_visits, 1, MPI_INT, 0, procid, cur_comm);
		}

		/* checking */
        for (int i = 0; i < num_vertices(partitioner.orig_g); ++i) {
            congestion[i].recalc_occ = 0; 
        }

		for (int i = procid*pow(2, current_level); i < nets_to_route.size(); i += initial_num_procs) {
			for (int j = 0; j < pow(2, current_level) && i+j < nets_to_route.size(); ++j) {
				net_t *net = nets_to_route[i+j].second;

				check_route_tree(route_trees[net->local_id], *net, routed_sinks[net->local_id], partitioner.orig_g);
				recalculate_occ(route_trees[net->local_id], partitioner.orig_g, congestion);
			}
        }

		//int *recalc_occ = new int[num_vertices(partitioner.orig_g)];
		//for (int i = 0; i < num_vertices(partitioner.orig_g); ++i) {
			//recalc_occ[i] = congestion[i].recalc_occ;	
		//}
		//if (procid == 0) {
			//int *recv_recalc_occ = new int[num_vertices(partitioner.orig_g)];
			//MPI_Reduce(recalc_occ, recv_recalc_occ, num_vertices(partitioner.orig_g), MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);
			//for (int i = 0; i < num_vertices(partitioner.orig_g); ++i) {
				//congestion[i].recalc_occ = recv_recalc_occ[i];
			//}
			//delete [] recv_recalc_occ;
		//} else {
			//MPI_Reduce(recalc_occ, nullptr, num_vertices(partitioner.orig_g), MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);
		//}
		//delete [] recalc_occ;
		
		//MPI_Op recalc_sum_op;
		//MPI_Op_create((MPI_User_function *)recalc_sum, 1, &recalc_sum_op);	
		//if (procid == 0) {
			//MPI_Reduce(MPI_IN_PLACE, congestion, num_vertices(partitioner.orig_g), recalc_occ_dt, recalc_sum_op, 0, MPI_COMM_WORLD);
		//} else {
			//MPI_Reduce(congestion, nullptr, num_vertices(partitioner.orig_g), recalc_occ_dt, recalc_sum_op, 0, MPI_COMM_WORLD);
		//}
		//MPI_Op_free(&recalc_sum_op);

		//if (procid == 0) {
			//int num_recvs = 0;
			//for (int i = 0; i < num_vertices(*graphs[0]); ++i) {
				////int rr_node_pid = partitioner.result_pid_by_level[current_level][i];
				////int from_pid = rr_node_pid == -1 ? 0 : rr_node_pid;

				////if (from_pid != procid) {
				//for (int from_pid = 1; from_pid < num_procs; ++from_pid) {
					//int recalc_occ;
					//MPI_Recv(&recalc_occ, 1, MPI_INT, from_pid, i, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
					//zlog_level(delta_log, ROUTER_V3, "%d Recvd %d recalc_occ %d\n", from_pid, i, recalc_occ);
					//congestion[i].recalc_occ += recalc_occ;

					//++num_recvs;
				//}
			//}
			//printf("Num recvs: %d\n", num_recvs);
		//} else {
			//int num_sends = 0;

			//for (int i = 0; i < num_vertices(*graphs[0]); ++i) {
				////int pid = partitioner.result_pid_by_level[current_level][i];
				////if (pid == procid) {
					//zlog_level(delta_log, ROUTER_V3, "Sending %d recalc_occ %d\n", i, congestion[i].recalc_occ);
					//MPI_Send(&congestion[i].recalc_occ, 1, MPI_INT, 0, i, MPI_COMM_WORLD);

					//++num_sends;
				////}
			//}

			//printf("[%d] Num sends: %d\n", procid, num_sends);
		//}

		bool valid = true;
		for (int i = 0; i < num_vertices(partitioner.orig_g); ++i) {
			sprintf_rr_node(i, buffer);
			if (congestion[i].recalc_occ != congestion[i].occ) {
				zlog_error(delta_log, "Node %s occ mismatch, recalc: %d original: %d\n", buffer, congestion[i].recalc_occ, congestion[i].occ);
				valid = false;
			}
		}
		assert(valid);

        int overused_total_bb_rank = 0;
        int num_congested_nets = 0;

		for (int i = procid*pow(2, current_level); i < nets_to_route.size(); i += initial_num_procs) {
			for (int j = 0; j < pow(2, current_level) && i+j < nets_to_route.size(); ++j) {
				net_t *net = nets_to_route[i+j].second;

				vector<int> overused_rr_node;
				assert(route_trees[net->local_id].root_rt_nodes.size() == 1);
				get_overused_nodes(route_trees[net->local_id], route_trees[net->local_id].root_rt_nodes[0], partitioner.orig_g, congestion, overused_rr_node);
				if (!overused_rr_node.empty()) {
					zlog_level(delta_log, ROUTER_V1, "Net %d bb_rank %d has %lu overused nodes:\n", net->vpr_id, net->bb_area_rank, overused_rr_node.size());
					for (const auto &item : overused_rr_node) {
						zlog_level(delta_log, ROUTER_V1, "%d ", item);
					}
					zlog_level(delta_log, ROUTER_V1, "\n");
					overused_total_bb_rank += net->bb_area_rank;
					++num_congested_nets;
				}
			}
		}

        zlog_level(delta_log, ROUTER_V1, "Average overused net bb rank: %g\n", (float)overused_total_bb_rank/num_congested_nets);
        zlog_level(delta_log, ROUTER_V1, "Num congested nets: %d\n", num_congested_nets);

        iter_start = clock::now();

        auto analyze_timing_start = clock::now();

		//sync_net_delay(nets_to_route, procid, num_procs, initial_num_procs, current_level, cur_comm, net_timing);

		crit_path_delay = analyze_timing(net_timing);

		//float *all_crits = nullptr;

		//if (procid == 0) {


			//all_crits = all_delays;

			//idx = 0;
			//for (int i = 0; i < num_procs; ++i) {
				//idx = 0;
				//for (int j = i; j < nets_to_route.size(); j += num_procs) {
					//net_t *net = nets_to_route[j].second;
					//zlog_level(delta_log, ROUTER_V3, "Net %d send crits:\n", net->vpr_id);
					//for (int k = 1; k <= net->sinks.size(); ++k) {
						//all_crits[displs[i] + idx] = net_timing[net->vpr_id].timing_criticality[k];
						//zlog_level(delta_log, ROUTER_V3, "\t%g\n", all_crits[displs[i] + idx]);
						//++idx;
					//}
				//}
				//assert(idx == recvcounts[i]);
			//}
		//} else {
			//all_crits = nullptr;
		//}

		//float *crits = delays;

		//MPI_Scatterv(all_crits, recvcounts, displs, MPI_FLOAT, crits, recvcounts[procid], MPI_FLOAT, 0, cur_comm);

		//idx = 0;
		//for (int i = procid; i < nets_to_route.size(); i += num_procs) {
			//net_t *net = nets_to_route[i].second;
			//zlog_level(delta_log, ROUTER_V3, "Net %d recv crits:\n", net->vpr_id);
			//for (int j = 1; j <= net->sinks.size(); ++j) {
				//net_timing[net->vpr_id].timing_criticality[j] = crits[idx];
				//zlog_level(delta_log, ROUTER_V3, "\t%g\n", net_timing[net->vpr_id].timing_criticality[j]);
				//++idx;
			//}
		//}

        analyze_timing_time = clock::now()-analyze_timing_start;


		int m_routed = (feasible_routing(partitioner.orig_g, congestion) && !has_unroutable_sinks) ? 1 : 0;
		int reduced_routed;
		MPI_Allreduce(&m_routed, &reduced_routed, 1, MPI_INT, MPI_LAND, cur_comm);

        zlog_level(delta_log, ROUTER_V1, "m_routed: %d reduced_routed: %d\n", m_routed, reduced_routed);

        if (reduced_routed) {
            //dump_route(*current_traces_ptr, "route.txt");
			routed = true;	
        } else {
            unsigned long num_overused_nodes = 0;
			vector<int> overused_nodes_by_type(NUM_RR_TYPES, 0);
            for (int i = 0; i < num_vertices(partitioner.orig_g); ++i) {
                if (congestion[i].occ > get_vertex_props(partitioner.orig_g, i).capacity) {
					const auto &v_p = get_vertex_props(partitioner.orig_g, i);
					++overused_nodes_by_type[v_p.type];

                    ++num_overused_nodes;
                }
            }

			static const char *name_type[] = { "SOURCE", "SINK", "IPIN", "OPIN",
				"CHANX", "CHANY", "INTRA_CLUSTER_EDGE" };
            zlog_info(delta_log, "Num overused nodes: %lu/%d (%.2f)\n", num_overused_nodes, num_vertices(partitioner.orig_g), num_overused_nodes*100.0/num_vertices(partitioner.orig_g));
			for (int i = 0; i < overused_nodes_by_type.size(); ++i) {
				zlog_info(delta_log, "\t%s: %d (%g)\n", name_type[i], overused_nodes_by_type[i], overused_nodes_by_type[i]*100.0/num_overused_nodes);
			}

			int *all_overused_nodes_by_type = new int[num_procs*NUM_RR_TYPES];
			int *overused_nodes_by_type_send = new int[NUM_RR_TYPES];
			for (int i = 0; i < overused_nodes_by_type.size(); ++i) {
				overused_nodes_by_type_send[i] = overused_nodes_by_type[i];
			}

			MPI_Gather(overused_nodes_by_type_send, NUM_RR_TYPES, MPI_INT, all_overused_nodes_by_type, NUM_RR_TYPES, MPI_INT, 0, cur_comm);

			unsigned long *all_num_overused_nodes = new unsigned long[num_procs];
			MPI_Gather(&num_overused_nodes, 1, MPI_UNSIGNED_LONG, all_num_overused_nodes, 1, MPI_UNSIGNED_LONG, 0, cur_comm);

			if (procid == 0) {
				printf("Num overused nodes: ");
				for (int i = 0; i < num_procs; ++i) {
					printf("%lu/%d (%.2f) ", all_num_overused_nodes[i], num_vertices(partitioner.orig_g), all_num_overused_nodes[i]*100.0/num_vertices(partitioner.orig_g));
				}
				printf("\n");

				for (int i = 0; i < overused_nodes_by_type.size(); ++i) {
					printf("\t%s: ", name_type[i]);
					for (int j = 0; j < num_procs; ++j) {
						printf("%d (%g) ", all_overused_nodes_by_type[j*NUM_RR_TYPES+i], all_overused_nodes_by_type[j*NUM_RR_TYPES+i]*100.0/all_num_overused_nodes[j]);
					}
					printf("\n");
				}
			} 

			delete [] all_overused_nodes_by_type;
			delete [] overused_nodes_by_type_send;
			delete [] all_num_overused_nodes;

			int not_decreasing = (num_overused_nodes > prev_num_overused_nodes && iter > 10) ? 1 : 0;
			//int not_decreasing = current_level+1 < partitioner.result_pid_by_level.size(); [> testing <]
			int reduced_not_decreasing;
			MPI_Allreduce(&not_decreasing, &reduced_not_decreasing, 1, MPI_INT, MPI_LOR, cur_comm);

			prev_num_overused_nodes = num_overused_nodes;

			zlog_level(delta_log, ROUTER_V1, "not_decreasing: %d reduced_not_decreasing: %d\n", not_decreasing, reduced_not_decreasing);

			if (reduced_not_decreasing) {
				/* need to send route tree over */
				if (procid % 2 == 0) {
					/*receiver*/
					for (int i = (procid+1)*pow(2, current_level); i < nets_to_route.size(); i += initial_num_procs) {
						for (int j = 0; j < pow(2, current_level) && i+j < nets_to_route.size(); ++j) {
							net_t *net = nets_to_route[i+j].second;

							zlog_level(delta_log, ROUTER_V3, "Recving net index %d from %d\n", i+j, procid+1);

							recv_route_tree(net, partitioner.orig_g, routed_sinks, states, congestion, params.pres_fac, route_trees, net_timing, procid+1, cur_comm);
						}
					}
				} else {
					/*sender*/
					for (int i = procid*pow(2, current_level); i < nets_to_route.size(); i += initial_num_procs) {
						for (int j = 0; j < pow(2, current_level) && i+j < nets_to_route.size(); ++j) {
							net_t *net = nets_to_route[i+j].second;

							zlog_level(delta_log, ROUTER_V3, "Sending net index %d from %d\n", i+j, procid-1);

							send_route_tree(net, partitioner.orig_g, routed_sinks, route_trees, procid-1, cur_comm);
						}
					}
					
				}

				assert(num_procs % 2 == 0);
				num_procs /= 2;
				MPI_Comm new_comm;
				MPI_Comm_split(cur_comm, procid%2, procid, &new_comm);

				if (procid % 2 == 1) {
					/* early exit code */
					break;
				}

				cur_comm = new_comm;
				MPI_Comm_rank(cur_comm, &procid);

				++current_level;

				assert(current_level < partitioner.result_pid_by_level.size());

				printf("[%d] Transitioned to level %d at iteration %d\n", procid, current_level, iter);
				zlog_level(delta_log, ROUTER_V1, "Transitioned to level %d at iteration %d\n", current_level, iter);
				zlog_level(delta_log, ROUTER_V1, "New pid %d for initial pid %d\n", procid, initial_procid);

				for (int i = procid*pow(2, current_level); i < nets_to_route.size(); i += initial_num_procs) {
					for (int j = 0; j < pow(2, current_level) && i+j < nets_to_route.size(); ++j) {
						net_t *net = nets_to_route[i+j].second;

						unroutable_sinks[net->local_id].clear();

						//for (const auto &rs : routed_sinks[net->local_id]) {
							//assert(find(begin(fixed_sinks[net->local_id]), end(fixed_sinks[net->local_id]), rs) == end(fixed_sinks[net->local_id]));

							//zlog_level(delta_log, ROUTER_V3, "Fixing net %d sink %d\n", net->vpr_id, rs->id);

							//fixed_sinks[net->local_id].push_back(rs);
						//}
					}
				}

				prev_num_overused_nodes = std::numeric_limits<unsigned long>::max();
				has_unroutable_sinks = false;

				for (int i = 0; i < num_vertices(partitioner.orig_g); ++i) {
					congestion[i].recalc_occ = 0; 
				}

				for (int i = procid*pow(2, current_level); i < nets_to_route.size(); i += initial_num_procs) {
					for (int j = 0; j < pow(2, current_level) && i+j < nets_to_route.size(); ++j) {
						net_t *net = nets_to_route[i+j].second;

						if (!routed_sinks[net->local_id].empty()) {
							zlog_level(delta_log, ROUTER_V3, "Checking net index %d\n", i+j);

							check_route_tree(route_trees[net->local_id], *net, routed_sinks[net->local_id], partitioner.orig_g);
							recalculate_occ(route_trees[net->local_id], partitioner.orig_g, congestion);
						} else {
							zlog_level(delta_log, ROUTER_V3, "Not checking net index %d because of empty route tree\n", i+j);
						}
					}
				}

				bool valid = true;
				for (int i = 0; i < num_vertices(partitioner.orig_g); ++i) {
					sprintf_rr_node(i, buffer);
					if (congestion[i].recalc_occ != congestion[i].occ) {
						zlog_error(delta_log, "Node %s occ mismatch, recalc: %d original: %d\n", buffer, congestion[i].recalc_occ, congestion[i].occ);
						valid = false;
					}
				}
				assert(valid);
			} 

			//if ((float)total_num_interpartition_sinks/total_num_sinks_to_route > 0.1) {
				//++current_level;
				//assert(current_level < partitioner.result_pid_by_level.size());
			//}

            auto update_cost_start = clock::now();

            if (iter == 0) {
                params.pres_fac = opts->initial_pres_fac;
                //update_costs_mpi(partitioner.orig_g[0], partitioner.result_pid_by_level[current_level], procid, congestion, win, params.pres_fac, 0);
				update_costs(partitioner.orig_g, congestion, params.pres_fac, 0);
            } else {
                params.pres_fac *= opts->pres_fac_mult;

				/* Avoid overflow for high iteration counts, even if acc_cost is big */
                params.pres_fac = std::min(params.pres_fac, static_cast<float>(HUGE_POSITIVE_FLOAT / 1e5));

                //update_costs_mpi(partitioner.orig_g[0], partitioner.result_pid_by_level[current_level], procid, congestion, win, params.pres_fac, opts->acc_fac);
				update_costs(partitioner.orig_g, congestion, params.pres_fac, opts->acc_fac);
            }

            update_cost_time = clock::now()-update_cost_start;
        }

        iter_time += clock::now()-iter_start;

        //total_route_time += route_time;
        total_greedy_route_time += greedy_route_time;
        total_update_cost_time += update_cost_time;
        total_analyze_timing_time += analyze_timing_time;
        total_iter_time += iter_time;

		if (procid == 0) {
			printf("Iteration time: %g s.\n", duration_cast<nanoseconds>(iter_time).count() / 1e9);
			printf("\tRoute time: %g s.\n", duration_cast<nanoseconds>(greedy_route_time).count() / 1e9);
			printf("\tUpdate cost time: %g s.\n", duration_cast<nanoseconds>(update_cost_time).count() / 1e9);
			printf("\tAnalyze timing time: %g s.\n", duration_cast<nanoseconds>(analyze_timing_time).count() / 1e9);
			printf("Critical path: %g ns\n", crit_path_delay);
		}

        //clock::time_point greedy_earliest_end_time = *std::min_element(begin(greedy_end_time), begin(greedy_end_time)+num_procs);

        //printf("greedy wait time: ");
        //for (int i = 0; i < num_procs; ++i) {
			//printf("%g (%g), ", duration_cast<nanoseconds>(greedy_end_time[i]-greedy_earliest_end_time).count() / 1e9, 100.0*(greedy_end_time[i]-greedy_earliest_end_time)/greedy_route_time);
        //}
		//printf("\n");
    }

    //sprintf(buffer, "%s_run_%d_rr_graph_occ.txt", s_circuit_name, run);
    //dump_rr_graph_occ(congestion, num_vertices(g), buffer);
	if (initial_procid == 0) {
		if (routed) {
			printf("Routed in %d iterations. Total iteration time: %g\n", iter, duration_cast<nanoseconds>(total_iter_time).count() / 1e9);
			printf("\tTotal route time: %g s.\n", duration_cast<nanoseconds>(total_greedy_route_time).count() / 1e9);
			printf("\t\tTotal greedy route time: %g s.\n", duration_cast<nanoseconds>(total_greedy_route_time).count() / 1e9);
			printf("\tTotal update cost time: %g s.\n", duration_cast<nanoseconds>(total_update_cost_time).count() / 1e9);
			printf("\tTotal analyze timing time: %g s.\n", duration_cast<nanoseconds>(total_analyze_timing_time).count() / 1e9);

			printf("Final critical path: %g ns\n", crit_path_delay);
		} else {
			printf("Failed to route in %d iterations. Total iteration time: %g\n", opts->max_router_iterations, duration_cast<nanoseconds>(total_iter_time).count() / 1e9);
			printf("\tTotal route time: %g s.\n", duration_cast<nanoseconds>(total_greedy_route_time).count() / 1e9);
			printf("\t\tTotal greedy route time: %g s.\n", duration_cast<nanoseconds>(total_greedy_route_time).count() / 1e9);
			printf("\tTotal update cost time: %g s.\n", duration_cast<nanoseconds>(total_update_cost_time).count() / 1e9);
			printf("\tTotal analyze timing time: %g s.\n", duration_cast<nanoseconds>(total_analyze_timing_time).count() / 1e9);
		}
	} 

	MPI_Barrier(MPI_COMM_WORLD);

	MPI_Finalize();
	exit(0);

    //delete_graph(g);
    delete_net_timing(nets, global_nets, net_timing);	
    delete [] congestion;
    delete [] net_timing;
    delete [] states;

    return routed;
}

bool mpi_spatial_route_rma(t_router_opts *opts, struct s_det_routing_arch det_routing_arch, t_direct_inf *directs, int num_directs, t_segment_inf *segment_inf, t_timing_inf timing_inf)
{
    using std::chrono::duration_cast;
    using std::chrono::nanoseconds;
    using clock = std::chrono::high_resolution_clock;

	int num_procs, procid;
    MPI_Comm_size(MPI_COMM_WORLD, &num_procs);
    MPI_Comm_rank(MPI_COMM_WORLD, &procid);
	init_congestion_mpi_datatype();
	init_datatypes();

    init_logging();
    zlog_set_record("custom_output", delta_log_output);
    zlog_set_record("missing_edge", missing_edge_log_output);
    zlog_set_record("ss", ss_log_output);
    delta_log_files.resize(opts->max_router_iterations);
    for (int i = 0; i < opts->max_router_iterations; ++i) {
        delta_log_files[i].resize(num_procs, nullptr);
    }
    missing_edge_log_files.resize(opts->max_router_iterations);
    for (int i = 0; i < opts->max_router_iterations; ++i) {
        missing_edge_log_files[i].resize(num_procs, nullptr);
    }

    //fclose(fopen(LOG_PATH_PREFIX"iter__tid_.log", "w"));

    test_fast_graph();
    test_topo();
    test_fm();
    test_filter_graph();
    test_partition_graph();
    test_connected_components();

    vector<RRGraph *> graphs;
    rr_graph_partitioner partitioner;
    partitioner.init(opts, &det_routing_arch, directs, num_directs, segment_inf, &timing_inf);
    partitioner.partition(num_procs, graphs);
    //for (const auto &g : graphs) {
        //routability(*g);
    //}
    //
    //RRGraph combined_g;
    //add_vertex(combined_g, num_vertices(partitioner.orig_g));

    //for (const auto &g : graphs) {
        //for (const auto &e : get_edges(*g)) {
            //int from = get_source(*g, e);
            //int to = get_target(*g, e);
            //const auto &from_ver = get_vertex_props(*g, from);
            //const auto &to_ver = get_vertex_props(*g, to);

            //if (is_channel(from_ver) && is_channel(to_ver)) {
                //assert(!has_edge(combined_g, from, to));
                //add_edge(combined_g, from, to);
            //}
        //}
    //}

    //for (const auto &e : get_edges(partitioner.orig_g)) {
        //int from = get_source(partitioner.orig_g, e);
        //int to = get_target(partitioner.orig_g, e);
        //const auto &from_ver = get_vertex_props(partitioner.orig_g, from);
        //const auto &to_ver = get_vertex_props(partitioner.orig_g, to);
        //if (!is_channel(from_ver) || !is_channel(to_ver)) {
            //assert(!has_edge(combined_g, from, to));
            //add_edge(combined_g, from, to);
        //}
    //}

    //printf("Combined/Orig graph has %d/%d (%g) edges.\n", num_edges(combined_g), num_edges(partitioner.orig_g), 100.0*num_edges(combined_g)/num_edges(partitioner.orig_g));

    //goto lol;

    //extern int num_types;
    //extern struct s_type_descriptor *type_descriptors;
    //extern int nx, ny;
    //extern struct s_grid_tile **grid;

    //free_rr_graph();

    //int warnings;

    //build_rr_graph(GRAPH_UNIDIR, num_types, type_descriptors, nx, ny, grid,
            //opts->fixed_channel_width, NULL, det_routing_arch.switch_block_type,
            //det_routing_arch.Fs, det_routing_arch.num_segment,
            //det_routing_arch.num_switch, segment_inf,
            //det_routing_arch.global_route_switch,
            //det_routing_arch.delayless_switch, timing_inf,
            //det_routing_arch.wire_to_ipin_switch, opts->base_cost_type,
            //directs, num_directs, FALSE,
            //&warnings);

    //RRGraph channel_with_interior_g;
    //init_channel_only_graph(channel_with_interior_g);

    //dump_rr_graph(channel_with_interior_g, "/Volumes/DATA/rr_graph_channel_with_interior.txt");

    //RRGraph orig_g;
    //init_graph(orig_g);

    //dump_rr_graph(orig_g, "/Volumes/DATA/rr_graph.txt");

    //free_rr_graph();
    //for (int i = 0; i < det_routing_arch.num_segment; ++i) {
        //for (int j = 0; j < segment_inf[i].sb_len; ++j) {
            //if (j != 0 && j != segment_inf[i].sb_len-1) {
                //segment_inf[i].sb[j] = FALSE;
            //}
        //}
    //}
    //build_rr_graph(GRAPH_UNIDIR, num_types, type_descriptors, nx, ny, grid,
            //opts->fixed_channel_width, NULL, det_routing_arch.switch_block_type,
            //det_routing_arch.Fs, det_routing_arch.num_segment,
            //det_routing_arch.num_switch, segment_inf,
            //det_routing_arch.global_route_switch,
            //det_routing_arch.delayless_switch, timing_inf,
            //det_routing_arch.wire_to_ipin_switch, opts->base_cost_type,
            //directs, num_directs, FALSE,
            //&warnings);

    //RRGraph channel_without_interior_g;
    //init_channel_only_graph(channel_without_interior_g);

    //dump_rr_graph(channel_without_interior_g, "/Volumes/DATA/rr_graph_channel_without_interior.txt");

    //vector<vector<vector<int>>> all_partition_components(opts->num_threads);
    ////for (int x = 0; x < nx+1; ++x) {
        ////for (int y = 0; y < ny+1; ++y) {
            ////init_partitioned_graph_5(opts->num_threads, channel_with_interior_g, channel_without_interior_g, orig_g, graphs, x, y);
        ////}
    ////}
    ////init_partitioned_graph_4(opts->num_threads, channel_with_interior_g, channel_without_interior_g, orig_g, graphs, 14, 14);
    //init_partitioned_graph_4(opts->num_threads, channel_with_interior_g, channel_without_interior_g, orig_g, graphs, all_partition_components);

    //for (int i = 0; i < graphs.size(); ++i) {
        //char filename[256];

        //sprintf(filename, "/Volumes/DATA/rr_graph_%d.txt", i);
        //dump_rr_graph(*graphs[i], filename);

        //sprintf(filename, "/Volumes/DATA/rr_graph_%d_edges.txt", i);
        //dump_edges(*graphs[i], filename);
    //}

    vector<net_t> nets;
    vector<net_t> global_nets;
    init_nets(nets, global_nets, opts->bb_factor, opts->large_bb);	

    vector<net_t *> nets_ptr(nets.size());
    for (int i = 0; i < nets.size(); ++i) {
        nets_ptr[i] = &nets[i];
    }
    extern s_bb *route_bb;
    sort(nets_ptr.begin(), nets_ptr.end(), [] (const net_t *a, const net_t *b) -> bool {
            return get_bounding_box_area(route_bb[a->vpr_id]) > get_bounding_box_area(route_bb[b->vpr_id]);
            });
    int rank = 0;
    for (auto &net : nets_ptr) {
        net->bb_area_rank = rank++;
    }

	route_state_t *states;
	congestion_t *congestion;
	vector<route_tree_t> route_trees;
	t_net_timing *net_timing;
	MPI_Win win;
	init_route_structs_mpi(partitioner.orig_g, nets, global_nets, &states, &congestion, &win, route_trees, &net_timing);

    route_parameters_t params;
    params.criticality_exp = opts->criticality_exp;
    params.astar_fac = opts->astar_fac;
    params.max_criticality = opts->max_criticality;
    params.bend_cost = opts->bend_cost;
	params.pres_fac = opts->first_iter_pres_fac; /* Typically 0 -> ignore cong. */

    char buffer[256];

    vector<pair<box, net_t *>> nets_to_route;
    vector<pair<box, net_t *>> nets_to_partition;
    //vector<vector<int>> overlaps;
    vector<vector<pair<int, int>>> overlaps;
    vector<vector<int>> partitions;
    vector<bool> has_interpartition_overlap;
    
    for (auto &net : nets) {
        box b = bg::make_inverse<box>();

        bg::expand(b, point(net.source.x, net.source.y));
        for (const auto &sink : net.sinks) {
            bg::expand(b, point(sink.x, sink.y));
        }
        bg::subtract_value(b.min_corner(), 1+opts->bb_factor);
        bg::add_value(b.max_corner(), opts->bb_factor);

        nets_to_route.push_back(make_pair(b, &net));
    }
    std::sort(begin(nets_to_route), end(nets_to_route), [] (const pair<box, net_t *> &a, const pair<box, net_t *> &b) -> bool {
            return a.second->sinks.size() > b.second->sinks.size();
            });

    bool routed = false;

    //vector<clock::time_point> greedy_end_time(num_procs);
    clock::duration total_greedy_route_time = clock::duration::zero();
    clock::duration total_update_cost_time = clock::duration::zero();
    clock::duration total_analyze_timing_time = clock::duration::zero();
    clock::duration total_iter_time = clock::duration::zero();

    int iter;
    float crit_path_delay;
	int current_level = 0;

    for (iter = 0; iter < opts->max_router_iterations && !routed; ++iter) {
        clock::duration greedy_route_time = clock::duration::zero();
        clock::duration update_cost_time = clock::duration::zero();
        clock::duration partitioning_time = clock::duration::zero();
        clock::duration analyze_timing_time = clock::duration::zero();
        clock::duration iter_time = clock::duration::zero();

        //for (int i = 0; i < opts->num_threads; ++i) {
            //sprintf(buffer, LOG_PATH_PREFIX"iter_%d_tid_%d.log", iter, i);
            //fclose(fopen(buffer, "w"));
        //}

		MPI_Barrier(MPI_COMM_WORLD);
        
        auto iter_start = clock::now();

        //auto route_start = clock::now();

        //tbb::enumerable_thread_specific<state_t *> state_tls;
        //
		vector<perf_t> perfs(num_procs);
        vector<int> thread_num_nets_routed(num_procs, 0);
        vector<int> thread_num_nets_to_route(num_procs, 0);
        vector<int> thread_num_sinks_routed(num_procs, 0);
        vector<int> thread_num_sinks_to_route(num_procs, 0);
		vector<int> thread_num_interpartition_sinks(num_procs, 0);
        bool has_interpartition_sinks = false;
        vector<vector<interpartition_sink_t>> interpartition_sinks(nets.size()); 

		sprintf(buffer, "%d", iter);
		zlog_put_mdc("iter", buffer);

		sprintf(buffer, "%d", procid);
		zlog_put_mdc("tid", buffer);

		zlog_error(delta_log, "We're fucking here\n");

        zlog_info(delta_log, "Routing iteration: %d\n", iter);
        printf("Routing iteration: %d\n", iter);

        auto greedy_route_start = clock::now();

		perfs[procid].num_heap_pushes = 0;
		perfs[procid].num_heap_pops = 0;
		perfs[procid].num_neighbor_visits = 0;

		int i = procid;
		while (i < nets_to_route.size()) {
			net_t *net = nets_to_route[i].second;

            update_sink_criticalities(*net, net_timing[net->vpr_id], params);

			//printf("Routing net %d\n", net->vpr_id);

			//auto rip_up_start = clock::now();
			//if (greedy_rip_up_all) {
			route_tree_mark_all_nodes_to_be_ripped(route_trees[net->local_id], partitioner.orig_g);
			//} else {
			//route_tree_mark_congested_nodes_to_be_ripped(route_trees[net->local_id], *graphs[procid], congestion);
			//}
			route_tree_rip_up_marked_mpi_rma(route_trees[net->local_id], partitioner.orig_g, partitioner.result_pid_by_level[current_level], procid, congestion, win, params.pres_fac);

			//local_perf.total_rip_up_time += clock::now()-rip_up_start;

			//auto route_start = clock::now();

			vector<sink_t *> sinks;	
			get_sinks_to_route(net, route_trees[net->local_id], sinks);

			if (!sinks.empty()) {
				route_net_mpi_rma(partitioner.orig_g, partitioner.result_pid_by_level[current_level], procid, win, net->vpr_id, &net->source, sinks, params, states, congestion, route_trees[net->local_id], net_timing[net->vpr_id], interpartition_sinks[net->local_id], &perfs[procid]);

				if (!has_interpartition_sinks) {
					has_interpartition_sinks = !interpartition_sinks[net->local_id].empty();
				}

				thread_num_interpartition_sinks[procid] += interpartition_sinks[net->local_id].size();

				++thread_num_nets_routed[procid];
				++thread_num_nets_to_route[procid];

				thread_num_sinks_to_route[procid] += sinks.size();
				thread_num_sinks_routed[procid] += sinks.size();
			}

			//local_perf.total_route_time += clock::now()-rip_up_start;
			i += num_procs;
		}

		//greedy_end_time = clock::now();
		MPI_Barrier(MPI_COMM_WORLD);

        greedy_route_time = clock::now()-greedy_route_start;

		//if (has_interpartition_sinks) {
			//for (const auto &net_interpartition_sinks : interpartition_sinks) {
				//for (const auto &is : net_interpartition_sinks) {
					//vector<RRNode> added_nodes;
					//for (const auto &node : is.path) {
						//if (node.update_cost) {
							//added_nodes.push_back(node.rr_node_id);
						//}
					//}
					//update_one_cost(partitioner.orig_g, congestion, added_nodes.begin(), added_nodes.end(), 1, params.pres_fac, false, nullptr);
				//}
			//}
		//}

        //if (greedy_rip_up_all) {
            //next_greedy_rip_up_iter += greedy_rip_up_all_period;
            //++greedy_rip_up_all_period;
            //prev_num_overused_nodes = std::numeric_limits<unsigned long>::max();
        //}

        //route_time = clock::now()-route_start;

        iter_time = clock::now()-iter_start;

		int total_num_sinks_to_route;
		if (procid == 0) {
			int total_num_nets_to_route = thread_num_nets_to_route[0];
			int total_num_nets_routed = thread_num_nets_routed[0];
			int total_num_sinks_routed = thread_num_sinks_routed[0];
			total_num_sinks_to_route = thread_num_sinks_to_route[0];
			int total_num_interpartition_sinks = thread_num_interpartition_sinks[0];

			for (int i = 1; i < num_procs; ++i) {
				MPI_Recv(&thread_num_nets_to_route[i], 1, MPI_INT, i, i, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
				MPI_Recv(&thread_num_nets_routed[i], 1, MPI_INT, i, i, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
				MPI_Recv(&thread_num_sinks_to_route[i], 1, MPI_INT, i, i, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
				MPI_Recv(&thread_num_sinks_routed[i], 1, MPI_INT, i, i, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
				MPI_Recv(&thread_num_interpartition_sinks[i], 1, MPI_INT, i, i, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

				total_num_nets_to_route += thread_num_nets_to_route[i];
				total_num_nets_routed += thread_num_nets_routed[i];
				total_num_sinks_routed += thread_num_sinks_routed[i];
				total_num_sinks_to_route += thread_num_sinks_to_route[i];
				total_num_interpartition_sinks += thread_num_interpartition_sinks[i];
			}
			assert(total_num_nets_to_route == total_num_nets_routed);
			assert(total_num_nets_to_route == nets.size());
			assert(total_num_sinks_to_route == total_num_sinks_routed);

			printf("num nets routed: ");
			for (int i = 0; i < num_procs; ++i) {
				printf("%d/%d (%g), ", thread_num_nets_routed[i], thread_num_nets_to_route[i], thread_num_nets_routed[i]*100.0/thread_num_nets_to_route[i]);
			}
			printf("\n");

			printf("num sinks routed: ");
			for (int i = 0; i < num_procs; ++i) {
				printf("%d/%d (%g), ", thread_num_sinks_routed[i], thread_num_sinks_to_route[i], thread_num_sinks_routed[i]*100.0/thread_num_sinks_to_route[i]);
			}
			printf("\n");

			printf("Total num nets routed: %d/%d (%g)\n", total_num_nets_routed, total_num_nets_to_route, total_num_nets_routed*100.0/total_num_nets_to_route);
			printf("Total num sinks routed: %d/%d (%g)\n", total_num_sinks_routed, total_num_sinks_to_route, (total_num_sinks_routed)*100.0/total_num_sinks_to_route);
			printf("Total num interpartition sinks routed: %d/%d (%g)\n", total_num_interpartition_sinks, total_num_sinks_to_route, total_num_interpartition_sinks*100.0/total_num_sinks_to_route);
		} else {
			MPI_Send(&thread_num_nets_to_route[procid], 1, MPI_INT, 0, procid, MPI_COMM_WORLD);
			MPI_Send(&thread_num_nets_routed[procid], 1, MPI_INT, 0, procid, MPI_COMM_WORLD);
			MPI_Send(&thread_num_sinks_to_route[procid], 1, MPI_INT, 0, procid, MPI_COMM_WORLD);
			MPI_Send(&thread_num_sinks_routed[procid], 1, MPI_INT, 0, procid, MPI_COMM_WORLD);
			MPI_Send(&thread_num_interpartition_sinks[procid], 1, MPI_INT, 0, procid, MPI_COMM_WORLD);
		}

		sprintf(buffer, "/Volumes/DATA/occ_before_%d.txt", procid);
		FILE *occ = fopen(buffer, "w");
		for (const auto &v : get_vertices(partitioner.orig_g)) {
			fprintf(occ, "%d\n", congestion[v].occ);
		}
		fclose(occ);

		//MPI_Win_fence(0, win);
		for (const auto &rr_node : get_vertices(partitioner.orig_g)) {
			int rr_node_pid = partitioner.result_pid_by_level[current_level][rr_node];
			int from_pid = rr_node_pid == -1 ? 0 : rr_node_pid;

			if (from_pid != procid) {
				congestion[rr_node].occ = std::numeric_limits<int>::min();
				assert(MPI_Win_lock(MPI_LOCK_SHARED, from_pid, 0, win) == MPI_SUCCESS);
				assert(MPI_Get(&congestion[rr_node].occ, 1, get_occ_dt(),
							from_pid,
							get_occ_disp(rr_node), 1, get_occ_dt(),
							win) == MPI_SUCCESS);
				assert(MPI_Win_unlock(from_pid, win) == MPI_SUCCESS);
				assert(congestion[rr_node].occ != std::numeric_limits<int>::min());
			}
		}
		//MPI_Win_fence(0, win);

		sprintf(buffer, "/Volumes/DATA/occ_after_%d.txt", procid);
		occ = fopen(buffer, "w");
		for (const auto &v : get_vertices(partitioner.orig_g)) {
			fprintf(occ, "%d\n", congestion[v].occ);
		}
		fclose(occ);

		/* checking */
        for (int i = 0; i < num_vertices(partitioner.orig_g); ++i) {
            congestion[i].recalc_occ = 0; 
        }

        for (int i = 0; i < nets_to_route.size(); ++i) {
			if (i % num_procs == procid) {
				net_t *net = nets_to_route[i].second;

				check_route_tree(route_trees[net->local_id], *net, partitioner.orig_g);
				recalculate_occ(route_trees[net->local_id], partitioner.orig_g, congestion);
			}
        }

		int *recalc_occ = new int[num_vertices(partitioner.orig_g)];
		for (int i = 0; i < num_vertices(partitioner.orig_g); ++i) {
			recalc_occ[i] = congestion[i].recalc_occ;	
		}
		if (procid == 0) {
			int *recv_recalc_occ = new int[num_vertices(partitioner.orig_g)];
			MPI_Reduce(recalc_occ, recv_recalc_occ, num_vertices(partitioner.orig_g), MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);
			for (int i = 0; i < num_vertices(partitioner.orig_g); ++i) {
				congestion[i].recalc_occ = recv_recalc_occ[i];
			}
			delete [] recv_recalc_occ;
		} else {
			MPI_Reduce(recalc_occ, nullptr, num_vertices(partitioner.orig_g), MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);
		}
		delete [] recalc_occ;
		
		//MPI_Op recalc_sum_op;
		//MPI_Op_create((MPI_User_function *)recalc_sum, 1, &recalc_sum_op);	
		//if (procid == 0) {
			//MPI_Reduce(MPI_IN_PLACE, congestion, num_vertices(partitioner.orig_g), recalc_occ_dt, recalc_sum_op, 0, MPI_COMM_WORLD);
		//} else {
			//MPI_Reduce(congestion, nullptr, num_vertices(partitioner.orig_g), recalc_occ_dt, recalc_sum_op, 0, MPI_COMM_WORLD);
		//}
		//MPI_Op_free(&recalc_sum_op);

		//if (procid == 0) {
			//int num_recvs = 0;
			//for (int i = 0; i < num_vertices(*graphs[0]); ++i) {
				////int rr_node_pid = partitioner.result_pid_by_level[current_level][i];
				////int from_pid = rr_node_pid == -1 ? 0 : rr_node_pid;

				////if (from_pid != procid) {
				//for (int from_pid = 1; from_pid < num_procs; ++from_pid) {
					//int recalc_occ;
					//MPI_Recv(&recalc_occ, 1, MPI_INT, from_pid, i, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
					//zlog_level(delta_log, ROUTER_V3, "%d Recvd %d recalc_occ %d\n", from_pid, i, recalc_occ);
					//congestion[i].recalc_occ += recalc_occ;

					//++num_recvs;
				//}
			//}
			//printf("Num recvs: %d\n", num_recvs);
		//} else {
			//int num_sends = 0;

			//for (int i = 0; i < num_vertices(*graphs[0]); ++i) {
				////int pid = partitioner.result_pid_by_level[current_level][i];
				////if (pid == procid) {
					//zlog_level(delta_log, ROUTER_V3, "Sending %d recalc_occ %d\n", i, congestion[i].recalc_occ);
					//MPI_Send(&congestion[i].recalc_occ, 1, MPI_INT, 0, i, MPI_COMM_WORLD);

					//++num_sends;
				////}
			//}

			//printf("[%d] Num sends: %d\n", procid, num_sends);
		//}

		if (procid == 0) {
			bool valid = true;
//#undef sprintf_rr_node
			for (int i = 0; i < num_vertices(partitioner.orig_g); ++i) {
				sprintf_rr_node(i, buffer);
				if (congestion[i].recalc_occ != congestion[i].occ) {
					zlog_error(delta_log, "Node %s occ mismatch, recalc: %d original: %d\n", buffer, congestion[i].recalc_occ, congestion[i].occ);
					valid = false;
				}
			}
//#define sprintf_rr_node(...)
			assert(valid);
		}

        int overused_total_bb_rank = 0;
        int num_congested_nets = 0;

        for (int i = 0; i < nets_to_route.size(); ++i) {
			if (i % num_procs == procid) {
				net_t *net = nets_to_route[i].second;

				vector<int> overused_rr_node;
				assert(route_trees[net->local_id].root_rt_nodes.size() == 1);
				get_overused_nodes(route_trees[net->local_id], route_trees[net->local_id].root_rt_nodes[0], partitioner.orig_g, congestion, overused_rr_node);
				if (!overused_rr_node.empty()) {
					zlog_level(delta_log, ROUTER_V1, "Net %d bb_rank %d has %lu overused nodes:\n", net->vpr_id, net->bb_area_rank, overused_rr_node.size());
					for (const auto &item : overused_rr_node) {
						zlog_level(delta_log, ROUTER_V1, "%d ", item);
					}
					zlog_level(delta_log, ROUTER_V1, "\n");
					overused_total_bb_rank += net->bb_area_rank;
					++num_congested_nets;
				}
			}
        }

        zlog_level(delta_log, ROUTER_V1, "Average overused net bb rank: %g\n", (float)overused_total_bb_rank/num_congested_nets);
        zlog_level(delta_log, ROUTER_V1, "Num congested nets: %d\n", num_congested_nets);

        iter_start = clock::now();

        if (feasible_routing(partitioner.orig_g, congestion)) {
            //dump_route(*current_traces_ptr, "route.txt");
            routed = true;
        } else {
            unsigned long num_overused_nodes = 0;
			vector<int> overused_nodes_by_type(NUM_RR_TYPES, 0);
            for (int i = 0; i < num_vertices(partitioner.orig_g); ++i) {
                if (congestion[i].occ > get_vertex_props(partitioner.orig_g, i).capacity) {
					const auto &v_p = get_vertex_props(partitioner.orig_g, i);
					++overused_nodes_by_type[v_p.type];

                    ++num_overused_nodes;
                }
            }
            zlog_info(delta_log, "Num overused nodes: %lu/%d (%.2f)\n", num_overused_nodes, num_vertices(partitioner.orig_g), num_overused_nodes*100.0/num_vertices(partitioner.orig_g));

			if (procid == 0) {
				printf("Num overused nodes: %lu/%d (%.2f)\n", num_overused_nodes, num_vertices(partitioner.orig_g), num_overused_nodes*100.0/num_vertices(partitioner.orig_g));
				static const char *name_type[] = { "SOURCE", "SINK", "IPIN", "OPIN",
					"CHANX", "CHANY", "INTRA_CLUSTER_EDGE" };
				for (int i = 0; i < overused_nodes_by_type.size(); ++i) {
					printf("\t%s: %d (%g)\n", name_type[i], overused_nodes_by_type[i], overused_nodes_by_type[i]*100.0/num_overused_nodes);
				}
			}

			//if ((float)total_num_interpartition_sinks/total_num_sinks_to_route > 0.1) {
				//++current_level;
				//assert(current_level < partitioner.result_pid_by_level.size());
			//}

            auto update_cost_start = clock::now();

            if (iter == 0) {
                params.pres_fac = opts->initial_pres_fac;
                //update_costs_mpi(partitioner.orig_g[0], partitioner.result_pid_by_level[current_level], procid, congestion, win, params.pres_fac, 0);
				update_costs(partitioner.orig_g, congestion, params.pres_fac, 0);
            } else {
                params.pres_fac *= opts->pres_fac_mult;

				/* Avoid overflow for high iteration counts, even if acc_cost is big */
                params.pres_fac = std::min(params.pres_fac, static_cast<float>(HUGE_POSITIVE_FLOAT / 1e5));

                //update_costs_mpi(partitioner.orig_g[0], partitioner.result_pid_by_level[current_level], procid, congestion, win, params.pres_fac, opts->acc_fac);
				update_costs(partitioner.orig_g, congestion, params.pres_fac, opts->acc_fac);
            }

            update_cost_time = clock::now()-update_cost_start;
        }

        auto analyze_timing_start = clock::now();

		//MPI_Allgather(, , ,
				//, , ,
				//MPI_COMM_WORLD);
		int num_delays = 0;
		for (int i = procid; i < nets_to_route.size(); i += num_procs) {
			net_t *net = nets_to_route[i].second;
			num_delays += net->sinks.size();
		}

		int *recvcounts = new int[num_procs];
		int *displs = new int[num_procs];

		for (int i = 0; i < num_procs; ++i) {
			recvcounts[i] = 0;
		}
		for (int i = 0; i < nets_to_route.size(); ++i) {
			net_t *net = nets_to_route[i].second;
			int to_pid = i % num_procs;
			recvcounts[to_pid] += net->sinks.size();
		}
		for (int i = 0; i < num_procs; ++i) {
			displs[i] = 0;
			for (int j = 0; j < i; ++j) {
				displs[i] += recvcounts[j];
			}
		}

		float *delays = new float[num_delays];
		int idx = 0;
		for (int i = procid; i < nets_to_route.size(); i += num_procs) {
			net_t *net = nets_to_route[i].second;
			zlog_level(delta_log, ROUTER_V3, "Net %d send delays:\n", net->vpr_id);
			for (int j = 1; j <= net->sinks.size(); ++j) {
				delays[idx] = net_timing[net->vpr_id].delay[j];
				zlog_level(delta_log, ROUTER_V3, "\t%g\n", delays[idx]);
				++idx;
			}
		}

		assert(idx == num_delays);
		assert(idx == recvcounts[procid]);
		int temp_sum = 0;
		for (int i = 0; i < num_procs; ++i) {
			temp_sum += recvcounts[i];
		}

		float *all_delays;

		if (procid == 0) {
			assert(temp_sum == total_num_sinks_to_route);

			all_delays = new float[total_num_sinks_to_route];
		} else {
			all_delays = nullptr;
		}
		
		MPI_Gatherv(delays, recvcounts[procid], MPI_FLOAT, all_delays, recvcounts, displs, MPI_FLOAT, 0, MPI_COMM_WORLD);

		float *all_crits = nullptr;

		if (procid == 0) {
			for (int i = 0; i < num_procs; ++i) {
				idx = 0;
				for (int j = i; j < nets_to_route.size(); j += num_procs) {
					net_t *net = nets_to_route[j].second;
					zlog_level(delta_log, ROUTER_V3, "Net %d recv delays:\n", net->vpr_id);
					for (int k = 1; k <= net->sinks.size(); ++k) {
						net_timing[net->vpr_id].delay[k] = all_delays[displs[i] + idx];
						zlog_level(delta_log, ROUTER_V3, "\t%g\n", net_timing[net->vpr_id].delay[k]);
						++idx;
					}
				}
				assert(idx == recvcounts[i]);
			}

			crit_path_delay = analyze_timing(net_timing);

			all_crits = all_delays;

			idx = 0;
			for (int i = 0; i < num_procs; ++i) {
				idx = 0;
				for (int j = i; j < nets_to_route.size(); j += num_procs) {
					net_t *net = nets_to_route[j].second;
					zlog_level(delta_log, ROUTER_V3, "Net %d send crits:\n", net->vpr_id);
					for (int k = 1; k <= net->sinks.size(); ++k) {
						all_crits[displs[i] + idx] = net_timing[net->vpr_id].timing_criticality[k];
						zlog_level(delta_log, ROUTER_V3, "\t%g\n", all_crits[displs[i] + idx]);
						++idx;
					}
				}
				assert(idx == recvcounts[i]);
			}
		} else {
			all_crits = nullptr;
		}

		float *crits = delays;

		MPI_Scatterv(all_crits, recvcounts, displs, MPI_FLOAT, crits, recvcounts[procid], MPI_FLOAT, 0, MPI_COMM_WORLD);

		idx = 0;
		for (int i = procid; i < nets_to_route.size(); i += num_procs) {
			net_t *net = nets_to_route[i].second;
			zlog_level(delta_log, ROUTER_V3, "Net %d recv crits:\n", net->vpr_id);
			for (int j = 1; j <= net->sinks.size(); ++j) {
				net_timing[net->vpr_id].timing_criticality[j] = crits[idx];
				zlog_level(delta_log, ROUTER_V3, "\t%g\n", net_timing[net->vpr_id].timing_criticality[j]);
				++idx;
			}
		}

		delete [] delays;
		delete [] all_delays;

        analyze_timing_time = clock::now()-analyze_timing_start;

        iter_time += clock::now()-iter_start;

        //total_route_time += route_time;
        total_greedy_route_time += greedy_route_time;
        total_update_cost_time += update_cost_time;
        total_analyze_timing_time += analyze_timing_time;
        total_iter_time += iter_time;

		if (procid == 0) {
			printf("Iteration time: %g s.\n", duration_cast<nanoseconds>(iter_time).count() / 1e9);
			printf("\tRoute time: %g s.\n", duration_cast<nanoseconds>(greedy_route_time).count() / 1e9);
			printf("\tUpdate cost time: %g s.\n", duration_cast<nanoseconds>(update_cost_time).count() / 1e9);
			printf("\tAnalyze timing time: %g s.\n", duration_cast<nanoseconds>(analyze_timing_time).count() / 1e9);
			printf("Critical path: %g ns\n", crit_path_delay);

			printf("num_heap_pushes: ");
			for (int i = 0; i < num_procs; ++i) {
				printf("%lu ", perfs[i].num_heap_pushes);
			}
			printf("\n");

			printf("num_heap_pops: ");
			for (int i = 0; i < num_procs; ++i) {
				printf("%lu ", perfs[i].num_heap_pops);
			}
			printf("\n");

			printf("num_neighbor_visits: ");
			for (int i = 0; i < num_procs; ++i) {
				printf("%lu ", perfs[i].num_neighbor_visits);
			}
			printf("\n");

			unsigned long total_num_heap_pushes = 0;
			unsigned long total_num_heap_pops = 0;
			unsigned long total_num_neighbor_visits = 0;

			for (int i = 0; i < num_procs; ++i) {
				total_num_heap_pushes += perfs[i].num_heap_pushes;
				total_num_heap_pops += perfs[i].num_heap_pops;
				total_num_neighbor_visits += perfs[i].num_neighbor_visits;
			}

			printf("total_num_heap_pushes: %lu\n", total_num_heap_pushes);
			printf("total_num_heap_pops: %lu\n", total_num_heap_pops); 
			printf("total_num_neighbor_visits: %lu\n", total_num_neighbor_visits);
		}

        //clock::time_point greedy_earliest_end_time = *std::min_element(begin(greedy_end_time), begin(greedy_end_time)+num_procs);

        //printf("greedy wait time: ");
        //for (int i = 0; i < num_procs; ++i) {
			//printf("%g (%g), ", duration_cast<nanoseconds>(greedy_end_time[i]-greedy_earliest_end_time).count() / 1e9, 100.0*(greedy_end_time[i]-greedy_earliest_end_time)/greedy_route_time);
        //}
		//printf("\n");
    }

    //sprintf(buffer, "%s_run_%d_rr_graph_occ.txt", s_circuit_name, run);
    //dump_rr_graph_occ(congestion, num_vertices(g), buffer);

    if (routed) {
        printf("Routed in %d iterations. Total iteration time: %g\n", iter, duration_cast<nanoseconds>(total_iter_time).count() / 1e9);
            printf("\tTotal route time: %g s.\n", duration_cast<nanoseconds>(total_greedy_route_time).count() / 1e9);
                printf("\t\tTotal greedy route time: %g s.\n", duration_cast<nanoseconds>(total_greedy_route_time).count() / 1e9);
            printf("\tTotal update cost time: %g s.\n", duration_cast<nanoseconds>(total_update_cost_time).count() / 1e9);
            printf("\tTotal analyze timing time: %g s.\n", duration_cast<nanoseconds>(total_analyze_timing_time).count() / 1e9);

        printf("Final critical path: %g ns\n", crit_path_delay);
    } else {
        printf("Failed to route in %d iterations. Total iteration time: %g\n", opts->max_router_iterations, duration_cast<nanoseconds>(total_iter_time).count() / 1e9);
            printf("\tTotal route time: %g s.\n", duration_cast<nanoseconds>(total_greedy_route_time).count() / 1e9);
                printf("\t\tTotal greedy route time: %g s.\n", duration_cast<nanoseconds>(total_greedy_route_time).count() / 1e9);
            printf("\tTotal update cost time: %g s.\n", duration_cast<nanoseconds>(total_update_cost_time).count() / 1e9);
            printf("\tTotal analyze timing time: %g s.\n", duration_cast<nanoseconds>(total_analyze_timing_time).count() / 1e9);
    }

	MPI_Finalize();
	exit(0);

    //delete_graph(g);
    delete_net_timing(nets, global_nets, net_timing);	
    delete [] congestion;
    delete [] net_timing;
    delete [] states;

    return routed;
}

bool spatial_route(t_router_opts *opts, struct s_det_routing_arch det_routing_arch, t_direct_inf *directs, int num_directs, t_segment_inf *segment_inf, t_timing_inf timing_inf)
{
    using std::chrono::duration_cast;
    using std::chrono::nanoseconds;
    using clock = std::chrono::high_resolution_clock;

    init_logging();
    zlog_set_record("custom_output", delta_log_output);
    zlog_set_record("missing_edge", missing_edge_log_output);
    zlog_set_record("ss", ss_log_output);
    delta_log_files.resize(opts->max_router_iterations);
    for (int i = 0; i < opts->max_router_iterations; ++i) {
        delta_log_files[i].resize(opts->num_threads, nullptr);
    }
    missing_edge_log_files.resize(opts->max_router_iterations);
    for (int i = 0; i < opts->max_router_iterations; ++i) {
        missing_edge_log_files[i].resize(opts->num_threads, nullptr);
    }

    //fclose(fopen(LOG_PATH_PREFIX"iter__tid_.log", "w"));

    test_fast_graph();
    test_topo();
    test_fm();
    test_filter_graph();
    test_partition_graph();
    test_connected_components();

    vector<RRGraph *> graphs;
    rr_graph_partitioner partitioner;
    partitioner.init(opts, &det_routing_arch, directs, num_directs, segment_inf, &timing_inf);
    partitioner.partition(opts->num_threads, graphs);
    //for (const auto &g : graphs) {
        //routability(*g);
    //}
    //
    RRGraph combined_g;
    add_vertex(combined_g, num_vertices(partitioner.orig_g));

    for (const auto &g : graphs) {
        for (const auto &e : get_edges(*g)) {
            int from = get_source(*g, e);
            int to = get_target(*g, e);
            const auto &from_ver = get_vertex_props(*g, from);
            const auto &to_ver = get_vertex_props(*g, to);

            if (is_channel(from_ver) && is_channel(to_ver)) {
                assert(!has_edge(combined_g, from, to));
                add_edge(combined_g, from, to);
            }
        }
    }

    for (const auto &e : get_edges(partitioner.orig_g)) {
        int from = get_source(partitioner.orig_g, e);
        int to = get_target(partitioner.orig_g, e);
        const auto &from_ver = get_vertex_props(partitioner.orig_g, from);
        const auto &to_ver = get_vertex_props(partitioner.orig_g, to);
        if (!is_channel(from_ver) || !is_channel(to_ver)) {
            assert(!has_edge(combined_g, from, to));
            add_edge(combined_g, from, to);
        }
    }

    printf("Combined/Orig graph has %d/%d (%g) edges.\n", num_edges(combined_g), num_edges(partitioner.orig_g), 100.0*num_edges(combined_g)/num_edges(partitioner.orig_g));

    //goto lol;

    //extern int num_types;
    //extern struct s_type_descriptor *type_descriptors;
    //extern int nx, ny;
    //extern struct s_grid_tile **grid;

    //free_rr_graph();

    //int warnings;

    //build_rr_graph(GRAPH_UNIDIR, num_types, type_descriptors, nx, ny, grid,
            //opts->fixed_channel_width, NULL, det_routing_arch.switch_block_type,
            //det_routing_arch.Fs, det_routing_arch.num_segment,
            //det_routing_arch.num_switch, segment_inf,
            //det_routing_arch.global_route_switch,
            //det_routing_arch.delayless_switch, timing_inf,
            //det_routing_arch.wire_to_ipin_switch, opts->base_cost_type,
            //directs, num_directs, FALSE,
            //&warnings);

    //RRGraph channel_with_interior_g;
    //init_channel_only_graph(channel_with_interior_g);

    //dump_rr_graph(channel_with_interior_g, "/Volumes/DATA/rr_graph_channel_with_interior.txt");

    //RRGraph orig_g;
    //init_graph(orig_g);

    //dump_rr_graph(orig_g, "/Volumes/DATA/rr_graph.txt");

    //free_rr_graph();
    //for (int i = 0; i < det_routing_arch.num_segment; ++i) {
        //for (int j = 0; j < segment_inf[i].sb_len; ++j) {
            //if (j != 0 && j != segment_inf[i].sb_len-1) {
                //segment_inf[i].sb[j] = FALSE;
            //}
        //}
    //}
    //build_rr_graph(GRAPH_UNIDIR, num_types, type_descriptors, nx, ny, grid,
            //opts->fixed_channel_width, NULL, det_routing_arch.switch_block_type,
            //det_routing_arch.Fs, det_routing_arch.num_segment,
            //det_routing_arch.num_switch, segment_inf,
            //det_routing_arch.global_route_switch,
            //det_routing_arch.delayless_switch, timing_inf,
            //det_routing_arch.wire_to_ipin_switch, opts->base_cost_type,
            //directs, num_directs, FALSE,
            //&warnings);

    //RRGraph channel_without_interior_g;
    //init_channel_only_graph(channel_without_interior_g);

    //dump_rr_graph(channel_without_interior_g, "/Volumes/DATA/rr_graph_channel_without_interior.txt");

    //vector<vector<vector<int>>> all_partition_components(opts->num_threads);
    ////for (int x = 0; x < nx+1; ++x) {
        ////for (int y = 0; y < ny+1; ++y) {
            ////init_partitioned_graph_5(opts->num_threads, channel_with_interior_g, channel_without_interior_g, orig_g, graphs, x, y);
        ////}
    ////}
    ////init_partitioned_graph_4(opts->num_threads, channel_with_interior_g, channel_without_interior_g, orig_g, graphs, 14, 14);
    //init_partitioned_graph_4(opts->num_threads, channel_with_interior_g, channel_without_interior_g, orig_g, graphs, all_partition_components);

    //for (int i = 0; i < graphs.size(); ++i) {
        //char filename[256];

        //sprintf(filename, "/Volumes/DATA/rr_graph_%d.txt", i);
        //dump_rr_graph(*graphs[i], filename);

        //sprintf(filename, "/Volumes/DATA/rr_graph_%d_edges.txt", i);
        //dump_edges(*graphs[i], filename);
    //}

    vector<net_t> nets;
    vector<net_t> global_nets;
    init_nets(nets, global_nets, opts->bb_factor, opts->large_bb);	

    vector<net_t *> nets_ptr(nets.size());
    for (int i = 0; i < nets.size(); ++i) {
        nets_ptr[i] = &nets[i];
    }
    extern s_bb *route_bb;
    sort(nets_ptr.begin(), nets_ptr.end(), [] (const net_t *a, const net_t *b) -> bool {
            return get_bounding_box_area(route_bb[a->vpr_id]) > get_bounding_box_area(route_bb[b->vpr_id]);
            });
    int rank = 0;
    for (auto &net : nets_ptr) {
        net->bb_area_rank = rank++;
    }

	vector<route_state_t *> states;
	congestion_locked_t *congestion;
	vector<route_tree_t> route_trees;
	t_net_timing *net_timing;
	init_route_structs_locked(partitioner.orig_g, nets, global_nets, opts->num_threads, states, &congestion, route_trees, &net_timing);

    route_parameters_t params;
    params.criticality_exp = opts->criticality_exp;
    params.astar_fac = opts->astar_fac;
    params.max_criticality = opts->max_criticality;
    params.bend_cost = opts->bend_cost;
	params.pres_fac = opts->first_iter_pres_fac; /* Typically 0 -> ignore cong. */

    char buffer[256];

    vector<pair<box, net_t *>> nets_to_route;
    vector<pair<box, net_t *>> nets_to_partition;
    //vector<vector<int>> overlaps;
    vector<vector<pair<int, int>>> overlaps;
    vector<vector<int>> partitions;
    vector<bool> has_interpartition_overlap;
    
    for (auto &net : nets) {
        box b = bg::make_inverse<box>();

        bg::expand(b, point(net.source.x, net.source.y));
        for (const auto &sink : net.sinks) {
            bg::expand(b, point(sink.x, sink.y));
        }
        bg::subtract_value(b.min_corner(), 1+opts->bb_factor);
        bg::add_value(b.max_corner(), opts->bb_factor);

        nets_to_route.push_back(make_pair(b, &net));
    }
    std::sort(begin(nets_to_route), end(nets_to_route), [] (const pair<box, net_t *> &a, const pair<box, net_t *> &b) -> bool {
            return a.second->sinks.size() > b.second->sinks.size();
            });

    bool routed = false;

    vector<perf_t> perfs(opts->num_threads);
    vector<lock_perf_t> lock_perfs(opts->num_threads);
    vector<clock::time_point> greedy_end_time(opts->num_threads);
    vector<clock::time_point> partitioned_end_time(opts->num_threads);
    clock::duration total_greedy_route_time = clock::duration::zero();
    clock::duration total_partitioned_route_time = clock::duration::zero();
    clock::duration total_update_cost_time = clock::duration::zero();
    clock::duration total_partitioning_time = clock::duration::zero();
    clock::duration total_analyze_timing_time = clock::duration::zero();
    clock::duration total_iter_time = clock::duration::zero();

    int iter;
    float crit_path_delay;
	int current_level = 0;

    for (iter = 0; iter < opts->max_router_iterations && !routed; ++iter) {
        clock::duration greedy_route_time = clock::duration::zero();
        clock::duration partitioned_route_time = clock::duration::zero();
        clock::duration update_cost_time = clock::duration::zero();
        clock::duration partitioning_time = clock::duration::zero();
        clock::duration analyze_timing_time = clock::duration::zero();
        clock::duration iter_time = clock::duration::zero();

        //for (int i = 0; i < opts->num_threads; ++i) {
            //sprintf(buffer, LOG_PATH_PREFIX"iter_%d_tid_%d.log", iter, i);
            //fclose(fopen(buffer, "w"));
        //}

        zlog_info(delta_log, "Routing iteration: %d\n", iter);
        printf("Routing iteration: %d\n", iter);
        
        auto iter_start = clock::now();

        //auto route_start = clock::now();

        for (auto &net : nets) {
            update_sink_criticalities(net, net_timing[net.vpr_id], params);
        }

        //tbb::enumerable_thread_specific<state_t *> state_tls;
        //
        tbb::atomic<int> net_index = 0;
        vector<tbb::spin_mutex> debug_lock(opts->num_threads);
        vector<int> thread_num_nets_routed(opts->num_threads);
        vector<int> thread_num_nets_to_route(opts->num_threads);
        vector<int> thread_num_sinks_routed(opts->num_threads);
        vector<int> thread_num_sinks_to_route(opts->num_threads);
        vector<int> thread_bfs_num_sinks_routed(opts->num_threads);
		vector<int> thread_num_interpartition_sinks(opts->num_threads);
        vector<int> graph_used_by_net(nets.size(), -1);
        vector<vector<pseudo_net_t *>> partition_pseudo_nets_0(opts->num_threads);
        vector<vector<pseudo_net_t *>> partition_pseudo_nets_1(opts->num_threads);
        vector<vector<pseudo_net_t *>> *current_partition_pseudo_nets = &partition_pseudo_nets_0;
        vector<vector<pseudo_net_t *>> *new_partition_pseudo_nets = &partition_pseudo_nets_1;
        vector<tbb::spin_mutex> pseudo_nets_locks(opts->num_threads);
        vector<vector<pseudo_net_t *>> net_pseudo_nets(nets.size());
        tbb::atomic<bool> has_interpartition_sinks = false;
        vector<vector<interpartition_sink_t>> interpartition_sinks(nets.size()); 
        vector<int> net_next_pid(nets.size(), -1);
        vector<int> net_initial_pid(nets.size(), -1);

        auto greedy_route_start = clock::now();
		int num_threads = pow(2, partitioner.num_levels-current_level-1);

        tbb::parallel_for(tbb::blocked_range<int>(0, opts->num_threads, 1),
                [&] (const tbb::blocked_range<int> &range) {

                assert(range.end()-range.begin() == 1);

                int tid = range.begin();

                char local_buffer[256];

                sprintf(local_buffer, "%d", iter);
                zlog_put_mdc("iter", local_buffer);

                sprintf(local_buffer, "%d", tid);
                zlog_put_mdc("tid", local_buffer);

                //assert(debug_lock[tid].try_lock());
                
                perf_t local_perf;
                local_perf.num_heap_pushes = 0;
                local_perf.num_heap_pops = 0;
                local_perf.num_neighbor_visits = 0;

                lock_perf_t local_lock_perf;
                local_lock_perf.num_lock_tries = 0;
                local_lock_perf.num_lock_waits = 0;
                local_lock_perf.total_wait_time = clock::duration::zero();

                int local_num_nets_to_route = 0;
                int local_num_nets_routed = 0;
                int local_num_sinks_routed = 0;
                int local_num_sinks_to_route = 0;
                int local_bfs_num_sinks_routed = 0;
				int local_num_interpartition_sinks = 0;

				if (tid < num_threads) {
					int i = tid;
					//while ((i = net_index++) < nets_to_route.size()) {
					while (i < nets_to_route.size()) {
						net_t *net = nets_to_route[i].second;

						//auto rip_up_start = clock::now();
						//if (greedy_rip_up_all) {
						route_tree_mark_all_nodes_to_be_ripped(route_trees[net->local_id], *graphs[tid]);
						//} else {
						//route_tree_mark_congested_nodes_to_be_ripped(route_trees[net->local_id], *graphs[tid], congestion);
						//}
						route_tree_rip_up_marked(route_trees[net->local_id], *graphs[tid], congestion, params.pres_fac, true, &local_lock_perf);

						//local_perf.total_rip_up_time += clock::now()-rip_up_start;

						//auto route_start = clock::now();

						vector<sink_t *> sinks;	
						get_sinks_to_route(net, route_trees[net->local_id], sinks);

						if (!sinks.empty()) {
							route_net_with_high_interpartition_cost(partitioner.orig_g, partitioner.result_pid_by_level[current_level], tid, net->vpr_id, &net->source, sinks, params, states[tid], congestion, route_trees[net->local_id], net_timing[net->vpr_id], interpartition_sinks[net->local_id], true, &local_perf, &local_lock_perf);

							if (!has_interpartition_sinks) {
								has_interpartition_sinks = !interpartition_sinks[net->local_id].empty();
							}

							local_num_interpartition_sinks += interpartition_sinks[net->local_id].size();
							++local_num_nets_routed;
							++local_num_nets_to_route;
							local_num_sinks_to_route += sinks.size();

							local_num_sinks_routed += sinks.size();
						}

						//local_perf.total_route_time += clock::now()-rip_up_start;
						i += num_threads;
					}
				}

                greedy_end_time[tid] = clock::now();

                perfs[tid] = local_perf;
                lock_perfs[tid] = local_lock_perf; 
                thread_num_nets_routed[tid] = local_num_nets_routed;
                thread_num_nets_to_route[tid] = local_num_nets_to_route;
                thread_num_sinks_routed[tid] = local_num_sinks_routed;
                thread_num_sinks_to_route[tid] = local_num_sinks_to_route;
                thread_bfs_num_sinks_routed[tid] = local_bfs_num_sinks_routed;
                thread_num_interpartition_sinks[tid] = local_num_interpartition_sinks;

                //debug_lock[tid].unlock();
                });

        greedy_route_time = clock::now()-greedy_route_start;

        auto partitioned_route_start = clock::now();

		//if (has_interpartition_sinks) {
			//for (const auto &net_interpartition_sinks : interpartition_sinks) {
				//for (const auto &is : net_interpartition_sinks) {
					//vector<RRNode> added_nodes;
					//for (const auto &node : is.path) {
						//if (node.update_cost) {
							//added_nodes.push_back(node.rr_node_id);
						//}
					//}
					//update_one_cost(partitioner.orig_g, congestion, added_nodes.begin(), added_nodes.end(), 1, params.pres_fac, false, nullptr);
				//}
			//}
		//}

        partitioned_route_time = clock::now()-partitioned_route_start;

        //if (greedy_rip_up_all) {
            //next_greedy_rip_up_iter += greedy_rip_up_all_period;
            //++greedy_rip_up_all_period;
            //prev_num_overused_nodes = std::numeric_limits<unsigned long>::max();
        //}

        //route_time = clock::now()-route_start;

        iter_time = clock::now()-iter_start;

        int total_num_nets_routed = 0;
        int total_num_nets_to_route = 0;
        int total_num_sinks_routed = 0;
        int total_num_sinks_to_route = 0;
        int total_bfs_num_sinks_routed = 0;
		int total_num_interpartition_sinks = 0;
        for (int i = 0; i < opts->num_threads; ++i) {
            total_num_nets_to_route += thread_num_nets_to_route[i];
            total_num_nets_routed += thread_num_nets_routed[i];
            total_num_sinks_routed += thread_num_sinks_routed[i];
            total_num_sinks_to_route += thread_num_sinks_to_route[i];
            total_bfs_num_sinks_routed += thread_bfs_num_sinks_routed[i];
			total_num_interpartition_sinks += thread_num_interpartition_sinks[i];
        }

		printf("num nets routed: ");
		for (int i = 0; i < opts->num_threads; ++i) {
			printf("%d/%d (%g), ", thread_num_nets_routed[i], thread_num_nets_to_route[i], thread_num_nets_routed[i]*100.0/thread_num_nets_to_route[i]);
		}
		printf("\n");

		printf("num sinks routed: ");
		for (int i = 0; i < opts->num_threads; ++i) {
			printf("%d/%d (%g), ", thread_num_sinks_routed[i], thread_num_sinks_to_route[i], thread_num_sinks_routed[i]*100.0/thread_num_sinks_to_route[i]);
		}
		printf("\n");

        printf("Total num nets routed: %d/%d (%g)\n", total_num_nets_routed, total_num_nets_to_route, total_num_nets_routed*100.0/total_num_nets_to_route);
        printf("Total num interpartition sinks routed: %d/%d (%g)\n", total_num_interpartition_sinks, total_num_sinks_to_route, total_num_interpartition_sinks*100.0/total_num_sinks_to_route);
        printf("Total num sinks routed: %d/%d (%g)\n", total_bfs_num_sinks_routed+total_num_sinks_routed, total_num_sinks_to_route, (total_bfs_num_sinks_routed+total_num_sinks_routed)*100.0/total_num_sinks_to_route);

		/* checking */
        for (int i = 0; i < num_vertices(*graphs[0]); ++i) {
            congestion[i].cong.recalc_occ = 0; 
        }

        for (const auto &net : nets) {
            check_route_tree(route_trees[net.local_id], net, *graphs[0]);
            recalculate_occ(route_trees[net.local_id], *graphs[0], congestion);
        }

        bool valid = true;
        for (int i = 0; i < num_vertices(*graphs[0]); ++i) {
            sprintf_rr_node(i, buffer);
            if (congestion[i].cong.recalc_occ != congestion[i].cong.occ) {
                zlog_error(delta_log, "Node %s occ mismatch, recalc: %d original: %d\n", buffer, congestion[i].cong.recalc_occ, congestion[i].cong.occ);
                valid = false;
            }
        }
        assert(valid);

        int overused_total_bb_rank = 0;
        int num_congested_nets = 0;
        for (const auto &net : nets) {
            vector<int> overused_rr_node;
            assert(route_trees[net.local_id].root_rt_nodes.size() == 1);
            get_overused_nodes(route_trees[net.local_id], route_trees[net.local_id].root_rt_nodes[0], *graphs[0], congestion, overused_rr_node);
            if (!overused_rr_node.empty()) {
                zlog_level(delta_log, ROUTER_V1, "Net %d bb_rank %d has %lu overused nodes:\n", net.vpr_id, net.bb_area_rank, overused_rr_node.size());
                for (const auto &item : overused_rr_node) {
                    zlog_level(delta_log, ROUTER_V1, "%d ", item);
                }
                zlog_level(delta_log, ROUTER_V1, "\n");
                overused_total_bb_rank += net.bb_area_rank;
                ++num_congested_nets;
            }
        }

        zlog_level(delta_log, ROUTER_V1, "Average overused net bb rank: %g\n", (float)overused_total_bb_rank/num_congested_nets);
        zlog_level(delta_log, ROUTER_V1, "Num congested nets: %d\n", num_congested_nets);

        iter_start = clock::now();

        if (feasible_routing(partitioner.orig_g, congestion)) {
            //dump_route(*current_traces_ptr, "route.txt");
            routed = true;
        } else {
            unsigned long num_overused_nodes = 0;
            for (int i = 0; i < num_vertices(*graphs[0]); ++i) {
                if (congestion[i].cong.occ > get_vertex_props(*graphs[0], i).capacity) {
                    ++num_overused_nodes;
                }
            }
            zlog_info(delta_log, "Num overused nodes: %lu/%d (%.2f)\n", num_overused_nodes, num_vertices(*graphs[0]), num_overused_nodes*100.0/num_vertices(*graphs[0]));
            printf("Num overused nodes: %lu/%d (%.2f)\n", num_overused_nodes, num_vertices(*graphs[0]), num_overused_nodes*100.0/num_vertices(*graphs[0]));

			//if ((float)total_num_interpartition_sinks/total_num_sinks_to_route > 0.1) {
				//++current_level;
				//assert(current_level < partitioner.result_pid_by_level.size());
			//}

            auto update_cost_start = clock::now();

            if (iter == 0) {
                params.pres_fac = opts->initial_pres_fac;
                update_costs(*graphs[0], congestion, params.pres_fac, 0);
            } else {
                params.pres_fac *= opts->pres_fac_mult;

				/* Avoid overflow for high iteration counts, even if acc_cost is big */
                params.pres_fac = std::min(params.pres_fac, static_cast<float>(HUGE_POSITIVE_FLOAT / 1e5));

                update_costs(*graphs[0], congestion, params.pres_fac, opts->acc_fac);
            }

            update_cost_time = clock::now()-update_cost_start;
        }

        auto analyze_timing_start = clock::now();

        crit_path_delay = analyze_timing(net_timing);

        analyze_timing_time = clock::now()-analyze_timing_start;

        iter_time += clock::now()-iter_start;

        //total_route_time += route_time;
        total_greedy_route_time += greedy_route_time;
        total_partitioned_route_time += partitioned_route_time;
        total_update_cost_time += update_cost_time;
        total_partitioning_time += partitioning_time;
        total_analyze_timing_time += analyze_timing_time;
        total_iter_time += iter_time;

        printf("Iteration time: %g s.\n", duration_cast<nanoseconds>(iter_time).count() / 1e9);
            printf("\tRoute time: %g s.\n", duration_cast<nanoseconds>(greedy_route_time+partitioned_route_time).count() / 1e9);
                printf("\t\tGreedy route time: %g s (%g).\n", duration_cast<nanoseconds>(greedy_route_time).count() / 1e9, 100.0*duration_cast<nanoseconds>(greedy_route_time).count()/duration_cast<nanoseconds>(greedy_route_time+partitioned_route_time).count());
                printf("\t\tPartitioned route time: %g s (%g).\n", duration_cast<nanoseconds>(partitioned_route_time).count() / 1e9, 100.0*duration_cast<nanoseconds>(partitioned_route_time).count()/duration_cast<nanoseconds>(greedy_route_time+partitioned_route_time).count());
            printf("\tUpdate cost time: %g s.\n", duration_cast<nanoseconds>(update_cost_time).count() / 1e9);
            printf("\tPartitioning time: %g s.\n", duration_cast<nanoseconds>(partitioning_time).count() / 1e9);
            printf("\tAnalyze timing time: %g s.\n", duration_cast<nanoseconds>(analyze_timing_time).count() / 1e9);
        printf("Critical path: %g ns\n", crit_path_delay);

        unsigned long total_num_heap_pushes = 0;
        unsigned long total_num_heap_pops = 0;
        unsigned long total_num_neighbor_visits = 0;

        for (int i = 0; i < opts->num_threads; ++i) {
            total_num_heap_pushes += perfs[i].num_heap_pushes;
            total_num_heap_pops += perfs[i].num_heap_pops;
            total_num_neighbor_visits += perfs[i].num_neighbor_visits;
        }

        printf("num lock waits/tries: ");
        for (int i = 0; i < opts->num_threads; ++i) {
			printf("%lu/%lu (%g) ", lock_perfs[i].num_lock_waits, lock_perfs[i].num_lock_tries, (float)lock_perfs[i].num_lock_waits*100/lock_perfs[i].num_lock_tries); 
		}
		printf("\n");

		printf("total wait time: ");
		for (int i = 0; i < opts->num_threads; ++i) {
			printf("%g (%g) ", duration_cast<nanoseconds>(lock_perfs[i].total_wait_time).count() / 1e9, 100.0*lock_perfs[i].total_wait_time/(greedy_route_time+partitioned_route_time));
		}
		printf("\n");

		printf("num_heap_pushes: ");
		for (int i = 0; i < opts->num_threads; ++i) {
			printf("%lu ", perfs[i].num_heap_pushes);
		}
		printf("\n");

		printf("num_heap_pops: ");
		for (int i = 0; i < opts->num_threads; ++i) {
			printf("%lu ", perfs[i].num_heap_pops);
		}
		printf("\n");

		printf("num_neighbor_visits: ");
		for (int i = 0; i < opts->num_threads; ++i) {
			printf("%lu ", perfs[i].num_neighbor_visits);
		}
		printf("\n");

        printf("total_num_heap_pushes: %lu\n", total_num_heap_pushes);
        printf("total_num_heap_pops: %lu\n", total_num_heap_pops); 
        printf("total_num_neighbor_visits: %lu\n", total_num_neighbor_visits);

        clock::time_point greedy_earliest_end_time = *std::min_element(begin(greedy_end_time), begin(greedy_end_time)+num_threads);

        printf("greedy wait time: ");
        for (int i = 0; i < opts->num_threads; ++i) {
			printf("%g (%g), ", duration_cast<nanoseconds>(greedy_end_time[i]-greedy_earliest_end_time).count() / 1e9, 100.0*(greedy_end_time[i]-greedy_earliest_end_time)/greedy_route_time);
        }
		printf("\n");
    }

    //sprintf(buffer, "%s_run_%d_rr_graph_occ.txt", s_circuit_name, run);
    //dump_rr_graph_occ(congestion, num_vertices(g), buffer);

    if (routed) {
        printf("Routed in %d iterations. Total iteration time: %g\n", iter, duration_cast<nanoseconds>(total_iter_time).count() / 1e9);
            printf("\tTotal route time: %g s.\n", duration_cast<nanoseconds>(total_greedy_route_time+total_partitioned_route_time).count() / 1e9);
                printf("\t\tTotal greedy route time: %g s.\n", duration_cast<nanoseconds>(total_greedy_route_time).count() / 1e9);
                printf("\t\tTotal partitioned route time: %g s.\n", duration_cast<nanoseconds>(total_partitioned_route_time).count() / 1e9);
            printf("\tTotal update cost time: %g s.\n", duration_cast<nanoseconds>(total_update_cost_time).count() / 1e9);
            printf("\tTotal partitioning time: %g s.\n", duration_cast<nanoseconds>(total_partitioning_time).count() / 1e9);
            printf("\tTotal analyze timing time: %g s.\n", duration_cast<nanoseconds>(total_analyze_timing_time).count() / 1e9);

        printf("Final critical path: %g ns\n", crit_path_delay);
    } else {
        printf("Failed to route in %d iterations. Total iteration time: %g\n", opts->max_router_iterations, duration_cast<nanoseconds>(total_iter_time).count() / 1e9);
            printf("\tTotal route time: %g s.\n", duration_cast<nanoseconds>(total_greedy_route_time+total_partitioned_route_time).count() / 1e9);
                printf("\t\tTotal greedy route time: %g s.\n", duration_cast<nanoseconds>(total_greedy_route_time).count() / 1e9);
                printf("\t\tTotal partitioned route time: %g s.\n", duration_cast<nanoseconds>(total_partitioned_route_time).count() / 1e9);
            printf("\tTotal update cost time: %g s.\n", duration_cast<nanoseconds>(total_update_cost_time).count() / 1e9);
            printf("\tTotal partitioning time: %g s.\n", duration_cast<nanoseconds>(total_partitioning_time).count() / 1e9);
            printf("\tTotal analyze timing time: %g s.\n", duration_cast<nanoseconds>(total_analyze_timing_time).count() / 1e9);
    }

    //delete_graph(g);
    delete_net_timing(nets, global_nets, net_timing);	
    delete [] congestion;
    delete [] net_timing;
    for (int i = 0; i < opts->num_threads; ++i) {
        delete [] states[i];
    }

    return routed;
}

bool spatial_route_2(t_router_opts *opts, struct s_det_routing_arch det_routing_arch, t_direct_inf *directs, int num_directs, t_segment_inf *segment_inf, t_timing_inf timing_inf)
{
	//using std::chrono::duration_cast;
	//using std::chrono::nanoseconds;
	//using clock = std::chrono::high_resolution_clock;

	//init_logging();
	//zlog_set_record("custom_output", delta_log_output);
	//zlog_set_record("missing_edge", missing_edge_log_output);
	//zlog_set_record("ss", ss_log_output);
	//delta_log_files.resize(opts->max_router_iterations);
	//for (int i = 0; i < opts->max_router_iterations; ++i) {
		//delta_log_files[i].resize(opts->num_threads, nullptr);
	//}
	//missing_edge_log_files.resize(opts->max_router_iterations);
	//for (int i = 0; i < opts->max_router_iterations; ++i) {
		//missing_edge_log_files[i].resize(opts->num_threads, nullptr);
	//}

	////fclose(fopen(LOG_PATH_PREFIX"iter__tid_.log", "w"));

	//test_fast_graph();
	//test_topo();
	//test_fm();
	//test_filter_graph();
	//test_partition_graph();
	//test_connected_components();

	//vector<RRGraph *> graphs;
	//rr_graph_partitioner partitioner;
	//partitioner.init(opts, &det_routing_arch, directs, num_directs, segment_inf, &timing_inf);
	//partitioner.partition(opts->num_threads, graphs);
	////for (const auto &g : graphs) {
		////routability(*g);
	////}
	////
	//RRGraph combined_g;
	//add_vertex(combined_g, num_vertices(partitioner.orig_g));

	//for (const auto &g : graphs) {
		//for (const auto &e : get_edges(*g)) {
			//int from = get_source(*g, e);
			//int to = get_target(*g, e);
			//const auto &from_ver = get_vertex_props(*g, from);
			//const auto &to_ver = get_vertex_props(*g, to);

			//if (is_channel(from_ver) && is_channel(to_ver)) {
				//assert(!has_edge(combined_g, from, to));
				//add_edge(combined_g, from, to);
			//}
		//}
	//}

	//for (const auto &e : get_edges(partitioner.orig_g)) {
		//int from = get_source(partitioner.orig_g, e);
		//int to = get_target(partitioner.orig_g, e);
		//const auto &from_ver = get_vertex_props(partitioner.orig_g, from);
		//const auto &to_ver = get_vertex_props(partitioner.orig_g, to);
		//if (!is_channel(from_ver) || !is_channel(to_ver)) {
			//assert(!has_edge(combined_g, from, to));
			//add_edge(combined_g, from, to);
		//}
	//}

	//printf("Combined/Orig graph has %d/%d (%g) edges.\n", num_edges(combined_g), num_edges(partitioner.orig_g), 100.0*num_edges(combined_g)/num_edges(partitioner.orig_g));

	////goto lol;

	////extern int num_types;
	////extern struct s_type_descriptor *type_descriptors;
	////extern int nx, ny;
	////extern struct s_grid_tile **grid;

	////free_rr_graph();

	////int warnings;

	////build_rr_graph(GRAPH_UNIDIR, num_types, type_descriptors, nx, ny, grid,
			////opts->fixed_channel_width, NULL, det_routing_arch.switch_block_type,
			////det_routing_arch.Fs, det_routing_arch.num_segment,
			////det_routing_arch.num_switch, segment_inf,
			////det_routing_arch.global_route_switch,
			////det_routing_arch.delayless_switch, timing_inf,
			////det_routing_arch.wire_to_ipin_switch, opts->base_cost_type,
			////directs, num_directs, FALSE,
			////&warnings);

	////RRGraph channel_with_interior_g;
	////init_channel_only_graph(channel_with_interior_g);

	////dump_rr_graph(channel_with_interior_g, "/Volumes/DATA/rr_graph_channel_with_interior.txt");

	////RRGraph orig_g;
	////init_graph(orig_g);

	////dump_rr_graph(orig_g, "/Volumes/DATA/rr_graph.txt");

	////free_rr_graph();
	////for (int i = 0; i < det_routing_arch.num_segment; ++i) {
		////for (int j = 0; j < segment_inf[i].sb_len; ++j) {
			////if (j != 0 && j != segment_inf[i].sb_len-1) {
				////segment_inf[i].sb[j] = FALSE;
			////}
		////}
	////}
	////build_rr_graph(GRAPH_UNIDIR, num_types, type_descriptors, nx, ny, grid,
			////opts->fixed_channel_width, NULL, det_routing_arch.switch_block_type,
			////det_routing_arch.Fs, det_routing_arch.num_segment,
			////det_routing_arch.num_switch, segment_inf,
			////det_routing_arch.global_route_switch,
			////det_routing_arch.delayless_switch, timing_inf,
			////det_routing_arch.wire_to_ipin_switch, opts->base_cost_type,
			////directs, num_directs, FALSE,
			////&warnings);

	////RRGraph channel_without_interior_g;
	////init_channel_only_graph(channel_without_interior_g);

	////dump_rr_graph(channel_without_interior_g, "/Volumes/DATA/rr_graph_channel_without_interior.txt");

	////vector<vector<vector<int>>> all_partition_components(opts->num_threads);
	//////for (int x = 0; x < nx+1; ++x) {
		//////for (int y = 0; y < ny+1; ++y) {
			//////init_partitioned_graph_5(opts->num_threads, channel_with_interior_g, channel_without_interior_g, orig_g, graphs, x, y);
		//////}
	//////}
	//////init_partitioned_graph_4(opts->num_threads, channel_with_interior_g, channel_without_interior_g, orig_g, graphs, 14, 14);
	////init_partitioned_graph_4(opts->num_threads, channel_with_interior_g, channel_without_interior_g, orig_g, graphs, all_partition_components);

	////for (int i = 0; i < graphs.size(); ++i) {
		////char filename[256];

		////sprintf(filename, "/Volumes/DATA/rr_graph_%d.txt", i);
		////dump_rr_graph(*graphs[i], filename);

		////sprintf(filename, "/Volumes/DATA/rr_graph_%d_edges.txt", i);
		////dump_edges(*graphs[i], filename);
	////}

	//vector<net_t> nets;
	//vector<net_t> global_nets;
	//init_nets(nets, global_nets, opts->bb_factor);	

	//vector<net_t *> nets_ptr(nets.size());
	//for (int i = 0; i < nets.size(); ++i) {
		//nets_ptr[i] = &nets[i];
	//}
	//extern s_bb *route_bb;
	//sort(nets_ptr.begin(), nets_ptr.end(), [] (const net_t *a, const net_t *b) -> bool {
			//return get_bounding_box_area(route_bb[a->vpr_id]) > get_bounding_box_area(route_bb[b->vpr_id]);
			//});
	//int rank = 0;
	//for (auto &net : nets_ptr) {
		//net->bb_area_rank = rank++;
	//}

	//t_net_timing *net_timing = new t_net_timing[nets.size()+global_nets.size()];
	//init_net_timing(nets, global_nets, net_timing);

	//vector<route_tree_t> route_trees(nets.size());
	//for (int i = 0; i < nets.size(); ++i) {
		//route_tree_init(route_trees[i]);
	//}

	//vector<route_state_t *> states(opts->num_threads);

	//for (int i = 0; i < opts->num_threads; ++i) {
		//states[i] = new route_state_t[num_vertices(*graphs[0])];
		//for (int j = 0; j < num_vertices(*graphs[0]); ++j) {
			//states[i][j].rr_node = -1;
			//states[i][j].known_cost = std::numeric_limits<float>::max();
			//states[i][j].cost = std::numeric_limits<float>::max();
			//states[i][j].prev_edge = RRGraph::null_edge();
			//states[i][j].upstream_R = -1;
			//states[i][j].delay = std::numeric_limits<float>::max();
		//}
	//}

	//congestion_t *congestion = new congestion_t[num_vertices(*graphs[0])];
	//for (int i = 0; i < num_vertices(*graphs[0]); ++i) {
		//congestion[i].acc_cost = 1;
		//congestion[i].pres_cost = 1;
		//congestion[i].occ = 0;
	//}

	//congestion_t *temp_congestion = new congestion_t[num_vertices(*graphs[0])];
	//for (int i = 0; i < num_vertices(*graphs[0]); ++i) {
		//temp_congestion[i].acc_cost = 1;
		//temp_congestion[i].pres_cost = 1;
		//temp_congestion[i].occ = 0;
	//}

	//route_parameters_t params;
	//params.criticality_exp = opts->criticality_exp;
	//params.astar_fac = opts->astar_fac;
	//params.max_criticality = opts->max_criticality;
	//params.bend_cost = opts->bend_cost;
	//params.pres_fac = opts->first_iter_pres_fac; [> Typically 0 -> ignore cong. <]

	//char buffer[256];

	//vector<pair<box, net_t *>> nets_to_route;
	//vector<pair<box, net_t *>> nets_to_partition;
	////vector<vector<int>> overlaps;
	//vector<vector<pair<int, int>>> overlaps;
	//vector<vector<int>> partitions;
	//vector<bool> has_interpartition_overlap;
	
	//for (auto &net : nets) {
		//box b = bg::make_inverse<box>();

		//bg::expand(b, point(net.source.x, net.source.y));
		//for (const auto &sink : net.sinks) {
			//bg::expand(b, point(sink.x, sink.y));
		//}
		//bg::subtract_value(b.min_corner(), 1+opts->bb_factor);
		//bg::add_value(b.max_corner(), opts->bb_factor);

		//nets_to_route.push_back(make_pair(b, &net));
	//}
	//std::sort(begin(nets_to_route), end(nets_to_route), [] (const pair<box, net_t *> &a, const pair<box, net_t *> &b) -> bool {
			//return a.second->sinks.size() > b.second->sinks.size();
			//});

	//bool routed = false;

	//vector<perf_t> perfs(opts->num_threads);
	//vector<lock_perf_t> lock_perfs(opts->num_threads);
	//vector<clock::time_point> greedy_end_time(opts->num_threads);
	//vector<clock::time_point> partitioned_end_time(opts->num_threads);
	//clock::duration total_greedy_route_time = clock::duration::zero();
	//clock::duration total_partitioned_route_time = clock::duration::zero();
	//clock::duration total_update_cost_time = clock::duration::zero();
	//clock::duration total_partitioning_time = clock::duration::zero();
	//clock::duration total_analyze_timing_time = clock::duration::zero();
	//clock::duration total_iter_time = clock::duration::zero();

	//int iter;
	//float crit_path_delay;


	//for (iter = 0; iter < opts->max_router_iterations && !routed; ++iter) {
		//clock::duration greedy_route_time = clock::duration::zero();
		//clock::duration partitioned_route_time = clock::duration::zero();
		//clock::duration update_cost_time = clock::duration::zero();
		//clock::duration partitioning_time = clock::duration::zero();
		//clock::duration analyze_timing_time = clock::duration::zero();
		//clock::duration iter_time = clock::duration::zero();

		////for (int i = 0; i < opts->num_threads; ++i) {
			////sprintf(buffer, LOG_PATH_PREFIX"iter_%d_tid_%d.log", iter, i);
			////fclose(fopen(buffer, "w"));
		////}

		//zlog_info(delta_log, "Routing iteration: %d\n", iter);
		//printf("Routing iteration: %d\n", iter);
		
		//auto iter_start = clock::now();

		////auto route_start = clock::now();

		//for (auto &net : nets) {
			//update_sink_criticalities(net, net_timing[net.vpr_id], params);
		//}

		////tbb::enumerable_thread_specific<state_t *> state_tls;
		////
		//tbb::atomic<int> net_index = 0;
		//vector<tbb::spin_mutex> debug_lock(opts->num_threads);
		//vector<int> thread_num_nets_routed(opts->num_threads);
		//vector<int> thread_num_nets_to_route(opts->num_threads);
		//vector<int> thread_num_sinks_routed(opts->num_threads);
		//vector<int> thread_num_sinks_to_route(opts->num_threads);
		//vector<int> thread_bfs_num_sinks_routed(opts->num_threads);
		//vector<int> graph_used_by_net(nets.size(), -1);
		//vector<vector<pseudo_net_t *>> partition_pseudo_nets_0(opts->num_threads);
		//vector<vector<pseudo_net_t *>> partition_pseudo_nets_1(opts->num_threads);
		//vector<vector<pseudo_net_t *>> *current_partition_pseudo_nets = &partition_pseudo_nets_0;
		//vector<vector<pseudo_net_t *>> *new_partition_pseudo_nets = &partition_pseudo_nets_1;
		//vector<tbb::spin_mutex> pseudo_nets_locks(opts->num_threads);
		//vector<vector<pseudo_net_t *>> net_pseudo_nets(nets.size());
		//tbb::atomic<bool> has_unrouted_sinks = false;
		//vector<unrouted_t> unrouted(nets.size()); 
		//vector<int> net_next_pid(nets.size(), -1);
		//vector<int> net_initial_pid(nets.size(), -1);

		//auto greedy_route_start = clock::now();

		//tbb::parallel_for(tbb::blocked_range<int>(0, opts->num_threads, 1),
				//[&] (const tbb::blocked_range<int> &range) {

				//assert(range.end()-range.begin() == 1);

				//int tid = range.begin();

				//char local_buffer[256];

				//sprintf(local_buffer, "%d", iter);
				//zlog_put_mdc("iter", local_buffer);

				//sprintf(local_buffer, "%d", tid);
				//zlog_put_mdc("tid", local_buffer);

				////assert(debug_lock[tid].try_lock());
				
				//perf_t local_perf;
				//local_perf.num_heap_pushes = 0;
				//local_perf.num_heap_pops = 0;
				//local_perf.num_neighbor_visits = 0;

				//lock_perf_t local_lock_perf;
				//local_lock_perf.num_lock_tries = 0;
				//local_lock_perf.num_lock_waits = 0;
				//local_lock_perf.total_wait_time = clock::duration::zero();

				//int local_num_nets_to_route = 0;
				//int local_num_nets_routed = 0;
				//int local_num_sinks_routed = 0;
				//int local_num_sinks_to_route = 0;
				//int local_bfs_num_sinks_routed = 0;

				//int i = tid;
				////while ((i = net_index++) < nets_to_route.size()) {
				//while (i < nets_to_route.size()) {
					//net_t *net = nets_to_route[i].second;

					//net_initial_pid[net->local_id] = tid;
					//net_next_pid[net->local_id] = (tid+1) % opts->num_threads;

					////auto rip_up_start = clock::now();
					////if (greedy_rip_up_all) {
						//route_tree_mark_all_nodes_to_be_ripped(route_trees[net->local_id], *graphs[tid]);
					////} else {
						////route_tree_mark_congested_nodes_to_be_ripped(route_trees[net->local_id], *graphs[tid], congestion);
					////}
					//route_tree_rip_up_marked(route_trees[net->local_id], *graphs[tid], congestion, params.pres_fac, true, &local_lock_perf);

					////local_perf.total_rip_up_time += clock::now()-rip_up_start;

					////auto route_start = clock::now();

					//vector<sink_t *> sinks;	
					//sinks.reserve(net->sinks.size());
					//for (auto &sink : net->sinks) {
						//RouteTreeNode sink_rt_node = route_tree_get_rt_node(route_trees[net->local_id], sink.rr_node);
						//if (sink_rt_node == RouteTree::null_vertex()) {
							//sinks.push_back(&sink);
						//} else {
							//assert(!get_vertex_props(route_trees[net->local_id].graph, sink_rt_node).pending_rip_up);
						//}
					//}

					//if (!sinks.empty()) {
						////route_net_4(partitioner.orig_g, partitioner.result_pid, tid, net->vpr_id, &net->source, sinks, params, states[tid], congestion, route_trees[net->local_id], net_timing[net->vpr_id], unrouted[net->local_id], opts->num_threads, true, &local_perf, &local_lock_perf);

						//if (unrouted[net->local_id].unrouted_sinks.empty()) {
							//++local_num_nets_routed;
						//} else {
							//auto pnet = get_pseudo_net(unrouted[net->local_id], net, partitioner.orig_g, partitioner.result_pid, net_next_pid[net->local_id], opts->num_threads);
							//int to_pid = net_next_pid[net->local_id];
							//net_next_pid[net->local_id] = (net_next_pid[net->local_id]+1) % opts->num_threads;

							//assert(pnet);

							//if (pnet) {
								//pseudo_nets_locks[to_pid].lock();
								//(*current_partition_pseudo_nets)[to_pid].push_back(pnet);
								//pseudo_nets_locks[to_pid].unlock();

								//net_pseudo_nets[net->local_id].push_back(pnet);

								//zlog_level(delta_log, ROUTER_V3, "Generated pseudo net [pid = %d net = %d %lu sources %lu/%lu sinks]\n", to_pid, pnet->net->vpr_id, pnet->pseudo_sources.size(), pnet->pseudo_sinks.size(), pnet->net->sinks.size());

								//for (const auto &psink : pnet->pseudo_sinks) {
									//zlog_level(delta_log, ROUTER_V3, "\tSink: %d\n", psink.sink->rr_node);
								//}

							//} else {
								//zlog_level(delta_log, ROUTER_V3, "Failed to generate pnet for [pid = %d net = %d]\n", to_pid, net->vpr_id);
							//}

							//has_unrouted_sinks = true;
						//}

						//++local_num_nets_to_route;
						//local_num_sinks_to_route += sinks.size();

						//local_num_sinks_routed += sinks.size()-unrouted[net->local_id].unrouted_sinks.size();
					//}

					////local_perf.total_route_time += clock::now()-rip_up_start;
					//i += opts->num_threads;
				//}

				//greedy_end_time[tid] = clock::now();

				//perfs[tid] = local_perf;
				//lock_perfs[tid] = local_lock_perf; 
				//thread_num_nets_routed[tid] = local_num_nets_routed;
				//thread_num_nets_to_route[tid] = local_num_nets_to_route;
				//thread_num_sinks_routed[tid] = local_num_sinks_routed;
				//thread_num_sinks_to_route[tid] = local_num_sinks_to_route;
				//thread_bfs_num_sinks_routed[tid] = local_bfs_num_sinks_routed;

				////debug_lock[tid].unlock();
				//});

		//greedy_route_time = clock::now()-greedy_route_start;

		//auto partitioned_route_start = clock::now();

		//while (has_unrouted_sinks) {
			//has_unrouted_sinks = false;

			//for (int pid = 0; pid < opts->num_threads; ++pid) {
				//for (int inet = 0; inet < (*current_partition_pseudo_nets)[pid].size(); ++inet) {
					//pseudo_net_t *pnet = (*current_partition_pseudo_nets)[pid][inet];

					//assert(!pnet->pseudo_sources.empty());
					//assert(!pnet->pseudo_sinks.empty());

					//zlog_level(delta_log, ROUTER_V1, "-- Routing pseudo net [pid = %d, initial pid = %d net = %d, %lu sources, %lu sinks]\n", pid, net_initial_pid[pnet->net->local_id], pnet->net->vpr_id, pnet->pseudo_sources.size(), pnet->pseudo_sinks.size());
					//zlog_level(delta_log, ROUTER_V1, "-- Pseudo sources:  ");

					//route_tree_t rt;
					//route_tree_init(rt);
					//for (const auto &s : pnet->pseudo_sources) {
						//zlog_level(delta_log, ROUTER_V1, "%d ", s);
						//RouteTreeNode rt_node = route_tree_add_rr_node(rt, s.node, *graphs[pid]);
						//const auto &rr_node = get_vertex_props(*graphs[pid], s.node);
						//route_tree_set_node_properties(get_vertex_props(rt.graph, rt_node), true, RRGraph::null_edge(), rr_node.R, 0.5 * rr_node.R * rr_node.C);
						//route_tree_add_root(rt, s.node);
					//}
					//zlog_level(delta_log, ROUTER_V1, "\n");

					//unrouted_t local_unrouted;

					//vector<sink_t *> local_sinks;
					//for (const auto &psink : pnet->pseudo_sinks) {
						//if (psink.path.empty()) {
							//local_sinks.push_back(psink.sink);
						//}
					//}

					//route_net_3(partitioner.orig_g, partitioner.result_pid, pid, pnet->net->vpr_id, nullptr, local_sinks, params, states[pid], congestion, rt, net_timing[pnet->net->vpr_id], local_unrouted, opts->num_threads, true, nullptr, nullptr);

					//int num_paths_connecting_to_source = 0;
					//vector<pseudo_sink_t *> routed_pseudo_sinks;
					//for (auto &psink : pnet->pseudo_sinks) {
						//bool routed = psink.path.empty() && find_if(begin(local_unrouted.unrouted_sinks), end(local_unrouted.unrouted_sinks), [&psink] (const auto &sink) -> bool { return sink->rr_node == psink.sink->rr_node; }) == end(local_unrouted.unrouted_sinks);

						//if (routed) {
							//routed_pseudo_sinks.push_back(&psink);

							//zlog_level(delta_log, ROUTER_V3, "-- Merging route tree for sink %d: %d\n", psink.sink->id, psink.sink->rr_node);

							//const auto &path_to_sink = route_tree_get_path(rt, psink.sink->rr_node);

							//auto psource = find_if(begin(pnet->pseudo_sources), end(pnet->pseudo_sources), [&path_to_sink] (const auto &psource) -> bool { return psource.node == path_to_sink.back().rr_node_id; });

							//if (psource != end(pnet->pseudo_sources)) {
								//++num_paths_connecting_to_source;
							//}

							//if (psource != end(pnet->pseudo_sources)) {
								//assert(psink.pseudo_source == nullptr);
								//psink.pseudo_source = &(*psource);

								//const auto &bn = get_source(partitioner.orig_g, psource->prev_edge);

								////auto iter2 = find_if(begin(pnet->pseudo_sources), end(pnet->pseudo_sources), [&bn, &partitioner] (const auto &psource) -> bool { return get_source(partitioner.orig_g, psource.prev_edge) == bn; });
								////assert(iter2 != end(pnet->pseudo_sources));

								////int num_matches = 0;
								////for (const auto &psink : pnet->pseudo_sinks) {
									////for (const auto &bnode : psink.->boundary_nodes) {
										////if (bnode.rr_node == bn) {
											////++num_matches;
										////}
									////}
								////}
								////assert(num_matches > 0);

								//int counts = count_if(begin(pnet->pseudo_sources), end(pnet->pseudo_sources), [&partitioner, &bn] (const auto &psource) -> bool { return get_source(partitioner.orig_g, psource.prev_edge) == bn; });
								//assert(counts > 0);

								//zlog_level(delta_log, ROUTER_V3, "-- Path starts from pseudo source %d\n", psource->node);

								//[> it is possible that there is already a path to the boundary node in the route tree <]
								//RouteTreeNode bn_rt_node = route_tree_get_rt_node(route_trees[pnet->net->local_id], bn);
								//if (bn_rt_node == RouteTree::null_vertex()) {
									//unrouted_t dummy;
									//sink_t bn_sink;
									//bn_sink.id = -1;
									//bn_sink.rr_node = bn;
									//bn_sink.criticality_fac = 1;
									//bn_sink.current_bounding_box = psink.sink->current_bounding_box;

									//zlog_level(delta_log, ROUTER_V3, "-- Routing to boundary node: %d\n", bn);

									//route_net_3(partitioner.orig_g, partitioner.result_pid, net_initial_pid[pnet->net->local_id], pnet->net->vpr_id, nullptr, { &bn_sink }, params, states[pid], congestion, route_trees[pnet->net->local_id], net_timing[pnet->net->vpr_id], dummy, opts->num_threads, true, nullptr, nullptr);

									//assert(dummy.unrouted_sinks.empty());

									////const auto &bn_path = psource->bnode->path;
									////route_tree_add_path(route_trees[pnet->net->local_id], bn_path, partitioner.orig_g, nullptr, false);

									////vector<RRNode> added_nodes;
									////for (const auto &n : bn_path) {
										////if (n.update_cost) {
											////added_nodes.push_back(n.rr_node_id);
										////}
									////}
									////update_one_cost(*graphs[pid], congestion, added_nodes.begin(), added_nodes.end(), 1, params.pres_fac, false, nullptr);

									//assert(route_tree_get_rt_node(route_trees[pnet->net->local_id], bn) != RouteTree::null_vertex());

									////path_node_t bnode;

									////bnode.rr_node_id = bn;
									////bnode.prev_edge = RRGraph::null_edge();
									////bnode.update_cost = false;

									////path_to_sink.push_back(bnode);
								//}
							//}
							
							//psink.path = path_to_sink;
						//}
					//}

					//for (const auto &psink : routed_pseudo_sinks) {
						//if (psink->pseudo_source) {
							//RouteTreeNode ps_rt_node = route_tree_get_rt_node(route_trees[pnet->net->local_id], psink->pseudo_source->node);
							//if (ps_rt_node == RouteTree::null_vertex()) {
								//zlog_level(delta_log, ROUTER_V3, "-- Adding pseudo source %d to route tree\n", psink->pseudo_source->node);

								//const auto &bn = get_source(partitioner.orig_g, psink->pseudo_source->prev_edge);

								//ps_rt_node = route_tree_add_rr_node(route_trees[pnet->net->local_id], psink->pseudo_source->node, *graphs[pid]);
								//const auto &ps_rr_node_p = get_vertex_props(*graphs[pid], psink->pseudo_source->node);
								//route_tree_set_node_properties(get_vertex_props(route_trees[pnet->net->local_id].graph, ps_rt_node), true, RRGraph::null_edge(), ps_rr_node_p.R, 0.5 * ps_rr_node_p.R * ps_rr_node_p.C);
								//route_tree_add_edge_between_rr_node(route_trees[pnet->net->local_id], bn, psink->pseudo_source->node);

								//assert(ps_rr_node_p.type != SOURCE);

								//update_one_cost_internal(psink->pseudo_source->node, ps_rr_node_p, congestion[psink->pseudo_source->node], 1, params.pres_fac, true, nullptr);
							//}
						//}

						//zlog_level(delta_log, ROUTER_V3, "-- Adding path to sink\n");
						//route_tree_add_path(route_trees[pnet->net->local_id], psink->path, partitioner.orig_g, nullptr);
					//}

					//if (!local_unrouted.unrouted_sinks.empty()) {
						////auto new_pnet = get_pseudo_net(unrouted_sinks[pnet->net->local_id], pnet->net, partitioner.orig_g, partitioner.result_pid, net_next_pid[pnet->net->local_id], opts->num_threads);
						//update_pseudo_sources(*pnet, unrouted[pnet->net->local_id], partitioner.orig_g, partitioner.result_pid, opts->num_threads, net_next_pid[pnet->net->local_id]);

						//int to_pid = net_next_pid[pnet->net->local_id];
						//assert(to_pid != net_initial_pid[pnet->net->local_id]);
						//net_next_pid[pnet->net->local_id] = (net_next_pid[pnet->net->local_id]+1) % opts->num_threads;

						//(*new_partition_pseudo_nets)[to_pid].push_back(pnet);

						//zlog_level(delta_log, ROUTER_V3, "Updating pseudo net [pid = %d, net = %d, %lu sources, %lu sinks]\n", to_pid, pnet->net->vpr_id, pnet->pseudo_sources.size(), pnet->pseudo_sinks.size());
						//for (const auto &psink : pnet->pseudo_sinks) {
							//zlog_level(delta_log, ROUTER_V3, "\t%d\n", psink.sink->rr_node);
						//}

						//has_unrouted_sinks = true;
					//} else {
						//zlog_level(delta_log, ROUTER_V3, "Routed all sinks of net %d\n", pnet->net->vpr_id);
					//}

				//}
				//(*current_partition_pseudo_nets)[pid].clear();
			//}
			
			//std::swap(current_partition_pseudo_nets, new_partition_pseudo_nets);
		//}
		 
		//partitioned_route_time = clock::now()-partitioned_route_start;

		//int lol = 0;
		//for (const auto &pnets : net_pseudo_nets) {
			//if (!pnets.empty()) {
				//assert(pnets.size() == 1);
				//const auto &pnet = pnets[0];
				//for (const auto &psink : pnet->pseudo_sinks) {
					//if (psink.path.empty()) {
						//++lol;
					//} else {
					//}
				//}
			//}
		//}

		//printf("Num unrouted sinks: %d\n", lol);
		//assert(lol == 0);

		////if (greedy_rip_up_all) {
			////next_greedy_rip_up_iter += greedy_rip_up_all_period;
			////++greedy_rip_up_all_period;
			////prev_num_overused_nodes = std::numeric_limits<unsigned long>::max();
		////}

		////route_time = clock::now()-route_start;

		//iter_time = clock::now()-iter_start;

		//int total_num_nets_routed = 0;
		//int total_num_nets_to_route = 0;
		//int total_num_sinks_routed = 0;
		//int total_num_sinks_to_route = 0;
		//int total_bfs_num_sinks_routed = 0;
		//for (int i = 0; i < opts->num_threads; ++i) {
			//printf("Thread %d num nets routed: %d/%d (%g)\n", i, thread_num_nets_routed[i], thread_num_nets_to_route[i], thread_num_nets_routed[i]*100.0/thread_num_nets_to_route[i]);
			//printf("Thread %d num sinks routed: %d/%d (%g)\n", i, thread_num_sinks_routed[i], thread_num_sinks_to_route[i], thread_num_sinks_routed[i]*100.0/thread_num_sinks_to_route[i]);
			//total_num_nets_to_route += thread_num_nets_to_route[i];
			//total_num_nets_routed += thread_num_nets_routed[i];
			//total_num_sinks_routed += thread_num_sinks_routed[i];
			//total_num_sinks_to_route += thread_num_sinks_to_route[i];
			//total_bfs_num_sinks_routed += thread_bfs_num_sinks_routed[i];
		//}
		//printf("Total num nets routed: %d/%d (%g)\n", total_num_nets_routed, total_num_nets_to_route, total_num_nets_routed*100.0/total_num_nets_to_route);
		//printf("Total partition num sinks routed: %d/%d (%g)\n", total_num_sinks_routed, total_num_sinks_to_route, total_num_sinks_routed*100.0/total_num_sinks_to_route);
		//printf("Total BFS num sinks routed: %d/%d (%g)\n", total_bfs_num_sinks_routed, total_num_sinks_to_route, total_bfs_num_sinks_routed*100.0/total_num_sinks_to_route);

		//printf("Total num sinks routed: %d/%d (%g)\n", total_bfs_num_sinks_routed+total_num_sinks_routed, total_num_sinks_to_route, (total_bfs_num_sinks_routed+total_num_sinks_routed)*100.0/total_num_sinks_to_route);

		//[> checking <]
		//for (int i = 0; i < num_vertices(*graphs[0]); ++i) {
			//congestion[i].recalc_occ = 0; 
		//}

		//for (const auto &net : nets) {
			//check_route_tree(route_trees[net.local_id], net, *graphs[0]);
			//recalculate_occ(route_trees[net.local_id], *graphs[0], congestion);
		//}

		//bool valid = true;
		//for (int i = 0; i < num_vertices(*graphs[0]); ++i) {
			//sprintf_rr_node(i, buffer);
			//if (congestion[i].recalc_occ != congestion[i].occ) {
				//zlog_error(delta_log, "Node %s occ mismatch, recalc: %d original: %d\n", buffer, congestion[i].recalc_occ, congestion[i].occ);
				//valid = false;
			//}
		//}
		//assert(valid);

		//int overused_total_bb_rank = 0;
		//int num_congested_nets = 0;
		//for (const auto &net : nets) {
			//vector<int> overused_rr_node;
			//assert(route_trees[net.local_id].root_rt_nodes.size() == 1);
			//get_overused_nodes(route_trees[net.local_id], route_trees[net.local_id].root_rt_nodes[0], *graphs[0], congestion, overused_rr_node);
			//if (!overused_rr_node.empty()) {
				//zlog_level(delta_log, ROUTER_V1, "Net %d bb_rank %d has %lu overused nodes:\n", net.vpr_id, net.bb_area_rank, overused_rr_node.size());
				//for (const auto &item : overused_rr_node) {
					//zlog_level(delta_log, ROUTER_V1, "%d ", item);
				//}
				//zlog_level(delta_log, ROUTER_V1, "\n");
				//overused_total_bb_rank += net.bb_area_rank;
				//++num_congested_nets;
			//}
		//}

		//zlog_level(delta_log, ROUTER_V1, "Average overused net bb rank: %g\n", (float)overused_total_bb_rank/num_congested_nets);
		//zlog_level(delta_log, ROUTER_V1, "Num congested nets: %d\n", num_congested_nets);

		//iter_start = clock::now();

		//if (feasible_routing(partitioner.orig_g, congestion)) {
			////dump_route(*current_traces_ptr, "route.txt");
			//routed = true;
		//} else {
			//unsigned long num_overused_nodes = 0;
			//for (int i = 0; i < num_vertices(*graphs[0]); ++i) {
				//if (congestion[i].occ > get_vertex_props(*graphs[0], i).capacity) {
					//++num_overused_nodes;
				//}
			//}
			//zlog_info(delta_log, "Num overused nodes: %lu/%d (%.2f)\n", num_overused_nodes, num_vertices(*graphs[0]), num_overused_nodes*100.0/num_vertices(*graphs[0]));
			//printf("Num overused nodes: %lu/%d (%.2f)\n", num_overused_nodes, num_vertices(*graphs[0]), num_overused_nodes*100.0/num_vertices(*graphs[0]));

			//auto update_cost_start = clock::now();

			//if (iter == 0) {
				//params.pres_fac = opts->initial_pres_fac;
				//update_costs(*graphs[0], congestion, params.pres_fac, 0);
			//} else {
				//params.pres_fac *= opts->pres_fac_mult;

				//[> Avoid overflow for high iteration counts, even if acc_cost is big <]
				//params.pres_fac = std::min(params.pres_fac, static_cast<float>(HUGE_POSITIVE_FLOAT / 1e5));

				//update_costs(*graphs[0], congestion, params.pres_fac, opts->acc_fac);
			//}

			//update_cost_time = clock::now()-update_cost_start;
		//}

		//auto analyze_timing_start = clock::now();

		//crit_path_delay = analyze_timing(net_timing);

		//analyze_timing_time = clock::now()-analyze_timing_start;

		//iter_time += clock::now()-iter_start;

		////total_route_time += route_time;
		//total_greedy_route_time += greedy_route_time;
		//total_partitioned_route_time += partitioned_route_time;
		//total_update_cost_time += update_cost_time;
		//total_partitioning_time += partitioning_time;
		//total_analyze_timing_time += analyze_timing_time;
		//total_iter_time += iter_time;

		//printf("Iteration time: %g s.\n", duration_cast<nanoseconds>(iter_time).count() / 1e9);
			//printf("\tRoute time: %g s.\n", duration_cast<nanoseconds>(greedy_route_time+partitioned_route_time).count() / 1e9);
				//printf("\t\tGreedy route time: %g s (%g).\n", duration_cast<nanoseconds>(greedy_route_time).count() / 1e9, 100.0*duration_cast<nanoseconds>(greedy_route_time).count()/duration_cast<nanoseconds>(greedy_route_time+partitioned_route_time).count());
				//printf("\t\tPartitioned route time: %g s (%g).\n", duration_cast<nanoseconds>(partitioned_route_time).count() / 1e9, 100.0*duration_cast<nanoseconds>(partitioned_route_time).count()/duration_cast<nanoseconds>(greedy_route_time+partitioned_route_time).count());
			//printf("\tUpdate cost time: %g s.\n", duration_cast<nanoseconds>(update_cost_time).count() / 1e9);
			//printf("\tPartitioning time: %g s.\n", duration_cast<nanoseconds>(partitioning_time).count() / 1e9);
			//printf("\tAnalyze timing time: %g s.\n", duration_cast<nanoseconds>(analyze_timing_time).count() / 1e9);
		//printf("Critical path: %g ns\n", crit_path_delay);

		//unsigned long total_num_heap_pushes = 0;
		//unsigned long total_num_heap_pops = 0;
		//unsigned long total_num_neighbor_visits = 0;

		//for (int i = 0; i < opts->num_threads; ++i) {
			//total_num_heap_pushes += perfs[i].num_heap_pushes;
			//total_num_heap_pops += perfs[i].num_heap_pops;
			//total_num_neighbor_visits += perfs[i].num_neighbor_visits;
		//}

		//for (int i = 0; i < opts->num_threads; ++i) {
			//printf("Thread %d num lock waits/tries: %lu/%lu (%g)\n", i, lock_perfs[i].num_lock_waits, lock_perfs[i].num_lock_tries, (float)lock_perfs[i].num_lock_waits*100/lock_perfs[i].num_lock_tries);
			//printf("Thread %d total wait time: %g (%g)\n", i, duration_cast<nanoseconds>(lock_perfs[i].total_wait_time).count() / 1e9, 100.0*lock_perfs[i].total_wait_time/(greedy_route_time+partitioned_route_time));
			//printf("Thread %d num_heap_pushes: %lu\n", i, perfs[i].num_heap_pushes);
			//printf("Thread %d num_heap_pops: %lu\n", i, perfs[i].num_heap_pops);
			//printf("Thread %d num_neighbor_visits: %lu\n", i, perfs[i].num_neighbor_visits);
		//}

		//printf("total_num_heap_pushes: %lu\n", total_num_heap_pushes);
		//printf("total_num_heap_pops: %lu\n", total_num_heap_pops); 
		//printf("total_num_neighbor_visits: %lu\n", total_num_neighbor_visits);

		//clock::time_point greedy_earliest_end_time = *std::min_element(begin(greedy_end_time), end(greedy_end_time));
		//clock::time_point partitioned_earliest_end_time = *std::min_element(begin(partitioned_end_time), end(partitioned_end_time));

		//for (int i = 0; i < opts->num_threads; ++i) {
			//printf("Thread %d greedy wait time %g (%g)\n", i, duration_cast<nanoseconds>(greedy_end_time[i]-greedy_earliest_end_time).count() / 1e9, 100.0*(greedy_end_time[i]-greedy_earliest_end_time)/greedy_route_time);
		//}
		//for (int i = 0; i < opts->num_threads; ++i) {
			//if (partitioned_route_time > clock::duration::zero()) {
				//printf("Thread %d partitioned wait time %g (%g)\n", i, duration_cast<nanoseconds>(partitioned_end_time[i]-partitioned_earliest_end_time).count() / 1e9, 100.0*(partitioned_end_time[i]-partitioned_earliest_end_time)/partitioned_route_time);
			//}
		//}
	//}

	////sprintf(buffer, "%s_run_%d_rr_graph_occ.txt", s_circuit_name, run);
	////dump_rr_graph_occ(congestion, num_vertices(g), buffer);

	//if (routed) {
		//printf("Routed in %d iterations. Total iteration time: %g\n", iter, duration_cast<nanoseconds>(total_iter_time).count() / 1e9);
			//printf("\tTotal route time: %g s.\n", duration_cast<nanoseconds>(total_greedy_route_time+total_partitioned_route_time).count() / 1e9);
				//printf("\t\tTotal greedy route time: %g s.\n", duration_cast<nanoseconds>(total_greedy_route_time).count() / 1e9);
				//printf("\t\tTotal partitioned route time: %g s.\n", duration_cast<nanoseconds>(total_partitioned_route_time).count() / 1e9);
			//printf("\tTotal update cost time: %g s.\n", duration_cast<nanoseconds>(total_update_cost_time).count() / 1e9);
			//printf("\tTotal partitioning time: %g s.\n", duration_cast<nanoseconds>(total_partitioning_time).count() / 1e9);
			//printf("\tTotal analyze timing time: %g s.\n", duration_cast<nanoseconds>(total_analyze_timing_time).count() / 1e9);

		//printf("Final critical path: %g ns\n", crit_path_delay);
	//} else {
		//printf("Failed to route in %d iterations. Total iteration time: %g\n", opts->max_router_iterations, duration_cast<nanoseconds>(total_iter_time).count() / 1e9);
			//printf("\tTotal route time: %g s.\n", duration_cast<nanoseconds>(total_greedy_route_time+total_partitioned_route_time).count() / 1e9);
				//printf("\t\tTotal greedy route time: %g s.\n", duration_cast<nanoseconds>(total_greedy_route_time).count() / 1e9);
				//printf("\t\tTotal partitioned route time: %g s.\n", duration_cast<nanoseconds>(total_partitioned_route_time).count() / 1e9);
			//printf("\tTotal update cost time: %g s.\n", duration_cast<nanoseconds>(total_update_cost_time).count() / 1e9);
			//printf("\tTotal partitioning time: %g s.\n", duration_cast<nanoseconds>(total_partitioning_time).count() / 1e9);
			//printf("\tTotal analyze timing time: %g s.\n", duration_cast<nanoseconds>(total_analyze_timing_time).count() / 1e9);
	//}

	////delete_graph(g);
	//delete_net_timing(nets, global_nets, net_timing);	
	//delete [] congestion;
	//delete [] net_timing;
	//for (int i = 0; i < opts->num_threads; ++i) {
		//delete [] states[i];
	//}

	//return routed;
}

bool spatial_route_3(t_router_opts *opts, struct s_det_routing_arch det_routing_arch, t_direct_inf *directs, int num_directs, t_segment_inf *segment_inf, t_timing_inf timing_inf)
{
	//using std::chrono::duration_cast;
	//using std::chrono::nanoseconds;
	//using clock = std::chrono::high_resolution_clock;

	//init_logging();
	//zlog_set_record("custom_output", delta_log_output);
	//zlog_set_record("missing_edge", missing_edge_log_output);
	//delta_log_files.resize(opts->num_threads, nullptr);
	//missing_edge_log_files.resize(opts->num_threads, nullptr);

	////fclose(fopen(LOG_PATH_PREFIX"iter__tid_.log", "w"));

	//test_topo();
	//test_fm();
	//test_filter_graph();
	//test_partition_graph();
	//test_connected_components();

	//vector<RRGraph *> graphs;
	//rr_graph_partitioner partitioner;
	//partitioner.init(opts, &det_routing_arch, directs, num_directs, segment_inf, &timing_inf);
	//partitioner.partition(opts->num_threads, graphs);
	////for (const auto &g : graphs) {
		////routability(*g);
	////}
	////
	//RRGraph combined_g;
	//add_vertex(combined_g, num_vertices(partitioner.orig_g));

	//for (const auto &g : graphs) {
		//for (const auto &e : get_edges(*g)) {
			//int from = get_source(*g, e);
			//int to = get_target(*g, e);
			//const auto &from_ver = get_vertex_props(*g, from);
			//const auto &to_ver = get_vertex_props(*g, to);

			//if (is_channel(from_ver) && is_channel(to_ver)) {
				//assert(!has_edge(combined_g, from, to));
				//add_edge(combined_g, from, to);
			//}
		//}
	//}

	//for (const auto &e : get_edges(partitioner.orig_g)) {
		//int from = get_source(partitioner.orig_g, e);
		//int to = get_target(partitioner.orig_g, e);
		//const auto &from_ver = get_vertex_props(partitioner.orig_g, from);
		//const auto &to_ver = get_vertex_props(partitioner.orig_g, to);
		//if (!is_channel(from_ver) || !is_channel(to_ver)) {
			//assert(!has_edge(combined_g, from, to));
			//add_edge(combined_g, from, to);
		//}
	//}

	//printf("Combined/Orig graph has %d/%d (%g) edges.\n", num_edges(combined_g), num_edges(partitioner.orig_g), 100.0*num_edges(combined_g)/num_edges(partitioner.orig_g));

	////goto lol;

	////extern int num_types;
	////extern struct s_type_descriptor *type_descriptors;
	////extern int nx, ny;
	////extern struct s_grid_tile **grid;

	////free_rr_graph();

	////int warnings;

	////build_rr_graph(GRAPH_UNIDIR, num_types, type_descriptors, nx, ny, grid,
			////opts->fixed_channel_width, NULL, det_routing_arch.switch_block_type,
			////det_routing_arch.Fs, det_routing_arch.num_segment,
			////det_routing_arch.num_switch, segment_inf,
			////det_routing_arch.global_route_switch,
			////det_routing_arch.delayless_switch, timing_inf,
			////det_routing_arch.wire_to_ipin_switch, opts->base_cost_type,
			////directs, num_directs, FALSE,
			////&warnings);

	////RRGraph channel_with_interior_g;
	////init_channel_only_graph(channel_with_interior_g);

	////dump_rr_graph(channel_with_interior_g, "/Volumes/DATA/rr_graph_channel_with_interior.txt");

	////RRGraph orig_g;
	////init_graph(orig_g);

	////dump_rr_graph(orig_g, "/Volumes/DATA/rr_graph.txt");

	////free_rr_graph();
	////for (int i = 0; i < det_routing_arch.num_segment; ++i) {
		////for (int j = 0; j < segment_inf[i].sb_len; ++j) {
			////if (j != 0 && j != segment_inf[i].sb_len-1) {
				////segment_inf[i].sb[j] = FALSE;
			////}
		////}
	////}
	////build_rr_graph(GRAPH_UNIDIR, num_types, type_descriptors, nx, ny, grid,
			////opts->fixed_channel_width, NULL, det_routing_arch.switch_block_type,
			////det_routing_arch.Fs, det_routing_arch.num_segment,
			////det_routing_arch.num_switch, segment_inf,
			////det_routing_arch.global_route_switch,
			////det_routing_arch.delayless_switch, timing_inf,
			////det_routing_arch.wire_to_ipin_switch, opts->base_cost_type,
			////directs, num_directs, FALSE,
			////&warnings);

	////RRGraph channel_without_interior_g;
	////init_channel_only_graph(channel_without_interior_g);

	////dump_rr_graph(channel_without_interior_g, "/Volumes/DATA/rr_graph_channel_without_interior.txt");

	////vector<vector<vector<int>>> all_partition_components(opts->num_threads);
	//////for (int x = 0; x < nx+1; ++x) {
		//////for (int y = 0; y < ny+1; ++y) {
			//////init_partitioned_graph_5(opts->num_threads, channel_with_interior_g, channel_without_interior_g, orig_g, graphs, x, y);
		//////}
	//////}
	//////init_partitioned_graph_4(opts->num_threads, channel_with_interior_g, channel_without_interior_g, orig_g, graphs, 14, 14);
	////init_partitioned_graph_4(opts->num_threads, channel_with_interior_g, channel_without_interior_g, orig_g, graphs, all_partition_components);

	////for (int i = 0; i < graphs.size(); ++i) {
		////char filename[256];

		////sprintf(filename, "/Volumes/DATA/rr_graph_%d.txt", i);
		////dump_rr_graph(*graphs[i], filename);

		////sprintf(filename, "/Volumes/DATA/rr_graph_%d_edges.txt", i);
		////dump_edges(*graphs[i], filename);
	////}

	//vector<net_t> nets;
	//vector<net_t> global_nets;
	//init_nets(nets, global_nets, opts->bb_factor);	

	//t_net_timing *net_timing = new t_net_timing[nets.size()+global_nets.size()];
	//init_net_timing(nets, global_nets, net_timing);

	//vector<route_tree_t> route_trees(nets.size());
	//for (int i = 0; i < nets.size(); ++i) {
		//route_tree_init(route_trees[i]);
	//}

	//vector<route_state_t *> states(opts->num_threads);

	//for (int i = 0; i < opts->num_threads; ++i) {
		//states[i] = new route_state_t[num_vertices(*graphs[0])];
		//for (int j = 0; j < num_vertices(*graphs[0]); ++j) {
			//states[i][j].rr_node = -1;
			//states[i][j].known_cost = std::numeric_limits<float>::max();
			//states[i][j].cost = std::numeric_limits<float>::max();
			//states[i][j].prev_edge = RRGraph::null_edge(); 
			//states[i][j].upstream_R = -1;
			//states[i][j].delay = std::numeric_limits<float>::max();
		//}
	//}

	//congestion_t *congestion = new congestion_t[num_vertices(*graphs[0])];
	//for (int i = 0; i < num_vertices(*graphs[0]); ++i) {
		//congestion[i].acc_cost = 1;
		//congestion[i].pres_cost = 1;
		//congestion[i].occ = 0;
	//}

	//congestion_t *temp_congestion = new congestion_t[num_vertices(*graphs[0])];
	//for (int i = 0; i < num_vertices(*graphs[0]); ++i) {
		//temp_congestion[i].acc_cost = 1;
		//temp_congestion[i].pres_cost = 1;
		//temp_congestion[i].occ = 0;
	//}

	//route_parameters_t params;
	//params.criticality_exp = opts->criticality_exp;
	//params.astar_fac = opts->astar_fac;
	//params.max_criticality = opts->max_criticality;
	//params.bend_cost = opts->bend_cost;
	//params.pres_fac = opts->first_iter_pres_fac; [> Typically 0 -> ignore cong. <]

	//char buffer[256];

	//vector<pair<box, net_t *>> nets_to_route;
	//vector<pair<box, net_t *>> nets_to_partition;
	////vector<vector<int>> overlaps;
	//vector<vector<pair<int, int>>> overlaps;
	//vector<vector<int>> partitions;
	//vector<bool> has_interpartition_overlap;
	
	//for (auto &net : nets) {
		//box b = bg::make_inverse<box>();

		//bg::expand(b, point(net.source.x, net.source.y));
		//for (const auto &sink : net.sinks) {
			//bg::expand(b, point(sink.x, sink.y));
		//}
		//bg::subtract_value(b.min_corner(), 1+opts->bb_factor);
		//bg::add_value(b.max_corner(), opts->bb_factor);

		//nets_to_route.push_back(make_pair(b, &net));
	//}
	//std::sort(begin(nets_to_route), end(nets_to_route), [] (const pair<box, net_t *> &a, const pair<box, net_t *> &b) -> bool {
			//return a.second->sinks.size() > b.second->sinks.size();
			//});

	//bool routed = false;

	//vector<perf_t> perfs(opts->num_threads);
	//vector<lock_perf_t> lock_perfs(opts->num_threads);
	//vector<clock::time_point> greedy_end_time(opts->num_threads);
	//vector<clock::time_point> partitioned_end_time(opts->num_threads);
	//clock::duration total_greedy_route_time = clock::duration::zero();
	//clock::duration total_partitioned_route_time = clock::duration::zero();
	//clock::duration total_update_cost_time = clock::duration::zero();
	//clock::duration total_partitioning_time = clock::duration::zero();
	//clock::duration total_analyze_timing_time = clock::duration::zero();
	//clock::duration total_iter_time = clock::duration::zero();

	//int iter;
	//float crit_path_delay;


	//for (iter = 0; iter < opts->max_router_iterations && !routed; ++iter) {
		//clock::duration greedy_route_time = clock::duration::zero();
		//clock::duration partitioned_route_time = clock::duration::zero();
		//clock::duration update_cost_time = clock::duration::zero();
		//clock::duration partitioning_time = clock::duration::zero();
		//clock::duration analyze_timing_time = clock::duration::zero();
		//clock::duration iter_time = clock::duration::zero();

		////for (int i = 0; i < opts->num_threads; ++i) {
			////sprintf(buffer, LOG_PATH_PREFIX"iter_%d_tid_%d.log", iter, i);
			////fclose(fopen(buffer, "w"));
		////}

		//zlog_info(delta_log, "Routing iteration: %d\n", iter);
		//printf("Routing iteration: %d\n", iter);
		
		//auto iter_start = clock::now();

		////auto route_start = clock::now();

		//for (auto &net : nets) {
			//update_sink_criticalities(net, net_timing[net.vpr_id], params);
		//}

		////tbb::enumerable_thread_specific<state_t *> state_tls;
		////
		//tbb::atomic<int> net_index = 0;
		//vector<tbb::spin_mutex> debug_lock(opts->num_threads);
		//vector<int> thread_num_nets_routed(opts->num_threads);
		//vector<int> thread_num_nets_to_route(opts->num_threads);
		//vector<int> thread_num_sinks_routed(opts->num_threads);
		//vector<int> thread_num_sinks_to_route(opts->num_threads);
		//vector<int> thread_bfs_num_sinks_routed(opts->num_threads);
		//vector<int> graph_used_by_net(nets.size(), -1);

		//auto greedy_route_start = clock::now();

		//tbb::parallel_for(tbb::blocked_range<int>(0, opts->num_threads, 1),
				//[&] (const tbb::blocked_range<int> &range) {

				//assert(range.end()-range.begin() == 1);

				//int tid = range.begin();

				//char local_buffer[256];

				//sprintf(local_buffer, "%d", iter);
				//zlog_put_mdc("iter", local_buffer);

				//sprintf(local_buffer, "%d", tid);
				//zlog_put_mdc("tid", local_buffer);

				////assert(debug_lock[tid].try_lock());
				
				//perf_t local_perf;
				//local_perf.num_heap_pushes = 0;
				//local_perf.num_heap_pops = 0;
				//local_perf.num_neighbor_visits = 0;

				//lock_perf_t local_lock_perf;
				//local_lock_perf.num_lock_tries = 0;
				//local_lock_perf.num_lock_waits = 0;
				//local_lock_perf.total_wait_time = clock::duration::zero();

				//int local_num_nets_to_route = 0;
				//int local_num_nets_routed = 0;
				//int local_num_sinks_routed = 0;
				//int local_num_sinks_to_route = 0;
				//int local_bfs_num_sinks_routed = 0;

				//int i;
				//while ((i = net_index++) < nets_to_route.size()) {
					//net_t *net = nets_to_route[i].second;

					////auto rip_up_start = clock::now();
					////if (greedy_rip_up_all) {
						////route_tree_mark_all_nodes_to_be_ripped(route_trees[net->local_id], g);
					////} else {
						//route_tree_mark_congested_nodes_to_be_ripped(route_trees[net->local_id], *graphs[tid], congestion);
					////}
					//route_tree_rip_up_marked(route_trees[net->local_id], *graphs[tid], congestion, params.pres_fac, true, &local_lock_perf);

					////local_perf.total_rip_up_time += clock::now()-rip_up_start;

					////auto route_start = clock::now();

					//vector<sink_t *> sinks;	
					//sinks.reserve(net->sinks.size());
					//for (auto &sink : net->sinks) {
						//RouteTreeNode sink_rt_node = route_tree_get_rt_node(route_trees[net->local_id], sink.rr_node);
						//if (!valid(sink_rt_node))  {
							//sinks.push_back(&sink);
						//} else {
							//assert(!get_vertex_props(route_trees[net->local_id].graph, sink_rt_node).pending_rip_up);
						//}
					//}
					//if (!sinks.empty()) {
						//const auto &unrouted_sinks = route_net_2(*graphs[tid], net->vpr_id, &net->source, sinks, params, states[tid], congestion, route_trees[net->local_id], net_timing[net->vpr_id], true, &local_perf, &local_lock_perf);

						//if (unrouted_sinks.size() == 0) {
							//++local_num_nets_routed;
						//} else {
							//route_tree_t temp_route_tree;
							//route_tree_init(temp_route_tree);

							//zlog_level(delta_log, ROUTER_V1, "Routing net %d with original RR graph\n", net->vpr_id);

							//const auto &orig_unrouted_sinks = route_net_2(partitioner.orig_g, net->vpr_id, &net->source, unrouted_sinks, params, states[tid], temp_congestion, temp_route_tree, net_timing[net->vpr_id], true, nullptr, nullptr);

							//assert(orig_unrouted_sinks.size() == 0);

							//int num_missing_edges = 0;
							//int num_missing_interior_edges = 0;
							//zlog_debug(missing_edge_log, "Net %d routed in graph %d requires edges below: \n", net->vpr_id, tid);
							//for (const auto &rt_edge_id : get_edges(temp_route_tree.graph)) {
								//const auto &rr_edge = get_edge_props(temp_route_tree.graph, rt_edge_id).rr_edge;
								//int from = get_source(partitioner.orig_g, rr_edge);
								//int to = get_target(partitioner.orig_g, rr_edge);
								//if (!has_edge(*graphs[tid], from, to)) {
									//int found = -1;
									//int num_found = 0;
									//for (int i = 0; i < opts->num_threads; ++i) {
										//if (i != tid && has_edge(*graphs[i], from, to)) {
											//found = i;
											//++num_found;
										//}
									//}
									//assert(num_found == 1 || num_found == 0);

									//if (has_edge(partitioner.channel_with_interior_g, from, to) && !has_edge(partitioner.channel_without_interior_g, from, to)) {
										//++num_missing_interior_edges;
									//}
									//char buffer[256];

									//sprintf_rr_node(from, buffer);
									//zlog_debug(missing_edge_log, "\t Edge %d %s ->", rr_edge, buffer);

									//sprintf_rr_node(to, buffer);
									//zlog_debug(missing_edge_log, " %s is in graph %d\n", buffer, found);

									//++num_missing_edges;
								//}
							//}
							//zlog_debug(missing_edge_log, "Net %d num missing edges = %d/%d (%g) num missing interior edges = %d/%d (%g)\n\n", net->vpr_id, num_missing_edges, num_edges(temp_route_tree.graph), num_missing_edges*100.0/num_edges(temp_route_tree.graph), num_missing_interior_edges, num_missing_edges, num_missing_interior_edges*100.0/num_missing_edges);

							//if (false) {
							//struct bfs_router_t {
								//set<int> visited_nodes;
								//set<int> visited_edges;
								//vector<int> pred;

								//bfs_router_t(int num_vertices) :
									//pred(num_vertices, -1)
								//{
								//}

								//bool tree_edge(const RREdge &e, const RRGraph &g)
								//{
									//int to = get_target(g, e);
									//assert(pred[to] == -1);
									//pred[to] = get_source(g, e);
									//assert(visited_edges.find(e) == visited_edges.end());
									//visited_edges.insert(e);
									//return true;
								//}

								//void examine_edge(const RREdge &e, const RRGraph &g)
								//{
								//}

								//void discover_vertex(int v, const RRGraph &g)
								//{
									//assert(visited_nodes.find(v) == visited_nodes.end());
									//visited_nodes.insert(v);
								//}

								//void examine_vertex(int v, const RRGraph &g)
								//{
								//}
							//} router_visitor(num_vertices(combined_g));

							//vector<VertexColor> color(num_vertices(combined_g), VertexColor::WHITE);
							//bfs(combined_g, { static_cast<unsigned long>(net->source.rr_node) }, color, router_visitor);
							//assert(router_visitor.visited_edges.size()+1 == router_visitor.visited_nodes.size());
							//for (int i = 0; i < unrouted_sinks.size(); ++i) {
								//if (router_visitor.visited_nodes.find(unrouted_sinks[i]->rr_node) != router_visitor.visited_nodes.end()) {
									//++local_bfs_num_sinks_routed;
								//} else {
									////for (const auto &c : all_partition_components[tid]) {
									////assert(find(begin(c), end(c), net->source.rr_node) == end(c)
									////|| find(begin(c), end(c), sinks[i]->rr_node) == end(c));
									////}
								//}
							//}
							//} else {
								//set<int> routed;
								//for (int t = 0; t < opts->num_threads; ++t) {
									//struct bfs_router_t {
										//set<int> visited_nodes;
										//set<int> visited_edges;
										//vector<int> pred;

										//bfs_router_t(int num_vertices) :
											//pred(num_vertices, -1)
										//{
										//}

										//bool tree_edge(const RREdge &e, const RRGraph &g)
										//{
											//int to = get_target(g, e);
											//assert(pred[to] == -1);
											//pred[to] = get_source(g, e);
											//assert(visited_edges.find(e) == visited_edges.end());
											//visited_edges.insert(e);
											//return true;
										//}

										//void examine_edge(const RREdge &e, const RRGraph &g)
										//{
										//}

										//void discover_vertex(int v, const RRGraph &g)
										//{
											//assert(visited_nodes.find(v) == visited_nodes.end());
											//visited_nodes.insert(v);
										//}

										//void examine_vertex(int v, const RRGraph &g)
										//{
										//}
									//} router_visitor(num_vertices(*graphs[t]));

									//vector<VertexColor> color(num_vertices(*graphs[t]), VertexColor::WHITE);
									//bfs(*graphs[t], { static_cast<unsigned long>(net->source.rr_node) }, color, router_visitor);
									//assert(router_visitor.visited_edges.size()+1 == router_visitor.visited_nodes.size());
									//for (int i = 0; i < unrouted_sinks.size(); ++i) {
										//if (router_visitor.visited_nodes.find(unrouted_sinks[i]->rr_node) != router_visitor.visited_nodes.end()) {
											//if (routed.find(unrouted_sinks[i]->rr_node) == routed.end()) {
												//routed.insert(unrouted_sinks[i]->rr_node);
												//++local_bfs_num_sinks_routed;
											//} 
										//} else {
											////for (const auto &c : all_partition_components[tid]) {
											////assert(find(begin(c), end(c), net->source.rr_node) == end(c)
											////|| find(begin(c), end(c), sinks[i]->rr_node) == end(c));
											////}
										//}
									//}
									//set<int> confirm;
									//for (const auto &u : unrouted_sinks) {
										//confirm.insert(u->rr_node);
									//}
									//assert(confirm.size() >= routed.size());
								//}
							//}
							////assert(local_bfs_num_sinks_routed == sinks.size());
							////assert(local_bfs_num_sinks_routed == num_sinks_routed);

							//graph_used_by_net[net->vpr_id] = tid;
						//}

						//++local_num_nets_to_route;

						//local_num_sinks_to_route += sinks.size();
						//local_num_sinks_routed += sinks.size()-unrouted_sinks.size();
					//}

					////local_perf.total_route_time += clock::now()-rip_up_start;
				//}

				//greedy_end_time[tid] = clock::now();

				//perfs[tid] = local_perf;
				//lock_perfs[tid] = local_lock_perf; 
				//thread_num_nets_routed[tid] = local_num_nets_routed;
				//thread_num_nets_to_route[tid] = local_num_nets_to_route;
				//thread_num_sinks_routed[tid] = local_num_sinks_routed;
				//thread_num_sinks_to_route[tid] = local_num_sinks_to_route;
				//thread_bfs_num_sinks_routed[tid] = local_bfs_num_sinks_routed;

				////debug_lock[tid].unlock();
				//});

		//greedy_route_time = clock::now()-greedy_route_start;

		////if (greedy_rip_up_all) {
			////next_greedy_rip_up_iter += greedy_rip_up_all_period;
			////++greedy_rip_up_all_period;
			////prev_num_overused_nodes = std::numeric_limits<unsigned long>::max();
		////}

		////route_time = clock::now()-route_start;

		//iter_time = clock::now()-iter_start;

		//int total_num_nets_routed = 0;
		//int total_num_nets_to_route = 0;
		//int total_num_sinks_routed = 0;
		//int total_num_sinks_to_route = 0;
		//int total_bfs_num_sinks_routed = 0;
		//for (int i = 0; i < opts->num_threads; ++i) {
			//printf("Thread %d num nets routed: %d/%d (%g)\n", i, thread_num_nets_routed[i], thread_num_nets_to_route[i], thread_num_nets_routed[i]*100.0/thread_num_nets_to_route[i]);
			//printf("Thread %d num sinks routed: %d/%d (%g)\n", i, thread_num_sinks_routed[i], thread_num_sinks_to_route[i], thread_num_sinks_routed[i]*100.0/thread_num_sinks_to_route[i]);
			//total_num_nets_to_route += thread_num_nets_to_route[i];
			//total_num_nets_routed += thread_num_nets_routed[i];
			//total_num_sinks_routed += thread_num_sinks_routed[i];
			//total_num_sinks_to_route += thread_num_sinks_to_route[i];
			//total_bfs_num_sinks_routed += thread_bfs_num_sinks_routed[i];
		//}
		//printf("Total num nets routed: %d/%d (%g)\n", total_num_nets_routed, total_num_nets_to_route, total_num_nets_routed*100.0/total_num_nets_to_route);
		//printf("Total partition num sinks routed: %d/%d (%g)\n", total_num_sinks_routed, total_num_sinks_to_route, total_num_sinks_routed*100.0/total_num_sinks_to_route);
		//printf("Total BFS num sinks routed: %d/%d (%g)\n", total_bfs_num_sinks_routed, total_num_sinks_to_route, total_bfs_num_sinks_routed*100.0/total_num_sinks_to_route);

		//printf("Total num sinks routed: %d/%d (%g)\n", total_bfs_num_sinks_routed+total_num_sinks_routed, total_num_sinks_to_route, (total_bfs_num_sinks_routed+total_num_sinks_routed)*100.0/total_num_sinks_to_route);

		//[> checking <]
		//for (int i = 0; i < num_vertices(*graphs[0]); ++i) {
			//congestion[i].recalc_occ = 0; 
		//}

		//for (const auto &net : nets) {
			//check_route_tree(route_trees[net.local_id], net, *graphs[0]);
			//recalculate_occ(route_trees[net.local_id], *graphs[0], congestion);
		//}

		//bool valid = true;
		//for (int i = 0; i < num_vertices(*graphs[0]); ++i) {
			//sprintf_rr_node(i, buffer);
			//if (congestion[i].recalc_occ != congestion[i].occ) {
				//zlog_error(delta_log, "Node %s occ mismatch, recalc: %d original: %d\n", buffer, congestion[i].recalc_occ, congestion[i].occ);
				//valid = false;
			//}
		//}
		//assert(valid);

		//for (const auto &net : nets) {
			//vector<int> overused_rr_node;
			//get_overused_nodes(route_trees[net.local_id], route_trees[net.local_id].root_rt_node_id), *graphs[0], congestion, overused_rr_node);
			//if (!overused_rr_node.empty()) {
				//zlog_level(delta_log, ROUTER_V1, "Net %d bb_rank %d overused nodes:\n", net.vpr_id, net.bb_area_rank);
				//for (const auto &item : overused_rr_node) {
					//zlog_level(delta_log, ROUTER_V1, "%d ", item);
				//}
				//zlog_level(delta_log, ROUTER_V1, "\n");
			//}
		//}

		//iter_start = clock::now();

		//if (feasible_routing(*graphs[0], congestion)) {
			////dump_route(*current_traces_ptr, "route.txt");
			//routed = true;
		//} else {
			//unsigned long num_overused_nodes = 0;
			//for (int i = 0; i < num_vertices(*graphs[0]); ++i) {
				//if (congestion[i].occ > get_vertex_props(*graphs[0], i).capacity) {
					//++num_overused_nodes;
				//}
			//}
			//zlog_info(delta_log, "Num overused nodes: %lu/%d (%.2f)\n", num_overused_nodes, num_vertices(*graphs[0]), num_overused_nodes*100.0/num_vertices(*graphs[0]));
			//printf("Num overused nodes: %lu/%d (%.2f)\n", num_overused_nodes, num_vertices(*graphs[0]), num_overused_nodes*100.0/num_vertices(*graphs[0]));

			//auto update_cost_start = clock::now();

			//if (iter == 0) {
				//params.pres_fac = opts->initial_pres_fac;
				//update_costs(*graphs[0], congestion, params.pres_fac, 0);
			//} else {
				//params.pres_fac *= opts->pres_fac_mult;

				//[> Avoid overflow for high iteration counts, even if acc_cost is big <]
				//params.pres_fac = std::min(params.pres_fac, static_cast<float>(HUGE_POSITIVE_FLOAT / 1e5));

				//update_costs(*graphs[0], congestion, params.pres_fac, opts->acc_fac);
			//}

			//update_cost_time = clock::now()-update_cost_start;
		//}

		//auto analyze_timing_start = clock::now();

		//crit_path_delay = analyze_timing(net_timing);

		//analyze_timing_time = clock::now()-analyze_timing_start;

		//iter_time += clock::now()-iter_start;

		////total_route_time += route_time;
		//total_greedy_route_time += greedy_route_time;
		//total_partitioned_route_time += partitioned_route_time;
		//total_update_cost_time += update_cost_time;
		//total_partitioning_time += partitioning_time;
		//total_analyze_timing_time += analyze_timing_time;
		//total_iter_time += iter_time;

		//printf("Iteration time: %g s.\n", duration_cast<nanoseconds>(iter_time).count() / 1e9);
			//printf("\tRoute time: %g s.\n", duration_cast<nanoseconds>(greedy_route_time+partitioned_route_time).count() / 1e9);
				//printf("\t\tGreedy route time: %g s.\n", duration_cast<nanoseconds>(greedy_route_time).count() / 1e9);
				//printf("\t\tPartitioned route time: %g s.\n", duration_cast<nanoseconds>(partitioned_route_time).count() / 1e9);
			//printf("\tUpdate cost time: %g s.\n", duration_cast<nanoseconds>(update_cost_time).count() / 1e9);
			//printf("\tPartitioning time: %g s.\n", duration_cast<nanoseconds>(partitioning_time).count() / 1e9);
			//printf("\tAnalyze timing time: %g s.\n", duration_cast<nanoseconds>(analyze_timing_time).count() / 1e9);
		//printf("Critical path: %g ns\n", crit_path_delay);

		//unsigned long total_num_heap_pushes = 0;
		//unsigned long total_num_heap_pops = 0;
		//unsigned long total_num_neighbor_visits = 0;

		//for (int i = 0; i < opts->num_threads; ++i) {
			//total_num_heap_pushes += perfs[i].num_heap_pushes;
			//total_num_heap_pops += perfs[i].num_heap_pops;
			//total_num_neighbor_visits += perfs[i].num_neighbor_visits;
		//}

		//for (int i = 0; i < opts->num_threads; ++i) {
			//printf("Thread %d num lock waits/tries: %lu/%lu (%g)\n", i, lock_perfs[i].num_lock_waits, lock_perfs[i].num_lock_tries, (float)lock_perfs[i].num_lock_waits*100/lock_perfs[i].num_lock_tries);
			//printf("Thread %d total wait time: %g (%g)\n", i, duration_cast<nanoseconds>(lock_perfs[i].total_wait_time).count() / 1e9, 100.0*lock_perfs[i].total_wait_time/(greedy_route_time+partitioned_route_time));
			//printf("Thread %d num_heap_pushes: %lu\n", i, perfs[i].num_heap_pushes);
			//printf("Thread %d num_heap_pops: %lu\n", i, perfs[i].num_heap_pops);
			//printf("Thread %d num_neighbor_visits: %lu\n", i, perfs[i].num_neighbor_visits);
		//}

		//printf("total_num_heap_pushes: %lu\n", total_num_heap_pushes);
		//printf("total_num_heap_pops: %lu\n", total_num_heap_pops); 
		//printf("total_num_neighbor_visits: %lu\n", total_num_neighbor_visits);

		//clock::time_point greedy_earliest_end_time = *std::min_element(begin(greedy_end_time), end(greedy_end_time));
		//clock::time_point partitioned_earliest_end_time = *std::min_element(begin(partitioned_end_time), end(partitioned_end_time));

		//for (int i = 0; i < opts->num_threads; ++i) {
			//printf("Thread %d greedy wait time %g (%g)\n", i, duration_cast<nanoseconds>(greedy_end_time[i]-greedy_earliest_end_time).count() / 1e9, 100.0*(greedy_end_time[i]-greedy_earliest_end_time)/greedy_route_time);
		//}
		//for (int i = 0; i < opts->num_threads; ++i) {
			//if (partitioned_route_time > clock::duration::zero()) {
				//printf("Thread %d partitioned wait time %g (%g)\n", i, duration_cast<nanoseconds>(partitioned_end_time[i]-partitioned_earliest_end_time).count() / 1e9, 100.0*(partitioned_end_time[i]-partitioned_earliest_end_time)/partitioned_route_time);
			//}
		//}
	//}

	////sprintf(buffer, "%s_run_%d_rr_graph_occ.txt", s_circuit_name, run);
	////dump_rr_graph_occ(congestion, num_vertices(g), buffer);

	//if (routed) {
		//printf("Routed in %d iterations. Total iteration time: %g\n", iter, duration_cast<nanoseconds>(total_iter_time).count() / 1e9);
			//printf("\tTotal route time: %g s.\n", duration_cast<nanoseconds>(total_greedy_route_time+total_partitioned_route_time).count() / 1e9);
				//printf("\t\tTotal greedy route time: %g s.\n", duration_cast<nanoseconds>(total_greedy_route_time).count() / 1e9);
				//printf("\t\tTotal partitioned route time: %g s.\n", duration_cast<nanoseconds>(total_partitioned_route_time).count() / 1e9);
			//printf("\tTotal update cost time: %g s.\n", duration_cast<nanoseconds>(total_update_cost_time).count() / 1e9);
			//printf("\tTotal partitioning time: %g s.\n", duration_cast<nanoseconds>(total_partitioning_time).count() / 1e9);
			//printf("\tTotal analyze timing time: %g s.\n", duration_cast<nanoseconds>(total_analyze_timing_time).count() / 1e9);

		//printf("Final critical path: %g ns\n", crit_path_delay);
	//} else {
		//printf("Failed to route in %d iterations. Total iteration time: %g\n", opts->max_router_iterations, duration_cast<nanoseconds>(total_iter_time).count() / 1e9);
			//printf("\tTotal route time: %g s.\n", duration_cast<nanoseconds>(total_greedy_route_time+total_partitioned_route_time).count() / 1e9);
				//printf("\t\tTotal greedy route time: %g s.\n", duration_cast<nanoseconds>(total_greedy_route_time).count() / 1e9);
				//printf("\t\tTotal partitioned route time: %g s.\n", duration_cast<nanoseconds>(total_partitioned_route_time).count() / 1e9);
			//printf("\tTotal update cost time: %g s.\n", duration_cast<nanoseconds>(total_update_cost_time).count() / 1e9);
			//printf("\tTotal partitioning time: %g s.\n", duration_cast<nanoseconds>(total_partitioning_time).count() / 1e9);
			//printf("\tTotal analyze timing time: %g s.\n", duration_cast<nanoseconds>(total_analyze_timing_time).count() / 1e9);
	//}

	////delete_graph(g);
	//delete_net_timing(nets, global_nets, net_timing);	
	//delete [] congestion;
	//delete [] net_timing;
	//for (int i = 0; i < opts->num_threads; ++i) {
		//delete [] states[i];
	//}

	//return routed;
}

