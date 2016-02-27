#include "pch.h"

#include "vpr_types.h"
#include "path_delay.h"
#include "rr_graph.h"

#include "log.h"
#include "barrier.h"
#include "graph.h"
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
#include "func.h"
#include "partition.h"
#include "fm.h"

void get_stats(const RRGraph &g)
{
	for (const auto &v : get_vertices(g)) {
		const auto &ver = get_vertex(g, v);
		if (ver.properties.type == CHANX || ver.properties.type == CHANY) {
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

		fprintf(dump, "%s Num out edges: %d\n", buffer, num_out_edges(g, v));

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

bool starts_at(const RRNode &rr_node, int x, int y)
{
	bool starts;
	assert(rr_node.properties.type == CHANX || rr_node.properties.type == CHANY);

	if (rr_node.properties.inc_direction) {
		starts = rr_node.properties.xlow == x && rr_node.properties.ylow == y;
	} else {
		starts = rr_node.properties.xhigh == x && rr_node.properties.yhigh == y;
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
		const auto &v = get_vertex(g, i);
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
		//vwgt[i] = num_out_edges(g, get_vertex(g, i));
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

		//if (get_vertex(g, rr_node_id).properties.type != CHANX && get_vertex(g, rr_node_id).properties.type != CHANY) {assert(false);}
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
			//if (!((current.properties.type == CHANX || current.properties.type == CHANY) &&
			//(neighbor.properties.type == CHANX || neighbor.properties.type == CHANY))
			//|| current.properties.inc_direction == neighbor.properties.inc_direction) {
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

			//if (current.properties.type == CHANX || current.properties.type == CHANY) {
			//if ((horizontal && neighbor.properties.type == CHANX) || (!horizontal && neighbor.properties.type == CHANY)) {
			//if (neighbor.properties.inc_direction == inc_direction) {
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
		//if (valid_vertex(get_vertex(g, get_source(g, e))) && valid_vertex(get_vertex(g, get_target(g, e)))) {
			//visited_edges.push_back(e);
		//}
	//}

	//void examine_edge(int e, const RRGraph &g)
	//{
	//}

	//void discover_vertex(int v, const RRGraph &g)
	//{
		////if (valid_vertex(get_vertex(g, v))) {
			////visited_nodes.push_back(v);
		////}
	//}

	//void examine_vertex(int v, const RRGraph &g)
	//{
	//}
//};

struct connected_components_visitor_t {
	vector<int> visited_edges;
	vector<int> visited_nodes;

	template<typename Graph>
	bool tree_edge(int e, const Graph &g)
	{
		assert(find(begin(visited_edges), end(visited_edges), e) == end(visited_edges));
		visited_edges.push_back(e);
		return true;
	}

	template<typename Graph>
	void examine_edge(int e, const Graph &g)
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
		get_vertex(out, v).properties = get_vertex(in, v).properties;
	}

	for (const auto &e : get_edges(in)) {
		int from = get_source(in, e);
		int to = get_target(in, e);
		/* add reverse edge */
		add_edge(out, from, to).properties = get_edge(in, e).properties;
		add_edge(out, to, from).properties = get_edge(in, e).properties;
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
		//if (valid_vertex(get_vertex(g, v)) && color[v] == VertexColor::WHITE) {

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
	graph_t<int, int> g;

	add_vertex(g, 4);

	add_edge(g, 0, 1);
	add_edge(g, 0, 2);
	add_edge(g, 1, 3);

	/* no filter */
	auto fg1 = make_filtered_graph(g, [] (unsigned long v) -> bool { return true; }, [] (unsigned long e) -> bool { return true; });

	const auto &viter1 = get_vertices(fg1);
	vector<unsigned long> vs1(begin(viter1), end(viter1));

	assert(vs1.size() == 4);
	assert(vs1[0] == 0);
	assert(vs1[1] == 1);
	assert(vs1[2] == 2);
	assert(vs1[3] == 3);

	const auto &eiter1 = get_edges(fg1);
	vector<int> es1(begin(eiter1), end(eiter1));

	assert(es1.size() == 3);
	assert(get_source(fg1, es1[0]) == 0 && get_target(fg1, es1[0]) == 1);
	assert(get_source(fg1, es1[1]) == 0 && get_target(fg1, es1[1]) == 2);
	assert(get_source(fg1, es1[2]) == 1 && get_target(fg1, es1[2]) == 3);

	/* filtered */
	vector<bool> vf1 = { true, false, false, true };
	vector<bool> ef1 = { true, false, true };

	auto fg2 = make_filtered_graph(g, [&vf1] (unsigned long v) -> bool { return vf1[v]; }, [&ef1] (unsigned long e) -> bool { return ef1[e]; });

	const auto &viter2 = get_vertices(fg2);
	vector<unsigned long> vs2(begin(viter2), end(viter2));

	assert(vs2.size() == 2);
	assert(vs2[0] == 0);
	assert(vs2[1] == 3);

	const auto &eiter2 = get_edges(fg2);
	vector<int> es2(begin(eiter2), end(eiter2));

	assert(es2.size() == 2);
	assert(get_source(fg2, es2[0]) == 0 && get_target(fg2, es2[0]) == 1);
	assert(get_source(fg2, es2[1]) == 1 && get_target(fg2, es2[1]) == 3);
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
		auto &v = get_vertex(channel_g, i);
		v.properties.type = rr_node[i].type;
		v.properties.inc_direction = rr_node[i].direction == INC_DIRECTION;
		v.properties.xlow = rr_node[i].xlow;
		v.properties.ylow = rr_node[i].ylow;
		v.properties.xhigh = rr_node[i].xhigh;
		v.properties.yhigh = rr_node[i].yhigh;

		extern struct s_grid_tile **grid;
		auto type = &grid[v.properties.xlow][v.properties.ylow];

		//v.properties.real_xlow = rr_node[i].xlow;
		//v.properties.real_ylow = rr_node[i].ylow;
		//v.properties.real_xhigh = rr_node[i].xhigh;
		//v.properties.real_yhigh = rr_node[i].ylow + type->offset;
		v.properties.R = rr_node[i].R;
		v.properties.C = rr_node[i].C;
		v.properties.cost_index = rr_node[i].cost_index;
		v.properties.capacity = rr_node[i].capacity;

		char buffer[256];
		sprintf_rr_node(i, buffer);
		//zlog_debug(rr_log, "%s: real_xlow: %d real_xhigh: %d real_ylow: %d real_yhigh: %d\n", buffer, v.properties.real_xlow, v.properties.real_xhigh, v.properties.real_ylow, v.properties.real_yhigh);

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
			
			e.properties.buffered = switch_inf[si].buffered; 
			e.properties.switch_delay = switch_inf[si].Tdel; 
			e.properties.R = switch_inf[si].R; 
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
	vector<int> visited_edges;
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

	bool tree_edge(int e, const RRGraph &g)
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

	void examine_edge(int e, const RRGraph &g)
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

		const auto &ver = get_vertex(g, v);

		const auto &start = get_node_start(v);
		const auto &key = make_tuple(start.first, start.second, ver.properties.type, ver.properties.inc_direction);

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
	vector<int> visited_edges;
	vector<int> visited_nodes;
	set<tuple<int, int, t_rr_type, bool>> visited_channels;
	int num_nodes_discovered;

	visitor_t() :
		num_nodes_discovered(0)
	{
	}

	bool tree_edge(int e, const RRGraph &g)
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

	void examine_edge(int e, const RRGraph &g)
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
		const auto &ver = get_vertex(g, v);
		const auto &start = get_node_start(v);
		const auto &key = make_tuple(start.first, start.second, ver.properties.type, ver.properties.inc_direction);
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
		const auto &ver = get_vertex(g, v);
		const auto &start = get_node_start(v);
		const auto &key = make_tuple(start.first, start.second, ver.properties.type, ver.properties.inc_direction);

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
		const auto &ver = get_vertex(g, v);
		return (ver.properties.type == CHANX || ver.properties.type == CHANY) && pid[v] == p;
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

	bool operator()(unsigned long e) const
	{
		int from = get_source(g, e);
		int to = get_target(g, e);

		return vertex_valid(from) && vertex_valid(to);
	}
};

class rr_graph_partitioner {
	private:
		int num_tracks;

		RRGraph orig_g;
		RRGraph undirected_orig_g;
		RRGraph channel_with_interior_g;
		RRGraph channel_without_interior_g;

		vector<vector<vector<int>>> num_starting_tracks; // [type][x][y]
		vector<int> real_track_num; // [rr_node]

		int num_partitions;
		//vector<int> initial_partition;
		//int partition_index_offset;
		vector<vector<vector<vector<int>>>> num_tracks_in_partition; //[type][part][x][y]

		vector<int> result_pid;
		int total_num_chans;

	public:
		void init_graphs(t_router_opts *opts, struct s_det_routing_arch *det_routing_arch, t_direct_inf *directs, int num_directs, t_segment_inf *segment_inf, t_timing_inf *timing_inf)
		{
			extern int num_types;
			extern struct s_type_descriptor *type_descriptors;
			extern int nx, ny;
			extern struct s_grid_tile **grid;

			free_rr_graph();

			int warnings;

			build_rr_graph(GRAPH_UNIDIR, num_types, type_descriptors, nx, ny, grid,
					opts->fixed_channel_width, NULL, det_routing_arch->switch_block_type,
					det_routing_arch->Fs, det_routing_arch->num_segment,
					det_routing_arch->num_switch, segment_inf,
					det_routing_arch->global_route_switch,
					det_routing_arch->delayless_switch, *timing_inf,
					det_routing_arch->wire_to_ipin_switch, opts->base_cost_type,
					directs, num_directs, FALSE,
					&warnings);

			init_channel_only_graph(channel_with_interior_g);

			dump_rr_graph(channel_with_interior_g, "/Volumes/DATA/rr_graph_channel_with_interior.txt");

			init_graph(orig_g);

			dump_rr_graph(orig_g, "/Volumes/DATA/rr_graph.txt");

			undirected_orig_g = orig_g;

			for (const auto &e : get_edges(orig_g)) {
				int from = get_source(orig_g, e);
				int to = get_target(orig_g, e);
				add_edge(undirected_orig_g, to, from);
			}

			free_rr_graph();
			for (int i = 0; i < det_routing_arch->num_segment; ++i) {
				for (int j = 0; j < segment_inf[i].sb_len; ++j) {
					if (j != 0 && j != segment_inf[i].sb_len-1) {
						segment_inf[i].sb[j] = FALSE;
					}
				}
			}
			build_rr_graph(GRAPH_UNIDIR, num_types, type_descriptors, nx, ny, grid,
					opts->fixed_channel_width, NULL, det_routing_arch->switch_block_type,
					det_routing_arch->Fs, det_routing_arch->num_segment,
					det_routing_arch->num_switch, segment_inf,
					det_routing_arch->global_route_switch,
					det_routing_arch->delayless_switch, *timing_inf,
					det_routing_arch->wire_to_ipin_switch, opts->base_cost_type,
					directs, num_directs, FALSE,
					&warnings);

			init_channel_only_graph(channel_without_interior_g);

			dump_rr_graph(channel_without_interior_g, "/Volumes/DATA/rr_graph_channel_without_interior.txt");
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
							const auto &source = get_vertex(channel_without_interior_g, rr_node);

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
			total_num_chans = 0;
			for (const auto &v : get_vertices(orig_g)) {
				const auto &ver = get_vertex(orig_g, v);
				if (ver.properties.type == CHANX || ver.properties.type == CHANY) {
					++total_num_chans;
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
			load_real_track_numbers();
			alloc_num_tracks_in_partition();
			result_pid.resize(num_vertices(orig_g));
		}

		//void load_initial_partitions()
		//{
			//initial_partition.resize(num_vertices(orig_g));
			//for (const auto &v : get_vertices(orig_g)) {
				//const auto &ver = get_vertex(orig_g, v);
				//if (ver.properties.type == CHANX || ver.properties.type == CHANY) {
					//assert(real_track_num[v] >= 0 && real_track_num[v] < num_tracks);
					//initial_partition[v] = real_track_num[v] % 2;
					//int type = ver.properties.type == CHANX ? 0 : 1;
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
		//const auto &ver = get_vertex(orig_g, v);
		//if (ver.properties.type == CHANX || ver.properties.type == CHANY) {
		//assert(real_track_num[v] >= 0 && real_track_num[v] < num_tracks);
		//initial_partition[v] = real_track_num[v] % num_partitions;
		//int type = ver.properties.type == CHANX ? 0 : 1;
		//int x, y;
		//std::tie(x, y) = get_node_start(v);

		//++num_tracks_in_partition[type][initial_partition[v]][x][y];
		//}
		//}
		//}

		int type_index(int v) const
		{
			const auto &ver = get_vertex(orig_g, v);
			assert(ver.properties.type == CHANX || ver.properties.type == CHANY);
			return ver.properties.type == CHANX ? 0 : 1;
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
			imba = after-average;

			assert(after >= 0);

			return after == 0;
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
		//const auto &ver = get_vertex(orig_g, v);
		//return (ver.properties.type == CHANX || ver.properties.type == CHANY) && (initial_partition[v] == start || initial_partition[v] == end);
		//},
		//[&] (unsigned long e) { 
		//const auto &from = get_vertex(orig_g, get_source(orig_g, e));
		//const auto &to = get_vertex(orig_g, get_target(orig_g, e));
		//return (from.properties.type == CHANX || from.properties.type == CHANY) && (initial_partition[id(from)] == start || initial_partition[id(from)] == end) &&
		//(to.properties.type == CHANX || to.properties.type == CHANY) && (initial_partition[id(to)] == start || initial_partition[id(to)] == end);
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
			void grow_initial_partition(const Graph &g, vector<int> &initial_partition)
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

					bool tree_edge(int e, const Graph &g)
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

					void examine_edge(int e, const Graph &g)
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
						const auto &ver = get_vertex(g, v);
						const auto &start = get_node_start(v);
						const auto &key = make_tuple(start.first, start.second, ver.properties.type, ver.properties.inc_direction);
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

				int current_pid = 0;
				vector<int> num_nodes_in_partition(2, 0);
				vector<VertexColor> color(num_vertices(g), VertexColor::WHITE);
				int num_nodes = 0;
				set<int> all_visited_nodes;
				for (const auto &v : get_vertices(g)) {
					if (color[v] == VertexColor::WHITE) {
						visitor_t visitor;
						visitor.record_source(g, v);

						auto start = clock::now();
						bfs(g, { static_cast<unsigned long>(v) }, color, visitor);
						auto time = clock::now()-start;

						start = clock::now();
						assert(visitor.visited_nodes.size() == visitor.visited_edges.size()+1);
						for (const auto &v : visitor.visited_nodes) {
							initial_partition[v] = current_pid;
							color[v] = VertexColor::BLACK;
							assert(all_visited_nodes.find(v) == all_visited_nodes.end());
							all_visited_nodes.insert(v);
						}
						//for (const auto &v : get_vertices(g)) {
							//if (all_visited_nodes.find(v) == end(all_visited_nodes)) {
								//assert(color[v] == VertexColor::WHITE);
							//} else {
								//assert(color[v] == VertexColor::BLACK);
							//}
						//}
						auto update_time = clock::now()-start;

						num_nodes_in_partition[current_pid] += visitor.visited_nodes.size();
						current_pid = (current_pid+1) % 2;

						//printf("ver: %d time = %lld ms update time = %lld ms visited nodes = %lu\n", v, std::chrono::duration_cast<std::chrono::milliseconds>(time).count(), std::chrono::duration_cast<std::chrono::milliseconds>(update_time).count(), visitor.visited_nodes.size());
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
				fm<Graph, rr_graph_partitioner> fm;
				fm.init(undirected_g, initial_partition, *this);
				fm.run();

				const vector<int> &pid = fm.get_pid();

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

		void partition(int _num_partitions)
		{
			assert(std::pow(2, std::log2(_num_partitions)) == _num_partitions);

			num_partitions = _num_partitions;

			std::fill(begin(result_pid), end(result_pid), -1);

			vector<int> pid(num_vertices(orig_g), 0);
			
			recursive_bipartition(
					make_subgraph(orig_g, partition_vertex_predicate_t<typename decltype(orig_g)::base>(orig_g, pid, 0)), 
					make_subgraph(undirected_orig_g, partition_vertex_predicate_t<typename decltype(orig_g)::base>(undirected_orig_g, pid, 0)), 
					std::log2(num_partitions)-1, 0);

			int num_nodes_partitioned = 0;
			vector<int> num_nodes_in_partition(num_partitions, 0);
			for (const auto &p : result_pid) {
				if (p != -1) {
					assert(p >= 0 && p < num_partitions);
					++num_nodes_in_partition[p];
					++num_nodes_partitioned;
				}
			}
			assert(num_nodes_partitioned == total_num_chans);
			for (int i = 0; i < num_partitions; ++i) {
				printf("Num nodes in partition %d = %d\n", i, num_nodes_in_partition[i]);
			}
		}
};

void init_partitioned_graph_6(int num_partitions, const RRGraph &channel_with_interior_g, const RRGraph &channel_without_interior_g, const RRGraph &orig_g, vector<RRGraph *> &graphs, int source_x, int source_y)
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

	//while (ptc < type->num_class && type->class_inf[ptc].type != DRIVER) {
		//++ptc;
	//}

	//while (num_visited_edges < num_edges(orig_g) && ptc < type->num_class) {
		//int source_rr_node_id = get_rr_node_index(source_x, source_y, SOURCE, ptc, rr_node_indices);
		//
	set<int> all_visited_nodes;
	set<int> all_visited_edges;
	vector<VertexColor> color(num_vertices(channel_without_interior_g), VertexColor::WHITE);
	vector<int> partition(num_vertices(orig_g), -1);
	vector<vector<int>> partition_nodes;
	vector<map<tuple<int, int, t_rr_type, bool>, int>> partition_visited_channels;
	int part = 0;
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
		vector<int> starting_wires;
		for (int i = 0; i < width; ++i) {
			int source_rr_node_id = get_rr_node_index(source_x, source_y, pass == 0 ? CHANX : CHANY, i, rr_node_indices);
			const auto &source = get_vertex(channel_without_interior_g, source_rr_node_id);

			if (/*starts_at(source, source_x, source_y) && */color[source_rr_node_id] == VertexColor::WHITE) {
				starting_wires.push_back(source_rr_node_id);
			}
		}
		int inc = ceil((float)starting_wires.size()/num_partitions);
		//int inc = ceil((float)real_num_tracks/num_partitions);
		for (int i = 0, batch = 0; i < starting_wires.size(); ++batch, ++part) {
		//for (int i = 0, batch = 0; i < real_num_tracks; ++batch, ++part) {
			vector<unsigned long> sources;
			for (int j = 0; j < inc && i < starting_wires.size(); ++j, ++i) {
			//for (int j = 0; j < inc && i < real_num_tracks; ++j, ++i) {
				int source_rr_node_id = starting_wires[i];
				//int source_rr_node_id = get_rr_node_index(source_x, source_y, pass == 0 ? CHANX : CHANY, i, rr_node_indices); 
				assert(color[source_rr_node_id] == VertexColor::WHITE);
				sources.push_back(source_rr_node_id);
				char buffer[256];
				sprintf_rr_node(source_rr_node_id, buffer);
				printf("first stage batch %d source %d: %s\n", batch, j, buffer);
			}

			visitor_2_t visitor(num_vertices(channel_without_interior_g));
			bfs_checked(channel_without_interior_g, sources, color, visitor, all_visited_nodes, all_visited_edges);
			visitor.print_path(89280);

			//set<tuple<int, int, t_rr_type, bool>> visited_channels_new;
			//for (const auto &n : visitor.visited_nodes) {
				//const auto &start = get_node_start(n);
				//const auto &node = get_vertex(channel_without_interior_g, n);
				//const auto &key = make_tuple(start.first, start.second, node.properties.type, node.properties.inc_direction);
				//assert(visited_channels_new.find(key) == visited_channels_new.end());
				//visited_channels_new.insert(key);
			//}

			char filename[256];
			sprintf(filename, "/Volumes/DATA/visited_pass_%d_batch_%d.txt", pass, batch);
			FILE *file = fopen(filename, "w");
			for (const auto &v : visitor.visited_nodes) {
				const auto &start = get_node_start(v);
				fprintf(file, "%d %d\n", start.first, start.second);
			}
			fclose(file);

			partition_nodes.push_back(visitor.visited_nodes);
			for (const auto &v : visitor.visited_nodes) {
				partition[v] = part;
			}

			int num_unvisited_nodes = (nx+2)*(ny+2)*4 - visitor.visited_nodes.size();
			//for (int x = 0; x < nx+2; ++x) {
				//for (int y = 0; y < ny+2; ++y) {
					//if (visited_channels_new.find(make_tuple(x, y, CHANX, true)) == visited_channels_new.end()) {
						////printf("Did not visit (%d,%d) CHANX increasing\n", x, y);
						//++num_unvisited_nodes;
					//}
					//if (visited_channels_new.find(make_tuple(x, y, CHANX, false)) == visited_channels_new.end()) {
						////printf("Did not visit (%d,%d) CHANX decreasing\n", x, y);
						//++num_unvisited_nodes;
					//}
					//if (visited_channels_new.find(make_tuple(x, y, CHANY, true)) == visited_channels_new.end()) {
						////printf("Did not visit (%d,%d) CHANY increasing\n", x, y);
						//++num_unvisited_nodes;
					//}
					//if (visited_channels_new.find(make_tuple(x, y, CHANY, false)) == visited_channels_new.end()) {
						////printf("Did not visit (%d,%d) CHANY decreasing\n", x, y);
						//++num_unvisited_nodes;
					//}
				//}
			//}

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
		const auto &ver = get_vertex(orig_g, v);
		if (ver.properties.type == CHANX || ver.properties.type == CHANY) {
			++num_channels;
		}
	}
	printf("%d,%d num_all_visited_nodes = %lu num_channels = %d\n", source_x, source_y, all_visited_nodes.size(), num_channels);

	char filename[256];
	sprintf(filename, "/Volumes/DATA/unvisited.txt");
	FILE *file = fopen(filename, "w");
	int num_unvisited_nodes = 0;
	for (const auto &v : get_vertices(orig_g)) {
		const auto &ver = get_vertex(orig_g, v);
		if ((ver.properties.type == CHANX || ver.properties.type == CHANY) && all_visited_nodes.find(v) == all_visited_nodes.end()) {
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

	for (int i = 0; i < partition_nodes.size(); ++i) {
		for (const auto &v : partition_nodes[i]) {
			assert(partition[v] == i);
		}
	}

	for (const auto &e : get_edges(channel_with_interior_g)) {
		int from = get_source(channel_with_interior_g, e);
		int to = get_target(channel_with_interior_g, e);
		if (!has_edge(channel_without_interior_g, from, to) || all_visited_edges.find(get_edge(channel_without_interior_g, from, to)) ==
					all_visited_edges.end()) {
			if (!(partition[from] == -1 && partition[to] == -1)) {
				char buffer[256];
				//sprintf_rr_node
				assert(false);
			}
			if (partition[from] == -1 || partition[to] == -1) {
			} else if (partition[from] == partition[to]) {
				//partition_nodes.
			}
		}
	}

	//for (const auto &v : get_vertices(channel_without_interior_g)) {
		//if ((v.properties.type == CHANX || v.properties.type == CHANY)) {
			//if (all_visited_nodes.find(id(v)) == all_visited_nodes.end()) {
				//visitor_t visitor;
				//bfs_checked(channel_without_interior_g, { v }, color, visitor, all_visited_nodes);

				//int num_unvisited_nodes = (nx+2)*(ny+2)*4 - visitor.visited_nodes.size();

				//printf("source: %d\n", id(v));
				//printf("num_nodes_discovered = %d\n", visitor.num_nodes_discovered);
				//printf("num_visited_nodes = %lu\n", visitor.visited_nodes.size());
				//printf("num_unvisited_nodes = %d\n", num_unvisited_nodes);
				//printf("num_visited_edge = %lu\n", visitor.visited_edges.size());
			//} else {
				//assert(color[id(v)] == VertexColor::BLACK);
			//}
		//}
	//}

	return;

	vector<int> partitions;
	int num_chan_edges = 0;
	int num_misc_edges = -1;

	for (int i = 0; i < num_partitions; ++i) {
		RRGraph &new_g = *new RRGraph;
		add_vertex(new_g, num_vertices(channel_with_interior_g));

		for (int i = 0; i < num_vertices(new_g); ++i) {
			get_vertex(new_g, i).properties = get_vertex(orig_g, i).properties;
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

				RREdge &new_e = add_edge(new_g, from, to);
				const auto &edge = get_edge(channel_with_interior_g, e);
				new_e.properties = edge.properties;
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
				//if (get_vertex(orig_g, components[j][k]).properties.type != CHANX && 
//get_vertex(orig_g, components[j][k]).properties.type != CHANY) {
					//assert(false);
				//}
				fprintf(file, "\t%s Part: %d\n", buffer, partitions[components[j][k]]);
			}
		}
		fclose(file);

		int temp_num_misc_edges = 0;
		for (const auto &e : get_edges(orig_g)) {
			const auto &from = get_vertex(orig_g, get_source(orig_g, e));
			const auto &to = get_vertex(orig_g, get_target(orig_g, e));

			if ((from.properties.type == CHANX || from.properties.type == CHANY) && 
					(to.properties.type == CHANX || to.properties.type == CHANY)) {
				continue;
			}

			RREdge &new_e = add_edge(new_g, id(from), id(to));
			const auto &edge = get_edge(orig_g, e);
			new_e.properties = edge.properties;

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
		vector<int> visited_edges;
		vector<int> visited_nodes;
		set<tuple<int, int, t_rr_type, bool>> visited_channels;
		int num_nodes_discovered;

		visitor_t() :
			num_nodes_discovered(0)
		{
		}

		bool tree_edge(int e, const RRGraph &g)
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

		void examine_edge(int e, const RRGraph &g)
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
			const auto &ver = get_vertex(g, v);
			const auto &start = get_node_start(v);
			const auto &key = make_tuple(start.first, start.second, ver.properties.type, ver.properties.inc_direction);
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
			const auto &source = get_vertex(channel_without_interior_g, source_rr_node_id);

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
				const auto &node = get_vertex(channel_without_interior_g, n);
				const auto &key = make_tuple(start.first, start.second, node.properties.type, node.properties.inc_direction);
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
		const auto &ver = get_vertex(orig_g, v);
		if (ver.properties.type == CHANX || ver.properties.type == CHANY) {
			++num_channels;
		}
	}
	printf("%d,%d num_all_visited_nodes = %lu num_channels = %d\n", source_x, source_y, all_visited_nodes.size(), num_channels);

	char filename[256];
	sprintf(filename, "/Volumes/DATA/unvisited.txt");
	FILE *file = fopen(filename, "w");
	int num_unvisited_nodes = 0;
	for (const auto &v : get_vertices(orig_g)) {
		const auto &ver = get_vertex(orig_g, v);
		if ((ver.properties.type == CHANX || ver.properties.type == CHANY) && all_visited_nodes.find(v) == all_visited_nodes.end()) {
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
		const auto &ver = get_vertex(channel_without_interior_g, v);
		if ((ver.properties.type == CHANX || ver.properties.type == CHANY)) {
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
			get_vertex(new_g, i).properties = get_vertex(orig_g, i).properties;
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

				RREdge &new_e = add_edge(new_g, from, to);
				const auto &edge = get_edge(channel_with_interior_g, e);
				new_e.properties = edge.properties;
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
				//if (get_vertex(orig_g, components[j][k]).properties.type != CHANX && 
//get_vertex(orig_g, components[j][k]).properties.type != CHANY) {
					//assert(false);
				//}
				fprintf(file, "\t%s Part: %d\n", buffer, partitions[components[j][k]]);
			}
		}
		fclose(file);

		int temp_num_misc_edges = 0;
		for (const auto &e : get_edges(orig_g)) {
			const auto &from = get_vertex(orig_g, get_source(orig_g, e));
			const auto &to = get_vertex(orig_g, get_target(orig_g, e));

			if ((from.properties.type == CHANX || from.properties.type == CHANY) && 
					(to.properties.type == CHANX || to.properties.type == CHANY)) {
				continue;
			}

			RREdge &new_e = add_edge(new_g, id(from), id(to));
			const auto &edge = get_edge(orig_g, e);
			new_e.properties = edge.properties;

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

	vector<int> partition;
	vector<int> weights(num_vertices(channel_without_interior_g));
	for (const auto &v : get_vertices(channel_without_interior_g)) {
		const auto &ver = get_vertex(channel_without_interior_g, v);
		if (ver.properties.type == CHANX || ver.properties.type == CHANY) {
			weights[v] = 1;
		} else {
			weights[v] = 0;
		}
	}
	partition_graph(channel_without_interior_g, num_partitions, weights, 1.001, partition);
	vector<int> num_nodes(num_partitions, 0);
	vector<int> num_non_chan_nodes(num_partitions, 0);
	for (int i = 0; i < num_vertices(channel_without_interior_g); ++i) {
		const auto &v = get_vertex(channel_without_interior_g, i);
		if (v.properties.type != CHANX && v.properties.type != CHANY) {
			++num_non_chan_nodes[partition[i]];
		}
		++num_nodes[partition[i]];
	}
	int total_num_nodes = 0;
	vector<int> num_chan_nodes(num_partitions);
	for (int i = 0; i < num_partitions; ++i) {
		printf("part %d num nodes %d num_chan_nodes: %d num_non_chan_nodes: %d\n", i, num_nodes[i], num_nodes[i] - num_non_chan_nodes[i], num_non_chan_nodes[i]);
		total_num_nodes += num_nodes[i];
		num_chan_nodes[i] = num_nodes[i] - num_non_chan_nodes[i];
	}
	assert(total_num_nodes == num_vertices(channel_without_interior_g));

	for (int i = 0; i < num_vertices(orig_g); ++i) {
		const auto &ver = get_vertex(orig_g, i);
		if (ver.properties.type != CHANX && ver.properties.type != CHANY) {
			partition[i] = -1;
		}
	}

	int num_chan_edges = 0;
	int num_non_chan_edges = -1;

	

	for (int i = 0; i < num_partitions; ++i) {
		RRGraph &new_g = *new RRGraph;
		add_vertex(new_g, num_vertices(channel_with_interior_g));

		for (int i = 0; i < num_vertices(new_g); ++i) {
			get_vertex(new_g, i).properties = get_vertex(orig_g, i).properties;
		}

		int num_intra_partiton_interior_edges = 0;
		int num_interior_edges = 0;
		for (const auto &e : get_edges(channel_with_interior_g)) {
			int from = get_source(channel_with_interior_g, e);
			int to = get_target(channel_with_interior_g, e);

			if (!has_edge(channel_without_interior_g, from, to)) {
				++num_interior_edges;
			}

			if (partition[from] == i && partition[from] == partition[to]) {
				if (!has_edge(channel_without_interior_g, from, to)) {
					++num_intra_partiton_interior_edges;
				}

				RREdge &new_e = add_edge(new_g, from, to);
				const auto &edge = get_edge(channel_with_interior_g, e);
				new_e.properties = edge.properties;
			}
		}

		num_chan_edges += num_edges(new_g);

		printf("graph %d num intra partition interior edges/num interior edges = %d/%d (%g)\n", i, num_intra_partiton_interior_edges, num_interior_edges, (float)num_intra_partiton_interior_edges/num_interior_edges*100);

		int local_num_non_chan_edges = 0;
		for (const auto &e : get_edges(orig_g)) {
			const auto &from = get_vertex(orig_g, get_source(orig_g, e));
			const auto &to = get_vertex(orig_g, get_target(orig_g, e));

			if ((from.properties.type == CHANX || from.properties.type == CHANY) && 
					(to.properties.type == CHANX || to.properties.type == CHANY)) {
				continue;
			}

			RREdge &new_e = add_edge(new_g, id(from), id(to));
			const auto &edge = get_edge(orig_g, e);
			new_e.properties = edge.properties;

			++local_num_non_chan_edges;
		}

		if (num_non_chan_edges == -1) {
			num_non_chan_edges = local_num_non_chan_edges;
		}

		assert(num_edges(new_g) <= num_edges(orig_g));

		graphs.push_back(&new_g);

		vector<vector<int>> &components = all_partition_components[i];
		//connected_components(new_g, [&partition, &i] (const RRNode &v) -> bool { return partition[id(v)] == i && (v.properties.type == CHANX || v.properties.type == CHANY); }, components);
		//connected_components(new_g, [&partition, &i] (const RRNode &v) -> bool { return partition[id(v)] == i || (v.properties.type != CHANX && v.properties.type != CHANY); }, components);
		RRGraph temp_g;
		connected_components_convert_to_undirected(new_g, temp_g);

		int expected_total_num_nodes_in_components = 0;	
		vector<bool> vf(num_vertices(temp_g));
		for (const auto &v : get_vertices(temp_g)) {
			const auto &ver = get_vertex(temp_g, v);
			if (ver.properties.type != CHANX && ver.properties.type != CHANY) {
				assert(partition[v] == -1);
			}
			vf[v] = partition[v] == i || (ver.properties.type != CHANX && ver.properties.type != CHANY);
			if (vf[v]) {
				++expected_total_num_nodes_in_components;
			}
		}
		vector<bool> ef(num_edges(temp_g));
		for (const auto &e : get_edges(temp_g)) {
			int from = get_source(temp_g, e);
			int to = get_target(temp_g, e);
			ef[e] = vf[from] && vf[to];
		}

		//connected_components(make_filtered_graph(temp_g, &vf, &ef), components);
		int total_num_nodes_in_components = 0;
		for (const auto &comp : components) {
			total_num_nodes_in_components += comp.size();
		}
		//assert(total_num_nodes_in_components == num_chan_nodes[i]);
		assert(total_num_nodes_in_components == expected_total_num_nodes_in_components);
		printf("graph %d num connected components = %lu average num nodes = %g\n", i, components.size(), (float)total_num_nodes_in_components/components.size());

		//assert(total_num_nodes_in_components == num_vertices(orig_g));

		char filename[256];
		char buffer[256];
		sprintf(filename, "/Volumes/DATA/components/graph_%d_components.txt", i);
		FILE *file = fopen(filename, "w");
		for (int j = 0; j < components.size(); ++j) {
			fprintf(file, "Component %d Size %lu\n", j, components[j].size());
			for (int k = 0; k < components[j].size(); ++k) {
				sprintf_rr_node(components[j][k], buffer);
				//if (get_vertex(orig_g, components[j][k]).properties.type != CHANX && 
//get_vertex(orig_g, components[j][k]).properties.type != CHANY) {
					//assert(false);
				//}
				fprintf(file, "\t%s Part: %d\n", buffer, partition[components[j][k]]);
			}
		}
		fclose(file);

	}

	assert(num_non_chan_edges + num_chan_edges <= num_edges(orig_g));
	printf("num_chan_edges_after_partition: %d total_chan_edges: %d percentage: %g\n", num_chan_edges, num_edges(orig_g)-num_non_chan_edges, num_chan_edges*100.0/(num_edges(orig_g)-num_non_chan_edges));
	//assert(all_of(begin(edge_valid), end(edge_valid), [] (bool val) -> bool { return !val; }));
	//assert(num_edges_spanned == num_edges(g));
}

void init_partitioned_graph_3(int num_partitions, RRGraph &g, vector<RRGraph *> &graphs)
{
	extern s_rr_node *rr_node;
	extern int num_rr_nodes;
	extern t_rr_indexed_data *rr_indexed_data;
	extern struct s_switch_inf *switch_inf;

	set<int> duplicate_edges;
	RRGraph orig_g;

	add_vertex(g, num_rr_nodes);
	add_vertex(orig_g, num_rr_nodes);

	for (int i = 0; i < num_vertices(g); ++i) {
		auto &v = get_vertex(g, i);
		v.properties.type = rr_node[i].type;
		v.properties.inc_direction = rr_node[i].direction == INC_DIRECTION;
		v.properties.xlow = rr_node[i].xlow;
		v.properties.ylow = rr_node[i].ylow;
		v.properties.xhigh = rr_node[i].xhigh;
		v.properties.yhigh = rr_node[i].yhigh;

		extern struct s_grid_tile **grid;
		auto type = &grid[v.properties.xlow][v.properties.ylow];

		//v.properties.real_xlow = rr_node[i].xlow;
		//v.properties.real_ylow = rr_node[i].ylow;
		//v.properties.real_xhigh = rr_node[i].xhigh;
		//v.properties.real_yhigh = rr_node[i].ylow + type->offset;
		v.properties.R = rr_node[i].R;
		v.properties.C = rr_node[i].C;
		v.properties.cost_index = rr_node[i].cost_index;
		v.properties.capacity = rr_node[i].capacity;

		char buffer[256];
		sprintf_rr_node(i, buffer);
		//zlog_debug(rr_log, "%s: real_xlow: %d real_xhigh: %d real_ylow: %d real_yhigh: %d\n", buffer, v.properties.real_xlow, v.properties.real_xhigh, v.properties.real_ylow, v.properties.real_yhigh);

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
			
			e.properties.buffered = switch_inf[si].buffered; 
			e.properties.switch_delay = switch_inf[si].Tdel; 
			e.properties.R = switch_inf[si].R; 

			auto &orig_e = add_edge(orig_g, i, neighbor_id);
			orig_e.properties = e.properties;

			auto &e2 = add_edge(g, neighbor_id, i);

			e2.properties = e.properties; 

			duplicate_edges.insert(id(e2));
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
				RREdge &new_e = add_edge(new_g, from, to);
				const auto &edge = get_edge(g, e);
				new_e.properties = edge.properties;
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
		auto &v = get_vertex(g, i);
		v.properties.type = rr_node[i].type;
		v.properties.inc_direction = rr_node[i].direction == INC_DIRECTION;
		v.properties.xlow = rr_node[i].xlow;
		v.properties.ylow = rr_node[i].ylow;
		v.properties.xhigh = rr_node[i].xhigh;
		v.properties.yhigh = rr_node[i].yhigh;

		extern struct s_grid_tile **grid;
		auto type = &grid[v.properties.xlow][v.properties.ylow];

		//v.properties.real_xlow = rr_node[i].xlow;
		//v.properties.real_ylow = rr_node[i].ylow;
		//v.properties.real_xhigh = rr_node[i].xhigh;
		//v.properties.real_yhigh = rr_node[i].ylow + type->offset;
		v.properties.R = rr_node[i].R;
		v.properties.C = rr_node[i].C;
		v.properties.cost_index = rr_node[i].cost_index;
		v.properties.capacity = rr_node[i].capacity;

		char buffer[256];
		sprintf_rr_node(i, buffer);
		//zlog_debug(rr_log, "%s: real_xlow: %d real_xhigh: %d real_ylow: %d real_yhigh: %d\n", buffer, v.properties.real_xlow, v.properties.real_xhigh, v.properties.real_ylow, v.properties.real_yhigh);

		for (int j = 0; j < rr_node[i].num_edges; ++j) {
			int neighbor = rr_node[i].edges[j];

			//if ((rr_node[i].type == CHANX || rr_node[i].type == CHANY) &&
					//(rr_node[neighbor].type == CHANX || rr_node[neighbor].type == CHANY) &&
					//get_track_domain(rr_node[neighbor].ptc_num, num_partitions) != get_track_domain(rr_node[i].ptc_num, num_partitions)) {
				//continue;
			//}

			auto &e = add_edge(g, i, neighbor);

			int si = rr_node[i].switches[j];
			
			e.properties.buffered = switch_inf[si].buffered; 
			e.properties.switch_delay = switch_inf[si].Tdel; 
			e.properties.R = switch_inf[si].R; 
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

			const auto &source = get_vertex(g, source_rr_node_id);
			if (!starts_at(source, source_x, source_y)) {
				continue;
			}

			vector<int> visited_edges;
			//bfs(g, source, pass == 0, source.properties.inc_direction, visited, visited_edges);

			num_visited_edges += visited_edges.size();	
			//while (++ptc < type->num_class && type->class_inf[ptc].type != DRIVER) {
			//}

			printf("visited nodes size: %lu\n", visited_edges.size());

			RRGraph new_g;
			add_vertex(new_g, num_vertices(g));

			for (const auto &e_id : visited_edges) {
				const RREdge &e = get_edge(g, e_id);
				RREdge &new_e = add_edge(new_g, get_source(g, e_id), get_target(g, e_id));
				new_e.properties = e.properties;
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
	extern s_rr_node *rr_node;
	extern int num_rr_nodes;
	extern t_rr_indexed_data *rr_indexed_data;
	extern struct s_switch_inf *switch_inf;

	add_vertex(g, num_rr_nodes);
	for (int i = 0; i < num_vertices(g); ++i) {
		auto &v = get_vertex(g, i);
		v.properties.type = rr_node[i].type;
		v.properties.inc_direction = rr_node[i].direction == INC_DIRECTION;
		v.properties.xlow = rr_node[i].xlow;
		v.properties.ylow = rr_node[i].ylow;
		v.properties.xhigh = rr_node[i].xhigh;
		v.properties.yhigh = rr_node[i].yhigh;

		extern struct s_grid_tile **grid;
		auto type = &grid[v.properties.xlow][v.properties.ylow];

		//v.properties.real_xlow = rr_node[i].xlow;
		//v.properties.real_ylow = rr_node[i].ylow;
		//v.properties.real_xhigh = rr_node[i].xhigh;
		//v.properties.real_yhigh = rr_node[i].ylow + type->offset;
		v.properties.R = rr_node[i].R;
		v.properties.C = rr_node[i].C;
		v.properties.cost_index = rr_node[i].cost_index;
		v.properties.capacity = rr_node[i].capacity;

		char buffer[256];
		sprintf_rr_node(i, buffer);
		//zlog_debug(rr_log, "%s: real_xlow: %d real_xhigh: %d real_ylow: %d real_yhigh: %d\n", buffer, v.properties.real_xlow, v.properties.real_xhigh, v.properties.real_ylow, v.properties.real_yhigh);

		for (int j = 0; j < rr_node[i].num_edges; ++j) {
			int neighbor = rr_node[i].edges[j];

			//if ((rr_node[i].type == CHANX || rr_node[i].type == CHANY) &&
					//(rr_node[neighbor].type == CHANX || rr_node[neighbor].type == CHANY) &&
					//get_track_domain(rr_node[neighbor].ptc_num, num_partitions) != get_track_domain(rr_node[i].ptc_num, num_partitions)) {
				//continue;
			//}

			auto &e = add_edge(g, i, neighbor);

			int si = rr_node[i].switches[j];
			
			e.properties.buffered = switch_inf[si].buffered; 
			e.properties.switch_delay = switch_inf[si].Tdel; 
			e.properties.R = switch_inf[si].R; 
		}
	}
	printf("RR graph num vertices: %d\n", num_vertices(g));
	printf("RR graph num edges: %d\n", num_edges(g));

	dump_rr_graph(g, "rr_graph.txt");
	
	vector<bool> edge_valid(num_edges(g), true);
	vector<vector<int>> trees;
	int num_edges_spanned = 0;
	while (num_edges_spanned < num_edges(g)) {
		trees.push_back(vector<int>());

		auto &edges = trees.back();
		custom_minimum_spanning_tree(g, [&] (int e) -> bool {
				const auto &from = get_vertex(g, get_source(g, e));
				const auto &to = get_vertex(g, get_target(g, e));
				return edge_valid[e] && 
				(from.properties.type == CHANX || from.properties.type == CHANY) &&
				 (to.properties.type == CHANX || to.properties.type == CHANY);

				}, edges);
		num_edges_spanned += edges.size();
		if (edges.size() == 0) {
			break;
		}
		printf("mst edges size: %lu\n", edges.size());
		for (const auto &e : edges) {
			edge_valid[e] = false;
		}
	}
	//assert(all_of(begin(edge_valid), end(edge_valid), [] (bool val) -> bool { return !val; }));
	//assert(num_edges_spanned == num_edges(g));

	for (int i = 0; i < trees.size(); ++i) {
		RRGraph new_g;
		add_vertex(new_g, num_vertices(g));

		for (auto e : trees[i]) {
			RREdge &new_e = add_edge(new_g, get_source(g, e), get_target(g, e));
			const auto &edge = get_edge(g, e);
			new_e.properties = edge.properties;
		}
		char filename[256];
		sprintf(filename, "rr_graph_%d.txt", i);
		dump_rr_graph(new_g, filename);
	}

	int low = 0;
	int high = trees.size();
	int inc = ceil((float)trees.size()/(num_partitions*2));
	vector<int> used(trees.size(), 0);
	while (low < high) {
		RRGraph &new_g = *new RRGraph;
		add_vertex(new_g, num_vertices(g));
		for (int i = 0; i < num_vertices(g); ++i) {
			get_vertex(new_g, i).properties = get_vertex(g, i).properties;
		}
		for (int i = low; i < std::min(low+inc, high); ++i) {
			assert(used[i] == 0);
			++used[i];
			for (auto e : trees[i]) {
				RREdge &new_e = add_edge(new_g, get_source(g, e), get_target(g, e));
				const auto &edge = get_edge(g, e);
				new_e.properties = edge.properties;
			}
		}
		low += inc;
		for (int i = std::max(high-inc, low+1); i < high; ++i) {
			assert(used[i] == 0);
			++used[i];
			for (auto e : trees[i]) {
				RREdge &new_e = add_edge(new_g, get_source(g, e), get_target(g, e));
				const auto &edge = get_edge(g, e);
				new_e.properties = edge.properties;
			}
		}
		high -= inc;
		graphs.push_back(&new_g);
		assert(num_edges(new_g) > 0);
	}

	assert(all_of(begin(used), end(used), [] (int val) -> bool { return val == 1; }));
	assert(graphs.size() == num_partitions);

	for (const auto e : get_edges(g)) {
		if (get_vertex(g, get_source(g, e)).properties.type == SOURCE || get_vertex(g, get_target(g, e)).properties.type == SINK) {
			for (auto &new_g : graphs) {
				int from = get_source(g, e);
				int to = get_target(g, e);

				if (!has_edge(*new_g, from, to)) {
					RREdge &new_e = add_edge(*new_g, from, to);
					const auto &edge = get_edge(g, e);
					new_e.properties = edge.properties;
				}
			}
		}
	}

	//for (int i = 0; i < graphs.size(); ++i) {
		//char filename[256];
		//sprintf(filename, "rr_graph_%d.txt", i);
		//dump_rr_graph(*graphs[i], filename);
	//}
}

#define LOG_PATH_PREFIX "/Volumes/DATA/"

vector<FILE *> delta_log_files;
vector<FILE *> missing_edge_log_files;

static void concurrent_log_impl(zlog_msg_t *msg, vector<FILE *> &log_files, int tid)
{
	assert(tid >= 0 && tid < log_files.size());
	FILE *file = log_files[tid];
	if (!file) {
		char filename[256];
		sprintf(filename, "%s%s", LOG_PATH_PREFIX, msg->path);

		file = fopen(filename, "w");
		if (!file) {
			perror(nullptr);
			assert(false);
		}

		log_files[tid] = file;
	}
	fprintf(file, "%s", msg->buf);
	fflush(file);
}

static int delta_log_output(zlog_msg_t *msg)
{
	int tid;
	int iter;

	if (sscanf(msg->path, "iter_%d_tid_%d.log", &iter, &tid) != 2) {
		return 0;
	}

	concurrent_log_impl(msg, delta_log_files, tid);

	return 0;
}

static int missing_edge_log_output(zlog_msg_t *msg)
{
	int tid;
	int iter;

	if (sscanf(msg->path, "missing_edge_iter_%d_tid_%d.log", &iter, &tid) != 2) {
		return 0;
	}

	concurrent_log_impl(msg, missing_edge_log_files, tid);

	return 0;
}

bool spatial_route(t_router_opts *opts, struct s_det_routing_arch det_routing_arch, t_direct_inf *directs, int num_directs, t_segment_inf *segment_inf, t_timing_inf timing_inf)
{
	using std::chrono::duration_cast;
	using std::chrono::nanoseconds;
	using clock = std::chrono::high_resolution_clock;

	init_logging();
	zlog_set_record("custom_output", delta_log_output);
	zlog_set_record("missing_edge", missing_edge_log_output);
	delta_log_files.resize(opts->num_threads, nullptr);
	missing_edge_log_files.resize(opts->num_threads, nullptr);

	//fclose(fopen(LOG_PATH_PREFIX"iter__tid_.log", "w"));

	test_fm();
	test_filter_graph();
	test_partition_graph();
	test_connected_components();

	rr_graph_partitioner partitioner;
	partitioner.init(opts, &det_routing_arch, directs, num_directs, segment_inf, &timing_inf);
	partitioner.partition(opts->num_threads);
	return;

	extern int num_types;
	extern struct s_type_descriptor *type_descriptors;
	extern int nx, ny;
	extern struct s_grid_tile **grid;

	free_rr_graph();

	int warnings;

	build_rr_graph(GRAPH_UNIDIR, num_types, type_descriptors, nx, ny, grid,
			opts->fixed_channel_width, NULL, det_routing_arch.switch_block_type,
			det_routing_arch.Fs, det_routing_arch.num_segment,
			det_routing_arch.num_switch, segment_inf,
			det_routing_arch.global_route_switch,
			det_routing_arch.delayless_switch, timing_inf,
			det_routing_arch.wire_to_ipin_switch, opts->base_cost_type,
			directs, num_directs, FALSE,
			&warnings);

	RRGraph channel_with_interior_g;
	init_channel_only_graph(channel_with_interior_g);

	dump_rr_graph(channel_with_interior_g, "/Volumes/DATA/rr_graph_channel_with_interior.txt");

	RRGraph orig_g;
	init_graph(orig_g);

	dump_rr_graph(orig_g, "/Volumes/DATA/rr_graph.txt");

	free_rr_graph();
	for (int i = 0; i < det_routing_arch.num_segment; ++i) {
		for (int j = 0; j < segment_inf[i].sb_len; ++j) {
			if (j != 0 && j != segment_inf[i].sb_len-1) {
				segment_inf[i].sb[j] = FALSE;
			}
		}
	}
	build_rr_graph(GRAPH_UNIDIR, num_types, type_descriptors, nx, ny, grid,
			opts->fixed_channel_width, NULL, det_routing_arch.switch_block_type,
			det_routing_arch.Fs, det_routing_arch.num_segment,
			det_routing_arch.num_switch, segment_inf,
			det_routing_arch.global_route_switch,
			det_routing_arch.delayless_switch, timing_inf,
			det_routing_arch.wire_to_ipin_switch, opts->base_cost_type,
			directs, num_directs, FALSE,
			&warnings);

	RRGraph channel_without_interior_g;
	init_channel_only_graph(channel_without_interior_g);

	dump_rr_graph(channel_without_interior_g, "/Volumes/DATA/rr_graph_channel_without_interior.txt");

	vector<RRGraph *> graphs;
	vector<vector<vector<int>>> all_partition_components(opts->num_threads);
	//for (int x = 0; x < nx+1; ++x) {
		//for (int y = 0; y < ny+1; ++y) {
			//init_partitioned_graph_5(opts->num_threads, channel_with_interior_g, channel_without_interior_g, orig_g, graphs, x, y);
		//}
	//}
	//init_partitioned_graph_4(opts->num_threads, channel_with_interior_g, channel_without_interior_g, orig_g, graphs, 14, 14);
	init_partitioned_graph_4(opts->num_threads, channel_with_interior_g, channel_without_interior_g, orig_g, graphs, all_partition_components);

	for (int i = 0; i < graphs.size(); ++i) {
		char filename[256];

		sprintf(filename, "/Volumes/DATA/rr_graph_%d.txt", i);
		dump_rr_graph(*graphs[i], filename);

		sprintf(filename, "/Volumes/DATA/rr_graph_%d_edges.txt", i);
		dump_edges(*graphs[i], filename);
	}

	vector<net_t> nets;
	vector<net_t> global_nets;
	init_nets(nets, global_nets, opts->bb_factor);	

	t_net_timing *net_timing = new t_net_timing[nets.size()+global_nets.size()];
	init_net_timing(nets, global_nets, net_timing);

	vector<route_tree_t> route_trees(nets.size());
	for (int i = 0; i < nets.size(); ++i) {
		route_tree_init(route_trees[i]);
	}

	vector<route_state_t *> states(opts->num_threads);

	for (int i = 0; i < opts->num_threads; ++i) {
		states[i] = new route_state_t[num_vertices(*graphs[0])];
		for (int j = 0; j < num_vertices(*graphs[0]); ++j) {
			states[i][j].rr_node = -1;
			states[i][j].known_cost = std::numeric_limits<float>::max();
			states[i][j].cost = std::numeric_limits<float>::max();
			states[i][j].prev_edge = -1;
			states[i][j].upstream_R = -1;
			states[i][j].delay = std::numeric_limits<float>::max();
		}
	}

	congestion_t *congestion = new congestion_t[num_vertices(*graphs[0])];
	for (int i = 0; i < num_vertices(*graphs[0]); ++i) {
		congestion[i].acc_cost = 1;
		congestion[i].pres_cost = 1;
		congestion[i].occ = 0;
	}

	congestion_t *temp_congestion = new congestion_t[num_vertices(*graphs[0])];
	for (int i = 0; i < num_vertices(*graphs[0]); ++i) {
		temp_congestion[i].acc_cost = 1;
		temp_congestion[i].pres_cost = 1;
		temp_congestion[i].occ = 0;
	}

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
		vector<int> graph_used_by_net(nets.size(), -1);

		auto greedy_route_start = clock::now();

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

				int i;
				while ((i = net_index++) < nets_to_route.size()) {
					net_t *net = nets_to_route[i].second;

					//auto rip_up_start = clock::now();
					//if (greedy_rip_up_all) {
						//route_tree_mark_all_nodes_to_be_ripped(route_trees[net->local_id], g);
					//} else {
						route_tree_mark_congested_nodes_to_be_ripped(route_trees[net->local_id], *graphs[tid], congestion);
					//}
					route_tree_rip_up_marked(route_trees[net->local_id], *graphs[tid], congestion, params.pres_fac, true, &local_lock_perf);

					//local_perf.total_rip_up_time += clock::now()-rip_up_start;

					//auto route_start = clock::now();

					vector<sink_t *> sinks;	
					sinks.reserve(net->sinks.size());
					for (auto &sink : net->sinks) {
						const RouteTreeNode *sink_rt_node = route_tree_get_rt_node(route_trees[net->local_id], sink.rr_node);
						if (!sink_rt_node)  {
							sinks.push_back(&sink);
						} else {
							assert(!sink_rt_node->properties.pending_rip_up);
						}
					}
					if (!sinks.empty()) {
						int num_sinks_routed = route_net_2(*graphs[tid], net->vpr_id, &net->source, sinks, params, states[tid], congestion, route_trees[net->local_id], net_timing[net->vpr_id], true, &local_perf, &local_lock_perf);

						if (num_sinks_routed == sinks.size()) {
							++local_num_nets_routed;
						} else {
							route_tree_t temp_route_tree;
							route_tree_init(temp_route_tree);

							zlog_level(delta_log, ROUTER_V1, "Routing net %d with original RR graph\n", net->vpr_id);

							int orig_num_sinks_routed = route_net_2(orig_g, net->vpr_id, &net->source, sinks, params, states[tid], temp_congestion, temp_route_tree, net_timing[net->vpr_id], true, nullptr, nullptr);

							assert(orig_num_sinks_routed == sinks.size());

							int num_missing_edges = 0;
							int num_missing_interior_edges = 0;
							zlog_debug(missing_edge_log, "Net %d routed in graph %d requires edges below: \n", tid, net->vpr_id);
							for (const auto &rt_edge_id : get_edges(temp_route_tree.graph)) {
								int rr_edge = get_edge(temp_route_tree.graph, rt_edge_id).properties.rr_edge;
								int from = get_source(orig_g, rr_edge);
								int to = get_target(orig_g, rr_edge);
								if (!has_edge(*graphs[tid], from, to)) {
									int found = -1;
									int num_found = 0;
									for (int i = 0; i < opts->num_threads; ++i) {
										if (i != tid && has_edge(*graphs[i], from, to)) {
											found = i;
											++num_found;
										}
									}
									assert(num_found == 1 || num_found == 0);

									if (has_edge(channel_with_interior_g, from, to) && !has_edge(channel_without_interior_g, from, to)) {
										++num_missing_interior_edges;
									}
									char buffer[256];

									sprintf_rr_node(from, buffer);
									zlog_debug(missing_edge_log, "\t Edge %d %s ->", rr_edge, buffer);

									sprintf_rr_node(to, buffer);
									zlog_debug(missing_edge_log, " %s is in graph %d\n", buffer, found);

									++num_missing_edges;
								}
							}
							zlog_debug(missing_edge_log, "Net %d num missing edges = %d/%d (%g) num missing interior edges = %d/%d (%g)\n\n", net->vpr_id, num_missing_edges, num_edges(temp_route_tree.graph), num_missing_edges*100.0/num_edges(temp_route_tree.graph), num_missing_interior_edges, num_missing_edges, num_missing_interior_edges*100.0/num_missing_edges);

							struct bfs_router_t {
								set<int> visited_nodes;
								set<int> visited_edges;
								vector<int> pred;

								bfs_router_t(int num_vertices) :
									pred(num_vertices, -1)
								{
								}

								bool tree_edge(int e, const RRGraph &g)
								{
									int to = get_target(g, e);
									assert(pred[to] == -1);
									pred[to] = get_source(g, e);
									assert(visited_edges.find(e) == visited_edges.end());
									visited_edges.insert(e);
									return true;
								}

								void examine_edge(int e, const RRGraph &g)
								{
								}

								void discover_vertex(int v, const RRGraph &g)
								{
									assert(visited_nodes.find(v) == visited_nodes.end());
									visited_nodes.insert(v);
								}

								void examine_vertex(int v, const RRGraph &g)
								{
								}
							} router_visitor(num_vertices(*graphs[tid]));

							vector<VertexColor> color(num_vertices(*graphs[tid]), VertexColor::WHITE);
							bfs(*graphs[tid], { static_cast<unsigned long>(net->source.rr_node) }, color, router_visitor);
							assert(router_visitor.visited_edges.size()+1 == router_visitor.visited_nodes.size());
							int lmao = router_visitor.visited_nodes.size();
							for (int i = 0; i < sinks.size(); ++i) {
								if (router_visitor.visited_nodes.find(sinks[i]->rr_node) != router_visitor.visited_nodes.end()) {
									++local_bfs_num_sinks_routed;
								} else {
									//for (const auto &c : all_partition_components[tid]) {
										//assert(find(begin(c), end(c), net->source.rr_node) == end(c)
												//|| find(begin(c), end(c), sinks[i]->rr_node) == end(c));
									//}
								}
							}
							//assert(local_bfs_num_sinks_routed == sinks.size());
							//assert(local_bfs_num_sinks_routed == num_sinks_routed);

							graph_used_by_net[net->vpr_id] = tid;
						}

						++local_num_nets_to_route;

						local_num_sinks_to_route += sinks.size();
						local_num_sinks_routed += num_sinks_routed;
					}

					//local_perf.total_route_time += clock::now()-rip_up_start;
				}

				greedy_end_time[tid] = clock::now();

				perfs[tid] = local_perf;
				lock_perfs[tid] = local_lock_perf; 
				thread_num_nets_routed[tid] = local_num_nets_routed;
				thread_num_nets_to_route[tid] = local_num_nets_to_route;
				thread_num_sinks_routed[tid] = local_num_sinks_routed;
				thread_num_sinks_to_route[tid] = local_num_sinks_to_route;
				thread_bfs_num_sinks_routed[tid] = local_bfs_num_sinks_routed;

				//debug_lock[tid].unlock();
				});

		greedy_route_time = clock::now()-greedy_route_start;

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
		for (int i = 0; i < opts->num_threads; ++i) {
			printf("Thread %d num nets routed: %d/%d (%g)\n", i, thread_num_nets_routed[i], thread_num_nets_to_route[i], thread_num_nets_routed[i]*100.0/thread_num_nets_to_route[i]);
			printf("Thread %d num sinks routed: %d/%d (%g)\n", i, thread_num_sinks_routed[i], thread_num_sinks_to_route[i], thread_num_sinks_routed[i]*100.0/thread_num_sinks_to_route[i]);
			total_num_nets_to_route += thread_num_nets_to_route[i];
			total_num_nets_routed += thread_num_nets_routed[i];
			total_num_sinks_routed += thread_num_sinks_routed[i];
			total_num_sinks_to_route += thread_num_sinks_to_route[i];
			total_bfs_num_sinks_routed += thread_bfs_num_sinks_routed[i];
		}
		printf("Total num nets routed: %d/%d (%g)\n", total_num_nets_routed, total_num_nets_to_route, total_num_nets_routed*100.0/total_num_nets_to_route);
		printf("Total num sinks routed: %d/%d (%g)\n", total_num_sinks_routed, total_num_sinks_to_route, total_num_sinks_routed*100.0/total_num_sinks_to_route);
		printf("Total BFS num sinks routed: %d/%d (%g)\n", total_bfs_num_sinks_routed, total_num_sinks_to_route, total_bfs_num_sinks_routed*100.0/total_num_sinks_to_route);

		/* checking */
		for (int i = 0; i < num_vertices(*graphs[0]); ++i) {
			congestion[i].recalc_occ = 0; 
		}

		for (const auto &net : nets) {
			check_route_tree(route_trees[net.local_id], net, *graphs[0]);
			recalculate_occ(route_trees[net.local_id], *graphs[0], congestion);
		}

		bool valid = true;
		for (int i = 0; i < num_vertices(*graphs[0]); ++i) {
			sprintf_rr_node(i, buffer);
			if (congestion[i].recalc_occ != congestion[i].occ) {
				zlog_error(delta_log, "Node %s occ mismatch, recalc: %d original: %d\n", buffer, congestion[i].recalc_occ, congestion[i].occ);
				valid = false;
			}
		}
		assert(valid);

		for (const auto &net : nets) {
			vector<int> overused_rr_node;
			get_overused_nodes(route_trees[net.local_id], get_vertex(route_trees[net.local_id].graph, route_trees[net.local_id].root_rt_node_id), *graphs[0], congestion, overused_rr_node);
			if (!overused_rr_node.empty()) {
				zlog_level(delta_log, ROUTER_V1, "Net %d bb_rank %d overused nodes:\n", net.vpr_id, net.bb_area_rank);
				for (const auto &item : overused_rr_node) {
					zlog_level(delta_log, ROUTER_V1, "%d ", item);
				}
				zlog_level(delta_log, ROUTER_V1, "\n");
			}
		}

		iter_start = clock::now();

		if (feasible_routing(*graphs[0], congestion)) {
			//dump_route(*current_traces_ptr, "route.txt");
			routed = true;
		} else {
			unsigned long num_overused_nodes = 0;
			for (int i = 0; i < num_vertices(*graphs[0]); ++i) {
				if (congestion[i].occ > get_vertex(*graphs[0], i).properties.capacity) {
					++num_overused_nodes;
				}
			}
			zlog_info(delta_log, "Num overused nodes: %lu/%d (%.2f)\n", num_overused_nodes, num_vertices(*graphs[0]), num_overused_nodes*100.0/num_vertices(*graphs[0]));
			printf("Num overused nodes: %lu/%d (%.2f)\n", num_overused_nodes, num_vertices(*graphs[0]), num_overused_nodes*100.0/num_vertices(*graphs[0]));

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
				printf("\t\tGreedy route time: %g s.\n", duration_cast<nanoseconds>(greedy_route_time).count() / 1e9);
				printf("\t\tPartitioned route time: %g s.\n", duration_cast<nanoseconds>(partitioned_route_time).count() / 1e9);
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

		for (int i = 0; i < opts->num_threads; ++i) {
			printf("Thread %d num lock waits/tries: %lu/%lu (%g)\n", i, lock_perfs[i].num_lock_waits, lock_perfs[i].num_lock_tries, (float)lock_perfs[i].num_lock_waits*100/lock_perfs[i].num_lock_tries);
			printf("Thread %d total wait time: %g (%g)\n", i, duration_cast<nanoseconds>(lock_perfs[i].total_wait_time).count() / 1e9, 100.0*lock_perfs[i].total_wait_time/(greedy_route_time+partitioned_route_time));
			printf("Thread %d num_heap_pushes: %lu\n", i, perfs[i].num_heap_pushes);
			printf("Thread %d num_heap_pops: %lu\n", i, perfs[i].num_heap_pops);
			printf("Thread %d num_neighbor_visits: %lu\n", i, perfs[i].num_neighbor_visits);
		}

		printf("total_num_heap_pushes: %lu\n", total_num_heap_pushes);
		printf("total_num_heap_pops: %lu\n", total_num_heap_pops); 
		printf("total_num_neighbor_visits: %lu\n", total_num_neighbor_visits);

		clock::time_point greedy_earliest_end_time = *std::min_element(begin(greedy_end_time), end(greedy_end_time));
		clock::time_point partitioned_earliest_end_time = *std::min_element(begin(partitioned_end_time), end(partitioned_end_time));

		for (int i = 0; i < opts->num_threads; ++i) {
			printf("Thread %d greedy wait time %g (%g)\n", i, duration_cast<nanoseconds>(greedy_end_time[i]-greedy_earliest_end_time).count() / 1e9, 100.0*(greedy_end_time[i]-greedy_earliest_end_time)/greedy_route_time);
		}
		for (int i = 0; i < opts->num_threads; ++i) {
			if (partitioned_route_time > clock::duration::zero()) {
				printf("Thread %d partitioned wait time %g (%g)\n", i, duration_cast<nanoseconds>(partitioned_end_time[i]-partitioned_earliest_end_time).count() / 1e9, 100.0*(partitioned_end_time[i]-partitioned_earliest_end_time)/partitioned_route_time);
			}
		}
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
