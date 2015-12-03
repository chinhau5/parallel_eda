#ifndef GRAPH_H
#define GRAPH_H

#include <vector>
#include <algorithm>

using namespace std;

template<typename Properties>
struct edge_t {
	int m_id;
	int a;
	int b;
	Properties properties;	
};

template<typename VertexProperties, typename EdgeProperties>
struct vertex_t {
	int m_id;
	vector<int> edges;
	VertexProperties properties;	
};

template<typename VertexProperties, typename EdgeProperties>
struct graph_t {
	vector<vertex_t<VertexProperties, EdgeProperties>> vertices;
	vector<edge_t<EdgeProperties>> edges;
};

template<typename VertexProperties, typename EdgeProperties>
int id(const vertex_t<VertexProperties, EdgeProperties> &v)
{
	return v.m_id;
}

template<typename EdgeProperties>
int id(const edge_t<EdgeProperties> &e)
{
	return e.m_id;
}

template<typename VertexProperties, typename EdgeProperties>
int num_vertices(const graph_t<VertexProperties, EdgeProperties> &g)
{
	return g.vertices.size();
}

template<typename VertexProperties, typename EdgeProperties>
int num_edges(const graph_t<VertexProperties, EdgeProperties> &g)
{
	return g.edges.size();
}

template<typename VertexProperties, typename EdgeProperties>
int num_out_edges(const graph_t<VertexProperties, EdgeProperties> &g, const vertex_t<VertexProperties, EdgeProperties> &v)
{
	return v.edges.size();
}

template<typename VertexProperties, typename EdgeProperties>
vertex_t<VertexProperties, EdgeProperties> &get_vertex(graph_t<VertexProperties, EdgeProperties> &g, int v)
{
	assert(v < g.vertices.size());
	return g.vertices[v];
}

template<typename VertexProperties, typename EdgeProperties>
const vertex_t<VertexProperties, EdgeProperties> &get_vertex(const graph_t<VertexProperties, EdgeProperties> &g, int v)
{
	assert(v < g.vertices.size());
	return g.vertices[v];
}

template<typename VertexProperties, typename EdgeProperties>
void add_vertex(graph_t<VertexProperties, EdgeProperties> &g, int n = 1)
{
	int start = g.vertices.size();
	g.vertices.resize(g.vertices.size()+n);
	for (int i = start; i < g.vertices.size(); ++i) {
		g.vertices[i].m_id = i;
	}
}

template<typename VertexProperties, typename EdgeProperties>
edge_t<EdgeProperties> &add_edge(graph_t<VertexProperties, EdgeProperties> &g, vertex_t<VertexProperties, EdgeProperties> &v_a, vertex_t<VertexProperties, EdgeProperties> &v_b)
{
	edge_t<EdgeProperties> e;
	e.m_id = g.edges.size();
	e.a = v_a.m_id;
	e.b = v_b.m_id;

	assert(find_if(v_a.edges.begin(), v_a.edges.end(), [&g, &v_b] (int e) -> bool { return g.edges[e].b == v_b.m_id; }) == v_a.edges.end());

	g.edges.push_back(e);
	v_a.edges.push_back(e.m_id);

	return g.edges[e.m_id];
}

template<typename VertexProperties, typename EdgeProperties>
void remove_edge(graph_t<VertexProperties, EdgeProperties> &g, const edge_t<EdgeProperties> &e_del)
{
	for (auto &v : g.vertices) {
		for (auto &e : v.edges) {
			if (e > e_del.m_id) {
				--e;
			}
		}
	}

	for (auto &e : g.edges) {
		if (e.m_id > e_del.m_id) {
			--e.m_id;
		}
	}

	auto iter = find_if(g.vertices[e_del.a].edges.begin(), g.vertices[e_del.a].edges.end(), [&g, &e_del] (int e) -> bool { return g.edges[e].b == e_del.b; });
	assert(iter != g.vertices[e_del.a].edges.end());

	g.vertices[e_del.a].edges.erase(iter);
	g.edges.erase(g.edges.begin() + e_del.m_id);
}

template<typename VertexProperties, typename EdgeProperties>
void remove_edge(graph_t<VertexProperties, EdgeProperties> &g, int a, int b)
{
	auto iter = find_if(g.vertices[a].edges.begin(), g.vertices[a].edges.end(), [&g, &b] (int e) -> bool { return g.edges[e].b == b; });
	assert(iter != g.vertices[a].edges.end());

	for (auto &v : g.vertices) {
		for (auto &e : v.edges) {
			if (e > *iter) {
				--e;
			}
		}
	}

	for (auto &e : g.edges) {
		if (e.m_id > *iter) {
			--e.m_id;
		}
	}

	g.edges.erase(*iter);
	g.vertices[a].edges.erase(iter);
}

template<typename VertexProperties, typename EdgeProperties>
edge_t<EdgeProperties> &get_edge(graph_t<VertexProperties, EdgeProperties> &g, int e)
{
	assert(e < g.edges.size());
	return g.edges[e];
}

template<typename VertexProperties, typename EdgeProperties>
const edge_t<EdgeProperties> &get_edge(const graph_t<VertexProperties, EdgeProperties> &g, int e)
{
	assert(e < g.edges.size());
	return g.edges[e];
}

template<typename VertexProperties, typename EdgeProperties>
const edge_t<EdgeProperties> *get_edge(const graph_t<VertexProperties, EdgeProperties> &g, const vertex_t<VertexProperties, EdgeProperties> &v_a, const vertex_t<VertexProperties, EdgeProperties> &v_b)
{
	auto iter = find_if(v_a.edges.begin(), v_a.edges.end(), [&g, &v_b] (int e) -> bool { return g.edges[e].b == v_b.m_id; });

	const edge_t<EdgeProperties> *res;
	if (iter != v_a.edges.end()) {
		res = &g.edges[*iter];
	} else {
		res = nullptr;
	}
	return res;
}

template<typename VertexProperties, typename EdgeProperties>
const vertex_t<VertexProperties, EdgeProperties> &get_source(const graph_t<VertexProperties, EdgeProperties> &g, const edge_t<EdgeProperties> &e)
{
	assert(e.a < g.vertices.size());
	return g.vertices[e.a];
}

template<typename VertexProperties, typename EdgeProperties>
const vertex_t<VertexProperties, EdgeProperties> &get_target(const graph_t<VertexProperties, EdgeProperties> &g, const edge_t<EdgeProperties> &e)
{
	assert(e.b < g.vertices.size());
	return g.vertices[e.b];
}

template<typename VertexProperties, typename EdgeProperties, typename Func>
void for_all_edges(const graph_t<VertexProperties, EdgeProperties> &g, const Func &f)
{
	for (const auto &e : g.edges) {
		f(e);
	}
}

template<typename VertexProperties, typename EdgeProperties, typename Func>
void for_all_edges(graph_t<VertexProperties, EdgeProperties> &g, const Func &f)
{
	for (auto &e : g.edges) {
		f(e);
	}
}

template<typename VertexProperties, typename EdgeProperties, typename Func>
void for_all_vertices_breakable(graph_t<VertexProperties, EdgeProperties> &g, const Func &f)
{
	for (auto &v : g.vertices) {
		if (!f(v)) {
			break;
		}
	}
}

template<typename VertexProperties, typename EdgeProperties, typename Func>
void for_all_vertices_breakable(const graph_t<VertexProperties, EdgeProperties> &g, const Func &f)
{
	for (const auto &v : g.vertices) {
		if (!f(v)) {
			break;
		}
	}
}

template<typename VertexProperties, typename EdgeProperties, typename Func>
void for_all_vertices(graph_t<VertexProperties, EdgeProperties> &g, const Func &f)
{
	for (auto &v : g.vertices) {
		f(v);
	}
}

template<typename VertexProperties, typename EdgeProperties, typename Func>
void for_all_vertices(const graph_t<VertexProperties, EdgeProperties> &g, const Func &f)
{
	for (auto &v : g.vertices) {
		f(v);
	}
}

template<typename VertexProperties, typename EdgeProperties, typename Func>
void for_all_out_edges(const graph_t<VertexProperties, EdgeProperties> &g, const vertex_t<VertexProperties, EdgeProperties> &v, const Func &f)
{
	assert(v.m_id < g.vertices.size());
	for (const auto &n : g.vertices[v.m_id].edges) {
		f((g).edges[n]);
	}
}

template<typename VertexProperties, typename EdgeProperties, typename Func>
void clear_vertices(graph_t<VertexProperties, EdgeProperties> &g)
{
	g.vertices.clear();
	g.edges.clear();
}

template<typename VertexProperties, typename EdgeProperties>
void clear_edges(graph_t<VertexProperties, EdgeProperties> &g)
{
	for_all_vertices(g, [] (vertex_t<VertexProperties, EdgeProperties> &v) -> void {
			v.edges.clear();
			});
	g.edges.clear();
}

template<typename VertexProperties, typename EdgeProperties, typename VertexLabeler, typename EdgeLabeler, typename VertexFilter>
void write_graph(const graph_t<VertexProperties, EdgeProperties> &g, const char *filename, const VertexLabeler &vertex_label, const EdgeLabeler &edge_label, const VertexFilter &vertex_filter)
{
	FILE *file = fopen(filename, "w");
	fprintf(file, "digraph g {\n");
	for_all_vertices(g, [&g, &vertex_label, &vertex_filter, &file] (const vertex_t<VertexProperties, EdgeProperties> &v) -> void {
		if (!vertex_filter(v)) {
			fprintf(file, "%d [label=\"%s\"", id(v), vertex_label(v).c_str());
			fprintf(file, "];\n");
		}
	});
	for_all_edges(g, [&g, &edge_label, &vertex_filter, &file] (const edge_t<EdgeProperties> &e) -> void {
		const auto &from = get_source(g, e);
		const auto &to = get_target(g, e);
		if (!vertex_filter(from) && !vertex_filter(to)) {
			fprintf(file, "%d -> %d [label=\"%s\"];\n", id(from), id(to), edge_label(e).c_str());
		}
	});
	fprintf(file, "}\n");
	fclose(file);
}


#endif
