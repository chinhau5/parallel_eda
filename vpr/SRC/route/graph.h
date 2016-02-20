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
	typedef vector<vertex_t<VertexProperties, EdgeProperties>> Vertices;
	typedef vector<edge_t<EdgeProperties>> Edges;
	Vertices vertices;
	Edges edges;

	struct out_edges_iterator : public std::iterator<
					  typename std::forward_iterator_tag,
					  edge_t<EdgeProperties>,
					  ptrdiff_t,
					  edge_t<EdgeProperties> *,
					  edge_t<EdgeProperties> & 
					  > {
		graph_t &g;
		vector<int>::const_iterator iter;

		out_edges_iterator(graph_t &g, const vector<int>::const_iterator &iter) 
			: g(g), iter(iter)
	   	{
		}

		out_edges_iterator &operator++()
		{
			++iter;
			return *this;
		}

		typename out_edges_iterator::reference operator*() const
		{
			return g.edges[*iter];
		}

		typename out_edges_iterator::pointer operator->() const
		{
			return &g.edges[*iter];
		}

		bool operator==(const out_edges_iterator &other) const
		{
			return other.iter == iter;
		}

		bool operator!=(const out_edges_iterator &other) const
		{
			return other.iter != iter;
		}
	};

	struct out_edges_const_iterator : public std::iterator<
					  typename std::forward_iterator_tag,
					  edge_t<EdgeProperties>,
					  ptrdiff_t,
					  const edge_t<EdgeProperties> *,
					  const edge_t<EdgeProperties> & 
					  > {
		const graph_t &g;
		vector<int>::const_iterator iter;

		out_edges_const_iterator(const graph_t &g, const vector<int>::const_iterator &iter) 
			: g(g), iter(iter)
	   	{
		}

		out_edges_const_iterator &operator++()
		{
			++iter;
			return *this;
		}

		typename out_edges_const_iterator::reference operator*() const
		{
			return g.edges[*iter];
		}

		typename out_edges_const_iterator::pointer operator->() const
		{
			return &g.edges[*iter];
		}

		bool operator==(const out_edges_const_iterator &other) const
		{
			return other.iter == iter;
		}

		bool operator!=(const out_edges_const_iterator &other) const
		{
			return other.iter != iter;
		}
	};

	template<typename Value, typename Base>
	struct iterator : public std::iterator<
					  typename std::forward_iterator_tag,
					  Value,
					  ptrdiff_t,
					  Value *,
					  Value & 
					  > {
		Base iter;

		iterator(const Base &iter) 
			: iter(iter)
	   	{
			//printf("iterator constructor\n");
		}

		//iterator(const iterator &other) 
			//: iter(other.iter)
		   //{
			//printf("iterator copy constructor\n");
		//}

		//void operator=(const iterator &other) 
		   //{
			//iter = other.iter;
			//printf("iterator assignment\n");
		//}

		iterator &operator++()
		{
			++iter;
			return *this;
		}

		typename iterator::reference operator*() const
		{
			return *iter;
		}

		typename iterator::pointer operator->() const
		{
			return &(*iter);
		}

		bool operator==(const iterator &other) const
		{
			return other.iter == iter;
		}

		bool operator!=(const iterator &other) const
		{
			return other.iter != iter;
		}
	};

	typedef iterator<edge_t<EdgeProperties>, typename Edges::iterator> edge_iterator;
	typedef iterator<const edge_t<EdgeProperties>, typename Edges::const_iterator> const_edge_iterator;

	typedef iterator<vertex_t<VertexProperties, EdgeProperties>, typename Vertices::iterator> vertex_iterator;
	typedef iterator<const vertex_t<VertexProperties, EdgeProperties>, typename Vertices::const_iterator> const_vertex_iterator;
	//typedef typename Edges::iterator edge_iterator;
	//typedef typename Edges::const_iterator edge_const_iterator;

	//typedef typename Vertices::iterator vertex_iterator;
	//typedef typename Vertices::const_iterator vertex_const_iterator;
};

//template<typename VertexProperties, typename EdgeProperties>
//bool operator!=(typename graph_t<VertexProperties, EdgeProperties>::const_vertex_iterator &a,
		//typename graph_t<VertexProperties, EdgeProperties>::const_vertex_iterator &b)
//{
	//return a.iter != b.iter;
//}

//template<typename VertexProperties, typename EdgeProperties>
//bool operator!=(typename graph_t<VertexProperties, EdgeProperties>::vertex_iterator &a,
		//typename graph_t<VertexProperties, EdgeProperties>::vertex_iterator &b)
//{
	//return a.iter != b.iter;
//}

//template<typename VertexProperties, typename EdgeProperties, typename Value, typename Base>
//bool operator!=(const typename graph_t<VertexProperties, EdgeProperties>::template iterator<Value, Base> &a,
		//const typename graph_t<VertexProperties, EdgeProperties>::template iterator<Value, Base> &b)
//{
	//return a.iter != b.iter;
//}

template<typename VertexProperties, typename EdgeProperties, typename Value, typename Base>
bool operator==(const typename graph_t<VertexProperties, EdgeProperties>::template iterator<Value, Base> &a,
		const typename graph_t<VertexProperties, EdgeProperties>::template iterator<Value, Base> &b)
{
	return a.iter != b.iter;
}

template<typename VertexProperties, typename EdgeProperties, typename Value, typename Base>
bool operator==(typename graph_t<VertexProperties, EdgeProperties>::template iterator<Value, Base> &a,
		typename graph_t<VertexProperties, EdgeProperties>::template iterator<Value, Base> &b)
{
	return a.iter != b.iter;
}


template<typename Iterator>
struct adapter_t {
	Iterator b;
	Iterator e;

	adapter_t(Iterator b, Iterator e)
		: b(b), e(e)
	{
	}

	Iterator begin() const
	{
		return b;
	}

	Iterator end() const
	{
		return e;
	}
};

template<typename VertexProperties, typename EdgeProperties>
adapter_t<typename graph_t<VertexProperties, EdgeProperties>::Edges::const_iterator>
get_edges(const graph_t<VertexProperties, EdgeProperties> &g) 
{
	return
		adapter_t<typename graph_t<VertexProperties, EdgeProperties>::Edges::const_iterator>
		(g.edges.cbegin(), g.edges.cend());
}

template<typename VertexProperties, typename EdgeProperties>
adapter_t<typename graph_t<VertexProperties, EdgeProperties>::Edges::iterator>
get_edges(graph_t<VertexProperties, EdgeProperties> &g) 
{
	return
		adapter_t<typename graph_t<VertexProperties, EdgeProperties>::Edges::iterator>
		(g.edges.begin(), g.edges.end());
}

template<typename VertexProperties, typename EdgeProperties>
adapter_t<typename graph_t<VertexProperties, EdgeProperties>::Vertices::const_iterator>
get_vertices(const graph_t<VertexProperties, EdgeProperties> &g) 
{
	return
		adapter_t<typename graph_t<VertexProperties, EdgeProperties>::Vertices::const_iterator>
		(g.vertices.cbegin(), g.vertices.cend());
}

template<typename VertexProperties, typename EdgeProperties>
adapter_t<typename graph_t<VertexProperties, EdgeProperties>::Vertices::iterator>
get_vertices(graph_t<VertexProperties, EdgeProperties> &g) 
{
	return
		adapter_t<typename graph_t<VertexProperties, EdgeProperties>::Vertices::iterator>
		(g.vertices.begin(), g.vertices.end());
}

template<typename VertexProperties, typename EdgeProperties>
adapter_t<typename graph_t<VertexProperties, EdgeProperties>::out_edges_iterator>
get_out_edges(graph_t<VertexProperties, EdgeProperties> &g, int v)
{
	const auto &edges = g.vertices[v].edges;

	return adapter_t<typename graph_t<VertexProperties, EdgeProperties>::out_edges_iterator>
		(typename graph_t<VertexProperties, EdgeProperties>::out_edges_iterator(g, edges.begin()),
		 typename graph_t<VertexProperties, EdgeProperties>::out_edges_iterator(g, edges.end()));
}

template<typename VertexProperties, typename EdgeProperties>
adapter_t<typename graph_t<VertexProperties, EdgeProperties>::out_edges_const_iterator>
get_out_edges(const graph_t<VertexProperties, EdgeProperties> &g, int v)
{
	const auto &edges = g.vertices[v].edges;

	return adapter_t<typename graph_t<VertexProperties, EdgeProperties>::out_edges_const_iterator>
		(typename graph_t<VertexProperties, EdgeProperties>::out_edges_const_iterator(g, edges.begin()),
		 typename graph_t<VertexProperties, EdgeProperties>::out_edges_const_iterator(g, edges.end()));
}

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
int num_out_edges(const graph_t<VertexProperties, EdgeProperties> &g, int v)
{
	return g.vertices[v].edges.size();
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
edge_t<EdgeProperties> &add_edge(graph_t<VertexProperties, EdgeProperties> &g, int a, int b)
{
	assert(a < num_vertices(g));
	assert(b < num_vertices(g));

	edge_t<EdgeProperties> e;
	e.m_id = g.edges.size();
	e.a = a;
	e.b = b;

	auto &v_a = get_vertex(g, e.a);

	assert(find_if(v_a.edges.begin(), v_a.edges.end(), [&g, b] (int e) -> bool { return g.edges[e].b == b; }) == v_a.edges.end());

	g.edges.push_back(e);
	v_a.edges.push_back(e.m_id);

	return g.edges[e.m_id];
}

//template<typename VertexProperties, typename EdgeProperties>
//edge_t<EdgeProperties> &add_edge(graph_t<VertexProperties, EdgeProperties> &g, const vertex_t<VertexProperties, EdgeProperties> &v_a, const vertex_t<VertexProperties, EdgeProperties> &v_b)
//{
	//edge_t<EdgeProperties> e;
	//e.m_id = g.edges.size();
	//e.a = v_a.m_id;
	//e.b = v_b.m_id;

	//auto &real_v_a = get_vertex(g, e.a);

	//assert(find_if(real_v_a.edges.begin(), real_v_a.edges.end(), [&g, &v_b] (int e) -> bool { return g.edges[e].b == v_b.m_id; }) == real_v_a.edges.end());

	//g.edges.push_back(e);
	//real_v_a.edges.push_back(e.m_id);

	//return g.edges[e.m_id];
//}

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

//template<typename VertexProperties, typename EdgeProperties>
//void remove_edge(graph_t<VertexProperties, EdgeProperties> &g, int a, int b)
//{
	//auto iter = find_if(g.vertices[a].edges.begin(), g.vertices[a].edges.end(), [&g, &b] (int e) -> bool { return g.edges[e].b == b; });
	//assert(iter != g.vertices[a].edges.end());

	//for (auto &v : g.vertices) {
		//for (auto &e : v.edges) {
			//if (e > *iter) {
				//--e;
			//}
		//}
	//}

	//for (auto &e : g.edges) {
		//if (e.m_id > *iter) {
			//--e.m_id;
		//}
	//}

	//g.edges.erase(*iter);
	//g.vertices[a].edges.erase(iter);
//}

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
const edge_t<EdgeProperties> &get_edge(const graph_t<VertexProperties, EdgeProperties> &g, int a, int b)
{
	const auto &v_a = get_vertex(g, a);
	auto iter = find_if(v_a.edges.begin(), v_a.edges.end(), [&g, &b] (int e) -> bool { return g.edges[e].b == b; });
	assert(iter != v_a.edges.end());
	return g.edges[*iter];
}

template<typename VertexProperties, typename EdgeProperties>
edge_t<EdgeProperties> &get_edge(graph_t<VertexProperties, EdgeProperties> &g, int a, int b)
{
	const auto &v_a = get_vertex(g, a);
	auto iter = find_if(v_a.edges.begin(), v_a.edges.end(), [&g, &b] (int e) -> bool { return g.edges[e].b == b; });
	assert(iter != v_a.edges.end());
	return g.edges[*iter];
}

template<typename VertexProperties, typename EdgeProperties>
bool has_edge(const graph_t<VertexProperties, EdgeProperties> &g, int a, int b)
{
	const auto &v_a = get_vertex(g, a);
	auto iter = find_if(v_a.edges.begin(), v_a.edges.end(), [&g, &b] (int e) -> bool { return g.edges[e].b == b; });
	return iter != v_a.edges.end();
}

//template<typename VertexProperties, typename EdgeProperties>
//const edge_t<EdgeProperties> *get_edge(const graph_t<VertexProperties, EdgeProperties> &g, const vertex_t<VertexProperties, EdgeProperties> &v_a, const vertex_t<VertexProperties, EdgeProperties> &v_b)
//{
	//auto iter = find_if(v_a.edges.begin(), v_a.edges.end(), [&g, &v_b] (int e) -> bool { return g.edges[e].b == v_b.m_id; });

	//const edge_t<EdgeProperties> *res;
	//if (iter != v_a.edges.end()) {
		//res = &g.edges[*iter];
	//} else {
		//res = nullptr;
	//}
	//return res;
//}

template<typename VertexProperties, typename EdgeProperties>
const vertex_t<VertexProperties, EdgeProperties> &get_source(const graph_t<VertexProperties, EdgeProperties> &g, const edge_t<EdgeProperties> &e)
{
	assert(e.a < g.vertices.size());
	return g.vertices[e.a];
}

template<typename VertexProperties, typename EdgeProperties>
vertex_t<VertexProperties, EdgeProperties> &get_source(graph_t<VertexProperties, EdgeProperties> &g, const edge_t<EdgeProperties> &e)
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

template<typename VertexProperties, typename EdgeProperties>
vertex_t<VertexProperties, EdgeProperties> &get_target(graph_t<VertexProperties, EdgeProperties> &g, const edge_t<EdgeProperties> &e)
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
void write_graph(const graph_t<VertexProperties, EdgeProperties> &g, const char *filename, const VertexLabeler &vertex_attrs, const EdgeLabeler &edge_label, const VertexFilter &vertex_filter)
{
	FILE *file = fopen(filename, "w");
	fprintf(file, "digraph g {\n");
	for_all_vertices(g, [&g, &vertex_attrs, &vertex_filter, &file] (const vertex_t<VertexProperties, EdgeProperties> &v) -> void {
		if (!vertex_filter(v)) {
			fprintf(file, "%d [%s];\n", id(v), vertex_attrs(v).c_str());
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
