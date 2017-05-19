#ifndef CACHE_GRAPH_H
#define CACHE_GRAPH_H
#include <vector>
#include <algorithm>
#include <boost/iterator/counting_iterator.hpp>
#include <boost/iterator/iterator_facade.hpp>
#include <boost/range/counting_range.hpp>

//using namespace std;

template<typename Properties>
struct cache_edge_t {
	int a;
	int b;
	int properties;	

	bool operator==(const cache_edge_t &other) const
	{
		return properties == properties;
	}

	bool operator<(const cache_edge_t &other) const
	{
		return properties < other.properties;
	}
};

template<typename Properties>
struct cache_edge_storage_t {

};

template<typename VertexProperties, typename EdgeProperties>
struct cache_vertex_t {
	//int m_id;
	std::vector<cache_edge_t<EdgeProperties>> edges;
	std::vector<EdgeProperties> edge_props;
	VertexProperties properties;
};

//template<typename EdgeProperties>
//struct edge_to_target_t {
	//int operator()(const cache_edge_t<EdgeProperties> &e) const
	//{
		//return e.b;
	//}
//};

template<typename VertexProperties, typename EdgeProperties>
struct cache_graph_t {
	typedef std::vector<cache_vertex_t<VertexProperties, EdgeProperties>> Vertices;
	//typedef vector<cache_edge_t<EdgeProperties>> Edges;
	Vertices vertices;
	//Edges edges;

	using base = cache_graph_t<VertexProperties, EdgeProperties>;

	typedef VertexProperties vertex_properties;
	typedef EdgeProperties edge_properties;

	base &base_graph()
	{
		return *this;
	}

	const base &base_graph() const
	{
		return *this;
	}

	static constexpr int null_vertex()
	{
		return -1;
	}

	static cache_edge_t<EdgeProperties> null_edge()
	{
		return { -1, -1, -1 };
	}

	template<typename Graph, typename Value>
	class base_edge_iterator : public boost::iterator_facade<base_edge_iterator<Graph, Value>, Value, boost::forward_traversal_tag> {
		public:
			base_edge_iterator(Graph &_g, bool begin)
				: g(_g)
			{
				if (begin) {
					current_v = 0;
					while (current_v < num_vertices(g) && num_out_edges(g, current_v) == 0) {
						++current_v;
					}

					if (current_v < num_vertices(g)) {
						current_e = 0;
						end_e = num_out_edges(g, current_v);
					} else {
						assert(current_v == num_vertices(g));
						current_e = -1;
						end_e = -1;
					}
				} else {
					current_v = num_vertices(g);
					current_e = -1;
					end_e = -1;
				}
			}

		private:
			friend class boost::iterator_core_access;

			void increment()
			{
				if (current_e < end_e) {
					++current_e;

					if (current_e == end_e) {
						++current_v;

						while (current_v < num_vertices(g) && num_out_edges(g, current_v) == 0) {
							++current_v;
						}

						if (current_v < num_vertices(g)) {
							current_e = 0;
							end_e = num_out_edges(g, current_v);
						} else {
							assert(current_v == num_vertices(g));
							current_e = -1;
							end_e = -1;
						}
					}
				}
			}

			bool equal(const base_edge_iterator &other) const
			{
				return current_v == other.current_v && current_e == other.current_e;
			}

			Value &dereference() const
			{
				return g.vertices[current_v].edges[current_e];
			}
			
			Graph &g;
			int current_v;
			int current_e;
			int end_e;
	};

	typedef base_edge_iterator<const cache_graph_t<VertexProperties, EdgeProperties>, const cache_edge_t<EdgeProperties>> edge_iterator;

	//class const_edge_iterator : public boost::iterator_facade<const_edge_iterator, const cache_edge_t<EdgeProperties>, boost::forward_traversal_tag> {
		//public:
			//const_edge_iterator(const cache_graph_t<VertexProperties, EdgeProperties> &_g, int _current_v, int _end_v)
				//: g(_g), current_v(_current_v), end_v(_end_v)
			//{
				//if (current_v < end_v) {
					//while (current_v < end_v && num_out_edges(g, current_v) == 0) {
						//++current_v;
					//}

					//if (current_v < end_v) {
						//current_e = 0;
						//end_e = g.vertices[current_v].edges.size();
					//} else {
						//current_e = 0;
						//end_e = 0;
					//}
				//} else {
					//current_e = 0;
					//end_e = 0;
				//}
			//}

		//private:
			//friend class boost::iterator_core_access;

			//void increment()
			//{
				//if (current_e < end_e) {
					//++current_e;

					//if (current_e == end_e) {
						//++current_v;

						//while (current_v < end_v && num_out_edges(g, current_v) == 0) {
							//++current_v;
						//}

						//if (current_v < end_v) {
							//current_e = 0;
							//end_e = num_out_edges(g, current_v);
						//}
					//}
				//}
			//}

			//bool equal(const const_edge_iterator &other) const
			//{
				//bool eq = current_v == other.current_v 
				//&& (current_v == end_v || current_e == other.current_e);
				//return eq;
			//}

			//const cache_edge_t<EdgeProperties> &dereference() const
			//{
				//return g.vertices[current_v].edges[current_e];
			//}
			
			//const cache_graph_t<VertexProperties, EdgeProperties> &g;
			//int current_v;
			//int end_v;
			//int current_e;
			//int end_e;
	//};

	//class edge_iterator : public boost::iterator_facade<edge_iterator, cache_edge_t<EdgeProperties>, boost::forward_traversal_tag> {
		//public:
			//edge_iterator(cache_graph_t<VertexProperties, EdgeProperties> &_g, int _current_v, int _end_v)
				//: g(_g), current_v(_current_v), end_v(_end_v)
			//{
				//if (current_v < end_v) {
					//while (current_v < end_v && num_out_edges(g, current_v) == 0) {
						//++current_v;
					//}

					//if (current_v < end_v) {
						//current_e = 0;
						//end_e = g.vertices[current_v].edges.size();
					//} else {
						//current_e = 0;
						//end_e = 0;
					//}
				//} else {
					//current_e = 0;
					//end_e = 0;
				//}
			//}

		//private:
			//friend class boost::iterator_core_access;

			//void increment()
			//{
				//if (current_e < end_e) {
					//++current_e;

					//if (current_e == end_e) {
						//++current_v;

						//while (current_v < end_v && num_out_edges(g, current_v) == 0) {
							//++current_v;
						//}

						//if (current_v < end_v) {
							//current_e = 0;
							//end_e = num_out_edges(g, current_v);
						//}
					//}
				//}
			//}

			//bool equal(const edge_iterator &other) const
			//{
				//bool eq = current_v == other.current_v 
				//&& (current_v == end_v || current_e == other.current_e);
				//return eq;
			//}

			//cache_edge_t<EdgeProperties> &dereference() const
			//{
				//return g.vertices[current_v].edges[current_e];
			//}
			
			//cache_graph_t<VertexProperties, EdgeProperties> &g;
			//int current_v;
			//int end_v;
			//int current_e;
			//int end_e;
	//};

	//struct out_edge_iterator : public std::iterator<
					  //typename std::forward_iterator_tag,
					  //cache_edge_t<EdgeProperties>,
					  //ptrdiff_t,
					  //cache_edge_t<EdgeProperties> *,
					  //cache_edge_t<EdgeProperties> & 
					  //> {
		//cache_graph_t &g;
		//vector<int>::const_iterator iter;

		//out_edge_iterator(cache_graph_t &g, const vector<int>::const_iterator &iter) 
			//: g(g), iter(iter)
		   //{
		//}

		//out_edge_iterator &operator++()
		//{
			//++iter;
			//return *this;
		//}

		//typename out_edge_iterator::reference operator*() const
		//{
			//return g.edges[*iter];
		//}

		//typename out_edge_iterator::pointer operator->() const
		//{
			//return &g.edges[*iter];
		//}

		//bool operator==(const out_edge_iterator &other) const
		//{
			//return other.iter == iter;
		//}

		//bool operator!=(const out_edge_iterator &other) const
		//{
			//return other.iter != iter;
		//}
	//};

	//struct out_edges_const_iterator : public std::iterator<
					  //typename std::forward_iterator_tag,
					  //cache_edge_t<EdgeProperties>,
					  //ptrdiff_t,
					  //const cache_edge_t<EdgeProperties> *,
					  //const cache_edge_t<EdgeProperties> & 
					  //> {
		//const cache_graph_t &g;
		//vector<int>::const_iterator iter;

		//out_edges_const_iterator(const cache_graph_t &g, const vector<int>::const_iterator &iter) 
			//: g(g), iter(iter)
		   //{
		//}

		//out_edges_const_iterator &operator++()
		//{
			//++iter;
			//return *this;
		//}

		//typename out_edges_const_iterator::reference operator*() const
		//{
			//return g.edges[*iter];
		//}

		//typename out_edges_const_iterator::pointer operator->() const
		//{
			//return &g.edges[*iter];
		//}

		//bool operator==(const out_edges_const_iterator &other) const
		//{
			//return other.iter == iter;
		//}

		//bool operator!=(const out_edges_const_iterator &other) const
		//{
			//return other.iter != iter;
		//}
	//};

	//template<typename Value, typename Base>
	//struct iterator : public std::iterator<
					  //typename std::forward_iterator_tag,
					  //Value,
					  //ptrdiff_t,
					  //Value,
					  //Value 
					  //> {
		//Base iter;

		//iterator(const Base &iter) 
			//: iter(iter)
		   //{
			////printf("iterator constructor\n");
		//}

		////iterator(const iterator &other) 
			////: iter(other.iter)
		   ////{
			////printf("iterator copy constructor\n");
		////}

		////void operator=(const iterator &other) 
		   ////{
			////iter = other.iter;
			////printf("iterator assignment\n");
		////}

		//iterator &operator++()
		//{
			//++iter;
			//return *this;
		//}

		//typename iterator::reference operator*() const
		//{
			//return iter;
		//}

		////typename iterator::pointer operator->() const
		////{
			////return &(*iter);
		////}

		//bool operator==(const iterator &other) const
		//{
			//return other.iter == iter;
		//}

		//bool operator!=(const iterator &other) const
		//{
			//return other.iter != iter;
		//}
	//};

	typedef boost::counting_iterator<unsigned long> vertex_iterator;

	//typedef boost::transform_iterator<edge_to_target_t<EdgeProperties>, typename vector<cache_edge_t<EdgeProperties>>::const_iterator> out_edge_iterator;
	typedef typename std::vector<cache_edge_t<EdgeProperties>>::const_iterator out_edge_iterator;
};

template<typename VertexProperties, typename EdgeProperties>
boost::iterator_range<typename cache_graph_t<VertexProperties, EdgeProperties>::edge_iterator>
get_edges(const cache_graph_t<VertexProperties, EdgeProperties> &g) 
{
	return boost::make_iterator_range(
			typename cache_graph_t<VertexProperties, EdgeProperties>::edge_iterator(g, true), 
			typename cache_graph_t<VertexProperties, EdgeProperties>::edge_iterator(g, false) 
			);
}

template<typename VertexProperties, typename EdgeProperties>
boost::iterator_range<typename cache_graph_t<VertexProperties, EdgeProperties>::vertex_iterator>
get_vertices(const cache_graph_t<VertexProperties, EdgeProperties> &g) 
{
	return boost::counting_range(0ul, g.vertices.size());
}

//template<typename VertexProperties, typename EdgeProperties>
//boost::iterator_range<typename cache_graph_t<VertexProperties, EdgeProperties>::out_edge_iterator>
//get_out_edges(const cache_graph_t<VertexProperties, EdgeProperties> &g, int v)
//{
	//assert(v < g.vertices.size());
	
	//const auto &edges = g.vertices[v].edges;

	//return boost::make_iterator_range(
			//make_transform_iterator(edges.begin()), make_transform_iterator(edges.end())); 
//}

template<typename VertexProperties, typename EdgeProperties>
boost::iterator_range<typename cache_graph_t<VertexProperties, EdgeProperties>::out_edge_iterator>
get_out_edges(const cache_graph_t<VertexProperties, EdgeProperties> &g, int v)
{
	assert(v < g.vertices.size());
	
	const auto &edges = g.vertices[v].edges;

	return boost::make_iterator_range(begin(edges), end(edges));
}

//get_neighbors(const cache_graph_t<VertexProperties, EdgeProperties> &g, int v)
//{
//}

template<typename EdgeProperties>
bool valid(const cache_edge_t<EdgeProperties> &e)
{
	return e.a != -1 && e.b != -1;
}

//template<typename EdgeProperties>
//cache_edge_t<EdgeProperties> invalid_edge()
//{
	//return { -1, -1, nullptr };
//}

//template<typename VertexProperties, typename EdgeProperties>
//int id(const cache_vertex_t<VertexProperties, EdgeProperties> &v)
//{
	//return v.m_id;
//}

//template<typename EdgeProperties>
//int id(const cache_edge_t<EdgeProperties> &e)
//{
	//return e.m_id;
//}

template<typename VertexProperties, typename EdgeProperties>
int num_vertices(const cache_graph_t<VertexProperties, EdgeProperties> &g)
{
	return g.vertices.size();
}

template<typename VertexProperties, typename EdgeProperties>
int num_edges(const cache_graph_t<VertexProperties, EdgeProperties> &g)
{
	int n = 0;
	for (const auto &e : get_edges(g)) {
		++n;
	}
	return n;
}

template<typename VertexProperties, typename EdgeProperties>
int num_out_edges(const cache_graph_t<VertexProperties, EdgeProperties> &g, int v)
{
	return g.vertices[v].edges.size();
}

template<typename VertexProperties, typename EdgeProperties>
VertexProperties &get_vertex_props(cache_graph_t<VertexProperties, EdgeProperties> &g, int v)
{
	assert(v < g.vertices.size());
	return g.vertices[v].properties;
}

template<typename VertexProperties, typename EdgeProperties>
const VertexProperties &get_vertex_props(const cache_graph_t<VertexProperties, EdgeProperties> &g, int v)
{
	assert(v < g.vertices.size());
	return g.vertices[v].properties;
}

template<typename VertexProperties, typename EdgeProperties>
void add_vertex(cache_graph_t<VertexProperties, EdgeProperties> &g, int n = 1)
{
	//int start = g.vertices.size();
	g.vertices.resize(g.vertices.size()+n);
	//for (int i = start; i < g.vertices.size(); ++i) {
		//g.vertices[i].m_id = i;
	//}
}

template<typename VertexProperties, typename EdgeProperties>
const cache_edge_t<EdgeProperties> &add_edge(cache_graph_t<VertexProperties, EdgeProperties> &g, int a, int b)
{
	assert(a < num_vertices(g));
	assert(b < num_vertices(g));

	auto &v_a = g.vertices[a];

	assert(find_if(begin(v_a.edges), end(v_a.edges), [&g, b] (const cache_edge_t<EdgeProperties> &e) -> bool { return e.b == b; }) == end(v_a.edges));

	v_a.edge_props.emplace_back();
	v_a.edges.emplace_back(cache_edge_t<EdgeProperties>{ a, b, (int)(v_a.edge_props.size()-1) });

	return v_a.edges.back();
}

//template<typename VertexProperties, typename EdgeProperties>
//cache_edge_t<EdgeProperties> &add_edge(cache_graph_t<VertexProperties, EdgeProperties> &g, const cache_vertex_t<VertexProperties, EdgeProperties> &v_a, const cache_vertex_t<VertexProperties, EdgeProperties> &v_b)
//{
	//cache_edge_t<EdgeProperties> e;
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
void remove_edge(cache_graph_t<VertexProperties, EdgeProperties> &g, const cache_edge_t<EdgeProperties> &edge)
{
	auto &v_a = g.vertices[edge.a];

	auto iter = find_if(begin(v_a.edges), end(v_a.edges), [&g, &edge] (const cache_edge_t<EdgeProperties> &other) -> bool { return other.b == edge.b; });
	assert(iter != end(v_a.edges));

	v_a.edge_props.erase(begin(v_a.edge_props) + std::distance(begin(v_a.edges), iter));
	v_a.edges.erase(iter);
}

template<typename VertexProperties, typename EdgeProperties>
const cache_edge_t<EdgeProperties> &get_edge(const cache_graph_t<VertexProperties, EdgeProperties> &g, int a, int b)
{
	const auto &v_a = g.vertices[a];
	auto iter = find_if(begin(v_a.edges), end(v_a.edges), [&g, &b] (const cache_edge_t<EdgeProperties> &e) -> bool { return e.b == b; });
	assert(iter != v_a.edges.end());
	return *iter;
}

template<typename VertexProperties, typename EdgeProperties>
const cache_edge_t<EdgeProperties> &get_edge_by_index(const cache_graph_t<VertexProperties, EdgeProperties> &g, int n, int i)
{
	const auto &v = g.vertices[n];
	assert(i < v.edges.size());
	return v.edges[i];
}

template<typename VertexProperties, typename EdgeProperties>
bool has_edge(const cache_graph_t<VertexProperties, EdgeProperties> &g, int a, int b)
{
	const auto &v_a = g.vertices[a];
	auto iter = find_if(begin(v_a.edges), end(v_a.edges), [&g, &b] (const cache_edge_t<EdgeProperties> &e) -> bool { return e.b == b; });
	return iter != v_a.edges.end();
}

template<typename VertexProperties, typename EdgeProperties>
int get_source(const cache_graph_t<VertexProperties, EdgeProperties> &g, const cache_edge_t<EdgeProperties> &e)
{
	return e.a;
}

template<typename VertexProperties, typename EdgeProperties>
int get_target(const cache_graph_t<VertexProperties, EdgeProperties> &g, const cache_edge_t<EdgeProperties> &e)
{
	return e.b;
}

template<typename VertexProperties, typename EdgeProperties>
const EdgeProperties &get_edge_props(const cache_graph_t<VertexProperties, EdgeProperties> &g, const cache_edge_t<EdgeProperties> &e)
{
	return g.vertices[e.a].edge_props[e.properties];
}

template<typename VertexProperties, typename EdgeProperties>
EdgeProperties &get_edge_props(cache_graph_t<VertexProperties, EdgeProperties> &g, const cache_edge_t<EdgeProperties> &e)
{
	return g.vertices[e.a].edge_props[e.properties];
}

template<typename VertexProperties, typename EdgeProperties>
void clear_edges(cache_graph_t<VertexProperties, EdgeProperties> &g)
{
	for (const auto &v : get_vertices(g)) {
		g.vertices[v].edges.clear();
	}
}

template<typename VertexProperties, typename EdgeProperties>
void clear_vertices(cache_graph_t<VertexProperties, EdgeProperties> &g)
{
	//clear_edges(g);
	g.vertices.clear();
}

#endif

