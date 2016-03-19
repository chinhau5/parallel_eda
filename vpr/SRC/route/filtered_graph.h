#ifndef FILTERED_GRAPH_H
#define FILTERED_GRAPH_H

#include "graph.h"

template<typename Graph, typename VertexPredicate, typename EdgePredicate>
struct filtered_graph_t {
	using base = Graph;

	const base &g;
	VertexPredicate valid_vertex;
	EdgePredicate valid_edge;

	filtered_graph_t(const base &g, const VertexPredicate &valid_vertex, const EdgePredicate &valid_edge)
		: g(g), valid_vertex(valid_vertex), valid_edge(valid_edge)
	{
	}

	const base &base_graph() const
	{
		return g;
	}

	//const VertexPredicate &vertex_predicate() const
	//{
		//return valid_vertex;
	//}
};

template<typename Graph, typename VertexPredicate, typename EdgePredicate>
filtered_graph_t<Graph, VertexPredicate, EdgePredicate> make_filtered_graph(const Graph &g, const VertexPredicate &valid_vertex, const EdgePredicate &valid_edge)
{
	return filtered_graph_t<Graph, VertexPredicate, EdgePredicate>(g, valid_vertex, valid_edge);
}

template<typename Graph, typename VertexPredicate, typename EdgePredicate, typename VertexPredicate2, typename EdgePredicate2>
filtered_graph_t<Graph, VertexPredicate2, EdgePredicate2> make_filtered_graph(const filtered_graph_t<Graph, VertexPredicate, EdgePredicate> &fg, const VertexPredicate2 &valid_vertex, const EdgePredicate2 &valid_edge)
{
	return make_filtered_graph(fg.g, valid_vertex, valid_edge);
}

template<typename Graph, typename VertexPredicate, typename EdgePredicate>
//boost::iterator_range<boost::filter_iterator<graph_filter_t, typename graph_t<VertexProperties, EdgeProperties>::vertex_iterator>>
auto
get_vertices(const filtered_graph_t<Graph, VertexPredicate, EdgePredicate> &fg) 
{
	//graph_filter_t filter(fg.valid_vertex);

	using vertex_iterator = typename Graph::vertex_iterator;

	return boost::make_iterator_range(
			boost::make_filter_iterator(fg.valid_vertex, vertex_iterator(0ul), vertex_iterator(fg.g.vertices.size())),
			boost::make_filter_iterator(fg.valid_vertex, vertex_iterator(fg.g.vertices.size()), vertex_iterator(fg.g.vertices.size()))
				);
}

template<typename Graph, typename VertexPredicate, typename EdgePredicate>
//boost::iterator_range<boost::filter_iterator<graph_filter_t, typename graph_t<VertexProperties, EdgeProperties>::edge_iterator>>
auto 
get_edges(const filtered_graph_t<Graph, VertexPredicate, EdgePredicate> &fg) 
{
	//graph_filter_t filter(fg.valid_edge);

	using edge_iterator = typename Graph::const_edge_iterator;

	const auto &edges = get_edges(fg.g);

	return boost::make_iterator_range(
			boost::make_filter_iterator(fg.valid_edge, edge_iterator(begin(edges)), edge_iterator(end(edges))),
			boost::make_filter_iterator(fg.valid_edge, edge_iterator(end(edges)), edge_iterator(end(edges)))
			);
}

template<typename Graph, typename VertexPredicate, typename EdgePredicate>
//boost::iterator_range<boost::filter_iterator<graph_filter_t, typename graph_t<VertexProperties, EdgeProperties>::out_edges_iterator>>
auto 
get_out_edges(const filtered_graph_t<Graph, VertexPredicate, EdgePredicate> &fg, int v)
{
	assert(v < fg.g.vertices.size());
	
	const auto &edges = fg.g.vertices[v].edges;

	//graph_filter_t filter(fg.valid_edge);

	return boost::make_iterator_range(
			boost::make_filter_iterator(fg.valid_edge, edges.begin(), edges.end()),
			boost::make_filter_iterator(fg.valid_edge, edges.end(), edges.end())
				);
}

template<typename Graph, typename VertexPredicate, typename EdgePredicate>
int num_vertices(const filtered_graph_t<Graph, VertexPredicate, EdgePredicate> &fg)
{
	return num_vertices(fg.g);
}

template<typename Graph, typename VertexPredicate, typename EdgePredicate>
int num_edges(const filtered_graph_t<Graph, VertexPredicate, EdgePredicate> &fg)
{
	return num_edges(fg.g);
}

template<typename Graph, typename VertexPredicate, typename EdgePredicate, typename Edge>
int get_source(const filtered_graph_t<Graph, VertexPredicate, EdgePredicate> &fg, const Edge &e)
{
	return e.a;
}

template<typename Graph, typename VertexPredicate, typename EdgePredicate, typename Edge>
int get_target(const filtered_graph_t<Graph, VertexPredicate, EdgePredicate> &fg, const Edge &e)
{
	return e.b;
}

template<typename Graph, typename VertexPredicate, typename EdgePredicate>
auto &get_vertex_props(filtered_graph_t<Graph, VertexPredicate, EdgePredicate> &fg, int v)
{
	assert(v < fg.g.vertices.size() && fg.valid_vertex(v));
	return fg.g.vertices[v].properties;
}

template<typename Graph, typename VertexPredicate, typename EdgePredicate>
const auto &get_vertex_props(const filtered_graph_t<Graph, VertexPredicate, EdgePredicate> &fg, int v)
{
	assert(v < fg.g.vertices.size() && fg.valid_vertex(v));
	return fg.g.vertices[v].properties;
}

#endif
