#ifndef FILTERED_GRAPH_H
#define FILTERED_GRAPH_H

#include "graph.h"

template<typename VertexProperties, typename EdgeProperties, typename VertexPredicate, typename EdgePredicate>
struct filtered_graph_t {
	using base = graph_t<VertexProperties, EdgeProperties>;

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

template<typename VertexProperties, typename EdgeProperties, typename VertexPredicate, typename EdgePredicate>
filtered_graph_t<VertexProperties, EdgeProperties, VertexPredicate, EdgePredicate> make_filtered_graph(const graph_t<VertexProperties, EdgeProperties> &g, const VertexPredicate &valid_vertex, const EdgePredicate &valid_edge)
{
	return filtered_graph_t<VertexProperties, EdgeProperties, VertexPredicate, EdgePredicate>(g, valid_vertex, valid_edge);
}

template<typename VertexProperties, typename EdgeProperties, typename VertexPredicate, typename EdgePredicate, typename VertexPredicate2, typename EdgePredicate2>
filtered_graph_t<VertexProperties, EdgeProperties, VertexPredicate2, EdgePredicate2> make_filtered_graph(const filtered_graph_t<VertexProperties, EdgeProperties, VertexPredicate, EdgePredicate> &fg, const VertexPredicate2 &valid_vertex, const EdgePredicate2 &valid_edge)
{
	return make_filtered_graph(fg.g, valid_vertex, valid_edge);
}

template<typename VertexProperties, typename EdgeProperties, typename VertexPredicate, typename EdgePredicate>
//boost::iterator_range<boost::filter_iterator<graph_filter_t, typename graph_t<VertexProperties, EdgeProperties>::vertex_iterator>>
auto
get_vertices(const filtered_graph_t<VertexProperties, EdgeProperties, VertexPredicate, EdgePredicate> &fg) 
{
	//graph_filter_t filter(fg.valid_vertex);

	using vertex_iterator = typename graph_t<VertexProperties, EdgeProperties>::vertex_iterator;

	return boost::make_iterator_range(
			boost::make_filter_iterator(fg.valid_vertex, vertex_iterator(0ul), vertex_iterator(fg.g.vertices.size())),
			boost::make_filter_iterator(fg.valid_vertex, vertex_iterator(fg.g.vertices.size()), vertex_iterator(fg.g.vertices.size()))
				);
}

template<typename VertexProperties, typename EdgeProperties, typename VertexPredicate, typename EdgePredicate>
//boost::iterator_range<boost::filter_iterator<graph_filter_t, typename graph_t<VertexProperties, EdgeProperties>::edge_iterator>>
auto 
get_edges(const filtered_graph_t<VertexProperties, EdgeProperties, VertexPredicate, EdgePredicate> &fg) 
{
	//graph_filter_t filter(fg.valid_edge);

	using edge_iterator = typename graph_t<VertexProperties, EdgeProperties>::edge_iterator;

	return boost::make_iterator_range(
			boost::make_filter_iterator(fg.valid_edge, edge_iterator(0ul), edge_iterator(fg.g.edges.size())),
			boost::make_filter_iterator(fg.valid_edge, edge_iterator(fg.g.edges.size()), edge_iterator(fg.g.edges.size()))
			);
}

template<typename VertexProperties, typename EdgeProperties, typename VertexPredicate, typename EdgePredicate>
//boost::iterator_range<boost::filter_iterator<graph_filter_t, typename graph_t<VertexProperties, EdgeProperties>::out_edges_iterator>>
auto 
get_out_edges(const filtered_graph_t<VertexProperties, EdgeProperties, VertexPredicate, EdgePredicate> &fg, int v)
{
	assert(v < fg.g.vertices.size());
	
	const auto &edges = fg.g.vertices[v].edges;

	//graph_filter_t filter(fg.valid_edge);

	return boost::make_iterator_range(
			boost::make_filter_iterator(fg.valid_edge, edges.begin(), edges.end()),
			boost::make_filter_iterator(fg.valid_edge, edges.end(), edges.end())
				);
}

template<typename VertexProperties, typename EdgeProperties, typename VertexPredicate, typename EdgePredicate>
int num_vertices(const filtered_graph_t<VertexProperties, EdgeProperties, VertexPredicate, EdgePredicate> &fg)
{
	return fg.g.vertices.size();
}

template<typename VertexProperties, typename EdgeProperties, typename VertexPredicate, typename EdgePredicate>
int num_edges(const filtered_graph_t<VertexProperties, EdgeProperties, VertexPredicate, EdgePredicate> &fg)
{
	return fg.g.edges.size();
}

template<typename VertexProperties, typename EdgeProperties, typename VertexPredicate, typename EdgePredicate>
int get_source(const filtered_graph_t<VertexProperties, EdgeProperties, VertexPredicate, EdgePredicate> &fg, int e)
{
	assert(e < fg.g.edges.size() && fg.valid_edge(e));
	return fg.g.edges[e].a;
}

template<typename VertexProperties, typename EdgeProperties, typename VertexPredicate, typename EdgePredicate>
int get_source(filtered_graph_t<VertexProperties, EdgeProperties, VertexPredicate, EdgePredicate> &fg, int e)
{
	assert(e < fg.g.edges.size() && fg.valid_edge(e));
	return fg.g.edges[e].a;
}

template<typename VertexProperties, typename EdgeProperties, typename VertexPredicate, typename EdgePredicate>
int get_target(const filtered_graph_t<VertexProperties, EdgeProperties, VertexPredicate, EdgePredicate> &fg, int e)
{
	assert(e < fg.g.edges.size() && fg.valid_edge(e));
	return fg.g.edges[e].b;
}

template<typename VertexProperties, typename EdgeProperties, typename VertexPredicate, typename EdgePredicate>
int get_target(filtered_graph_t<VertexProperties, EdgeProperties, VertexPredicate, EdgePredicate> &fg, int e)
{
	assert(e < fg.g.edges.size() && fg.valid_edge(e));
	return fg.g.edges[e].b;
}

template<typename VertexProperties, typename EdgeProperties, typename VertexPredicate, typename EdgePredicate>
vertex_t<VertexProperties, EdgeProperties> &get_vertex(filtered_graph_t<VertexProperties, EdgeProperties, VertexPredicate, EdgePredicate> &fg, int v)
{
	assert(v < fg.g.vertices.size() && fg.valid_vertex(v));
	return fg.g.vertices[v];
}

template<typename VertexProperties, typename EdgeProperties, typename VertexPredicate, typename EdgePredicate>
const vertex_t<VertexProperties, EdgeProperties> &get_vertex(const filtered_graph_t<VertexProperties, EdgeProperties, VertexPredicate, EdgePredicate> &fg, int v)
{
	assert(v < fg.g.vertices.size() && fg.valid_vertex(v));
	return fg.g.vertices[v];
}

#endif
