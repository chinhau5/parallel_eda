#ifndef ROUTE_TREE_H
#define ROUTE_TREE_H

#include "log.h"
#include "utility.h"
#include "geometry.h"
#include "new_rr_graph.h"
#include "route.h"
#include <boost/range/iterator_range.hpp>
#include <boost/range/adaptor/filtered.hpp>
#include <mpi.h>

typedef struct rt_edge_property_t {
	RREdge rr_edge;
} rt_edge_property_t;

typedef cache_edge_t<rt_edge_property_t> RouteTreeEdge;

typedef struct rt_node_property_t {
	RRNode rr_node;
	bool valid;
	bool pending_rip_up;
	bool ripped_up;
	//bool branch_point;
	RouteTreeEdge rt_edge_to_parent;
	int num_iterations_fixed;
	//int saved_num_out_edges; [> only valid for SOURCE <]

	//int owner;

	/* public properties */
	bool reexpand;
	//RREdge rr_edge_to_parent;
	float upstream_R;	
	float delay;
	float downstream_C;
	/*float upstream_R_from_route_state;*/
} rt_node_property_t;

typedef cache_graph_t<rt_node_property_t, rt_edge_property_t> RouteTree;
typedef int RouteTreeNode;

typedef struct route_tree_t {
	typedef std::pair<segment, int> rtree_value;

	const RRGraph *rrg;
	RouteTree graph;
	std::vector<int> root_rt_nodes;
	int root_rt_node_id;
	int num_nodes;
	std::map<int, int> rr_node_to_rt_node;
	//map<int, vector<int>> sink_rr_node_to_path;
	std::map<RRNode, std::shared_ptr<std::vector<path_node_t>>> rr_node_to_path;
	//map<int, RouteTreeNode> path_branch_point;
	//bgi::rtree<rtree_value, bgi::rstar<64>> point_tree;
	//box scheduler_bounding_box;
	//std::map<int, std::vector<int>> sink_edges;

	//template<typename Graph, typename Value, typename Base>
	//struct iterator : public std::iterator<
					  //typename std::forward_iterator_tag,
					  //Value,
					  //ptrdiff_t,
					  //Value,
					  //Value 
					  //> { 

		//const Graph &g;
		//Base c;
		//Base e;

		//iterator(const Graph &g, const Base &c, const Base &e) 
			//: g(g), c(c), e(e)
		   //{
		//}

		//iterator &operator++()
		//{
			//do {
				//if (c != e) {
					//++c;
				//}
			//} while (c != e && get_vertex(g, c).properties.valid);
		
			//return *this;
		//}

		//typename iterator::reference operator*() const
		//{
			//return c;
		//}

		//bool operator==(const iterator &other) const
		//{
			//return other.c == c;
		//}

		//bool operator!=(const iterator &other) const
		//{
			//return other.c != c;
		//}
	//};

	//typedef iterator<RouteTree, int, int> vertex_iterator;
	typedef typename RouteTree::vertex_iterator vertex_iterator;
	//typedef iterator<const RouteTreeNode, typename RouteTree::vertex_iterator> vertex_const_iterator;
	typedef typename RouteTree::out_edge_iterator branch_iterator;
	//typedef typename RouteTree::out_edges_const_iterator branch_const_iterator;
} route_tree_t;

void route_tree_init(route_tree_t &rt, const RRGraph *rrg);

bool route_tree_empty(const route_tree_t &rt);

int route_tree_num_nodes(const route_tree_t &rt);

//void route_tree_is_only_path_internal(const RRGraph &g, const RRNode &node);

RouteTreeNode route_tree_add_rr_node(route_tree_t &rt, RRNode rr_node);

//int route_tree_num_nodes(const route_tree_t &rt);

const RouteTreeEdge &route_tree_add_edge_between_rr_node(route_tree_t &rt, RRNode rr_node_a, RRNode rr_node_b);

void route_tree_remove_edge(route_tree_t &rt, const RouteTreeEdge &rt_edge);

void route_tree_remove_node(route_tree_t &rt, RRNode rr_node);

bool route_tree_has_edge(const route_tree_t &rt, RRNode a, RRNode b);

RouteTreeNode route_tree_get_rt_node(const route_tree_t &rt, RRNode rr_node);

struct valid_rt_node {
	const RouteTree &g;
	valid_rt_node(const RouteTree &g) :
		g(g) {}

	bool operator()(unsigned long rt_node) const
	{
		return get_vertex_props(g, rt_node).valid; 
	}
};


boost::iterator_range<boost::filter_iterator<valid_rt_node, route_tree_t::vertex_iterator>>
route_tree_get_nodes(const route_tree_t &rt);

boost::iterator_range<route_tree_t::branch_iterator>
route_tree_get_branches(const route_tree_t &rt, RouteTreeNode rt_node);

void route_tree_mark_paths_to_be_ripped(route_tree_t &rt, const RRGraph &g, const std::vector<int> &pid, int this_pid, const std::vector<RRNode> &rr_nodes);

void route_tree_mark_all_nodes_to_be_ripped(route_tree_t &rt, const RRGraph &g);

void route_tree_mark_paths_to_be_ripped(route_tree_t &rt, const RRGraph &g, const std::vector<RRNode> &rr_nodes);

void route_tree_rip_up_marked(route_tree_t &rt, const RRGraph &g, congestion_locked_t *congestion, float pres_fac, bool lock, lock_perf_t *lock_perf);

void route_tree_rip_up_marked(route_tree_t &rt, const RRGraph &g, congestion_t *congestion, float pres_fac);

void route_tree_rip_up_marked(route_tree_t &rt, const RRGraph &g, congestion_local_t *congestion, float pres_fac);

void route_tree_rip_up_marked_mpi_rma(route_tree_t &rt, const RRGraph &g, const std::vector<int> &pid, int this_pid, congestion_t *congestion, MPI_Win win, float pres_fac);

void route_tree_rip_up_marked_mpi_send_recv(route_tree_t &rt, const RRGraph &g, congestion_t *congestion, float pres_fac, std::queue<RRNode> &cost_update_q);

template<typename Callback>
void route_tree_rip_up_marked_mpi_collective(route_tree_t &rt, const RRGraph &g, congestion_t *congestion, float pres_fac, const Callback &callback)
{
	char buffer[256];
	ongoing_transaction_t trans;
	trans.data = std::make_shared<std::vector<node_update_t>>();

	for (const auto &rt_node : route_tree_get_nodes(rt)) {
		auto &rt_node_p = get_vertex_props(rt.graph, rt_node);

		RRNode rr_node = rt_node_p.rr_node;
		const auto &rr_node_p = get_vertex_props(g, rr_node);
		sprintf_rr_node(rt_node_p.rr_node, buffer);

		/*const auto &bp = rt.path_branch_point.find(rr_node);*/
		/*if (bp != rt.path_branch_point.end()) {*/
			/*rt.path_branch_point.erase(bp);*/
		/*}*/

		if (rt_node_p.pending_rip_up) {
			zlog_level(delta_log, ROUTER_V2, "Ripping up node %s from route tree. Occ: %d Cap: %d\n", buffer, congestion[rt_node_p.rr_node].occ, rr_node_p.capacity);

			/*if (rr_node_p.type == SOURCE) {*/
				/*[>assert(rt_node.saved_num_out_edges > 0);<]*/
				/*[>update_one_cost_internal_mpi_send(rr_node, g, congestion, -num_out_edges(rt.graph, rt_node), pres_fac, this_pid, num_procs, comm); <]*/
				/*update_one_cost_internal(rr_node, g, congestion, -num_out_edges(rt.graph, rt_node), pres_fac); */

				/*d.delta = -num_out_edges(rt.graph, rt_node);*/
			/*} else {*/
				/*[>update_one_cost_internal_mpi_send(rr_node, g, congestion, -1, pres_fac, this_pid, num_procs, comm); <]*/
				/*update_one_cost_internal(rr_node, g, congestion, -1, pres_fac); */

				/*d.delta = -1;*/
			/*}*/

			callback(rr_node);

			const auto &edge = rt_node_p.rt_edge_to_parent;
			if (valid(edge)) {
				/*RouteTreeNode parent_rt_node = get_source(rt.graph, edge);*/
				/*const auto &parent_rt_node_p = get_vertex_props(rt.graph, parent_rt_node);*/
				/*const auto &parent_rr_node_p = get_vertex_props(g, parent_rt_node_p.rr_node);*/
				/* since reconnection back to SOURCE always causes cost to be updated, we need to
				 * update the cost when ripping up also.
				 * if parent is pending rip up, cost will be updated that time. so dont handle it here */
				/*if (parent_rr_node_p.type == SOURCE && !parent_rt_node_p.pending_rip_up && !parent_rt_node_p.ripped_up) {*/
					/*[>update_one_cost_internal_mpi_send(parent_rt_node_p.rr_node, g, congestion, -1, pres_fac, this_pid, num_procs, comm); <]*/
					/*update_one_cost_internal(parent_rt_node_p.rr_node, g, congestion, -1, pres_fac); */

					/*d.rr_node = parent_rt_node_p.rr_node;*/
					/*d.delta = -1;*/

					/*trans.data->push_back(d);*/
				/*}*/

				route_tree_remove_edge(rt, edge);

				assert(!valid(edge));
			} 

			rt_node_p.pending_rip_up = false;
			rt_node_p.ripped_up = true;

			route_tree_remove_node(rt, rr_node);
		} else {
			zlog_level(delta_log, ROUTER_V2, "NOT ripping up node %s from route tree\n", buffer);
			/* invalid assertion because we might have ripped this up from another virtual net */
			/*assert(rt_node_p.ripped_up == false);*/
		}
	}
}

//void route_tree_rip_up_segment(route_tree_t &rt, int sink_rr_node, RRGraph &g, float pres_fac);

//void route_tree_rip_up_segment_2(route_tree_t &rt, int sink_rr_node, RRGraph &g, float pres_fac);

void route_tree_add_path(route_tree_t &rt, const std::shared_ptr<std::vector<path_node_t>> &path, const RRGraph &g, const route_state_t *state = nullptr, bool add_branch_point = true);

//void route_tree_add_path(route_tree_t &rt, const RRGraph &g, const route_state_t *state, const vector<int> &rr_nodes, int vpr_net_id);

//void route_tree_add_path(route_tree_t &rt, const RRGraph &g, const route_state_t *state, int sink_rr_node, int vpr_net_id, vector<int> &added_rr_nodes);

void route_tree_clear(route_tree_t &rt);

void route_tree_set_root(route_tree_t &rt, RRNode rr_node);

void route_tree_add_root(route_tree_t &rt, RRNode rr_node);

//bool route_tree_is_only_path(const route_tree_t &rt, int sink_rr_node, const RRGraph &g);

//pair<int, const RRNode *> route_tree_get_connection(const RRNode &current_rr_node, const RRGraph &g, const route_state_t *state, bool end);

void route_tree_set_node_properties(route_tree_t &rt, const RouteTreeNode &rt_node, bool reexpand, float upstream_R, float delay);

//RouteTreeNode &route_tree_add_or_get_rr_node(route_tree_t &rt, int rr_node_id, const RRGraph &g, const route_state_t *state, bool &update_cost, bool &stop_traceback);

//RouteTreeNode *route_tree_checked_add_rr_node(route_tree_t &rt, const RRNode &rr_node, const route_state_t *state);

//void route_tree_add_to_heap_internal(const route_tree_t &rt, RouteTreeNode rt_node, const RRGraph &g, const RRNode &target, float criticality_fac, float astar_fac, const bounding_box_t &current_bounding_box, std::priority_queue<route_state_t> &heap, perf_t *perf);

void route_tree_add_to_heap(const route_tree_t &rt, const RRGraph &g, RRNode target, float criticality_fac, float astar_fac, const bounding_box_t &current_bounding_box, std::priority_queue<route_state_t> &heap, perf_t *perf);

void route_tree_multi_root_add_to_heap(const route_tree_t &rt, const RRGraph &g, RRNode target, float criticality_fac, float astar_fac, const bounding_box_t &current_bounding_box, std::priority_queue<route_state_t> &heap, perf_t *perf);

template<typename Congestion>
bool route_tree_node_check_and_mark_congested_for_rip_up(route_tree_t &rt, RouteTreeNode rt_node, const RRGraph &g, const Congestion *congestion)
{
	auto &rt_node_p = get_vertex_props(rt.graph, rt_node);
	const auto &rr_node_p = get_vertex_props(g, rt_node_p.rr_node);

	if (valid(rt_node_p.rt_edge_to_parent)) {
		rt_node_p.pending_rip_up = get_vertex_props(rt.graph, get_source(rt.graph, rt_node_p.rt_edge_to_parent)).pending_rip_up;
	} else {
		rt_node_p.pending_rip_up = false;
	}

	rt_node_p.pending_rip_up |= (get_occ(congestion, rt_node_p.rr_node) > rr_node_p.capacity);
	/*rt_node_p.pending_rip_up = true;*/

	//if (rt_node_p.pending_rip_up) {
		//bg::expand(rt.scheduler_bounding_box, segment(point(rr_node_p.xlow, rr_node_p.ylow), point(rr_node_p.xhigh, rr_node_p.yhigh)));
	//}

	rt_node_p.ripped_up = false;

	return rt_node_p.pending_rip_up;
}

template<typename Congestion>
int route_tree_mark_congested_nodes_to_be_ripped_internal(route_tree_t &rt, const RRGraph &g, const Congestion *congestion, RouteTreeNode rt_node)
{
	const auto &rt_node_p = get_vertex_props(rt.graph, rt_node);
	
	assert(rt_node_p.valid);

	int num_marked = route_tree_node_check_and_mark_congested_for_rip_up(rt, rt_node, g, congestion) ? 1 : 0;

	for (auto &branch : route_tree_get_branches(rt, rt_node)) {
		const auto &child = get_target(rt.graph, branch);

		num_marked += route_tree_mark_congested_nodes_to_be_ripped_internal(rt, g, congestion, child);
	}

	return num_marked;
}

template<typename Congestion>
int route_tree_mark_congested_nodes_to_be_ripped(route_tree_t &rt, const RRGraph &g, const Congestion *congestion)
{
	//rt.scheduler_bounding_box = bg::make_inverse<box>();
	
	int num_marked = 0;

	assert(rt.root_rt_nodes.size() == 1 || rt.root_rt_nodes.size() == 0);

	for (const auto &root_rt_node : rt.root_rt_nodes) {
		num_marked += route_tree_mark_congested_nodes_to_be_ripped_internal(rt, g, congestion, root_rt_node);
	}

	return num_marked;
}

template<typename Congestion>
bool route_tree_node_is_congested(route_tree_t &rt, RouteTreeNode rt_node, const RRGraph &g, const Congestion *congestion)
{
	auto &rt_node_p = get_vertex_props(rt.graph, rt_node);
	const auto &rr_node_p = get_vertex_props(g, rt_node_p.rr_node);

	return (get_occ(congestion, rt_node_p.rr_node) > rr_node_p.capacity);
}

template<typename Congestion>
bool route_tree_is_congested_internal(route_tree_t &rt, RouteTreeNode rt_node, const RRGraph &g, const Congestion *congestion)
{
	bool is_congested = route_tree_node_is_congested(rt, rt_node, g, congestion);

	if (!is_congested) {
		auto bis = route_tree_get_branches(rt, rt_node);

		for (auto bi = std::begin(bis); !is_congested && bi != std::end(bis); ++bi) {
			const auto &child = get_target(rt.graph, *bi);

			is_congested = route_tree_is_congested_internal(rt, child, g, congestion);
		}
	}

	return is_congested;
}

template<typename Congestion>
bool route_tree_is_congested(route_tree_t &rt, const RRGraph &g, const Congestion *congestion)
{
	int num_marked = 0;

	assert(rt.root_rt_nodes.size() == 1 || rt.root_rt_nodes.size() == 0);

	bool is_congested = false;
	for (int i = 0; i < rt.root_rt_nodes.size() && !is_congested; ++i) {
		is_congested = route_tree_is_congested_internal(rt, rt.root_rt_nodes[i], g, congestion);
	}

	return is_congested;
}

//std::shared_ptr<vector<path_node_t>> route_tree_get_path(const route_tree_t &rt, RRNode to_node);

//RouteTreeNode route_tree_get_nearest_node(route_tree_t &rt, const point &p, const RRGraph &g, int *num_iters);

#endif
