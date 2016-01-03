#ifndef ROUTE_TREE_H
#define ROUTE_TREE_H

#include "log.h"
#include "utility.h"
#include "geometry.h"

typedef struct rt_node_property_t {
	int rr_node;
	bool valid;
	bool pending_rip_up;
	bool ripped_up;
	int rt_edge_to_parent;
	int num_iterations_fixed;
	int saved_num_out_edges; /* only valid for SOURCE */

	int owner;

	bool reexpand;
	const RREdge *rr_edge_to_parent;
	float upstream_R;	
	float delay;

	float downstream_C;
	/*float upstream_R_from_route_state;*/
} rt_node_property_t;

typedef struct rt_edge_property_t {
	const RREdge *rr_edge;
} rt_edge_property_t;

typedef graph_t<rt_node_property_t, rt_edge_property_t> RouteTree;
typedef vertex_t<struct rt_node_property_t, struct rt_edge_property_t> RouteTreeNode;
typedef edge_t<rt_edge_property_t> RouteTreeEdge;

typedef struct route_tree_t {
	typedef std::pair<segment, int> rtree_value;

	RouteTree graph;
	int root_rt_node_id;
	int num_nodes;
	std::map<int, int> rr_node_to_rt_node;
	/*map<int, vector<int>> sink_rr_node_to_path;*/
	bgi::rtree<rtree_value, bgi::rstar<64>> point_tree;
	box scheduler_bounding_box;
	map<int, vector<int>> sink_edges;

	template<typename Value, typename Base>
	struct iterator : public std::iterator<
					  typename std::forward_iterator_tag,
					  Value,
					  ptrdiff_t,
					  Value *,
					  Value & 
					  > { 
		Base c;
		Base e;

		iterator(const Base &c, const Base &e) 
			: c(c), e(e)
	   	{
		}

		iterator &operator++()
		{
			do {
				if (c != e) {
					++c;
				}
			} while (c != e && !c->properties.valid);
		
			return *this;
		}

		typename iterator::reference operator*() const
		{
			return *c;
		}

		bool operator==(const iterator &other) const
		{
			return other.c == c;
		}

		bool operator!=(const iterator &other) const
		{
			return other.c != c;
		}
	};

	typedef iterator<RouteTreeNode, typename RouteTree::vertex_iterator> vertex_iterator;
	typedef iterator<const RouteTreeNode, typename RouteTree::vertex_const_iterator> vertex_const_iterator;
	typedef typename RouteTree::out_edges_iterator branch_iterator;
	typedef typename RouteTree::out_edges_const_iterator branch_const_iterator;
} route_tree_t;

void route_tree_init(route_tree_t &rt);

bool route_tree_empty(const route_tree_t &rt);

void route_tree_is_only_path_internal(const RRGraph &g, const RRNode &node);

RouteTreeNode *route_tree_add_rr_node(route_tree_t &rt, const RRNode &rr_node);

int route_tree_num_nodes(const route_tree_t &rt);

RouteTreeEdge &route_tree_add_edge_between_rr_node(route_tree_t &rt, int rr_node_a, int rr_node_b);

RouteTreeNode *route_tree_get_rt_node(route_tree_t &rt, int rr_node);

const RouteTreeNode *route_tree_get_rt_node(const route_tree_t &rt, int rr_node);

adapter_t<route_tree_t::vertex_iterator>
route_tree_get_nodes(route_tree_t &rt);

adapter_t<route_tree_t::vertex_const_iterator>
route_tree_get_nodes(const route_tree_t &rt);

adapter_t<route_tree_t::branch_iterator>
route_tree_get_branches(route_tree_t &rt, RouteTreeNode &rt_node);

adapter_t<route_tree_t::branch_const_iterator>
route_tree_get_branches(const route_tree_t &rt, const RouteTreeNode &rt_node);

bool route_tree_mark_nodes_to_be_ripped(route_tree_t &rt, const RRGraph &g, int num_iterations_fixed_threshold);

void route_tree_rip_up_marked(route_tree_t &rt, RRGraph &g, float pres_fac);

void route_tree_rip_up_segment(route_tree_t &rt, int sink_rr_node, RRGraph &g, float pres_fac);

void route_tree_rip_up_segment_2(route_tree_t &rt, int sink_rr_node, RRGraph &g, float pres_fac);

void route_tree_add_path(route_tree_t &rt, const RRGraph &g, const route_state_t *state, const vector<int> &rr_nodes, int vpr_net_id);

void route_tree_add_path(route_tree_t &rt, const RRGraph &g, const route_state_t *state, int sink_rr_node, int vpr_net_id, vector<int> &added_rr_nodes);

void route_tree_clear(route_tree_t &rt);

void route_tree_set_root(route_tree_t &rt, int rr_node);

bool route_tree_is_only_path(const route_tree_t &rt, int sink_rr_node, const RRGraph &g);

pair<const RREdge *, const RRNode *> route_tree_get_connection(const RRNode &current_rr_node, const RRGraph &g, const route_state_t *state, bool end);

void route_tree_set_node_properties(RouteTreeNode &rt_node, bool reexpand, const RREdge *prev_edge, float upstream_R, float delay);

RouteTreeNode *route_tree_add_or_get_rr_node(route_tree_t &rt, int rr_node_id, const RRGraph &g, const route_state_t *state, bool &update_cost, bool &stop_traceback);

RouteTreeNode *route_tree_checked_add_rr_node(route_tree_t &rt, const RRNode &rr_node, const route_state_t *state);

void route_tree_add_to_heap_internal(const route_tree_t &rt, const RouteTreeNode *rt_node, const RRGraph &g, const RRNode &target, float criticality_fac, float astar_fac, const bounding_box_t &current_bounding_box, std::priority_queue<route_state_t> &heap, perf_t *perf);

void route_tree_add_to_heap(const route_tree_t &rt, const RRGraph &g, const RRNode &target, float criticality_fac, float astar_fac, const bounding_box_t &current_bounding_box, std::priority_queue<route_state_t> &heap, perf_t *perf);

RouteTreeNode *route_tree_get_nearest_node(route_tree_t &rt, const point &p, const RRGraph &g, int *num_iters);

#endif
