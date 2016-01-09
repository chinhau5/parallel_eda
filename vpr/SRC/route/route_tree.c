#include "route.h"
#include "route_tree.h"
#include "log.h"

extern zlog_category_t *delta_log;

void sprintf_rr_node(int rr_node, char *buffer);

void route_tree_init(route_tree_t &rt)
{
	rt.root_rt_node_id = -1;
	rt.num_nodes = 0;
}

void route_tree_clear(route_tree_t &rt)
{
	rt.root_rt_node_id = -1;
	clear_edges(rt.graph);
	for (const auto &item : rt.rr_node_to_rt_node) {
		RouteTreeNode &rt_node = get_vertex(rt.graph, item.second);
		rt_node.properties.valid = false;
	}
	rt.sink_edges.clear();
}

bool route_tree_empty(const route_tree_t &rt)
{
	int num_rt_nodes = 0;
	for (const auto &n : route_tree_get_nodes(rt)) {
		++num_rt_nodes;
	}
	if (rt.root_rt_node_id == -1) {
		assert(num_rt_nodes == 0);
	} else {
		assert(num_rt_nodes > 0);
	}
	return rt.root_rt_node_id == -1;
}

RouteTreeNode *route_tree_add_rr_node(route_tree_t &rt, const RRNode &rr_node)
{
	int rr_node_id = id(rr_node);

	auto iter = rt.rr_node_to_rt_node.find(rr_node_id);

	RouteTreeNode *rt_node = nullptr;

	if (iter == rt.rr_node_to_rt_node.end()) {
		add_vertex(rt.graph);
		int rt_node_id = num_vertices(rt.graph)-1;
		rt.rr_node_to_rt_node[rr_node_id] = rt_node_id;
		rt_node = &get_vertex(rt.graph, rt_node_id);
		rt_node->properties.rr_node = rr_node_id;
	} else {
		rt_node = &get_vertex(rt.graph, iter->second);
		assert(rt_node->properties.rr_node == rr_node_id);
		if (rt_node->properties.valid) {
			rt_node = nullptr;
		}
	}

	if (rt_node) {
		assert(!rt_node->properties.valid);
		rt_node->properties.valid = true;
		rt_node->properties.pending_rip_up = false;
		rt_node->properties.rt_edge_to_parent = -1;
		rt_node->properties.num_iterations_fixed = 0;

		segment seg(point(rr_node.properties.xlow, rr_node.properties.ylow), point(rr_node.properties.xhigh, rr_node.properties.yhigh));

		if (rr_node.properties.type != IPIN && rr_node.properties.type != SINK) {
			assert(rt.point_tree.count(make_pair(seg, id(*rt_node))) == 0);

			rt.point_tree.insert(make_pair(seg, id(*rt_node)));
		}

		++rt.num_nodes;

		/*assert(rt.num_nodes == rt.point_tree.size());*/
	}

	return rt_node;
}

void route_tree_remove_node(route_tree_t &rt, const RRNode &rr_node)
{
	RouteTreeNode *rt_node = route_tree_get_rt_node(rt, id(rr_node));

	assert(rt_node->properties.valid);
	rt_node->properties.valid = false;

	if (id(*rt_node) == rt.root_rt_node_id) {
		rt.root_rt_node_id = -1;
	}

	if (rr_node.properties.type != IPIN && rr_node.properties.type != SINK) {
		assert(rt.point_tree.remove(make_pair(
						segment(
							point(rr_node.properties.xlow, rr_node.properties.ylow),
							point(rr_node.properties.xhigh, rr_node.properties.yhigh)
							),
						id(*rt_node)
						)
					)
			  );
	}

	--rt.num_nodes;
	
	/*assert(rt.num_nodes == rt.point_tree.size());*/
}

int route_tree_num_nodes(const route_tree_t &rt)
{
	return rt.num_nodes;
}

RouteTreeEdge &route_tree_add_edge_between_rr_node(route_tree_t &rt, int rr_node_a, int rr_node_b)
{
	/*int rt_node_a = route_tree_add_rr_node(rt, rr_node_a);*/
	/*int rt_node_b = route_tree_add_rr_node(rt, rr_node_b);*/
	/*assert(rt_node_a < num_vertices(rt.graph));*/
	/*assert(rt_node_b < num_vertices(rt.graph));*/
	RouteTreeNode *rt_node_a = route_tree_get_rt_node(rt, rr_node_a);
	RouteTreeNode *rt_node_b = route_tree_get_rt_node(rt, rr_node_b);

	assert(rt_node_a);
	assert(rt_node_b);

	if (get_edge(rt.graph, *rt_node_a, *rt_node_b)) {
		char buffer_a[256];
		char buffer_b[256];
		sprintf_rr_node(rt_node_a->properties.rr_node, buffer_a);
		sprintf_rr_node(rt_node_b->properties.rr_node, buffer_b);
		zlog_error(delta_log, "Existing edge between %s and %s\n", buffer_a, buffer_b);
		assert(false);
	}

	auto &edge = add_edge(rt.graph, *rt_node_a, *rt_node_b);
	/* a route tree node can only have one driver */
	assert(rt_node_b->properties.rt_edge_to_parent == -1);
	rt_node_b->properties.rt_edge_to_parent = id(edge);

	return edge;
}

void route_tree_remove_edge(route_tree_t &rt, const RouteTreeEdge &edge)
{
	char s_from[256];
	char s_to[256];
	sprintf_rr_node(get_source(rt.graph, edge).properties.rr_node, s_from);
	sprintf_rr_node(get_target(rt.graph, edge).properties.rr_node, s_to);
	zlog_level(delta_log, ROUTER_V2, "Removing edge %s -> %s\n", s_from, s_to);

	for (auto &rt_node: route_tree_get_nodes(rt)) {
		if (rt_node.properties.rt_edge_to_parent > id(edge)) {
			--rt_node.properties.rt_edge_to_parent;
			assert(rt_node.properties.rt_edge_to_parent >= 0);
		}
	}

	auto &target = get_target(rt.graph, edge);
	target.properties.rt_edge_to_parent = -1;
	remove_edge(rt.graph, edge);
}

RouteTreeNode *route_tree_get_rt_node(route_tree_t &rt, int rr_node)
{
	RouteTreeNode *res;

	auto iter = rt.rr_node_to_rt_node.find(rr_node);

	if (iter == rt.rr_node_to_rt_node.end()) {
		res = nullptr;
	} else {
		RouteTreeNode &v = get_vertex(rt.graph, iter->second);
		if (v.properties.valid) {
			res = &v;
		} else {
			res = nullptr;
		}
	}

	return res;
}

const RouteTreeNode *route_tree_get_rt_node(const route_tree_t &rt, int rr_node)
{
	const RouteTreeNode *res;

	auto iter = rt.rr_node_to_rt_node.find(rr_node);
	
	if (iter == rt.rr_node_to_rt_node.end()) {
		res = nullptr;
	} else {	
		const RouteTreeNode &v = get_vertex(rt.graph, iter->second);
		if (v.properties.valid) {
			res = &v;
		} else {
			res = nullptr;
		}
	}

	return res;
}

adapter_t<route_tree_t::vertex_iterator>
route_tree_get_nodes(route_tree_t &rt)
{
	const auto &vertices = get_vertices(rt.graph);
 	auto b = begin(vertices);
	auto e = end(vertices);
	while (b != e && !b->properties.valid) {
		++b;
	}
	return adapter_t<route_tree_t::vertex_iterator>(
			route_tree_t::vertex_iterator(b, e), route_tree_t::vertex_iterator(e, e));
}

adapter_t<route_tree_t::vertex_const_iterator>
route_tree_get_nodes(const route_tree_t &rt)
{
	const auto &vertices = get_vertices(rt.graph);
 	auto b = begin(vertices);
	auto e = end(vertices);
	while (b != e && !b->properties.valid) {
		++b;
	}
	return adapter_t<route_tree_t::vertex_const_iterator>(
			route_tree_t::vertex_const_iterator(b, e), route_tree_t::vertex_const_iterator(e, e));
}

adapter_t<route_tree_t::branch_iterator>
route_tree_get_branches(route_tree_t &rt, RouteTreeNode &rt_node)
{
	const auto &out_edges = get_out_edges(rt.graph, rt_node);
	return adapter_t<route_tree_t::branch_iterator>(
			route_tree_t::branch_iterator(begin(out_edges)), route_tree_t::branch_iterator(end(out_edges)));
}

adapter_t<route_tree_t::branch_const_iterator>
route_tree_get_branches(const route_tree_t &rt, const RouteTreeNode &rt_node)
{
	const auto &out_edges = get_out_edges(rt.graph, rt_node);
	return adapter_t<route_tree_t::branch_const_iterator>(
			route_tree_t::branch_const_iterator(begin(out_edges)), route_tree_t::branch_const_iterator(end(out_edges)));
}

RouteTreeNode *route_tree_get_nearest_node(route_tree_t &rt, const point &p, const RRGraph &g, int *num_iters)
{
	if (!rt.point_tree.size()) {
		return nullptr;
	}

	RouteTreeNode *res = nullptr;

	/*if (num_nodes != rt.point_tree.size()) {*/
		/*vector<route_tree_t::rtree_value> in_rtree;*/
		/*assert(false);*/
	/*}*/

	extern zlog_category_t *delta_log;
	zlog_level(delta_log, ROUTER_V3, "Route tree has %d nodes\n", route_tree_num_nodes(rt));

	if (num_iters) {
		*num_iters = 0;
	}

	auto it = rt.point_tree.qbegin(bgi::nearest(p, rt.point_tree.size()));
	while (it != rt.point_tree.qend() && !res) {
		RouteTreeNode *rt_node = &get_vertex(rt.graph, it->second);
		const auto &rr_node = get_vertex(g, rt_node->properties.rr_node);

		/*auto bounding_box = bg::make_inverse<box>();*/

		/*bg::expand(bounding_box, p);*/

		/*int area = bg::area(bounding_box);*/

		/*char buffer[256];*/
		/*sprintf_rr_node(id(rr_node), buffer);*/
		/*zlog_level(delta_log, ROUTER_V3, "Current nearest to (%d,%d): %s BB: %d-%d %d-%d Area: %d Pending rip up: %d\n", p.get<0>(), p.get<1>(), buffer, bounding_box.min_corner().get<0>(), bounding_box.max_corner().get<0>(), bounding_box.max_corner().get<1>(), bounding_box.max_corner().get<1>(), area, rt_node->properties.pending_rip_up ? 1 : 0);*/

		assert(rr_node.properties.type != IPIN && rr_node.properties.type != SINK);

		if (!rt_node->properties.pending_rip_up) {
			res = rt_node;
		}

		++it;

		if (num_iters) {
			++(*num_iters);
		}
	}

	return res;
}

/*void route_tree_rip_up_marked_internal(route_tree_t &rt, RouteTreeNode &rt_node, RRGraph &g, float pres_fac)*/
/*{*/
/*}*/

void route_tree_rip_up_marked(route_tree_t &rt, RRGraph &g, float pres_fac, bool lock)
{
	char buffer[256];
	for (auto &rt_node : route_tree_get_nodes(rt)) {
		auto &rr_node = get_vertex(g, rt_node.properties.rr_node);
		sprintf_rr_node(id(rr_node), buffer);

		if (rt_node.properties.pending_rip_up) {
			zlog_level(delta_log, ROUTER_V2, "Ripping up node %s from route tree. Occ: %d Cap: %d\n", buffer, rr_node.properties.occ, rr_node.properties.capacity);

			if (rr_node.properties.type == SOURCE) {
				/*assert(rt_node.properties.saved_num_out_edges > 0);*/
				update_one_cost_internal(rr_node, -num_out_edges(rt.graph, rt_node), pres_fac, lock); 
			} else {
				update_one_cost_internal(rr_node, -1, pres_fac, lock); 
			}
			route_tree_remove_node(rt, rr_node);
			rt_node.properties.pending_rip_up = false;
			rt_node.properties.ripped_up = true;

			if (rt_node.properties.rt_edge_to_parent != -1) {
				auto &edge = get_edge(rt.graph, rt_node.properties.rt_edge_to_parent);
				const auto &parent_rt_node = get_source(rt.graph, edge);
				auto &parent_rr_node = get_vertex(g, parent_rt_node.properties.rr_node);
				/* since reconnection back to SOURCE always causes cost to be updated, we need to
				 * update the cost when ripping up also.
				 * if parent is pending rip up, cost will be updated that time. so dont handle it here */
				if (parent_rr_node.properties.type == SOURCE && !parent_rt_node.properties.pending_rip_up && !parent_rt_node.properties.ripped_up) {
					update_one_cost_internal(parent_rr_node, -1, pres_fac, lock); 
				}
				route_tree_remove_edge(rt, edge);

				assert(rt_node.properties.rt_edge_to_parent == -1);
			} 
		} else {
			zlog_level(delta_log, ROUTER_V2, "NOT ripping up node %s from route tree\n", buffer);
			/* invalid assertion because we might have ripped this up from another virtual net */
			/*assert(rt_node.properties.ripped_up == false);*/
		}
	}
}

void route_tree_rip_up_marked_2(route_tree_t &rt, RRGraph &g, float pres_fac)
{
	char buffer[256];
	for (auto &rt_node : route_tree_get_nodes(rt)) {
		auto &rr_node = get_vertex(g, rt_node.properties.rr_node);
		sprintf_rr_node(id(rr_node), buffer);
		if (rt_node.properties.pending_rip_up) {
			zlog_level(delta_log, ROUTER_V2, "Ripping up node %s from route tree\n", buffer);

			if (rr_node.properties.type == SOURCE) {
				if (rr_node.properties.capacity > 1) {
					assert(false);
					for (const auto &e : route_tree_get_branches(rt, rt_node)) {
						/*auto &target = get_target(rt.graph, e);*/
						update_one_cost_internal(rr_node, -1, pres_fac, false); 
					}
				} else {
					update_one_cost_internal(rr_node, -1, pres_fac, false); 
				}
			} else {
				update_one_cost_internal(rr_node, -1, pres_fac, false); 
			}
			route_tree_remove_node(rt, rr_node);
			rt_node.properties.pending_rip_up = false;

			if (rt_node.properties.rt_edge_to_parent != -1) {
				auto &edge = get_edge(rt.graph, rt_node.properties.rt_edge_to_parent);
				route_tree_remove_edge(rt, edge);
			} 
		} else {
			zlog_level(delta_log, ROUTER_V2, "NOT ripping up node %s from route tree\n", buffer);
		}
	}
}

bool route_tree_node_check_and_mark_for_rip_up(route_tree_t &rt, RouteTreeNode &rt_node, const RRGraph &g, int num_iterations_fixed_threshold)
{
	const auto &rr_node = get_vertex(g, rt_node.properties.rr_node);

	/* just a fast hack to handle the case for titan benchmark where SOURCE has capacity > 1 */
	if (rr_node.properties.type == SOURCE) {
		rt_node.properties.saved_num_out_edges = num_out_edges(rt.graph, rt_node);
	}

	++rt_node.properties.num_iterations_fixed;

	if (rt_node.properties.rt_edge_to_parent != -1) {
		rt_node.properties.pending_rip_up = get_source(rt.graph, get_edge(rt.graph, rt_node.properties.rt_edge_to_parent)).properties.pending_rip_up;
	} else {
		rt_node.properties.pending_rip_up = false;
	}

	rt_node.properties.pending_rip_up |= (rr_node.properties.occ > rr_node.properties.capacity) || (rt_node.properties.num_iterations_fixed > num_iterations_fixed_threshold);
	/*rt_node.properties.pending_rip_up = true;*/

	if (rt_node.properties.pending_rip_up) {
		rt_node.properties.num_iterations_fixed = 0;

		bg::expand(rt.scheduler_bounding_box, segment(point(rr_node.properties.xlow, rr_node.properties.ylow), point(rr_node.properties.xhigh, rr_node.properties.yhigh)));
	}

	rt_node.properties.ripped_up = false;

	return rt_node.properties.pending_rip_up;
}

bool route_tree_mark_nodes_to_be_ripped_internal(route_tree_t &rt, const RRGraph &g, RouteTreeNode &rt_node, int num_iterations_fixed_threshold)
{
	assert(rt_node.properties.valid);

	bool marked = route_tree_node_check_and_mark_for_rip_up(rt, rt_node, g, num_iterations_fixed_threshold);

	for (auto &branch : route_tree_get_branches(rt, rt_node)) {
		auto &child = get_target(rt.graph, branch);

		marked |= route_tree_mark_nodes_to_be_ripped_internal(rt, g, child, num_iterations_fixed_threshold);
	}

	return marked;
}

bool route_tree_mark_nodes_to_be_ripped(route_tree_t &rt, const RRGraph &g, int num_iterations_fixed_threshold)
{
	rt.scheduler_bounding_box = bg::make_inverse<box>();

	if (rt.root_rt_node_id == -1) {
		return false;
	}

	auto &root_rt_node = get_vertex(rt.graph, rt.root_rt_node_id);
	return route_tree_mark_nodes_to_be_ripped_internal(rt, g, root_rt_node, num_iterations_fixed_threshold);
}

void route_tree_rip_up_segment_2(route_tree_t &rt, int sink_rr_node, RRGraph &g, float pres_fac)
{
	RouteTreeNode *current = route_tree_get_rt_node(rt, sink_rr_node);
	if (!current) {
		return;
	}

	assert(num_out_edges(rt.graph, *current) == 0);

	bool stop;
	do {
		char buffer[256];
		sprintf_rr_node(current->properties.rr_node, buffer);

		zlog_level(delta_log, ROUTER_V2, "Ripping up node %s from route tree\n", buffer);

		update_one_cost_internal(get_vertex(g, current->properties.rr_node), -1, pres_fac, false);
		assert(current->properties.valid);
		current->properties.valid = false;
		current->properties.owner = -1;

		if (get_vertex(g, current->properties.rr_node).properties.type == SOURCE) {
			zlog_level(delta_log, ROUTER_V2, "Ripping up the source node. The route tree is now empty\n");
			rt.root_rt_node_id = -1;
			bool empty = true;
			for_all_vertices(rt.graph, [&empty] (const RouteTreeNode &rt_node) -> void {
					if (rt_node.properties.valid) {
					empty = false;
					}
					});
			assert(empty);
		}

		RouteTreeNode *parent;
		if (current->properties.rt_edge_to_parent != -1) {
			auto &edge = get_edge(rt.graph, current->properties.rt_edge_to_parent);
			parent = &get_source(rt.graph, edge);
			for_all_vertices(rt.graph, [&current] (RouteTreeNode &node) -> void {
					if (node.properties.rt_edge_to_parent > current->properties.rt_edge_to_parent) {
					--node.properties.rt_edge_to_parent;
					}
					});
			remove_edge(rt.graph, edge);
		} else {
			parent = nullptr;
		}

		current = parent;
	} while (current && num_out_edges(rt.graph, *current) == 0);
}

void route_tree_rip_up_segment(route_tree_t &rt, int sink_rr_node, RRGraph &g, float pres_fac)
{
	auto iter = rt.sink_edges.find(sink_rr_node);
	char s_src[256];
	char s_dst[256];

	if (iter == rt.sink_edges.end()) {
		sprintf_rr_node(sink_rr_node, s_src);
		zlog_level(delta_log, ROUTER_V2, "No existing branch for sink %s (sink_edges size %lu)\n", s_src, rt.sink_edges.size());
		return;
	}

	RouteTreeNode *src = nullptr;
	RouteTreeNode *dst = nullptr;

	bool should_remove_src = true;

	for (auto ei = iter->second.begin(); ei != iter->second.end() && should_remove_src; ) {
		auto &edge = get_edge(rt.graph, *ei);
		src = &get_source(rt.graph, edge);
		dst = &get_target(rt.graph, edge);
		sprintf_rr_node(src->properties.rr_node, s_src);
		sprintf_rr_node(dst->properties.rr_node, s_dst);
		zlog_level(delta_log, ROUTER_V3, "Current edge: %s -> %s, %s out edges: %d\n", s_src, s_dst, s_src, num_out_edges(rt.graph, *src));

		assert(num_out_edges(rt.graph, *dst) == 0);

		zlog_level(delta_log, ROUTER_V2, "Ripping up node %s from route tree\n", s_dst);
		update_one_cost_internal(get_vertex(g, dst->properties.rr_node), -1, pres_fac, false);
		assert(dst->properties.valid);
		dst->properties.valid = false;
		dst->properties.owner = -1;

		for (auto &other : rt.sink_edges) {
			for (auto &other_e : other.second) {
				if (other_e > *ei) {
					--other_e;
					assert(other_e >= 0);
				}
			}
		}
		remove_edge(rt.graph, edge);
		ei = iter->second.erase(ei);

		should_remove_src = num_out_edges(rt.graph, *src) == 0;

		if (!should_remove_src) {
			zlog_level(delta_log, ROUTER_V2, "Stopping rip up because %s has %d out edges (> 0)\n", s_src, num_out_edges(rt.graph, *src));
		}
	}

	if (should_remove_src) {
		zlog_level(delta_log, ROUTER_V2, "Ripping up node %s from route tree\n", s_src);
		update_one_cost_internal(get_vertex(g, src->properties.rr_node), -1, pres_fac, false);
		assert(src->properties.valid);
		src->properties.valid = false;
		src->properties.owner = -1;
		if (get_vertex(g, src->properties.rr_node).properties.type == SOURCE) {
			zlog_level(delta_log, ROUTER_V2, "Ripping up the source node. The route tree is now empty\n");
			rt.root_rt_node_id = -1;
			bool empty = true;
			for_all_vertices(rt.graph, [&empty] (const RouteTreeNode &rt_node) -> void {
					if (rt_node.properties.valid) {
					empty = false;
					}
					});
			assert(empty);
		}
	} else {
		zlog_level(delta_log, ROUTER_V2, "Not ripping up %s\n", s_src);
	}

	/*if (src) {*/
		/*sprintf_rr_node(src->properties.rr_node, buffer);*/
		/*[>if (num_out_edges(rt.graph, *src) <= 1) {<]*/
			/*[>zlog_level(delta_log, ROUTER_V2, "Ripping up node %s from route tree because num_out_edges == %d (<= 1)\n", buffer, num_out_edges(rt.graph, *src));<]*/
			/*zlog_level(delta_log, ROUTER_V2, "Ripping up node %s from route tree\n", buffer);*/
			/*update_one_cost_internal(get_vertex(g, src->properties.rr_node), -1, pres_fac);*/
			/*src->properties.valid = false;*/
		/*[>} else {<]*/
			/*[>zlog_level(delta_log, ROUTER_V2, "Not ripping up node %s from route tree\n", buffer);<]*/
		/*[>}<]*/
	/*}*/

	rt.sink_edges.erase(iter);
}

void route_tree_add_path(route_tree_t &rt, const RRGraph &g, const route_state_t *state, const vector<int> &rr_nodes, int vpr_net_id)
{
	char p[256];
	char c[256];

	assert(get_vertex(g, rr_nodes[0]).properties.type == SINK);
	/*assert(rt.sink_edges.find(rr_nodes[0]) == rt.sink_edges.end());*/
	/*vector<int> edges;*/
	/*rt.sinks.push_back(rr_nodes[0]); vector<edge*/

	for (int i = 0; i < rr_nodes.size()-1; ++i) {
		/* current */
		int current_rr_node_id = rr_nodes[i];
		const RRNode &current_rr_node = get_vertex(g, current_rr_node_id);
		RouteTreeNode *current_rt_node = route_tree_checked_add_rr_node(rt, current_rr_node, state);
		route_tree_set_node_properties(*current_rt_node, current_rr_node.properties.type != IPIN && current_rr_node.properties.type != SINK, state[current_rr_node_id].prev_edge, state[current_rr_node_id].upstream_R, state[current_rr_node_id].delay);

		/* parent */
		const RRNode &parent_rr_node = get_source(g, *state[current_rr_node_id].prev_edge);
		int parent_rr_node_id = id(parent_rr_node);
		assert(parent_rr_node_id == rr_nodes[i+1]);

		RouteTreeNode *parent_rt_node = route_tree_checked_add_rr_node(rt, parent_rr_node, state);
		auto connection = route_tree_get_connection(parent_rr_node, g, state, i+1 == rr_nodes.size()-1);
		route_tree_set_node_properties(*parent_rt_node, parent_rr_node.properties.type != IPIN && parent_rr_node.properties.type != SINK, connection.first, state[parent_rr_node_id].upstream_R, state[parent_rr_node_id].delay);

		/* setting up edges */
		auto &e = route_tree_add_edge_between_rr_node(rt, parent_rr_node_id, current_rr_node_id);
		e.properties.rr_edge = state[current_rr_node_id].prev_edge;
		/*edges.push_back(id(e));*/

		if (connection.first) {
			assert(connection.second);
			auto &e = route_tree_add_edge_between_rr_node(rt, id(*connection.second), parent_rr_node_id);
			e.properties.rr_edge = connection.first;
			/*edges.push_back(id(e));*/
		}

		/* printing */
		sprintf_rr_node(parent_rr_node_id, p);
		sprintf_rr_node(current_rr_node_id, c);
		zlog_level(delta_log, ROUTER_V2, "Net %d route tree edge: %s -> %s\n", vpr_net_id, p, c);

		if (connection.second) {
			sprintf_rr_node(id(*connection.second), p);
			sprintf_rr_node(parent_rr_node_id, c);
			zlog_level(delta_log, ROUTER_V2, "Net %d connection route tree edge: %s -> %s\n", vpr_net_id, p, c);
		}
	}

	/*rt.sink_edges.insert(make_pair(rr_nodes[0], edges));*/
}

void route_tree_add_path(route_tree_t &rt, const RRGraph &g, const route_state_t *state, int sink_rr_node, int vpr_net_id, vector<int> &rr_nodes_with_pending_cost_update)
{
	assert(get_vertex(g, sink_rr_node).properties.type == SINK);
	/*assert(rt.sink_edges.find(sink_rr_node) == rt.sink_edges.end());*/

	int current_rr_node_id = sink_rr_node;
	bool update_cost;
	bool stop_traceback;
	char p[256];
	char c[256];

	/*vector<int> edges;*/

	const RRNode &current_rr_node = get_vertex(g, current_rr_node_id);
	RouteTreeNode *current_rt_node = route_tree_add_or_get_rr_node(rt, current_rr_node_id, g, state, update_cost, stop_traceback);
	assert(current_rt_node && !stop_traceback);
	if (update_cost) {
		rr_nodes_with_pending_cost_update.push_back(current_rr_node_id);
	}
	current_rt_node->properties.owner = sink_rr_node;
	route_tree_set_node_properties(*current_rt_node, current_rr_node.properties.type != IPIN && current_rr_node.properties.type != SINK, state[current_rr_node_id].prev_edge, state[current_rr_node_id].upstream_R, state[current_rr_node_id].delay);

	while (state[current_rr_node_id].prev_edge && !stop_traceback) {
		/* parent */
		const RRNode &parent_rr_node = get_source(g, *state[current_rr_node_id].prev_edge);
		int parent_rr_node_id = id(parent_rr_node);

		RouteTreeNode *parent_rt_node = route_tree_add_or_get_rr_node(rt, parent_rr_node_id, g, state, update_cost, stop_traceback);
		if (update_cost) {
			rr_nodes_with_pending_cost_update.push_back(parent_rr_node_id);
		}
		parent_rt_node->properties.owner = sink_rr_node;
		route_tree_set_node_properties(*parent_rt_node, parent_rr_node.properties.type != IPIN && parent_rr_node.properties.type != SINK, nullptr, state[parent_rr_node_id].upstream_R, state[parent_rr_node_id].delay);

		/* setting up edges */
		auto &e = route_tree_add_edge_between_rr_node(rt, parent_rr_node_id, current_rr_node_id);
		e.properties.rr_edge = state[current_rr_node_id].prev_edge;

		/*edges.push_back(id(e));*/

		/* printing */
		sprintf_rr_node(parent_rr_node_id, p);
		sprintf_rr_node(current_rr_node_id, c);
		zlog_level(delta_log, ROUTER_V2, "Net %d route tree edge: %s -> %s\n", vpr_net_id, p, c);

		current_rr_node_id = parent_rr_node_id;
	}

	/*rt.sink_edges.insert(make_pair(sink_rr_node, edges));*/

	/*assert(rt.sink_rr_node_to_start_rt_node.find(sink_rr_node) != rt.sink_rr_node_to_start_rt_node.end());*/
	/*RouteTreeNode &new_branch_rt_node = route_tree_get_rt_node(rt, new_branch);*/
	/*rt.sink_rr_node_to_start_rt_node[sink_rr_node] = id(new_branch_rt_node);*/

	/*assert(rt.sink_rr_node_to_path.find(sink_rr_node) == rt.sink_rr_node_to_path.end());*/
	/*const auto &iter = rt.sink_rr_node_to_path.insert(make_pair(sink_rr_node, path)).first;*/

	/*return iter->second;*/
}

void route_tree_set_root(route_tree_t &rt, int rr_node)
{
	auto iter = rt.rr_node_to_rt_node.find(rr_node);
	assert(iter != rt.rr_node_to_rt_node.end());
	assert(rt.root_rt_node_id == -1);
	rt.root_rt_node_id = iter->second;
}

bool route_tree_is_only_path(const route_tree_t &rt, int sink_rr_node, const RRGraph &g)
{
	/*const RouteTreeNode *rt_node = route_tree_get_rt_node(rt, sink_rr_node);*/
	/*assert(rt_node);*/
	/*assert(rt_node->properties.parent >= 0);*/
	/*const RouteTreeNode *parent_rt_node = &get_vertex(rt.graph, rt_node->properties.parent);*/
	/*bool is_only_path = true;*/
	/*while ([>is_only_path &&<] parent_rt_node) {*/
		/*char buffer[256];*/
		/*sprintf_rr_node(rt_node->properties.rr_node, buffer);*/
		/*[>printf("current: %s\n", buffer);<]*/
		/*zlog_level(delta_log, ROUTER_V3, "current: %s\n", buffer);*/
		/*for_all_out_edges(g, get_vertex(g, parent_rt_node->properties.rr_node), [&g, &rt_node, &is_only_path, &buffer] (const RREdge &e) {*/
					/*auto &sibling = get_target(g, e);*/
					/*if (id(sibling) != rt_node->properties.rr_node && sibling.properties.type != IPIN && sibling.properties.type != SINK) {*/
						/*sprintf_rr_node(id(sibling), buffer);*/
						/*zlog_level(delta_log, ROUTER_V3, "sibling: %s occ/cap: %d/%d\n", buffer, sibling.properties.occ, sibling.properties.capacity);*/
						/*is_only_path = is_only_path & (sibling.properties.occ > sibling.properties.capacity);*/
					/*}*/
				/*});*/
		/*rt_node = parent_rt_node;*/
		/*if (rt_node && rt_node->properties.parent >= 0) {*/
			/*parent_rt_node = &get_vertex(rt.graph, rt_node->properties.parent);*/
		/*} else {*/
			/*parent_rt_node = nullptr;*/
		/*}*/
	/*}*/
	/*return is_only_path;*/
}

pair<const RREdge *, const RRNode *> route_tree_get_connection(const RRNode &current_rr_node, const RRGraph &g, const route_state_t *state, bool end)
{
	pair<const RREdge *, const RRNode *> connection;
	if (end) {
		if (current_rr_node.properties.type == SOURCE) {
			connection.first = nullptr;
			connection.second = nullptr;
		} else {
			connection.first = state[id(current_rr_node)].prev_edge;
			connection.second = &get_source(g, *connection.first);
			int connection_rr_node = id(*connection.second);
			if (state[connection_rr_node].prev_edge != nullptr) {
				char p[256];
				char c[256];
				sprintf_rr_node(connection_rr_node, c);
				sprintf_rr_node(id(get_source(g, *state[connection_rr_node].prev_edge)), p);
				zlog_warn(delta_log, "Warning: Connection %s is not from existing route tree because we found a shorter path to it via %s\n", c, p);
			}
		}
	} else {
		connection.first = nullptr;
		connection.second = nullptr;
	}
	return connection;
}

void route_tree_set_node_properties(RouteTreeNode &rt_node, bool reexpand, const RREdge *prev_edge, float upstream_R, float delay)
{
	rt_node.properties.reexpand = reexpand;
	rt_node.properties.rr_edge_to_parent = prev_edge;
	rt_node.properties.upstream_R = upstream_R;
	rt_node.properties.delay = delay;
}

RouteTreeNode *route_tree_add_or_get_rr_node(route_tree_t &rt, int rr_node_id, const RRGraph &g, const route_state_t *state, bool &update_cost, bool &stop_traceback)
{
	const auto &rr_node = get_vertex(g, rr_node_id);
	RouteTreeNode *rt_node = route_tree_add_rr_node(rt, rr_node);
	update_cost = false;
	stop_traceback = false;

	if (!rt_node) {
		rt_node = route_tree_get_rt_node(rt, rr_node_id);
		/* not added because there is an existing node */
		assert(rt_node);

		char buffer[256];
		sprintf_rr_node(rr_node_id, buffer);

		if (rt_node->properties.rr_edge_to_parent != nullptr) {
			char parent[256];
			sprintf_rr_node(id(get_source(g, *rt_node->properties.rr_edge_to_parent)), parent);
			zlog_error(delta_log, "Error: Existing route tree node %s has non-null rr_edge_to_parent that connects to %s\n", buffer, parent);
			assert(false);
		}

		if (state[rr_node_id].prev_edge != rt_node->properties.rr_edge_to_parent) {
			char s_state[256];
			char s_rt[256];
			sprintf_rr_node(id(get_source(g, *state[rr_node_id].prev_edge)), s_state);
			/*sprintf_rr_node(id(get_source(g, *rt_node->properties.rr_edge_to_parent)), s_rt);*/
			zlog_warn(delta_log, "Warning: Existing route tree node %s does not have a matching route state. (state.prev_edge: %s rt_node.properties.rr_edge_to_parent: ) because we have found a shorter path to that node\n", buffer, s_state);

			/*route_tree_add_edge_between_rr_node(rt, id(get_source(g, *state[rr_node_id].prev_edge)), rr_node_id);*/

			/*write_graph(rt.graph, "existing_route_tree_node.dot",*/
					/*[] (const RouteTreeNode &rt_node) -> string {*/
					/*char buffer[256];*/
					/*sprintf_rr_node(rt_node.properties.rr_node, buffer);*/
					/*return buffer;*/
					/*},*/
					/*[] (const RouteTreeEdge &rt_edge) -> string {*/
					/*return "";*/
					/*},*/
					/*[] (const RouteTreeNode &rt_node) -> bool {*/
					/*return false;*/
					/*});*/
			stop_traceback = true;
		} else {
			zlog_level(delta_log, ROUTER_V2, "Reconnecting to existing route tree node %s ", buffer);
			if (rr_node.properties.type == SOURCE) {
				if (rr_node.properties.capacity > 1) {
					zlog_level(delta_log, ROUTER_V2, "with capacity %d (> 1)", rr_node.properties.capacity);
				} else {
					zlog_level(delta_log, ROUTER_V2, "with capacity %d", rr_node.properties.capacity);
				}
				update_cost = true;
			} else {
				assert(rr_node.properties.type == CHANX || rr_node.properties.type == CHANY || rr_node.properties.type == OPIN);
			}
			zlog_level(delta_log, ROUTER_V2, "\n");
		}
	} else {
		update_cost = true;
	}

	return rt_node;
}

RouteTreeNode *route_tree_checked_add_rr_node(route_tree_t &rt, const RRNode &rr_node, const route_state_t *state)
{
	RouteTreeNode *rt_node = route_tree_add_rr_node(rt, rr_node);
	int rr_node_id = id(rr_node);
	if (!rt_node && (rt_node = route_tree_get_rt_node(rt, rr_node_id))) {
		if (rt_node->properties.upstream_R != state[rr_node_id].upstream_R ||
				rt_node->properties.delay != state[rr_node_id].delay) {
			zlog_error(delta_log, "Error: Adding a new RT node with different properties [old_upstream_R: %g new_upstream_R: %g] [old_delay: %g new_delay: %g]\n", rt_node->properties.upstream_R, state[rr_node_id].upstream_R, rt_node->properties.delay, state[rr_node_id].delay); 
			assert(false);
		} 
	}

	return rt_node;
}
void route_tree_add_to_heap_internal(const route_tree_t &rt, const RouteTreeNode *rt_node, const RRGraph &g, const RRNode &target, float criticality_fac, float astar_fac, const bounding_box_t &current_bounding_box, std::priority_queue<route_state_t> &heap, perf_t *perf)
{
	char buffer[256];
	sprintf_rr_node(rt_node->properties.rr_node, buffer);

	if (!rt_node->properties.valid) {
		zlog_level(delta_log, ROUTER_V2, "Invalid route tree node %s\n", buffer);
		assert(false);
	}

	const auto &prop = get_vertex(g, rt_node->properties.rr_node).properties;

	if (rt_node->properties.reexpand) {
		if (prop.xhigh < current_bounding_box.xmin
				|| prop.xlow > current_bounding_box.xmax
				|| prop.yhigh < current_bounding_box.ymin
				|| prop.ylow > current_bounding_box.ymax) {
			zlog_level(delta_log, ROUTER_V2, "Existing route tree node %s outside of current bounding box\n", buffer);
		} else {
			route_state_t item;

			item.rr_node = rt_node->properties.rr_node;
			item.known_cost = criticality_fac * rt_node->properties.delay;
			float expected_cost = get_timing_driven_expected_cost(get_vertex(g, rt_node->properties.rr_node), target, criticality_fac, rt_node->properties.upstream_R);
			item.cost = item.known_cost + astar_fac * expected_cost;
			item.prev_edge = nullptr;
			item.upstream_R = rt_node->properties.upstream_R;
			item.delay = rt_node->properties.delay;

			zlog_level(delta_log, ROUTER_V2, "Adding route tree node %s prev: %d cost: %g known_cost: %g delay: %g crit_fac: %g expected_cost: %g astar_fac: %g upstream_R: %g to heap\n", buffer, rt_node->properties.rt_edge_to_parent != -1 ? get_source(rt.graph, get_edge(rt.graph, rt_node->properties.rt_edge_to_parent)).properties.rr_node : -1, item.cost, item.known_cost, item.delay, criticality_fac, expected_cost, astar_fac, item.upstream_R);

			if (perf) {
				++perf->num_heap_pushes;
			}
			heap.push(item);
		}
	} else {
		zlog_level(delta_log, ROUTER_V2, "Not expanding route tree node %s\n", buffer);
	}

	/*for_all_out_edges(rt.graph, *rt_node, [&rt, &heap, &criticality_fac, &astar_fac, &g, &target, &current_bounding_box, &perf] (const RouteTreeEdge &e) -> void {*/
	for (const auto &branch : route_tree_get_branches(rt, *rt_node)) {
		const auto &neighbor = get_target(rt.graph, branch);

		route_tree_add_to_heap_internal(rt, &neighbor, g, target, criticality_fac, astar_fac, current_bounding_box, heap, perf);
	}
}

void route_tree_add_to_heap(const route_tree_t &rt, const RRGraph &g, const RRNode &target, float criticality_fac, float astar_fac, const bounding_box_t &current_bounding_box, std::priority_queue<route_state_t> &heap, perf_t *perf)
{
	route_tree_add_to_heap_internal(rt, &get_vertex(rt.graph, rt.root_rt_node_id), g, target, criticality_fac, astar_fac, current_bounding_box, heap, perf);
}

/*const vector<int> &route_tree_add_path(route_tree_t &rt, const RRGraph &g, const route_state_t *state, pair<vector<int>::const_iterator, vector<int>::const_iterator> &path)*/
/*{*/
	/*int new_branch = -1;*/

	/*auto iter = path.first*/

	/*for (; iter != path.second; ++iter) {*/
		/*int current_rr_node = *iter;*/
		/*int prev_rr_node = *(iter + 1);*/

		/*auto &prev_rt_node = route_tree_add_rr_node(rt, prev_rr_node);*/
		/*prev_rt_node.properties.reexpand = get_vertex(g, prev_rr_node).properties.type != IPIN;*/
		/*prev_rt_node.properties.rr_node = prev_rr_node;*/
		/*prev_rt_node.properties.prev_edge = nullptr;*/
		/*prev_rt_node.properties.upstream_R = state[prev_rr_node].upstream_R;*/
		/*prev_rt_node.properties.delay = state[prev_rr_node].delay;*/

		/*auto &current_rt_node = route_tree_add_rr_node(rt, current_rr_node);*/
		/*current_rt_node.properties.reexpand = get_vertex(g, current_rr_node).properties.type != IPIN;*/
		/*current_rt_node.properties.rr_node = current_rr_node;*/
		/*current_rt_node.properties.prev_edge = state[current_rr_node].prev_edge;*/
		/*current_rt_node.properties.upstream_R = state[current_rr_node].upstream_R;*/
		/*current_rt_node.properties.delay = state[current_rr_node].delay;*/

		/*char p[256];*/
		/*char c[256];*/
		/*sprintf_rr_node(prev_rr_node, p);*/
		/*sprintf_rr_node(current_rr_node, c);*/
		/*zlog_level(delta_log, ROUTER_V2, "Route tree add path: %s -> %s\n", p, c);*/

		/*auto &e = route_tree_add_edge_between_rr_node(rt, prev_rr_node, current_rr_node);*/
		/*e.properties.rr_edge = state[current_rr_node].prev_edge;*/

		/*path.push_back(current_rr_node);*/

		/*new_branch = current_rr_node;*/
	/*}*/

	/*[>assert(rt.sink_rr_node_to_start_rt_node.find(sink_rr_node) != rt.sink_rr_node_to_start_rt_node.end());<]*/
	/*RouteTreeNode &new_branch_rt_node = route_tree_get_rt_node(rt, new_branch);*/
	/*[>rt.sink_rr_node_to_start_rt_node[sink_rr_node] = id(new_branch_rt_node);<]*/

	/*assert(rt.sink_rr_node_to_path.find(sink_rr_node) == rt.sink_rr_node_to_path.end());*/
	/*const auto &iter = rt.sink_rr_node_to_path.insert(make_pair(sink_rr_node, path)).first;*/

	/*return iter->second;*/
/*}*/

/*void route_tree_is_only_path_internal(const RRGraph &g, const RRNode &node)*/
/*{*/

/*}*/

/*RouteTree &route_tree_get_sink_root(route_tree_t &rt, int sink_rr_node)*/
/*{*/
	/*assert(rt.sink_rr_node_to_start_rt_node.find(sink_rr_node) != rt.sink_rr_node_to_start_rt_node.end());*/
	/*return get_vertex(rt.graph, rt.sink_rr_node_to_start_rt_node[sink_rr_node]);*/
/*}*/

/*const vector<int> *route_tree_get_path_to_sink(const route_tree_t &rt, int sink_rr_node)*/
/*{*/
	/*const auto &iter = rt.sink_rr_node_to_path.find(sink_rr_node);*/
	/*const vector<int> *res = nullptr;*/
	/*if (iter != rt.sink_rr_node_to_path.end()) {*/
		/*res = &iter->second;*/
	/*}*/
	/*return res;*/
/*}*/

/*void route_tree_remove_path(route_tree_t &rt, int sink_rr_node)*/
/*{*/
	/*auto iter = rt.sink_rr_node_to_path.find(sink_rr_node);*/
	/*assert(iter != rt.sink_rr_node_to_path.end());*/
	/*vector<int> &path = iter->second;*/
	 /*[>important to sort because <]*/
	/*[>sort(path.begin(), path.end(), std::greater<int>());<]*/
	/*for (const auto &edge : path) {*/
		/*remove_edge(rt.graph, get_edge(rt.graph, edge));*/
		/*[> TODO: need to invalidate "path" after removing each edge <]*/
	/*}*/

	/*rt.sink_rr_node_to_path.erase(iter);*/
/*}*/

/*RouteTreeEdge &route_tree_add_edge(route_tree_t &rt, int rt_node_a, int rt_node_b)*/
/*{*/
	/*[>int rt_node_a = route_tree_add_rr_node(rt, rr_node_a);<]*/
	/*[>int rt_node_b = route_tree_add_rr_node(rt, rr_node_b);<]*/
	/*[>assert(rt_node_a < num_vertices(rt.graph));<]*/
	/*[>assert(rt_node_b < num_vertices(rt.graph));<]*/
	/*return add_edge(rt.graph, rt_node_a, rt_node_b);*/
/*}*/

void route_tree_merge(route_tree_t &merged_rt, const route_tree_t &other_rt)
{
	/*for_all_edges(other_rt.graph, [&merged_rt, &other_rt] (const RouteTreeEdge &e) -> void {*/
			/*const RouteTreeNode &other_source_rt_node = get_source(other_rt.graph, e);*/
			/*const RouteTreeNode &other_target_rt_node = get_target(other_rt.graph, e);*/
			/*int source_rr_node = other_source_rt_node.properties.rr_node;*/
			/*int target_rr_node = other_target_rt_node.properties.rr_node;*/

			/*RouteTreeNode &merged_source_rt_node = route_tree_add_rr_node(merged_rt, source_rr_node);*/
			/*merged_source_rt_node.properties = other_source_rt_node.properties;*/
			/*RouteTreeNode &merged_target_rt_node = route_tree_add_rr_node(merged_rt, target_rr_node);*/
			/*merged_target_rt_node.properties = other_target_rt_node.properties;*/

			/*route_tree_add_edge_between_rr_node(merged_rt, source_rr_node, target_rr_node);*/
			/*}*/
			/*);*/
}
