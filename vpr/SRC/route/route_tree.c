#include "route.h"
#include "router.h"
#include "route_tree.h"
#include "congestion.h"
#include "log.h"
#include "utility.h"

void route_tree_init(route_tree_t &rt)
{
	rt.root_rt_node_id = RouteTree::null_vertex();
	rt.num_nodes = 0;
}

bool route_tree_empty(const route_tree_t &rt)
{
	int num_rt_nodes = 0;
	for (const auto &n : route_tree_get_nodes(rt)) {
		++num_rt_nodes;
	}
	/*if (rt.root_rt_node_id == -1) {*/
		/*assert(num_rt_nodes == 0);*/
	/*} else {*/
		/*assert(num_rt_nodes > 0);*/
	/*}*/
	/*return rt.root_rt_node_id == -1;*/
	return num_rt_nodes == 0;
}

boost::iterator_range<boost::filter_iterator<valid_rt_node, route_tree_t::vertex_iterator>>
route_tree_get_nodes(const route_tree_t &rt)
{
	valid_rt_node validator(rt.graph);
	const auto &nodes = get_vertices(rt.graph);

	return boost::make_iterator_range(
			boost::make_filter_iterator<valid_rt_node>(validator, begin(nodes), end(nodes)),
			boost::make_filter_iterator<valid_rt_node>(validator, end(nodes), end(nodes))
			);
	/*return get_vertices(rt.graph) | boost::adaptors::filtered([&rt] (unsigned long rt_node) -> bool { return get_vertex_props(rt.graph, rt_node).valid; }); */
	/*return route_tree_get_nodes_impl<const route_tree_t, route_tree_t::vertex_iterator>(rt);*/
}


boost::iterator_range<route_tree_t::branch_iterator>
route_tree_get_branches(const route_tree_t &rt, int rt_node)
{
	return get_out_edges(rt.graph, rt_node);
}

RouteTreeNode route_tree_add_rr_node(route_tree_t &rt, RRNode rr_node, const RRGraph &g)
{
	auto iter = rt.rr_node_to_rt_node.find(rr_node);

	RouteTreeNode rt_node;
	rt_node_property_t *rt_node_p;

	if (iter == rt.rr_node_to_rt_node.end()) {
		add_vertex(rt.graph);
		rt_node = num_vertices(rt.graph)-1;
		rt.rr_node_to_rt_node[rr_node] = rt_node;
		rt_node_p = &get_vertex_props(rt.graph, rt_node);
		rt_node_p->rr_node = rr_node;
	} else {
		rt_node_p = &get_vertex_props(rt.graph, iter->second);
		assert(rt_node_p->rr_node == rr_node);
		if (rt_node_p->valid) {
			rt_node = RouteTree::null_vertex();
		}
	}

	if (rt_node != RouteTree::null_vertex()) {
		assert(!rt_node_p->valid);
		rt_node_p->valid = true;
		rt_node_p->pending_rip_up = false;
		rt_node_p->rt_edge_to_parent = RouteTree::null_edge();
		rt_node_p->num_iterations_fixed = 0;

		auto &rr_node_p = get_vertex_props(g, rr_node);

		segment seg(point(rr_node_p.xlow, rr_node_p.ylow), point(rr_node_p.xhigh, rr_node_p.yhigh));

		if (rr_node_p.type != IPIN && rr_node_p.type != SINK) {
			assert(rt.point_tree.count(make_pair(seg, rt_node)) == 0);

			rt.point_tree.insert(make_pair(seg, rt_node));
		}

		++rt.num_nodes;

		/*assert(rt.num_nodes == rt.point_tree.size());*/
	}

	return rt_node;
}

const RouteTreeEdge &route_tree_add_edge_between_rr_node(route_tree_t &rt, RRNode rr_node_a, RRNode rr_node_b)
{
	/*int rt_node_a = route_tree_add_rr_node(rt, rr_node_a);*/
	/*int rt_node_b = route_tree_add_rr_node(rt, rr_node_b);*/
	/*assert(rt_node_a < num_vertices(rt.graph));*/
	/*assert(rt_node_b < num_vertices(rt.graph));*/
	RouteTreeNode rt_node_a = route_tree_get_rt_node(rt, rr_node_a);
	RouteTreeNode rt_node_b = route_tree_get_rt_node(rt, rr_node_b);

	assert(rt_node_a != -1);
	assert(rt_node_b != -1);

	if (has_edge(rt.graph, rt_node_a, rt_node_b)) {
		char buffer_a[256];
		char buffer_b[256];
		sprintf_rr_node(get_vertex_props(rt.graph, rt_node_a).rr_node, buffer_a);
		sprintf_rr_node(get_vertex_props(rt.graph, rt_node_b).rr_node, buffer_b);
		zlog_error(delta_log, "Existing edge between %s and %s\n", buffer_a, buffer_b);
		assert(false);
	}

	const auto &edge = add_edge(rt.graph, rt_node_a, rt_node_b);
	/* a route tree node can only have one driver */
	auto &rt_node_b_p = get_vertex_props(rt.graph, rt_node_b);
	assert(!valid(rt_node_b_p.rt_edge_to_parent));
	rt_node_b_p.rt_edge_to_parent = edge;

	return edge;
}


RouteTreeNode route_tree_get_rt_node(const route_tree_t &rt, RRNode rr_node)
{
	RouteTreeNode res;

	auto iter = rt.rr_node_to_rt_node.find(rr_node);

	if (iter == rt.rr_node_to_rt_node.end()) {
		res = -1;
	} else {
		auto &v = get_vertex_props(rt.graph, iter->second);
		if (v.valid) {
			res = iter->second;
		} else {
			res = -1;
		}
	}

	return res;
}

void route_tree_set_node_properties(rt_node_property_t &rt_node_p, bool reexpand, const RREdge &prev_edge, float upstream_R, float delay)
{
	rt_node_p.reexpand = reexpand;
	rt_node_p.rr_edge_to_parent = prev_edge;
	rt_node_p.upstream_R = upstream_R;
	rt_node_p.delay = delay;
}

void route_tree_add_root(route_tree_t &rt, RRNode rr_node)
{
	auto iter = rt.rr_node_to_rt_node.find(rr_node);
	assert(iter != rt.rr_node_to_rt_node.end());
	assert(find(begin(rt.root_rt_nodes), end(rt.root_rt_nodes), iter->second) == end(rt.root_rt_nodes));
	rt.root_rt_nodes.push_back(iter->second);
}

void route_tree_set_root(route_tree_t &rt, RRNode rr_node)
{
	auto iter = rt.rr_node_to_rt_node.find(rr_node);
	assert(iter != rt.rr_node_to_rt_node.end());
	assert(rt.root_rt_node_id == RouteTree::null_vertex());
	rt.root_rt_node_id = iter->second;
}

void route_tree_add_to_heap_internal(const route_tree_t &rt, RouteTreeNode rt_node, const RRGraph &g, const RRNode &target, float criticality_fac, float astar_fac, const bounding_box_t &current_bounding_box, std::priority_queue<route_state_t> &heap, perf_t *perf)
{
	const auto &rt_node_p = get_vertex_props(rt.graph, rt_node);

	char buffer[256];
	sprintf_rr_node(rt_node_p.rr_node, buffer);

	if (!rt_node_p.valid) {
		zlog_level(delta_log, ROUTER_V2, "Invalid route tree node %s\n", buffer);
		assert(false);
	}

	const auto &rr_node_p = get_vertex_props(g, rt_node_p.rr_node);

	if (rt_node_p.reexpand) {
		if (rr_node_p.xhigh < current_bounding_box.xmin
				|| rr_node_p.xlow > current_bounding_box.xmax
				|| rr_node_p.yhigh < current_bounding_box.ymin
				|| rr_node_p.ylow > current_bounding_box.ymax) {
			zlog_level(delta_log, ROUTER_V2, "Existing route tree node %s outside of current bounding box\n", buffer);
		} else {
			route_state_t item;

			item.rr_node = rt_node_p.rr_node;
			item.known_cost = criticality_fac * rt_node_p.delay;
			float expected_cost = get_timing_driven_expected_cost(rr_node_p, get_vertex_props(g, target), criticality_fac, rt_node_p.upstream_R);
			item.cost = item.known_cost + astar_fac * expected_cost;
			item.prev_edge = RRGraph::null_edge();
			item.upstream_R = rt_node_p.upstream_R;
			item.delay = rt_node_p.delay;

			zlog_level(delta_log, ROUTER_V2, "Adding route tree node %s prev: %d cost: %g known_cost: %g delay: %g crit_fac: %g expected_cost: %g astar_fac: %g upstream_R: %g to heap\n", buffer, valid(rt_node_p.rt_edge_to_parent) ? get_source(rt.graph, rt_node_p.rt_edge_to_parent) : -1, item.cost, item.known_cost, item.delay, criticality_fac, expected_cost, astar_fac, item.upstream_R);

			if (perf) {
				++perf->num_heap_pushes;
			}
			heap.push(item);
		}
	} else {
		zlog_level(delta_log, ROUTER_V2, "Not expanding route tree node %s\n", buffer);
	}

	/*for_all_out_edges(rt.graph, *rt_node, [&rt, &heap, &criticality_fac, &astar_fac, &g, &target, &current_bounding_box, &perf] (const RouteTreeEdge &e) -> void {*/
	for (const auto &branch : route_tree_get_branches(rt, rt_node)) {
		const auto &neighbor = get_target(rt.graph, branch);

		/*zlog_level(delta_log, ROUTER_V2, "Neighbor of route tree node %s: %d\n", id(neighbor));*/

		route_tree_add_to_heap_internal(rt, neighbor, g, target, criticality_fac, astar_fac, current_bounding_box, heap, perf);
	}
}

void route_tree_multi_root_add_to_heap(const route_tree_t &rt, const RRGraph &g, RRNode target, float criticality_fac, float astar_fac, const bounding_box_t &current_bounding_box, std::priority_queue<route_state_t> &heap, perf_t *perf)
{
	for (const auto &root : rt.root_rt_nodes) {
		route_tree_add_to_heap_internal(rt, root, g, target, criticality_fac, astar_fac, current_bounding_box, heap, perf);
	}
}

void route_tree_add_to_heap(const route_tree_t &rt, const RRGraph &g, RRNode target, float criticality_fac, float astar_fac, const bounding_box_t &current_bounding_box, std::priority_queue<route_state_t> &heap, perf_t *perf)
{
	route_tree_add_to_heap_internal(rt, rt.root_rt_node_id, g, target, criticality_fac, astar_fac, current_bounding_box, heap, perf);
}

void route_tree_add_path(route_tree_t &rt, const vector<path_node_t> &path, const RRGraph &g, const route_state_t *state)
{
	assert(path.size() > 0);

	int current_rr_node_id = path[0].rr_node_id;

	RouteTreeNode current_rt_node = route_tree_add_rr_node(rt, current_rr_node_id, g);
	const auto &current_rr_node_p = get_vertex_props(g, current_rr_node_id);
	route_tree_set_node_properties(get_vertex_props(rt.graph, current_rt_node), current_rr_node_p.type != IPIN && current_rr_node_p.type != SINK, path[0].prev_edge, state[current_rr_node_id].upstream_R, state[current_rr_node_id].delay);

	for (int i = 1; i < path.size(); ++i) {
		int previous_rr_node_id = path[i].rr_node_id;

		if (i < path.size()-1) {
			RouteTreeNode previous_rt_node = route_tree_add_rr_node(rt, previous_rr_node_id, g);
			const auto &previous_rr_node_p = get_vertex_props(g, previous_rr_node_id);
			route_tree_set_node_properties(get_vertex_props(rt.graph, previous_rt_node), previous_rr_node_p.type != IPIN && previous_rr_node_p.type != SINK, path[i].prev_edge, state[previous_rr_node_id].upstream_R, state[previous_rr_node_id].delay);
		}

		route_tree_add_edge_between_rr_node(rt, previous_rr_node_id, current_rr_node_id);
		current_rr_node_id = previous_rr_node_id;
	}
}

bool route_tree_node_check_and_mark_congested_for_rip_up(route_tree_t &rt, RouteTreeNode rt_node, const RRGraph &g, const congestion_t &congestion)
{
	auto &rt_node_p = get_vertex_props(rt.graph, rt_node);
	const auto &rr_node_p = get_vertex_props(g, rt_node_p.rr_node);

	if (valid(rt_node_p.rt_edge_to_parent)) {
		rt_node_p.pending_rip_up = get_vertex_props(rt.graph, get_source(rt.graph, rt_node_p.rt_edge_to_parent)).pending_rip_up;
	} else {
		rt_node_p.pending_rip_up = false;
	}

	rt_node_p.pending_rip_up |= (congestion.occ > rr_node_p.capacity);
	/*rt_node_p.pending_rip_up = true;*/

	if (rt_node_p.pending_rip_up) {
		bg::expand(rt.scheduler_bounding_box, segment(point(rr_node_p.xlow, rr_node_p.ylow), point(rr_node_p.xhigh, rr_node_p.yhigh)));
	}

	rt_node_p.ripped_up = false;

	return rt_node_p.pending_rip_up;
}

bool route_tree_mark_congested_nodes_to_be_ripped_internal(route_tree_t &rt, const RRGraph &g, const congestion_t *congestion, RouteTreeNode rt_node)
{
	const auto &rt_node_p = get_vertex_props(rt.graph, rt_node);
	
	assert(rt_node_p.valid);

	bool marked = route_tree_node_check_and_mark_congested_for_rip_up(rt, rt_node, g, congestion[rt_node_p.rr_node]);

	for (auto &branch : route_tree_get_branches(rt, rt_node)) {
		auto child = get_target(rt.graph, branch);

		marked |= route_tree_mark_congested_nodes_to_be_ripped_internal(rt, g, congestion, child);
	}

	return marked;
}

bool route_tree_mark_congested_nodes_to_be_ripped(route_tree_t &rt, const RRGraph &g, const congestion_t *congestion)
{
	rt.scheduler_bounding_box = bg::make_inverse<box>();

	if (rt.root_rt_node_id == -1) {
		return false;
	}

	return route_tree_mark_congested_nodes_to_be_ripped_internal(rt, g, congestion, rt.root_rt_node_id);
}

void route_tree_remove_edge(route_tree_t &rt, const RouteTreeEdge &rt_edge)
{
	char s_from[256];
	char s_to[256];
	RRNode from = get_source(rt.graph, rt_edge);
	RRNode to = get_target(rt.graph, rt_edge);
	sprintf_rr_node(from, s_from);
	sprintf_rr_node(to, s_to);
	zlog_level(delta_log, ROUTER_V2, "Removing edge %s -> %s\n", s_from, s_to);

	auto &target = get_vertex_props(rt.graph, get_target(rt.graph, rt_edge));
	target.rt_edge_to_parent = RouteTree::null_edge();
	remove_edge(rt.graph, rt_edge);
}


void route_tree_remove_node(route_tree_t &rt, RRNode rr_node, const RRGraph &g)
{
	RouteTreeNode rt_node = route_tree_get_rt_node(rt, rr_node);
	auto &rt_node_p = get_vertex_props(rt.graph, rt_node);

	assert(rt_node_p.valid);
	rt_node_p.valid = false;

	if (rt_node == rt.root_rt_node_id) {
		rt.root_rt_node_id = RouteTree::null_vertex();
	}

	const auto &rr_node_p = get_vertex_props(g, rr_node);

	if (rr_node_p.type != IPIN && rr_node_p.type != SINK) {
		assert(rt.point_tree.remove(make_pair(
						segment(
							point(rr_node_p.xlow, rr_node_p.ylow),
							point(rr_node_p.xhigh, rr_node_p.yhigh)
							),
						rt_node
						)
					)
			  );
	}

	--rt.num_nodes;
	
	/*assert(rt.num_nodes == rt.point_tree.size());*/
}


void route_tree_rip_up_marked(route_tree_t &rt, const RRGraph &g, congestion_t *congestion, float pres_fac, bool lock, lock_perf_t *lock_perf)
{
	char buffer[256];
	for (auto rt_node : route_tree_get_nodes(rt)) {
		auto &rt_node_p = get_vertex_props(rt.graph, rt_node);

		RRNode rr_node = rt_node_p.rr_node;
		const auto &rr_node_p = get_vertex_props(g, rr_node);
		sprintf_rr_node(rt_node_p.rr_node, buffer);

		if (rt_node_p.pending_rip_up) {
			zlog_level(delta_log, ROUTER_V2, "Ripping up node %s from route tree. Occ: %d Cap: %d\n", buffer, congestion[rt_node_p.rr_node].occ, rr_node_p.capacity);

			if (rr_node_p.type == SOURCE) {
				/*assert(rt_node.saved_num_out_edges > 0);*/
				update_one_cost_internal(rr_node, rr_node_p, congestion[rr_node], -num_out_edges(rt.graph, rt_node), pres_fac, lock, lock_perf); 
			} else {
				update_one_cost_internal(rr_node, rr_node_p, congestion[rr_node], -1, pres_fac, lock, lock_perf); 
			}
			route_tree_remove_node(rt, rr_node, g);
			rt_node_p.pending_rip_up = false;
			rt_node_p.ripped_up = true;

			const auto &edge = rt_node_p.rt_edge_to_parent;
			if (valid(edge)) {
				RRNode parent_rr_node = get_source(rt.graph, edge);
				const auto &parent_rt_node_p = get_vertex_props(rt.graph, parent_rr_node);
				const auto &parent_rr_node_p = get_vertex_props(g, parent_rt_node_p.rr_node);
				/* since reconnection back to SOURCE always causes cost to be updated, we need to
				 * update the cost when ripping up also.
				 * if parent is pending rip up, cost will be updated that time. so dont handle it here */
				if (parent_rr_node_p.type == SOURCE && !parent_rt_node_p.pending_rip_up && !parent_rt_node_p.ripped_up) {
					update_one_cost_internal(parent_rr_node, parent_rr_node_p, congestion[parent_rr_node], -1, pres_fac, lock, lock_perf); 
				}

				route_tree_remove_edge(rt, edge);

				assert(!valid(edge));
			} 
		} else {
			zlog_level(delta_log, ROUTER_V2, "NOT ripping up node %s from route tree\n", buffer);
			/* invalid assertion because we might have ripped this up from another virtual net */
			/*assert(rt_node_p.ripped_up == false);*/
		}
	}
}

