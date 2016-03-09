#include "route.h"
#include "graph.h"
#include "route_tree.h"
#include "util.h"
#include "log.h"

void update_costs(const RRGraph &g, congestion_t *congestion, float pres_fac, float acc_fac)
{
	for (const auto &rr_node : get_vertices(g)) {
		int occ = congestion[rr_node].occ;
		int capacity = get_vertex_props(g, rr_node).capacity;
		if (occ > capacity) {
			congestion[rr_node].acc_cost += (occ - capacity) * acc_fac;
			congestion[rr_node].pres_cost = 1. + (occ + 1 - capacity) * pres_fac;
		} else if (occ == capacity) {
			/* If occ == capacity, we don't need to increase acc_cost, but a change    *
			 * in pres_fac could have made it necessary to recompute the cost anyway.  */
			congestion[rr_node].pres_cost = 1. + pres_fac;
		}
	}
}

void update_one_cost_internal(RRNode rr_node, const rr_node_property_t &rr_node_p, congestion_t &congestion, /*int net_id, */int delta, float pres_fac, bool lock, lock_perf_t *lock_perf)
{
	if (lock) {
		if (lock_perf) {
			++lock_perf->num_lock_tries;
		}

		if (!congestion.lock.try_lock()) {
			if (lock_perf) {
				++lock_perf->num_lock_waits;
			}
			using clock = std::chrono::high_resolution_clock;
			auto wait_start = clock::now();
			congestion.lock.lock();
			if (lock_perf) {
				lock_perf->total_wait_time += clock::now()-wait_start;
			}
		}
		/*rr_node_p.lock->lock();*/
	}
	
	congestion.occ += delta;

	assert(congestion.occ >= 0);

	if (congestion.occ < rr_node_p.capacity) {
		congestion.pres_cost = 1;
	} else {
		congestion.pres_cost = 1 + (congestion.occ + 1 - rr_node_p.capacity) * pres_fac;
	}

	/*if (delta > 0) {*/
		/*assert(std::find(begin(rr_node_p.users), begin(rr_node_p.users), net_id) == end(rr_node_p.users));*/
		/*rr_node_p.users.push_back(net_id);*/
	/*} else {*/
		/*assert(delta < 0);*/
		/*auto iter = std::find(begin(rr_node_p.users), begin(rr_node_p.users), net_id);*/
		/*assert(iter != end(rr_node_p.users));*/
		/*rr_node_p.users.erase(iter);*/
	/*}*/

	if (lock) {
		congestion.lock.unlock();
	}
		
	char buffer[256];
	sprintf_rr_node(rr_node, buffer);
	zlog_level(delta_log, ROUTER_V2, "Update cost of %s delta: %d new_occ: %d pres_fac: %g\n", buffer, delta, congestion.occ, pres_fac);
}

void update_one_cost(const RRGraph &g, congestion_t *congestion, const vector<RRNode>::const_iterator &rr_nodes_begin, const vector<RRNode>::const_iterator &rr_nodes_end, /*int net_id,*/ int delta, float pres_fac, bool lock, lock_perf_t *lock_perf)
{
	/*const RouteTreeNode *last = nullptr;*/
	for (auto iter = rr_nodes_begin; iter != rr_nodes_end; ++iter) {
		/* we don't update get_source because that's the link to existing route tree and the cost is handled by update_one_cost_internal */
		/*const RouteTreeNode &rt_node = get_target(rt.graph, get_edge(rt.graph, rt_edge_id));*/
		/*last = &get_target(rt.graph, get_edge(rt.graph, rt_edge_id));*/
		/*RRNode &rr_node = get_vertex(g, rt_node.rr_node);*/
		const auto &rr_node_p = get_vertex_props(g, *iter);
		update_one_cost_internal(*iter, rr_node_p, congestion[*iter], /*net_id,*/ delta, pres_fac, lock, lock_perf);
	}
	/*RRNode &rr_node = get_vertex(g, last->rr_node);*/
	/*update_one_cost_internal(rr_node, delta, pres_fac);*/
}

void update_one_cost(const RRGraph &g, congestion_t *congestion, route_tree_t &rt, RouteTreeNode rt_node, int delta, float pres_fac, bool lock)
{
	const auto &rt_node_p = get_vertex_props(rt.graph, rt_node);

	assert(rt_node_p.valid);

	int rr_node = rt_node_p.rr_node;

	const auto &rr_node_p = get_vertex_props(g, rr_node);

	update_one_cost_internal(rr_node, rr_node_p, congestion[rr_node], /*-1,*/ delta, pres_fac, lock, nullptr);

	for (const auto &e : get_out_edges(rt.graph, rt_node)) {
		update_one_cost(g, congestion, rt, get_target(rt.graph, e), delta, pres_fac, lock);
	}
}

