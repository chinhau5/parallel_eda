#include "route.h"
#include "graph.h"
#include "route_tree.h"
#include "util.h"
#include "log.h"

void update_costs(const RRGraph &g, congestion_t *congestion, float pres_fac, float acc_fac)
{
	for_all_vertices(g, [&congestion, &pres_fac, &acc_fac] (const RRNode &v) -> void {
			int rr_node = id(v);
			int occ = congestion[rr_node].occ;
			int capacity = v.properties.capacity;
			if (occ > capacity) {
				congestion[rr_node].acc_cost += (occ - capacity) * acc_fac;
				congestion[rr_node].pres_cost = 1. + (occ + 1 - capacity) * pres_fac;
			} else if (occ == capacity) {
				/* If occ == capacity, we don't need to increase acc_cost, but a change    *
				 * in pres_fac could have made it necessary to recompute the cost anyway.  */
				congestion[rr_node].pres_cost = 1. + pres_fac;
			}
			});
}

void update_one_cost_internal(const RRNode &rr_node, congestion_t &congestion, /*int net_id, */int delta, float pres_fac, bool lock, lock_perf_t *lock_perf)
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
		/*rr_node.properties.lock->lock();*/
	}
	
	congestion.occ += delta;

	assert(congestion.occ >= 0);

	if (congestion.occ < rr_node.properties.capacity) {
		congestion.pres_cost = 1;
	} else {
		congestion.pres_cost = 1 + (congestion.occ + 1 - rr_node.properties.capacity) * pres_fac;
	}

	/*if (delta > 0) {*/
		/*assert(std::find(begin(rr_node.properties.users), begin(rr_node.properties.users), net_id) == end(rr_node.properties.users));*/
		/*rr_node.properties.users.push_back(net_id);*/
	/*} else {*/
		/*assert(delta < 0);*/
		/*auto iter = std::find(begin(rr_node.properties.users), begin(rr_node.properties.users), net_id);*/
		/*assert(iter != end(rr_node.properties.users));*/
		/*rr_node.properties.users.erase(iter);*/
	/*}*/

	if (lock) {
		congestion.lock.unlock();
	}
		
	char buffer[256];
	sprintf_rr_node(id(rr_node), buffer);
	zlog_level(delta_log, ROUTER_V2, "Update cost of %s delta: %d new_occ: %d pres_fac: %g\n", buffer, delta, congestion.occ, pres_fac);
}

void update_one_cost(const RRGraph &g, congestion_t *congestion, const vector<int>::const_iterator &rr_nodes_begin, const vector<int>::const_iterator &rr_nodes_end, /*int net_id,*/ int delta, float pres_fac, bool lock, lock_perf_t *lock_perf)
{
	/*const RouteTreeNode *last = nullptr;*/
	for (auto iter = rr_nodes_begin; iter != rr_nodes_end; ++iter) {
		/* we don't update get_source because that's the link to existing route tree and the cost is handled by update_one_cost_internal */
		/*const RouteTreeNode &rt_node = get_target(rt.graph, get_edge(rt.graph, rt_edge_id));*/
		/*last = &get_target(rt.graph, get_edge(rt.graph, rt_edge_id));*/
		/*RRNode &rr_node = get_vertex(g, rt_node.properties.rr_node);*/
		const RRNode &rr_node = get_vertex(g, *iter);
		update_one_cost_internal(rr_node, congestion[*iter], /*net_id,*/ delta, pres_fac, lock, lock_perf);
	}
	/*RRNode &rr_node = get_vertex(g, last->properties.rr_node);*/
	/*update_one_cost_internal(rr_node, delta, pres_fac);*/
}

void update_one_cost(const RRGraph &g, congestion_t *congestion, route_tree_t &rt, const RouteTreeNode &node, int delta, float pres_fac, bool lock)
{
	assert(node.properties.valid);

	const RRNode &rr_node = get_vertex(g, node.properties.rr_node);

	update_one_cost_internal(rr_node, congestion[id(rr_node)], /*-1,*/ delta, pres_fac, lock, nullptr);

	for_all_out_edges(rt.graph, node, [&g, &rt, &congestion, &delta, &pres_fac, &lock] (const RouteTreeEdge &e) -> void {
			update_one_cost(g, congestion, rt, get_vertex(rt.graph, get_target(rt.graph, id(e))), delta, pres_fac, lock);
			});
}

