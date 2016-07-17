#ifndef EXPAND_H
#define EXPAND_H

float get_delay(const rr_edge_property_t &e, const rr_node_property_t &v, float unbuffered_upstream_R);

/* fine grained locking based expand neighbor */
template<typename ShouldExpandFunc>
void expand_neighbors_with_fine_grain_lock(const RRGraph &g, int current, const route_state_t *state, const congestion_locked_t *congestion, const rr_node_property_t &target, float criticality_fac, float astar_fac, std::priority_queue<route_state_t> &heap, const ShouldExpandFunc &should_expand, bool lock, perf_t *perf, lock_perf_t *lock_perf)
{
	for (const auto &e : get_out_edges(g, current)) {
		int neighbor_id = get_target(g, e);
		const auto &neighbor_p = get_vertex_props(g, neighbor_id);

		char buffer[256];
		sprintf_rr_node(neighbor_id, buffer);
		zlog_level(delta_log, ROUTER_V3, "\tNeighbor: %s ", buffer);

		if (perf) {
			++perf->num_neighbor_visits;
		}

		if (!should_expand(neighbor_id)) {
			continue;
		}

		route_state_t item;

		item.rr_node = neighbor_id;
		item.prev_edge = e;

		const auto &e_p = get_edge_props(g, e);

		const route_state_t *current_state = &state[current];

		float upstream_R = e_p.R + neighbor_p.R;
		if (!e_p.buffered) {
			upstream_R += current_state->upstream_R;
		}
		item.upstream_R = upstream_R;

		/*if (lock) {*/
		/*[>if (lock_perf) {<]*/
		/*[>++lock_perf->num_lock_tries;<]*/
		/*[>}<]*/
		/*[>if (!neighbor_p.lock->try_lock()) {<]*/
		/*[>if (lock_perf) {<]*/
		/*[>++lock_perf->num_lock_waits;<]*/
		/*[>}<]*/
		/*[>neighbor_p.lock->lock();<]*/
		/*[>} <]*/
		/*neighbor_p.lock->lock();*/
		/*}*/
		/*if (lock) {*/
		/*neighbor_p.lock->unlock();*/
		/*}*/
		//float delay = get_delay(e_p, neighbor_p, unbuffered_upstream_R);
		float delay;
		if (e_p.buffered) {
			delay = e_p.switch_delay + neighbor_p.C * (e_p.R + 0.5 * neighbor_p.R);
		} else {
			delay = e_p.switch_delay + neighbor_p.C * (current_state->upstream_R + e_p.R + 0.5 * neighbor_p.R);
		}

		item.delay = current_state->delay + delay;

		float congestion_cost = get_congestion_cost(congestion[item.rr_node].cong, neighbor_p.cost_index);

		float known_cost = criticality_fac * delay + (1 - criticality_fac) * congestion_cost;

		item.known_cost = current_state->known_cost + known_cost;

		float expected_cost = get_timing_driven_expected_cost(get_vertex_props(g, item.rr_node), target, criticality_fac, upstream_R);

		item.cost = item.known_cost + astar_fac * expected_cost;

		heap.push(item);

		if (perf) {
			++perf->num_heap_pushes;
		}

		zlog_level(delta_log, ROUTER_V3, " [cost: %g known_cost: %g][occ/cap: %d/%d pres: %g acc: %g][edge_delay: %g edge_R: %g node_R: %g node_C: %g] \n",
				item.cost, item.known_cost, 
				congestion[item.rr_node].cong.occ, neighbor_p.capacity, congestion[item.rr_node].cong.pres_cost, congestion[item.rr_node].cong.acc_cost,
				e_p.switch_delay, e_p.R, neighbor_p.R, neighbor_p.C);
	}
}

#endif
