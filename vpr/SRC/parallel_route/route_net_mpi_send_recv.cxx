#include "pch.h"

#include "vpr_types.h"
#include "path_delay.h"

#include "route.h"
#include "route_tree.h"
#include "congestion.h"
#include "trace.h"
#include "router.h"
#include "filtered_graph.h"
#include "log.h"

std::shared_ptr<vector<path_node_t>> get_path(int sink_rr_node_id, const route_state_t *state, const route_tree_t &rt, const RRGraph &g, int vpr_net_id);
float get_timing_driven_expected_cost(const rr_node_property_t &current, const rr_node_property_t &target, float criticality_fac, float R_upstream);
float get_delay(const rr_edge_property_t &e, const rr_node_property_t &v, float unbuffered_upstream_R);

void broadcast_pending_cost_updates(queue<RRNode> &cost_update_q, int delta, int this_pid, int num_procs, MPI_Comm comm, vector<ongoing_transaction_t> &transactions)
{
	if (cost_update_q.empty()) {
		return;
	}

	auto data = make_shared<vector<node_update_t>>();

	while (!cost_update_q.empty()) {
		node_update_t d;

		d.rr_node = cost_update_q.front();
		d.delta = delta;

		data->push_back(d);

		cost_update_q.pop();
	}

	assert(data->size() > 0 && data->size() < 4096);
	for (int i = 0; i < num_procs; ++i) {
		if (i != this_pid) {
			zlog_level(delta_log, ROUTER_V3, "Sent a path of length %lu from %d to %d\n", data->size(), this_pid, i);

			ongoing_transaction_t trans;

			trans.data = data;
			assert(MPI_Isend(trans.data->data(), trans.data->size()*2, MPI_INT, i, 3399, comm, &trans.req) == MPI_SUCCESS);

			transactions.push_back(trans);
			assert(MPI_Request_free(&trans.req) == MPI_SUCCESS);
			//assert(MPI_ISend(trans.data->data(), trans.data->size()*2, MPI_INT, i, 0, comm) == MPI_SUCCESS);
		}
	}
}

void sync(congestion_t *congestion, const RRGraph &g, float pres_fac, int this_pid, int num_procs, MPI_Comm comm, mpi_perf_t *mpi_perf)
{
	int flag = 0;

    using clock = std::chrono::high_resolution_clock;

	for (int pid = 0; pid < num_procs; ++pid) {
		if (pid != this_pid) {
			MPI_Status status;
			auto probe_start = clock::now();
			assert(MPI_Iprobe(pid, 0, comm, &flag, &status) == MPI_SUCCESS);
			if (mpi_perf) {
				mpi_perf->total_probe_time += clock::now()-probe_start;
			}
			//assert(MPI_Probe(pid, 0, comm, &status) == MPI_SUCCESS);

			while (flag) {
				int count = -1;
				assert(MPI_Get_count(&status, MPI_INT, &count) == MPI_SUCCESS);

				//if (count <= 0) {
					//zlog_warn(delta_log, "MPI update count from %d is %d\n", status.MPI_SOURCE, count);
					//break;
				//}

				//if (count > 0) {
				
				assert(count > 0 && count % 2 == 0);

				node_update_t *d = new node_update_t[count/2];

				assert(MPI_Recv(d, count, MPI_INT, status.MPI_SOURCE, status.MPI_TAG, comm, &status) == MPI_SUCCESS);

				int new_count;
				assert(MPI_Get_count(&status, MPI_INT, &new_count) == MPI_SUCCESS);
				assert(count == new_count);

				zlog_level(delta_log, ROUTER_V3, "Received a path of length %d from %d\n", count/2, status.MPI_SOURCE);

				for (int i = 0; i < count/2; ++i) {
					update_one_cost_internal(d[i].rr_node, g, congestion, d[i].delta, pres_fac);
					zlog_level(delta_log, ROUTER_V3, "MPI update, source %d node %d delta %d\n", status.MPI_SOURCE, d[i].rr_node, d[i].delta);
				}

				delete [] d;
				//}

				probe_start = clock::now();
				assert(MPI_Iprobe(pid, 0, comm, &flag, &status) == MPI_SUCCESS);
				if (mpi_perf) {
					mpi_perf->total_probe_time += clock::now()-probe_start;
				}
			}
		}
	}
}


template<typename ShouldExpandFunc>
void expand_neighbors_mpi_recv(const RRGraph &g, RRNode current, const route_state_t *state, congestion_t *congestion, int this_pid, int num_procs, MPI_Comm comm, const rr_node_property_t &target, float criticality_fac, float astar_fac, float pres_fac, const ShouldExpandFunc &should_expand, std::priority_queue<route_state_t> &heap, perf_t *perf)
{
	for (const auto &e : get_out_edges(g, current)) {
		int neighbor = get_target(g, e);
		const auto &neighbor_p = get_vertex_props(g, neighbor);

		char buffer[256];
		sprintf_rr_node(neighbor, buffer);
		zlog_level(delta_log, ROUTER_V3, "\tNeighbor: %s ", buffer);

		if (perf) {
			++perf->num_neighbor_visits;
		}

		if (!should_expand(neighbor)) {
			continue;
		}

		route_state_t item;

		item.rr_node = neighbor;
		item.prev_edge = e;

		const auto &e_p = get_edge_props(g, e);

		const route_state_t *current_state = &state[current];

		float unbuffered_upstream_R = current_state->upstream_R;
		float upstream_R = e_p.R + neighbor_p.R;
		if (!e_p.buffered) {
			upstream_R += unbuffered_upstream_R;
		}
		item.upstream_R = upstream_R;

		//float congestion_cost = get_congestion_cost_mpi_recv(item.rr_node, g, congestion, comm, this_pid, num_procs, pres_fac);
		extern t_rr_indexed_data *rr_indexed_data;
		float congestion_cost = rr_indexed_data[neighbor_p.cost_index].base_cost * congestion[neighbor].acc_cost * congestion[neighbor].pres_cost;

		float delay = get_delay(e_p, neighbor_p, unbuffered_upstream_R);

		item.delay = current_state->delay + delay;

		float known_cost = criticality_fac * delay + (1 - criticality_fac) * congestion_cost;

		item.known_cost = current_state->known_cost + known_cost;

		float expected_cost = get_timing_driven_expected_cost(neighbor_p, target, criticality_fac, upstream_R);

		item.cost = item.known_cost + astar_fac * expected_cost;

		heap.push(item);

		if (perf) {
			++perf->num_heap_pushes;
		}

		zlog_level(delta_log, ROUTER_V3, " [cost: %g known_cost: %g][occ/cap: %d/%d pres: %g acc: %g][edge_delay: %g edge_R: %g node_R: %g node_C: %g] \n",
				item.cost, item.known_cost, 
				congestion[item.rr_node].occ, neighbor_p.capacity, congestion[item.rr_node].pres_cost, congestion[item.rr_node].acc_cost,
				e_p.switch_delay, e_p.R, neighbor_p.R, neighbor_p.C);
	}
}

void route_net_mpi_send_recv(const RRGraph &g, int this_pid, int vpr_id, const source_t *source, const vector<sink_t *> &sinks, const route_parameters_t &params, route_state_t *state, congestion_t *congestion, int num_procs, MPI_Comm comm, vector<ongoing_transaction_t> &transactions, route_tree_t &rt, t_net_timing &net_timing, vector<sink_t *> &routed_sinks, vector<sink_t *> &unrouted_sinks, perf_t *perf, mpi_perf_t *mpi_perf)
{
    using clock = std::chrono::high_resolution_clock;

	std::priority_queue<route_state_t> heap;

	vector<int> modified;

	zlog_level(delta_log, ROUTER_V1, "Routing net %d (%lu sinks)\n", vpr_id, sinks.size());

	vector<sink_t *> sorted_sinks = sinks;
	std::sort(begin(sorted_sinks), end(sorted_sinks), [] (const sink_t *a, const sink_t *b) -> bool {
			return a->criticality_fac > b->criticality_fac;
			});

	char buffer[256];

	bool empty_route_tree;

	if (route_tree_empty(rt)) {
		/* special case */
		sprintf_rr_node(source->rr_node, buffer);
		zlog_level(delta_log, ROUTER_V2, "Empty route tree. Setting root to %s\n", buffer);

		RouteTreeNode root_rt_node = route_tree_add_rr_node(rt, source->rr_node, g);
		const auto &source_rr_node_p = get_vertex_props(g, source->rr_node);
		route_tree_set_node_properties(rt, root_rt_node, true, source_rr_node_p.R, 0.5 * source_rr_node_p.R * source_rr_node_p.C);
		route_tree_add_root(rt, source->rr_node);

		//update_one_cost_internal(source->rr_node, g, congestion, 1, params.pres_fac);
		empty_route_tree = true;
	} else {
		empty_route_tree = false;

		if (rt.root_rt_node_id != -1) {
			const auto &rt_root_p = get_vertex_props(rt.graph, rt.root_rt_node_id);
			if (source && rt_root_p.rr_node != source->rr_node) {
				char root[256];
				char source_str[256];
				sprintf_rr_node(rt_root_p.rr_node, root);
				sprintf_rr_node(source->rr_node, source_str);
				zlog_error(delta_log, "Error: Non empty route tree node has root %s that is different from net source %s\n",
						root, source_str);
				assert(false);
			}
		}
	}

	queue<RRNode> cost_update_q;
	int sync_freq = (int)floor(sqrt(sorted_sinks.size()));

	int isink = 0;
	for (const auto &sink : sorted_sinks) {
		sprintf_rr_node(sink->rr_node, buffer);
		zlog_level(delta_log, ROUTER_V1, "Net %d current sink %d: %s criticality: %g BB: %d-%d %d-%d Prev BB: %d-%d %d-%d\n", vpr_id, sink->id, buffer, sink->criticality_fac, sink->current_bounding_box.xmin, sink->current_bounding_box.xmax, sink->current_bounding_box.ymin, sink->current_bounding_box.ymax, sink->previous_bounding_box.xmin, sink->previous_bounding_box.xmax, sink->previous_bounding_box.ymin, sink->previous_bounding_box.ymax);

		const auto &sink_rr_node = get_vertex_props(g, sink->rr_node);

		route_tree_multi_root_add_to_heap(rt, g, sink->rr_node, sink->criticality_fac, params.astar_fac, sink->current_bounding_box, heap, perf);

		bool found_sink = false;
		while (!heap.empty() && !found_sink) {
			auto item = heap.top();
			heap.pop();

			if (perf) {
				++perf->num_heap_pops;
			}

			sprintf_rr_node(item.rr_node, buffer);
			const auto &v = get_vertex_props(g, item.rr_node);
			zlog_level(delta_log, ROUTER_V3, "Current: %s occ/cap: %d/%d prev: %d new_cost: %g new_known: %g new_delay: %g old_cost: %g old_known: %g old_delay: %g\n", buffer, congestion[item.rr_node].occ, v.capacity, valid(item.prev_edge) ? get_source(g, item.prev_edge) : -1, item.cost, item.known_cost, item.delay, state[item.rr_node].cost, state[item.rr_node].known_cost, state[item.rr_node].delay);

			//assert(pid[item.rr_node] == -1 || pid[item.rr_node] == 0);

			if (item.rr_node == sink->rr_node) {
				state[item.rr_node] = item;
				modified.push_back(item.rr_node);
				found_sink = true;
			} else if (item.cost < state[item.rr_node].cost) {
				assert(item.known_cost < state[item.rr_node].known_cost);
				state[item.rr_node] = item;
				modified.push_back(item.rr_node);

				expand_neighbors_mpi_recv(g, item.rr_node, state, congestion, this_pid, num_procs, comm, sink_rr_node, sink->criticality_fac, params.astar_fac, params.astar_fac, [&g, &sink, &sink_rr_node, &v, &item, &this_pid] (const RRNode &n) -> bool {

					/*if (trace_has_node(prev_trace, id(n))) {*/
						/*zlog_level(delta_log, ROUTER_V3, " existing node route tree ");*/
					/*}*/
					const auto &prop = get_vertex_props(g, n);

					if (prop.xhigh < sink->current_bounding_box.xmin
							|| prop.xlow > sink->current_bounding_box.xmax
							|| prop.yhigh < sink->current_bounding_box.ymin
							|| prop.ylow > sink->current_bounding_box.ymax) {
					zlog_level(delta_log, ROUTER_V3, "outside of bounding box\n");
					return false;
					}

					if (prop.type == IPIN
							&& (prop.xhigh != sink_rr_node.xhigh 
								|| prop.yhigh != sink_rr_node.yhigh)) {
					zlog_level(delta_log, ROUTER_V3, "not target IPIN\n");
					return false;
					}

					/*if (net.sinks.size() >= HIGH_FANOUT_NET_LIM) {*/
					/*if (prop.xhigh < target_x - highfanout_rlim*/
							/*|| prop.xlow > target_x + highfanout_rlim*/
							/*|| prop.yhigh < target_y - highfanout_rlim*/
							/*|| prop.ylow > target_y + highfanout_rlim) {*/
						/*return false;*/
					/*}*/
					/*}*/
					//if (pid[n] != -1 && pid[n] != this_pid) {
						//zlog_level(delta_log, ROUTER_V3, " pid %d not in current partition %d\n", pid[n], this_pid);
						//return false;
					//}

					return true;
				}, heap, perf);
			} else {
				zlog_level(delta_log, ROUTER_V3, "\tNot expanding neighbor because known_cost %g > %g or cost %g > %g\n", item.known_cost, state[item.rr_node].known_cost, item.cost, state[item.rr_node].cost);
			}
		}

		//bool sync_costs = (isink % sync_freq) == 0;
		bool sync_costs = true;

		if (!found_sink) {
			assert(heap.empty());
			zlog_error(delta_log, "Error: Failed to find sink %d\n", sink->rr_node);

			assert(find(begin(unrouted_sinks), end(unrouted_sinks), sink) == end(unrouted_sinks));
			unrouted_sinks.push_back(sink);

			if (sink->id != -1) {
				net_timing.delay[sink->id+1] = 0;
			}
		} else {
			assert(find(begin(routed_sinks), end(routed_sinks), sink) == end(routed_sinks));

			routed_sinks.push_back(sink);

			const auto &path = get_path(sink->rr_node, state, rt, g, vpr_id);

			route_tree_add_path(rt, path, g, state);

			bool update_source_cost = false;
			if (empty_route_tree) {
				assert(path->back().rr_node_id == source->rr_node);
				update_source_cost = true;
				empty_route_tree = false;
			}

			vector<RRNode> added_nodes;
			for (const auto &n : *path) {
				if (valid(n.prev_edge) || update_source_cost) {
					if (update_source_cost) {
						assert(n.rr_node_id == source->rr_node);
					}
					added_nodes.push_back(n.rr_node_id);
					cost_update_q.push(n.rr_node_id);
				}
			}

			update_one_cost(g, congestion, added_nodes.begin(), added_nodes.end(), 1, params.pres_fac);

			if (sync_costs) {
				auto broadcast_start = clock::now();
				broadcast_pending_cost_updates(cost_update_q, 1, this_pid, num_procs, comm, transactions);
				if (mpi_perf) {
					mpi_perf->total_broadcast_time += clock::now()-broadcast_start;
				}
			}

			if (sink->id != -1) {
				net_timing.delay[sink->id+1] = get_vertex_props(rt.graph, route_tree_get_rt_node(rt, sink->rr_node)).delay;
			}

			zlog_level(delta_log, ROUTER_V2, "Routed sink %d\n", sink->rr_node);
		}

		if (sync_costs) {
			auto sync_start = clock::now();
			sync(congestion, g, params.pres_fac, this_pid, num_procs, comm, mpi_perf);
			if (mpi_perf) {
				mpi_perf->total_sync_time += clock::now()-sync_start;
			}
		}

		for (const auto &m : modified)  {
			state[m].known_cost = std::numeric_limits<float>::max();
			state[m].cost = std::numeric_limits<float>::max();
		}

		modified.clear();

		heap = std::priority_queue<route_state_t>();

		++isink;
	}

	auto broadcast_start = clock::now();
	broadcast_pending_cost_updates(cost_update_q, 1, this_pid, num_procs, comm, transactions);
	if (mpi_perf) {
		mpi_perf->total_broadcast_time += clock::now()-broadcast_start;
	}
	//zlog_level(delta_log, ROUTER_V1, "\n");

	//return unrouted_sinks;

	/*delete [] state;*/

	/*check_route_tree(rt, net, g);*/
}

