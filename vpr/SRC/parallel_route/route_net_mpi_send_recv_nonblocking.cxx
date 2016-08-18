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
#include "route_net_mpi_send_recv_reduced_comm.h"
#include "clock.h"

std::shared_ptr<vector<path_node_t>> get_path(int sink_rr_node_id, const route_state_t *state, const route_tree_t &rt, const RRGraph &g, int vpr_net_id);
float get_timing_driven_expected_cost(const rr_node_property_t &current, const rr_node_property_t &target, float criticality_fac, float R_upstream);
float get_delay(const rr_edge_property_t &e, const rr_node_property_t &v, float unbuffered_upstream_R);

int get_free_buffer(mpi_context_t *mpi, int size)
{
	int data_index;

	if (!mpi->free_send_data_index.empty()) {
		data_index = mpi->free_send_data_index.front();

		mpi->free_send_data_index.pop();

		if (mpi->pending_send_data_raw[data_index].first->size() < size) {
			mpi->pending_send_data_raw[data_index].first->resize(2*size);
		}
	} else {
		data_index = mpi->pending_send_data_raw.size();

		mpi->pending_send_data_raw.emplace_back(new vector<int>(size), 0);
	}

	return data_index;
}

int get_free_req(mpi_context_t *mpi, int data_index)
{
	int req_index;

	if (!mpi->free_send_req_index.empty()) {
		req_index = mpi->free_send_req_index.front();

		mpi->free_send_req_index.pop();

		mpi->pending_send_req_data_ref[req_index] = data_index;

		zlog_level(delta_log, ROUTER_V3, "Existing req, req_index %d data_index %d\n", req_index, data_index);
	} else {
		req_index = mpi->pending_send_req.size();

		mpi->pending_send_req.emplace_back(MPI_REQUEST_NULL);
		mpi->pending_send_req_data_ref.emplace_back(data_index);

		zlog_level(delta_log, ROUTER_V3, "New req, req_index %d data_index %d\n", req_index, data_index);
	}

	return req_index;
}

void broadcast(int data_index, int size, int tag, mpi_context_t *mpi)
{
	for (int i = 0; i < mpi->comm_size; ++i) {
		if (i != mpi->rank) {
			zlog_level(delta_log, ROUTER_V2, "Sent packet from %d to %d\n", mpi->rank, i);

			int req_index = get_free_req(mpi, data_index);
			int error;

			assert(mpi->pending_send_req[req_index] == MPI_REQUEST_NULL);

			if (data_index < 0) {
				error = MPI_Isend(nullptr, 0, MPI_INT, i, tag, mpi->comm, &mpi->pending_send_req[req_index]);
			} else {
				pair<vector<int> *, int> *data = &mpi->pending_send_data_raw[data_index];
				assert(data->second >= 0);
				++data->second;
				error = MPI_Isend(data->first->data(), size, MPI_INT, i, tag, mpi->comm, &mpi->pending_send_req[req_index]);
			}

			assert(error == MPI_SUCCESS);
		}
	}
}

void broadcast_trailer(mpi_context_t *mpi)
{
	zlog_level(delta_log, ROUTER_V2, "Sending trailer packet\n");

	broadcast(-1, 0, LAST_TAG, mpi);
}

void broadcast_rip_up_no_progress_reuse(int net_id, mpi_context_t *mpi)
{
	zlog_level(delta_log, ROUTER_V3, "MPI sent rip up");

	broadcast(-1, 0, RIP_UP_TAG + net_id, mpi);
}

void broadcast_pending_cost_updates_reduced_no_progress_reuse(const vector<RRNode> &added_nodes, int net_id, int delta, mpi_context_t *mpi)
{
	if (mpi->comm_size < 2) {
		return;
	}

	assert(!added_nodes.empty());

	int data_index = get_free_buffer(mpi, added_nodes.size()*2);
	vector<int> *buffer = mpi->pending_send_data_raw[data_index].first;

	assert(buffer->size() >= added_nodes.size()*2);
	assert(buffer->size() <= mpi->buffer_size);

	node_update_t *cost_updates = (node_update_t *)buffer->data();
	for (int i = 0; i < added_nodes.size(); ++i) {
		cost_updates[i].rr_node = added_nodes[i];
		cost_updates[i].delta = delta;
	}

	zlog_level(delta_log, ROUTER_V3, "MPI sent a path of length %lu\n", added_nodes.size());

	broadcast(data_index, added_nodes.size()*2, COST_UPDATE_TAG+net_id, mpi);
}

void progress_sends_waitall(mpi_context_t *mpi)
{
	if (mpi->comm_size < 2) {
		return;
	}

	assert(!mpi->pending_send_req.empty());

	vector<MPI_Status> statuses(mpi->pending_send_req.size());

	int error = MPI_Waitall(mpi->pending_send_req.size(), mpi->pending_send_req.data(), statuses.data());

	assert(error == MPI_SUCCESS);

	for (int i = 0; i < statuses.size(); ++i) {
		if (statuses[i].MPI_TAG == MPI_ANY_TAG && statuses[i].MPI_SOURCE == MPI_ANY_SOURCE && statuses[i].MPI_ERROR == MPI_SUCCESS) {
			/* empty status */
		} else {
			int data_index = mpi->pending_send_req_data_ref[i];
			
			if (data_index >= 0) {
				assert(mpi->pending_send_data_raw[data_index].second > 0);
				--mpi->pending_send_data_raw[data_index].second;
			}

			if (mpi->pending_send_data_raw[data_index].second == 0) {
				mpi->free_send_data_index.push(data_index);
			}

			mpi->free_send_req_index.push(i);
		}
	}
}

void progress_sends_testsome(mpi_context_t *mpi)
{
	if (mpi->comm_size < 2) {
		return;
	}

	assert(!mpi->pending_send_req.empty());

	if (mpi->completed_send_indices.size() < mpi->pending_send_req.size()) {
		mpi->completed_send_indices.resize(mpi->pending_send_req.size());
	}

	int num_completed;

	int error = MPI_Testsome(mpi->pending_send_req.size(), mpi->pending_send_req.data(),
					&num_completed, mpi->completed_send_indices.data(), MPI_STATUSES_IGNORE);

	assert(error == MPI_SUCCESS);
	assert(num_completed != MPI_UNDEFINED);

	for (int i = 0; i < num_completed; ++i) {
		int completed_index = mpi->completed_send_indices[i];

		assert(mpi->pending_send_req[completed_index] == MPI_REQUEST_NULL);

		int data_index = mpi->pending_send_req_data_ref[completed_index];

		if (data_index >= 0) {
			zlog_level(delta_log, ROUTER_V3, "Reducing ref count, req_index %d data_index %d new_ref_count %d\n", completed_index, data_index, mpi->pending_send_data_raw[data_index].second);

			assert(mpi->pending_send_data_raw[data_index].second > 0);
			--mpi->pending_send_data_raw[data_index].second;

			if (mpi->pending_send_data_raw[data_index].second == 0) {
				zlog_level(delta_log, ROUTER_V3, "Freeing data, req_index %d data_index %d\n", completed_index, data_index);

				mpi->free_send_data_index.push(data_index);
			}
		}

		zlog_level(delta_log, ROUTER_V3, "Freeing req, req_index %d\n");

		mpi->free_send_req_index.push(completed_index);
	}
}

void progress_sends_waitsome(mpi_context_t *mpi)
{
	assert(!mpi->pending_send_req.empty());

	int num_completed;
	vector<int> completed_indices(mpi->pending_send_req.size());
	vector<MPI_Status> completed_statuses(mpi->pending_send_req.size());

	int error = MPI_Waitsome(mpi->pending_send_req.size(), mpi->pending_send_req.data(),
					&num_completed, completed_indices.data(), completed_statuses.data());

	assert(error == MPI_SUCCESS);

	vector<std::shared_ptr<vector<node_update_t>>> uncompleted_send_data;
	vector<MPI_Request> uncompleted_send_req;

	int cur_completed = 0;
	for (int i = 0; i < mpi->pending_send_req.size(); ++i) {
		if (mpi->pending_send_req[i] != MPI_REQUEST_NULL) {
			uncompleted_send_req.push_back(mpi->pending_send_req[i]);
			uncompleted_send_data.push_back(mpi->pending_send_data[i]);
		} else {
			assert(completed_indices[cur_completed] == i);
			++cur_completed;
		}
	}

	mpi->pending_send_req = uncompleted_send_req;
	mpi->pending_send_data = uncompleted_send_data;
}

void handle_packet_2(int *data, const MPI_Status &status, const vector<net_t> &nets, congestion_t *congestion, const RRGraph &g, float pres_fac, vector<vector<RRNode>> &net_route_trees, mpi_context_t *mpi, mpi_perf_t *mpi_perf)
{
	//auto update_start = clock::now();

	int num_recvd;
	int error = MPI_Get_count(&status, MPI_INT, &num_recvd);

	assert(error == MPI_SUCCESS);

	//if (mpi_perf) {
		//mpi_perf->total_update_time += clock::now()-update_start;
	//}

	if (status.MPI_TAG == LAST_TAG) {
		assert(num_recvd == 0);

		zlog_level(delta_log, ROUTER_V3, "MPI received last update from %d\n", status.MPI_SOURCE);

		assert(!mpi->received_last_update[status.MPI_SOURCE]);
		mpi->received_last_update[status.MPI_SOURCE] = true;

		//if (mpi_perf) {
			//mpi_perf->total_update_time += clock::now()-update_start;
		//}
	} else {
		//#define RIP_UP_TAG LAST_TAG+1
		//#define COST_UPDATE_TAG RIP_UP_TAG+0x1000000 //support for up to 16M nets
		//
		assert(status.MPI_TAG >= RIP_UP_TAG);
		if (status.MPI_TAG >= COST_UPDATE_TAG) {
			assert(num_recvd % 2 == 0);

			int net_id = status.MPI_TAG - COST_UPDATE_TAG;
			assert(net_id >= 0 && net_id < net_route_trees.size());

			zlog_level(delta_log, ROUTER_V3, "MPI received a net %d path of length %d from %d\n", nets[net_id].vpr_id, num_recvd/2, status.MPI_SOURCE);

			//auto update_start = clock::now();

			const node_update_t *d = (node_update_t *)data;

			//update_start = clock::now();

			for (int j = 0; j < num_recvd/2; ++j) {
				update_one_cost_internal(d[j].rr_node, g, congestion, d[j].delta, pres_fac);
				//zlog_level(delta_log, ROUTER_V3, "MPI update, source %d node %d delta %d\n", status.MPI_SOURCE, d[j].rr_node, d[j].delta);

				if (d[j].delta > 0) {
					//assert(find(begin(net_route_trees[net_id]), end(net_route_trees[net_id]), d[j].rr_node) == end(net_route_trees[net_id]));
				
					net_route_trees[net_id].push_back(d[j].rr_node);
				} else {
					assert(d[j].delta < 0);
					//assert(find(begin(net_route_trees[net_id]), end(net_route_trees[net_id]), d[j].rr_node) != end(net_route_trees[net_id]));
					net_route_trees[net_id].erase(std::remove(begin(net_route_trees[net_id]), end(net_route_trees[net_id]), d[j].rr_node), end(net_route_trees[net_id]));
				}
			}

			//if (mpi_perf) {
				//mpi_perf->total_update_time += clock::now()-update_start;
			//}
		} else {
			assert(num_recvd == 0);

			int net_id = status.MPI_TAG - RIP_UP_TAG;
			assert(net_id >= 0 && net_id < net_route_trees.size());

			zlog_level(delta_log, ROUTER_V3, "MPI ripping up net %d route tree of size %lu from %d\n", nets[net_id].vpr_id, net_route_trees[net_id].size(), status.MPI_SOURCE);

			//auto update_start = clock::now();

			assert(!net_route_trees[net_id].empty());

			for (const auto &node : net_route_trees[net_id]) {
				update_one_cost_internal(node, g, congestion, -1, pres_fac);
				//zlog_level(delta_log, ROUTER_V3, "MPI rip up, source %d node %d delta -1\n", status.MPI_SOURCE, node);
				//assert(find(begin(net_route_trees[net_id]), end(net_route_trees[net_id]), node) != end(net_route_trees[net_id]));
			}

			net_route_trees[net_id].clear();

			//if (mpi_perf) {
				//mpi_perf->total_update_time += clock::now()-update_start;
			//}
		}
	}
}

void bootstrap_irecv(mpi_context_t *mpi)
{
	assert(mpi->comm_size == mpi->pending_recv_data_flat.size());
	assert(mpi->comm_size == mpi->pending_recv_req_flat.size());

	for (int pid = 0; pid < mpi->comm_size; ++pid) {
		if (pid != mpi->rank) {
			assert(!mpi->received_last_update[pid]);
			assert(mpi->pending_recv_req_flat[pid] == MPI_REQUEST_NULL);

			int error = MPI_Irecv(mpi->pending_recv_data_flat[pid].data(), mpi->pending_recv_data_flat[pid].size(), MPI_INT, pid, MPI_ANY_TAG, mpi->comm, &mpi->pending_recv_req_flat[pid]);

			assert(error == MPI_SUCCESS);
		}
	}
}

void sync_irecv(const vector<net_t> &nets, congestion_t *congestion, const RRGraph &g, float pres_fac, vector<vector<RRNode>> &net_route_trees, mpi_context_t *mpi, mpi_perf_t *mpi_perf)
{
	if (mpi->comm_size < 2) {
		return;
	}

	int num_completed;

	assert(mpi->comm_size == mpi->pending_recv_req_flat.size());
	assert(mpi->comm_size == mpi->completed_recv_indices.size());
	assert(mpi->comm_size == mpi->completed_recv_statuses.size());

	int error = MPI_Testsome(mpi->pending_recv_req_flat.size(), mpi->pending_recv_req_flat.data(),
			&num_completed, mpi->completed_recv_indices.data(), mpi->completed_recv_statuses.data());
	assert(error == MPI_SUCCESS);

	while (num_completed > 0) {
		for (int i = 0; i < num_completed; ++i) {
			int completed_rank = mpi->completed_recv_indices[i];

			assert(completed_rank != mpi->rank);
			assert(mpi->pending_recv_req_flat[completed_rank] == MPI_REQUEST_NULL);

			handle_packet_2(mpi->pending_recv_data_flat[completed_rank].data(), mpi->completed_recv_statuses[i], nets, congestion, g, pres_fac, net_route_trees, mpi, mpi_perf);

			if (!mpi->received_last_update[completed_rank]) {
				error = MPI_Irecv(mpi->pending_recv_data_flat[completed_rank].data(), mpi->pending_recv_data_flat[completed_rank].size(), MPI_INT, completed_rank, MPI_ANY_TAG, mpi->comm, &mpi->pending_recv_req_flat[completed_rank]);

				assert(error == MPI_SUCCESS);
			}
		}

		error = MPI_Testsome(mpi->pending_recv_req_flat.size(), mpi->pending_recv_req_flat.data(),
				&num_completed, mpi->completed_recv_indices.data(), mpi->completed_recv_statuses.data());
		assert(error == MPI_SUCCESS);
	}
}

template<typename ShouldExpandFunc>
void expand_neighbors_mpi_recv(const RRGraph &g, RRNode current, const route_state_t *state, congestion_t *congestion, const rr_node_property_t &target, float criticality_fac, float astar_fac, float pres_fac, const ShouldExpandFunc &should_expand, std::priority_queue<route_state_t> &heap, perf_t *perf)
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

void route_net_mpi_send_recv_nonblocking(const RRGraph &g, int vpr_id, int net_id, const source_t *source, const vector<sink_t *> &sinks, const route_parameters_t &params, const vector<net_t> &nets, route_state_t *state, congestion_t *congestion, route_tree_t &rt, t_net_timing &net_timing, vector<vector<RRNode>> &net_route_trees, vector<sink_t *> &routed_sinks, vector<sink_t *> &unrouted_sinks, mpi_context_t *mpi, perf_t *perf, mpi_perf_t *mpi_perf, bool sync_only_once, bool delayed_progress)
{
    using clock = myclock;

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
		route_tree_set_node_properties(rt, root_rt_node, true, RRGraph::null_edge(), source_rr_node_p.R, 0.5 * source_rr_node_p.R * source_rr_node_p.C);
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

	int sync_freq = (int)floor(sqrt(sorted_sinks.size()));

	vector<RRNode> all_added_nodes;

	int isink = 0;
	for (const auto &sink : sorted_sinks) {
		if (delayed_progress && isink > 0) {
			auto progress_start = clock::now();
			progress_sends_testsome(mpi);
			if (mpi_perf) {
				mpi_perf->total_send_testsome_time += clock::now()-progress_start;
			}
		}

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

				expand_neighbors_mpi_recv(g, item.rr_node, state, congestion, sink_rr_node, sink->criticality_fac, params.astar_fac, params.astar_fac, [&g, &sink, &sink_rr_node, &v, &item] (const RRNode &n) -> bool {

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
					//if (pid[n] != -1 && pid[n] != mpi.rank) {
						//zlog_level(delta_log, ROUTER_V3, " pid %d not in current partition %d\n", pid[n], mpi.rank);
						//return false;
					//}

					return true;
				}, heap, perf);
			} else {
				zlog_level(delta_log, ROUTER_V3, "\tNot expanding neighbor because known_cost %g > %g or cost %g > %g\n", item.known_cost, state[item.rr_node].known_cost, item.cost, state[item.rr_node].cost);
			}
		}

		//bool sync_costs = (isink % sync_freq) == 0;
		//bool sync_costs = true;

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
				if (valid(n.prev_edge) || n.rr_node_id == source->rr_node) {
					added_nodes.push_back(n.rr_node_id);
					all_added_nodes.push_back(n.rr_node_id);
				}
			}

			update_one_cost(g, congestion, added_nodes.begin(), added_nodes.end(), 1, params.pres_fac);

			if (!sync_only_once) {
				auto broadcast_start = clock::now();
				broadcast_pending_cost_updates_reduced_no_progress_reuse(added_nodes, net_id, 1, mpi);
				if (mpi_perf) {
					mpi_perf->total_broadcast_time += clock::now()-broadcast_start;
				}
			}

			if (sink->id != -1) {
				net_timing.delay[sink->id+1] = get_vertex_props(rt.graph, route_tree_get_rt_node(rt, sink->rr_node)).delay;
			}

			zlog_level(delta_log, ROUTER_V2, "Routed sink %d\n", sink->rr_node);
		}

		if (!sync_only_once) {
			auto sync_start = clock::now();
			sync_irecv(nets, congestion, g, params.pres_fac, net_route_trees, mpi, mpi_perf);
			if (mpi_perf) {
				mpi_perf->total_sync_time += clock::now()-sync_start;
			}
		}

		if (!delayed_progress) {
			auto progress_start = clock::now();
			progress_sends_testsome(mpi);
			if (mpi_perf) {
				mpi_perf->total_send_testsome_time += clock::now()-progress_start;
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

	if (sync_only_once) {
		auto broadcast_start = clock::now();
		broadcast_pending_cost_updates_reduced_no_progress_reuse(all_added_nodes, net_id, 1, mpi);
		if (mpi_perf) {
			mpi_perf->total_broadcast_time += clock::now()-broadcast_start;
		}

		auto sync_start = clock::now();
		sync_irecv(nets, congestion, g, params.pres_fac, net_route_trees, mpi, mpi_perf);
		if (mpi_perf) {
			mpi_perf->total_sync_time += clock::now()-sync_start;
		}

		auto progress_start = clock::now();
		progress_sends_testsome(mpi);
		if (mpi_perf) {
			mpi_perf->total_send_testsome_time += clock::now()-progress_start;
		}
	}

	zlog_level(delta_log, ROUTER_V2, "Routed net %d\n", vpr_id);
	//zlog_level(delta_log, ROUTER_V1, "\n");

	//return unrouted_sinks;

	/*delete [] state;*/

	/*check_route_tree(rt, net, g);*/
}

