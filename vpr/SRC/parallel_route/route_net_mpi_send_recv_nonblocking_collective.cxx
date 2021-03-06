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
#include "route_net_mpi_send_recv_nonblocking_collective.h"
#include "clock.h"

float get_timing_driven_expected_cost(const rr_node_property_t &current, const rr_node_property_t &target, float criticality_fac, float R_upstream);

void decode_header(unsigned int header, PacketID &packet_id, int &meta)
{
	packet_id = static_cast<PacketID>(header >> 28);
	meta = header & 0xFFFFFFF;
}

int encode_header(PacketID packet_id, unsigned int net_id)
{
	assert(sizeof(unsigned int) == 4);
	assert((net_id & 0xFFFFFFF) == net_id);

	return (static_cast<unsigned int>(packet_id) << 28) | net_id;
}

int get_free_buffer(mpi_context_t *mpi, int size)
{
	int data_index;

	if (!mpi->free_send_data_index.empty()) {
		data_index = mpi->free_send_data_index.front();

		mpi->free_send_data_index.pop();

		if (mpi->pending_send_data_nbc[data_index]->size() < size) {
			mpi->pending_send_data_nbc[data_index]->resize(2*size);
		}
	} else {
		data_index = mpi->pending_send_data_nbc.size();

		mpi->pending_send_data_nbc.emplace_back(new vector<int>(size));
	}

	return data_index;
}

int get_free_buffer_in_place(mpi_context_t *mpi, int size)
{
	return get_free_buffer(mpi, size*mpi->comm_size);
}

static void reuse_req(mpi_context_t *mpi, int req_index, int meta_data_index, int data_index)
{
	assert(mpi->pending_send_req[req_index] == MPI_REQUEST_NULL);

	mpi->pending_send_req_meta_data_ref[req_index] = meta_data_index;
	mpi->pending_send_req_data_ref[req_index] = data_index;
}

static int get_free_req(mpi_context_t *mpi, int meta_data_index, int data_index)
{
	assert((mpi->send_req_queue_tail % 2) == 0);

	int req_index = (mpi->send_req_queue_tail + 2) % mpi->pending_send_req.size();

	mpi->pending_send_req_meta_data_ref[req_index] = meta_data_index;
	mpi->pending_send_req_data_ref[req_index] = data_index;

	zlog_level(delta_log, ROUTER_V3, "Queuing req, req_index %d meta_data_index %d data_index %d\n", req_index, meta_data_index, data_index);

	++mpi->send_req_queue_size;

	return req_index;
}

void start_actual_broadcast(mpi_context_t *mpi, const broadcast_data_t &bcast_d)
{
	int *recvcounts = mpi->pending_send_data_nbc[bcast_d.meta_ref]->data();

	int total_size = 0;

	vector<int> displs(mpi->comm_size);

	zlog_level(delta_log, ROUTER_V3, "Rank recvcount: ");

	for (int i = 0; i < mpi->comm_size; ++i) {
		displs[i] = total_size;
		total_size += recvcounts[i];
		zlog_level(delta_log, ROUTER_V3, "%d ", recvcounts[i]);
	}

	zlog_level(delta_log, ROUTER_V3, "\n");

	zlog_level(delta_log, ROUTER_V3, "Rank displs: ");
	for (int i = 0; i < mpi->comm_size; ++i) {
		zlog_level(delta_log, ROUTER_V3, "%d ", displs[i]);
	}
	zlog_level(delta_log, ROUTER_V3, "\n");

	int *local_data = mpi->pending_send_data_nbc[bcast_d.data_ref]->data();

	int global_data_index = get_free_buffer(mpi, total_size);

	int *global_data = mpi->pending_send_data_nbc[global_data_index]->data();

	memcpy(&global_data[displs[mpi->rank]], local_data, sizeof(int)*recvcounts[mpi->rank]);

	zlog_level(delta_log, ROUTER_V3, "Starting actual broadcast of size %d meta_data_index %d global_data_index %d\n", total_size, bcast_d.meta_ref, global_data_index);

	assert(mpi->active_req == MPI_REQUEST_NULL);

	mpi->active_bcast_data = bcast_d;
	mpi->active_bcast_data.is_meta = false;
	mpi->active_bcast_data.data_ref = global_data_index;

	int error = MPI_Iallgatherv(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, global_data, recvcounts, displs.data(), MPI_INT, mpi->comm, &mpi->active_req);

	assert(error == MPI_SUCCESS);

	zlog_level(delta_log, ROUTER_V3, "Freeing data_index %d\n", bcast_d.data_ref);

	mpi->free_send_data_index.push(bcast_d.data_ref);
}

/* bcast_d cannot be a ref because start_actual_broadcast modifies active_bcast_data */
static void handle_packet(broadcast_data_t bcast_d, const vector<net_t> &nets, congestion_t *congestion, const RRGraph &g, float pres_fac, vector<vector<RRNode>> &net_route_trees, mpi_context_t *mpi)
{
	zlog_level(delta_log, ROUTER_V3, "Completing req meta_data_index %d data_index %d\n", bcast_d.meta_ref, bcast_d.data_ref);

	if (bcast_d.is_meta) {
		start_actual_broadcast(mpi, bcast_d);
	} else {
		zlog_level(delta_log, ROUTER_V3, "Non-meta\n");

		int *recvcounts = mpi->pending_send_data_nbc[bcast_d.meta_ref]->data();
		int *data = mpi->pending_send_data_nbc[bcast_d.data_ref]->data();

		int offset = 0;
		for (int i = 0; i < mpi->comm_size; ++i) {
			assert(recvcounts[i] > 0);

			int *local_data = &data[offset];
			offset += recvcounts[i];

			if (i == mpi->rank) {
				continue;
			}

			PacketID current_pid;
			int net_id;
			decode_header(local_data[0], current_pid, net_id);

			zlog_level(delta_log, ROUTER_V3, "Rank %d header %d pid %d net_id %d\n", i, local_data[0], static_cast<int>(current_pid), net_id);

			assert(static_cast<int>(current_pid) >= 0 && current_pid < PacketID::NUM_PACKET_IDS);
			assert(net_id >= 0 && net_id < nets.size());

			int delta; 
			int *changed_rr_nodes = &local_data[1];

			switch (current_pid) {
				case PacketID::ROUTE:
				case PacketID::RIP_UP:
					delta = current_pid == PacketID::ROUTE ? 1 : -1;

					if (delta > 0) {
						zlog_level(delta_log, ROUTER_V3, "Rank %d net %d route\n", i, nets[net_id].vpr_id);
					} else {
						assert(delta < 0);
						zlog_level(delta_log, ROUTER_V3, "Rank %d net %d rip up\n", i, nets[net_id].vpr_id);
					}

					for (int j = 0; j < recvcounts[i]-1; ++j) {
						update_one_cost_internal(changed_rr_nodes[j], g, congestion, delta, pres_fac);
						//zlog_level(delta_log, ROUTER_V3, "MPI update, source %d node %d delta %d\n", status.MPI_SOURCE, d[j].rr_node, d[j].delta);

						if (delta > 0) {
							assert(find(begin(net_route_trees[net_id]), end(net_route_trees[net_id]), changed_rr_nodes[j]) == end(net_route_trees[net_id]));

							net_route_trees[net_id].push_back(changed_rr_nodes[j]);
						} else {
							assert(delta < 0);
							assert(find(begin(net_route_trees[net_id]), end(net_route_trees[net_id]), changed_rr_nodes[j]) != end(net_route_trees[net_id]));

							net_route_trees[net_id].erase(std::remove(begin(net_route_trees[net_id]), end(net_route_trees[net_id]), changed_rr_nodes[j]), end(net_route_trees[net_id]));
						}
					}

					break;

				case PacketID::RIP_UP_ALL:
					zlog_level(delta_log, ROUTER_V3, "Rank %d rip up entire up net %d of size %lu\n", i, nets[net_id].vpr_id, net_route_trees[net_id].size());

					assert(recvcounts[i] == 1);

					assert(!net_route_trees[net_id].empty());

					for (const auto &node : net_route_trees[net_id]) {
						update_one_cost_internal(node, g, congestion, -1, pres_fac);
						//zlog_level(delta_log, ROUTER_V3, "MPI rip up, source %d node %d delta -1\n", status.MPI_SOURCE, node);
						//assert(find(begin(net_route_trees[net_id]), end(net_route_trees[net_id]), node) != end(net_route_trees[net_id]));
					}

					net_route_trees[net_id].clear();

					break;

				case PacketID::NO_OP:
					zlog_level(delta_log, ROUTER_V3, "Rank %d no op\n", i);

					assert(recvcounts[i] == 1);
					assert(net_id == 0);

					break;

				default:
					zlog_level(delta_log, ROUTER_V3, "Unknown packet ID %d from rank %d\n", static_cast<int>(current_pid), i);

					assert(false);
					break;
			}
		}

		mpi->free_send_data_index.push(bcast_d.meta_ref);
		mpi->free_send_data_index.push(bcast_d.data_ref);

		zlog_level(delta_log, ROUTER_V3, "Freeing meta_data_index %d data_index %d\n", bcast_d.meta_ref, bcast_d.data_ref);
	}
}

void start_meta_broadcast(mpi_context_t *mpi)
{
	if (mpi->active_req == MPI_REQUEST_NULL && !q_empty(&mpi->pending_bcast_data_q)) {
		/* can call broadcast immediately */
		const broadcast_data_t &bcast_d = q_front(&mpi->pending_bcast_data_q);

		assert(bcast_d.is_meta);

		int *meta_data = mpi->pending_send_data_nbc[bcast_d.meta_ref]->data();

		mpi->active_bcast_data = bcast_d;

		int error = MPI_Iallgather(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, meta_data, 1, MPI_INT, mpi->comm, &mpi->active_req);

		assert(error == MPI_SUCCESS);

		assert(q_pop(&mpi->pending_bcast_data_q));
	}
}

void progress_broadcast(const vector<net_t> &nets, congestion_t *congestion, const RRGraph &g, float pres_fac, vector<vector<RRNode>> &net_route_trees, mpi_context_t *mpi)
{
	if (mpi->comm_size < 2) {
		return;
	}

	zlog_level(delta_log, ROUTER_V3, "Progressing broadcast\n");

	int flag = 1;
	while (mpi->active_req != MPI_REQUEST_NULL && flag) {
		int error = MPI_Test(&mpi->active_req, &flag, MPI_STATUS_IGNORE);

		assert(error == MPI_SUCCESS);

		if (flag) {
			assert(mpi->active_req == MPI_REQUEST_NULL);

			/* might start actual broadcast and modifies active_req and active_bcast_data */
			handle_packet(mpi->active_bcast_data, nets, congestion, g, pres_fac, net_route_trees, mpi);

			start_meta_broadcast(mpi);
		}
	}
}

void progress_broadcast_blocking(const vector<net_t> &nets, congestion_t *congestion, const RRGraph &g, float pres_fac, vector<vector<RRNode>> &net_route_trees, mpi_context_t *mpi)
{
	if (mpi->comm_size < 2) {
		return;
	}

	assert(mpi->active_req != MPI_REQUEST_NULL);

	int error = MPI_Wait(&mpi->active_req, MPI_STATUS_IGNORE);

	assert(error == MPI_SUCCESS);

	assert(mpi->active_req == MPI_REQUEST_NULL);

	/* might start actual broadcast and modifies active_req and active_bcast_data */
	handle_packet(mpi->active_bcast_data, nets, congestion, g, pres_fac, net_route_trees, mpi);

	start_meta_broadcast(mpi);
}

void broadcast_dynamic(int data_index, int size, mpi_context_t *mpi)
{
	if (mpi->comm_size < 2) {
		return;
	}

	mpi->max_send_data_size = std::max(mpi->max_send_data_size, size);

	int meta_data_index = get_free_buffer(mpi, mpi->comm_size);

	int *meta_data = mpi->pending_send_data_nbc[meta_data_index]->data();
	meta_data[mpi->rank] = size;

	broadcast_data_t bcast_d;
	bcast_d.is_meta = true;
	bcast_d.meta_ref = meta_data_index;
	bcast_d.data_ref = data_index;

	assert(q_push(&mpi->pending_bcast_data_q, bcast_d));

	mpi->max_send_req_size = std::max(mpi->max_send_req_size, q_size(&mpi->pending_bcast_data_q));

	zlog_level(delta_log, ROUTER_V3, "Starting meta broadcast of size %d meta_data_index %d data_index %d\n", size, meta_data_index, data_index);

	start_meta_broadcast(mpi);
}

//void broadcast_static(int data_index, int size, mpi_context_t *mpi)
//{
	//mpi->max_send_data_size = std::max(mpi->max_send_data_size, size);

	//int *data = mpi->pending_send_data_nbc[data_index]->data();

	//int req_index = get_free_req(mpi, -1, data_index);

	//MPI_Iallgather(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, data, size, MPI_INT, mpi->comm, &mpi->pending_send_req[req_index]);
//}

void broadcast_rip_up_all_collective(int net_id, mpi_context_t *mpi)
{
	if (mpi->comm_size < 2) {
		return;
	}

	int data_index = get_free_buffer(mpi, 1); /* 1 byte for header */
	int *buffer = mpi->pending_send_data_nbc[data_index]->data();

	buffer[0] = encode_header(PacketID::RIP_UP_ALL, net_id);

	zlog_level(delta_log, ROUTER_V3, "Rip up all packet header %d pid %d net_id %d\n", buffer[0], static_cast<int>(PacketID::RIP_UP_ALL), net_id);

	broadcast_dynamic(data_index, 1, mpi);
}

template<typename ShouldExpandFunc>
void expand_neighbors_mpi_recv(const RRGraph &g, RRNode current, const route_state_t *state, congestion_t *congestion, const rr_node_property_t &target, float criticality_fac, float astar_fac, const ShouldExpandFunc &should_expand, std::priority_queue<route_state_t> &heap, perf_t *perf)
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

		float upstream_R = e_p.R + neighbor_p.R;
		if (!e_p.buffered) {
			upstream_R += current_state->upstream_R;
		}
		item.upstream_R = upstream_R;

		//float congestion_cost = get_congestion_cost_mpi_recv(item.rr_node, g, congestion, comm, this_pid, num_procs, pres_fac);
		extern t_rr_indexed_data *rr_indexed_data;
		float congestion_cost = rr_indexed_data[neighbor_p.cost_index].base_cost * congestion[neighbor].acc_cost * congestion[neighbor].pres_cost;

		float delay;
		if (e_p.buffered) {
			delay = e_p.switch_delay + neighbor_p.C * (e_p.R + 0.5 * neighbor_p.R);
		} else {
			delay = e_p.switch_delay + neighbor_p.C * (current_state->upstream_R + e_p.R + 0.5 * neighbor_p.R);
		}

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

RREdge get_previous_edge(int rr_node_id, const route_state_t *state, const route_tree_t &rt, const RRGraph &g);

template<typename NodeCallback>
void get_path(int sink_rr_node_id, const route_state_t *state, const route_tree_t &rt, const RRGraph &g, const NodeCallback &callback)
{
	int current_rr_node_id = sink_rr_node_id;
	RREdge previous_edge = get_previous_edge(current_rr_node_id, state, rt, g);

	while (valid(previous_edge)) {

		//char c[256];
		/* printing */
		//sprintf_rr_node(current_rr_node_id, c);
		//zlog_level(delta_log, ROUTER_V2, "Net %d get path: %s\n", vpr_net_id, c);
		//
		callback(previous_edge);

		current_rr_node_id = get_source(g, previous_edge);
		previous_edge = get_previous_edge(current_rr_node_id, state, rt, g);
	}

	//callback(current_rr_node_id, RRGraph::null_edge());
}

void route_net_mpi_send_recv_nonblocking_collective(const RRGraph &g, int vpr_id, int net_id, const source_t *source, const vector<sink_t *> &sinks, const t_router_opts *params, float pres_fac, const vector<net_t> &nets, const vector<vector<net_t *>> &partition_nets, int current_net_index, route_state_t *state, congestion_t *congestion, route_tree_t &rt, t_net_timing &net_timing, vector<vector<RRNode>> &net_route_trees, mpi_context_t *mpi, perf_t *perf, mpi_perf_t *mpi_perf, bool delayed_progress)
{
    using clock = myclock;

	std::priority_queue<route_state_t> heap;

	vector<int> modified;

	zlog_level(delta_log, ROUTER_V1, "Routing net %d local id %d (%lu sinks)\n", vpr_id, net_id, sinks.size());

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

		//update_one_cost_internal(source->rr_node, g, congestion, 1, pres_fac);
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

	//int sync_freq = (int)floor(sqrt(sorted_sinks.size()));
	
	//bool blocking_wait = false;
	//for (int i = 0; i < mpi->comm_size && !blocking_wait; ++i) {
		//if (i == mpi->rank) {
			//continue;
		//}

		//if (current_net_index < partition_nets[i].size() && current_net_index < partition_nets[mpi->rank].size()) {
			//bool overlap = box_overlap(partition_nets[i][current_net_index]->bounding_box, partition_nets[mpi->rank][current_net_index]->bounding_box);

			//zlog_level(delta_log, ROUTER_V3, "Net %d BB: %d-%d %d-%d Net %d BB: %d-%d %d-%d Overlap: %d\n",
					//partition_nets[i][current_net_index]->local_id, partition_nets[i][current_net_index]->bounding_box.xmin, partition_nets[i][current_net_index]->bounding_box.xmax, partition_nets[i][current_net_index]->bounding_box.ymin, partition_nets[i][current_net_index]->bounding_box.ymax,
					//partition_nets[mpi->rank][current_net_index]->local_id, partition_nets[mpi->rank][current_net_index]->bounding_box.xmin, partition_nets[mpi->rank][current_net_index]->bounding_box.xmax, partition_nets[mpi->rank][current_net_index]->bounding_box.ymin, partition_nets[mpi->rank][current_net_index]->bounding_box.ymax,
					//overlap ? 1 : 0
					//);

			//blocking_wait = blocking_wait || overlap;
		//}
	//}

	//zlog_level(delta_log, ROUTER_V3, "Net %d blocking wait %d\n", net_id, blocking_wait ? 1 : 0);

	int isink = 0;
	for (const auto &sink : sorted_sinks) {
		if (delayed_progress && isink > 0) {
			auto progress_start = clock::now();
			//if (q_size(&mpi->pending_bcast_data_q) > 16) { 
				//progress_broadcast_blocking(nets, congestion, g, pres_fac, net_route_trees, mpi);
			//} else {
				progress_broadcast(nets, congestion, g, pres_fac, net_route_trees, mpi);
			//}
			if (mpi_perf) {
				mpi_perf->total_send_testsome_time += clock::now()-progress_start;
			}
		}

		sprintf_rr_node(sink->rr_node, buffer);
		zlog_level(delta_log, ROUTER_V1, "Net %d current sink %d: %s criticality: %g BB: %d-%d %d-%d\n", vpr_id, sink->id, buffer, sink->criticality_fac, sink->current_bounding_box.xmin, sink->current_bounding_box.xmax, sink->current_bounding_box.ymin, sink->current_bounding_box.ymax);

		const auto &sink_rr_node = get_vertex_props(g, sink->rr_node);

		route_tree_multi_root_add_to_heap(rt, g, sink->rr_node, sink->criticality_fac, params->astar_fac, sink->current_bounding_box, heap, perf);

		int count = 0;

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

				expand_neighbors_mpi_recv(g, item.rr_node, state, congestion, sink_rr_node, sink->criticality_fac, params->astar_fac, [&g, &sink, &sink_rr_node, &v, &item] (const RRNode &n) -> bool {

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

			if (count % params->progress_freq == 0) {
				auto progress_start = clock::now();
				progress_broadcast(nets, congestion, g, pres_fac, net_route_trees, mpi);
				if (mpi_perf) {
					mpi_perf->total_send_testsome_time += clock::now()-progress_start;
				}
			}

			++count;
		}

		//bool sync_costs = (isink % sync_freq) == 0;
		//bool sync_costs = true;

		if (!found_sink) {
			assert(heap.empty());
			zlog_error(delta_log, "Error: Failed to find sink %d\n", sink->rr_node);

			if (sink->id != -1) {
				net_timing.delay[sink->id+1] = 0;
			}
		} else {
			int data_index = get_free_buffer(mpi, 128);

			sprintf_rr_node(sink->rr_node, buffer);
			zlog_level(delta_log, ROUTER_V3, "Adding %s to route tree\n", buffer);

			RouteTreeNode rt_node = route_tree_add_rr_node(rt, sink->rr_node, g);
			assert(rt_node != RouteTree::null_vertex());
			route_tree_set_node_properties(rt, rt_node, false, state[sink->rr_node].upstream_R, state[sink->rr_node].delay);
			update_one_cost_internal(sink->rr_node, g, congestion, 1, pres_fac);

			/* TODO: encode sink->rr_node into data_index */
			vector<int> &data = *mpi->pending_send_data_nbc[data_index];
			data[0] = encode_header(PacketID::ROUTE, net_id);
			data[1] = sink->rr_node;

			zlog_level(delta_log, ROUTER_V3, "Route packet header %d pid %d net_id %d\n", data[0], static_cast<int>(PacketID::ROUTE), net_id);

			RRNode child_rr_node = sink->rr_node;

			int num_nodes_added = 1;

			get_path(sink->rr_node, state, rt, g, [&] (const RREdge &edge) -> void
					{
					RRNode parent_rr_node = get_source(g, edge);

					const auto &parent_rr_node_p = get_vertex_props(g, parent_rr_node);

					bool bcast_node = false;

					sprintf_rr_node(parent_rr_node, buffer);
					zlog_level(delta_log, ROUTER_V3, "Adding %s to route tree\n", buffer);

					RouteTreeNode rt_node = route_tree_add_rr_node(rt, parent_rr_node, g);
					if (rt_node != RouteTree::null_vertex()) {
						route_tree_set_node_properties(rt, rt_node, parent_rr_node_p.type != IPIN && parent_rr_node_p.type != SINK, state[parent_rr_node].upstream_R, state[parent_rr_node].delay);

						update_one_cost_internal(parent_rr_node, g, congestion, 1, pres_fac);
						bcast_node = true;
					} else if (get_vertex_props(g, parent_rr_node).type == SOURCE) {
						update_one_cost_internal(parent_rr_node, g, congestion, 1, pres_fac);
						bcast_node = true;
					}

					assert(child_rr_node == get_target(g, edge));

					char buffer2[256];
					sprintf_rr_node(child_rr_node, buffer2);

					zlog_level(delta_log, ROUTER_V3, "Adding edge %s -> %s\n", buffer, buffer2);

					const RouteTreeEdge &rt_edge = route_tree_add_edge_between_rr_node(rt, parent_rr_node, child_rr_node); 
					auto &rt_edge_props = get_edge_props(rt.graph, rt_edge);
					rt_edge_props.rr_edge = edge;

					//if (parent_rr_node_p.type == OPIN && existing_opin == RRGraph::null_vertex()) {
						//existing_opin = parent_rr_node;
					//}

					child_rr_node = parent_rr_node;

					if (bcast_node) {
						if (data.size() <= 1+num_nodes_added) {
							data.resize(data.size()*2);
						}

						data[1+num_nodes_added] = parent_rr_node;

						++num_nodes_added;
					}
					});

			auto broadcast_start = clock::now();
			broadcast_dynamic(data_index, 1+num_nodes_added, mpi);
			if (mpi_perf) {
				mpi_perf->total_broadcast_time += clock::now()-broadcast_start;
			}

			--mpi->num_broadcasts_required;
			assert(mpi->num_broadcasts_required >= 0);

			if (sink->id != -1) {
				net_timing.delay[sink->id+1] = get_vertex_props(rt.graph, route_tree_get_rt_node(rt, sink->rr_node)).delay;
			}

			zlog_level(delta_log, ROUTER_V2, "Routed sink %d\n", sink->rr_node);
		}

		if (!delayed_progress) {
			auto progress_start = clock::now();
			//if (q_size(&mpi->pending_bcast_data_q) > 256) { 
				//progress_broadcast_blocking(nets, congestion, g, pres_fac, net_route_trees, mpi);
			//} else {
				progress_broadcast(nets, congestion, g, pres_fac, net_route_trees, mpi);
			//}
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

	zlog_level(delta_log, ROUTER_V2, "Routed net %d\n", vpr_id);
	//zlog_level(delta_log, ROUTER_V1, "\n");

	//return unrouted_sinks;

	/*delete [] state;*/

	/*check_route_tree(rt, net, g);*/
}

