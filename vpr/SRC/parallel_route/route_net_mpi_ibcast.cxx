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
#include "route_net_mpi_ibcast.h"
#include "clock.h"

float get_timing_driven_expected_cost(const rr_node_property_t &current, const rr_node_property_t &target, float criticality_fac, float R_upstream);

void read_header(unsigned int *packet, IbcastPacketID &packet_id, int &payload_size, int &net_id)
{
	packet_id = static_cast<IbcastPacketID>(packet[0] >> 28);
	payload_size = packet[0] & 0xFFFFFFF;
	net_id = packet[1];
}

void write_header(int *packet, IbcastPacketID packet_id, unsigned int payload_size, int net_id)
{
	assert(sizeof(unsigned int) == 4);
	assert((payload_size & 0xFFFFFFF) == payload_size);

	packet[0] = (static_cast<unsigned int>(packet_id) << 28) | payload_size;
	packet[1] = net_id;
}

int get_ibcast_buffer(mpi_context_t *mpi)
{
	int data_index;

	int size = mpi->max_ibcast_count[mpi->rank];

	if (!mpi->free_send_data_index.empty()) {
		data_index = mpi->free_send_data_index.front();
		mpi->free_send_data_index.pop();

		zlog_level(delta_log, ROUTER_V3, "Existing data, data_index %d\n", data_index);

		if (mpi->pending_send_data_nbc[data_index]->size() < size) {
			mpi->pending_send_data_nbc[data_index]->resize(2*size);
		}

		assert(mpi->pending_send_data_ref_count[data_index] == 0);
	} else {
		data_index = mpi->pending_send_data_nbc.size();

		zlog_level(delta_log, ROUTER_V3, "New data, data_index %d\n", data_index);

		assert(mpi->pending_send_data_nbc.size() == mpi->pending_send_data_ref_count.size());

		mpi->pending_send_data_nbc.emplace_back(new vector<int>(size));

		mpi->pending_send_data_ref_count.emplace_back(0);
	}

	return data_index;
}

static int get_free_req(mpi_context_t *mpi, int data_index, int offset, int count)
{
	int req_index;

	if (!mpi->free_send_req_index.empty()) {
		req_index = mpi->free_send_req_index.front();
		mpi->free_send_req_index.pop();

		mpi->pending_req_meta[req_index].data_index = data_index;
		mpi->pending_req_meta[req_index].data_offset = offset;
		mpi->pending_req_meta[req_index].data_count = count;

		++mpi->pending_send_data_ref_count[data_index];

		zlog_level(delta_log, ROUTER_V3, "Existing req, req_index %d data_index %d\n", req_index, data_index);
	} else {
		assert(mpi->pending_req_meta.size() == mpi->pending_send_req.size());

		req_index = mpi->pending_send_req.size();
		mpi->pending_send_req.emplace_back(MPI_REQUEST_NULL);

		mpi->pending_req_meta.emplace_back(request_meta_t{ data_index, offset, count });

		++mpi->pending_send_data_ref_count[data_index];

		zlog_level(delta_log, ROUTER_V3, "New req, req_index %d data_index %d\n", req_index, data_index);
	}

	return req_index;
}

static void broadcast_as_root(int data_index, int offset, int count, mpi_context_t *mpi)
{
	int req_index = get_free_req(mpi, data_index, offset, count);

	zlog_level(delta_log, ROUTER_V3, "Starting broadcast as root, req_index %d data_index %d offset %d\n", req_index, data_index, offset);

	mpi->max_send_req_size = std::max((unsigned long)mpi->max_send_req_size, mpi->pending_send_req.size());

	assert(mpi->pending_send_req[req_index] == MPI_REQUEST_NULL);
	assert(mpi->pending_req_meta[req_index].data_index == data_index);
	assert(mpi->pending_req_meta[req_index].data_offset == offset);
	assert(mpi->pending_req_meta[req_index].data_count == count);
	assert(mpi->pending_send_data_nbc[data_index]->size() >= offset+count);

	int *data = mpi->pending_send_data_nbc[data_index]->data();

	int error = MPI_Ibcast(&data[offset], count, MPI_INT, mpi->rank, mpi->ibcast_comm[mpi->rank], &mpi->pending_send_req[req_index]);

	assert(error == MPI_SUCCESS);

	++mpi->num_pending_reqs;
}

void broadcast_as_root_large(int data_index, int count, mpi_context_t *mpi)
{
	assert(mpi->comm_size > 1);

	zlog_level(delta_log, ROUTER_V3, "Starting large broadcast as root, data_index %d count %d\n", data_index, count);

	broadcast_as_root(data_index, 0, mpi->max_ibcast_count[mpi->rank], mpi);
	int remaining_count = count-mpi->max_ibcast_count[mpi->rank];
	if (remaining_count > 0) {
		broadcast_as_root(data_index, mpi->max_ibcast_count[mpi->rank], remaining_count, mpi);
	}
}

static void broadcast_as_non_root(int data_index, int offset, int count, int rank, mpi_context_t *mpi)
{
	assert(mpi->comm_size > 1);

	assert(rank < mpi->comm_size && rank != mpi->rank);

	int req_index = rank;

	zlog_level(delta_log, ROUTER_V3, "Starting broadcast as non-root, req_index %d data_index %d offset %d count %d rank %d\n", req_index, data_index, offset, count, rank);

	mpi->max_send_req_size = std::max((unsigned long)mpi->max_send_req_size, mpi->pending_send_req.size());

	assert(mpi->pending_send_req[req_index] == MPI_REQUEST_NULL);
	mpi->pending_req_meta[req_index].data_index = data_index;
	mpi->pending_req_meta[req_index].data_offset = offset;
	mpi->pending_req_meta[req_index].data_count = count;
	assert(mpi->pending_send_data_nbc[data_index]->size() >= offset+count);

	int *data = mpi->pending_send_data_nbc[data_index]->data();

	int error = MPI_Ibcast(&data[offset], count, MPI_INT, rank, mpi->ibcast_comm[rank], &mpi->pending_send_req[req_index]);

	assert(error == MPI_SUCCESS);

	++mpi->num_pending_reqs;
}

/* bcast_d cannot be a ref because start_actual_broadcast modifies active_bcast_data */
static void handle_packet(int completed_index, const vector<net_t> &nets, congestion_t *congestion, const RRGraph &g, float pres_fac, vector<vector<RRNode>> &net_route_trees, mpi_context_t *mpi)
{
	assert(completed_index < mpi->comm_size && completed_index != mpi->rank);

	int data_index = mpi->pending_req_meta[completed_index].data_index;
	int *data = mpi->pending_send_data_nbc[data_index]->data();

	zlog_level(delta_log, ROUTER_V3, "Received broadcast, req %d data_index %d\n", completed_index, data_index);

	IbcastPacketID pid;
	int count;
	int net_id;
	read_header((unsigned int *)data, pid, count, net_id);

	int received_count = mpi->pending_req_meta[completed_index].data_offset + mpi->pending_req_meta[completed_index].data_count;

	if (received_count < count + 2) {
		/* more fragments to receive */
		zlog_level(delta_log, ROUTER_V3, "Received only %d out of %d\n", received_count, count+2);

		if (mpi->pending_send_data_nbc[data_index]->size() < count+2) {
			mpi->pending_send_data_nbc[data_index]->resize((count+2)*2);
		}

		broadcast_as_non_root(data_index, received_count, count+2-received_count, completed_index, mpi);
		return;
	}

	zlog_level(delta_log, ROUTER_V3, "Header %d pid %d count %d net_id %d\n", data[0], static_cast<int>(pid), count, net_id);

	assert(static_cast<int>(pid) >= 0 && pid < IbcastPacketID::NUM_PACKET_IDS);
	assert(net_id >= 0 && net_id < nets.size());

	int delta; 
	int *changed_rr_nodes = &data[2];

	switch (pid) {
		case IbcastPacketID::ROUTE:
		case IbcastPacketID::RIP_UP:
			delta = pid == IbcastPacketID::ROUTE ? 1 : -1;

			if (delta > 0) {
				zlog_level(delta_log, ROUTER_V3, "Rank %d net %d route\n", completed_index, nets[net_id].vpr_id);
			} else {
				assert(delta < 0);
				zlog_level(delta_log, ROUTER_V3, "Rank %d net %d rip up\n", completed_index, nets[net_id].vpr_id);
			}

			for (int i = 0; i < count; ++i) {
				update_one_cost_internal(changed_rr_nodes[i], g, congestion, delta, pres_fac);
				//zlog_level(delta_log, ROUTER_V3, "MPI update, source %d node %d delta %d\n", status.MPI_SOURCE, d[i].rr_node, d[i].delta);

				if (delta > 0) {
					//assert(find(begin(net_route_trees[net_id]), end(net_route_trees[net_id]), changed_rr_nodes[i]) == end(net_route_trees[net_id]));

					net_route_trees[net_id].push_back(changed_rr_nodes[i]);
				} else {
					assert(delta < 0);
					//assert(find(begin(net_route_trees[net_id]), end(net_route_trees[net_id]), changed_rr_nodes[i]) != end(net_route_trees[net_id]));

					net_route_trees[net_id].erase(std::remove(begin(net_route_trees[net_id]), end(net_route_trees[net_id]), changed_rr_nodes[i]), end(net_route_trees[net_id]));
				}
			}

			break;

		case IbcastPacketID::RIP_UP_ALL:
			zlog_level(delta_log, ROUTER_V3, "Rank %d rip up entire up net %d of size %lu\n", completed_index, nets[net_id].vpr_id, net_route_trees[net_id].size());

			assert(count == 0);

			assert(!net_route_trees[net_id].empty());

			for (const auto &node : net_route_trees[net_id]) {
				update_one_cost_internal(node, g, congestion, -1, pres_fac);
				//zlog_level(delta_log, ROUTER_V3, "MPI rip up, source %d node %d delta -1\n", status.MPI_SOURCE, node);
				//assert(find(begin(net_route_trees[net_id]), end(net_route_trees[net_id]), node) != end(net_route_trees[net_id]));
			}

			net_route_trees[net_id].clear();

			break;

		case IbcastPacketID::TRAILER:
			zlog_level(delta_log, ROUTER_V3, "Rank %d trailer\n", completed_index);

			assert(count == 0);
			assert(net_id == 0);

			mpi->received_last_update[completed_index] = true;

			break;

		default:
			zlog_level(delta_log, ROUTER_V3, "Unknown packet ID %d from rank %d\n", static_cast<int>(pid), completed_index);

			assert(false);
			break;
	}

	if (!mpi->received_last_update[completed_index]) {
		broadcast_as_non_root(data_index, 0, mpi->max_ibcast_count[completed_index], completed_index, mpi);
	} else {
		--mpi->pending_send_data_ref_count[completed_index];
		assert(mpi->pending_send_data_ref_count[completed_index] == 0);
	}
}

void handle_completed_request(int completed_index, const vector<net_t> &nets, congestion_t *congestion, const RRGraph &g, float pres_fac, vector<vector<RRNode>> &net_route_trees, mpi_context_t *mpi)
{
	assert(mpi->pending_send_req[completed_index] == MPI_REQUEST_NULL);

	if (completed_index < mpi->comm_size) {
		handle_packet(completed_index, nets, congestion, g, pres_fac, net_route_trees, mpi);
	} else {
		zlog_level(delta_log, ROUTER_V3, "Completed broadcast, req index %d\n", completed_index);

		/* i am root */
		int data_index = mpi->pending_req_meta[completed_index].data_index;

		--mpi->pending_send_data_ref_count[data_index];

		zlog_level(delta_log, ROUTER_V3, "data_index %d new ref count %d\n", data_index, mpi->pending_send_data_ref_count[data_index]);

		assert(mpi->pending_send_data_ref_count[data_index] >= 0);

		if (mpi->pending_send_data_ref_count[data_index] == 0) {
			zlog_level(delta_log, ROUTER_V3, "Freeing data_index %d\n", data_index);
			mpi->free_send_data_index.push(data_index);
		}

		zlog_level(delta_log, ROUTER_V3, "Freeing req_index %d\n", completed_index);
		mpi->free_send_req_index.push(completed_index);
	}

	--mpi->num_pending_reqs;
	assert(mpi->num_pending_reqs >= 0);
}

void progress_ibcast(const vector<net_t> &nets, congestion_t *congestion, const RRGraph &g, float pres_fac, vector<vector<RRNode>> &net_route_trees, mpi_context_t *mpi)
{
	if (mpi->num_pending_reqs == 0) {
		return;
	}

	assert(!mpi->pending_send_req.empty());

	zlog_level(delta_log, ROUTER_V3, "Progressing broadcast, num reqs %d\n", mpi->pending_send_req.size());

	if (mpi->completed_indices.size() < mpi->pending_send_req.size()) {
		mpi->completed_indices.resize(mpi->pending_send_req.size()*2);
	}

	int num_completed; 
	int error = MPI_Testsome(mpi->pending_send_req.size(), mpi->pending_send_req.data(), 
			&num_completed, mpi->completed_indices.data(), MPI_STATUSES_IGNORE);

	assert(error == MPI_SUCCESS);
	assert(num_completed != MPI_UNDEFINED);

	for (int i = 0; i < num_completed; ++i) {
		int completed_index = mpi->completed_indices[i];

		handle_completed_request(completed_index, nets, congestion, g, pres_fac, net_route_trees, mpi);
	}
}

bool all_done(mpi_context_t *mpi)
{
	bool all_received_last_update = true;

	for (int pid = 0; pid < mpi->comm_size && all_received_last_update; ++pid) {
		if (pid != mpi->rank && !mpi->received_last_update[pid]) {
			all_received_last_update = false;
		}
	}

	return all_received_last_update;
}

void progress_ibcast_blocking(const vector<net_t> &nets, congestion_t *congestion, const RRGraph &g, float pres_fac, vector<vector<RRNode>> &net_route_trees, mpi_context_t *mpi)
{
	if (mpi->num_pending_reqs == 0) {
		return;
	}

	do {
		vector<MPI_Request> prev_reqs = mpi->pending_send_req;

		zlog_level(delta_log, ROUTER_V3, "Progressing broadcast blocking, num reqs %d\n", mpi->pending_send_req.size());

		int error = MPI_Waitall(mpi->pending_send_req.size(), mpi->pending_send_req.data(), MPI_STATUSES_IGNORE);

		assert(error == MPI_SUCCESS);

		for (int i = 0; i < mpi->pending_send_req.size(); ++i) {
			assert(mpi->pending_send_req[i] == MPI_REQUEST_NULL);

			if (mpi->pending_send_req[i] == prev_reqs[i]) {
				continue;
			}

			assert(prev_reqs[i] != MPI_REQUEST_NULL);

			handle_completed_request(i, nets, congestion, g, pres_fac, net_route_trees, mpi);
		}
	} while (!all_done(mpi));
}

void bootstrap_ibcast(mpi_context_t *mpi)
{
	zlog_level(delta_log, ROUTER_V3, "Bootstrapping\n");

	for (int i = 0; i < mpi->comm_size; ++i) {
		if (i != mpi->rank) {
			++mpi->pending_send_data_ref_count[i];

			broadcast_as_non_root(i, 0, mpi->max_ibcast_count[i], i, mpi);
		} 
	}
}

//void broadcast_static(int data_index, int size, mpi_context_t *mpi)
//{
	//mpi->max_send_data_size = std::max(mpi->max_send_data_size, size);

	//int *data = mpi->pending_send_data_nbc[data_index]->data();

	//int req_index = get_free_req(mpi, -1, data_index);

	//MPI_Iallgather(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, data, size, MPI_INT, mpi->comm, &mpi->pending_send_req[req_index]);
//}

template<typename ShouldExpandFunc>
void expand_neighbors(const RRGraph &g, RRNode current, const route_state_t *state, congestion_t *congestion, const rr_node_property_t &target, float criticality_fac, float astar_fac, const ShouldExpandFunc &should_expand, std::priority_queue<route_state_t> &heap, perf_t *perf)
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

void route_net_mpi_ibcast(const RRGraph &g, int vpr_id, int net_id, const source_t *source, const vector<sink_t *> &sinks, const t_router_opts *params, float pres_fac, const vector<net_t> &nets, const vector<vector<net_t *>> &partition_nets, int current_net_index, route_state_t *state, congestion_t *congestion, route_tree_t &rt, t_net_timing &net_timing, vector<vector<RRNode>> &net_route_trees, mpi_context_t *mpi, perf_t *perf, mpi_perf_t *mpi_perf, bool delayed_progress)
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
		if (delayed_progress && isink > 0 && mpi->comm_size > 1) {
			auto progress_start = clock::now();
			//if (q_size(&mpi->pending_bcast_data_q) > 16) { 
				//progress_ibcast_blocking(nets, congestion, g, pres_fac, net_route_trees, mpi);
			//} else {
				progress_ibcast(nets, congestion, g, pres_fac, net_route_trees, mpi);
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

				expand_neighbors(g, item.rr_node, state, congestion, sink_rr_node, sink->criticality_fac, params->astar_fac, [&g, &sink, &sink_rr_node, &v, &item] (const RRNode &n) -> bool {

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
				progress_ibcast(nets, congestion, g, pres_fac, net_route_trees, mpi);
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
			sprintf_rr_node(sink->rr_node, buffer);
			zlog_level(delta_log, ROUTER_V3, "Adding %s to route tree\n", buffer);

			RouteTreeNode rt_node = route_tree_add_rr_node(rt, sink->rr_node, g);
			assert(rt_node != RouteTree::null_vertex());
			route_tree_set_node_properties(rt, rt_node, false, state[sink->rr_node].upstream_R, state[sink->rr_node].delay);
			update_one_cost_internal(sink->rr_node, g, congestion, 1, pres_fac);

			/* TODO: encode sink->rr_node into data_index */
			int data_index = get_ibcast_buffer(mpi);
			vector<int> &data = *mpi->pending_send_data_nbc[data_index];
			data[2] = sink->rr_node;

			zlog_level(delta_log, ROUTER_V3, "Route packet header %d pid %d net_id %d\n", data[0], static_cast<int>(IbcastPacketID::ROUTE), net_id);

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
						if (data.size() <= 2+num_nodes_added) {
							data.resize((2+num_nodes_added)*2);
						}

						data[2+num_nodes_added] = parent_rr_node;

						++num_nodes_added;
					}
					});

			write_header(&data[0], IbcastPacketID::ROUTE, num_nodes_added, net_id);

			if (mpi->comm_size > 1) {
				auto broadcast_start = clock::now();

				broadcast_as_root_large(data_index, 2+num_nodes_added, mpi);

				if (mpi_perf) {
					mpi_perf->total_broadcast_time += clock::now()-broadcast_start;
				}
			} else {
				assert(mpi->pending_send_data_ref_count[data_index] == 0);
				mpi->free_send_data_index.push(data_index);
			}

			if (sink->id != -1) {
				net_timing.delay[sink->id+1] = get_vertex_props(rt.graph, route_tree_get_rt_node(rt, sink->rr_node)).delay;
			}

			zlog_level(delta_log, ROUTER_V2, "Routed sink %d\n", sink->rr_node);
		}

		if (!delayed_progress && mpi->comm_size > 1) {
			auto progress_start = clock::now();
			//if (q_size(&mpi->pending_bcast_data_q) > 256) { 
				//progress_ibcast_blocking(nets, congestion, g, pres_fac, net_route_trees, mpi);
			//} else {
				progress_ibcast(nets, congestion, g, pres_fac, net_route_trees, mpi);
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

