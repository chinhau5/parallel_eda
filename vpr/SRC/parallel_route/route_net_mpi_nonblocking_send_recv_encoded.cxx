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
#include "route_net_mpi_nonblocking_send_recv_encoded.h"
#include "clock.h"
#include "path_codec.h"

float get_timing_driven_expected_cost(const rr_node_property_t &current, const rr_node_property_t &target, float criticality_fac, float R_upstream);

void read_header_sr(unsigned int *packet, NBSRPacketID &packet_id, int &net_id)
{
	packet_id = static_cast<NBSRPacketID>(packet[0] & 0xF);
	net_id = packet[0] >> 4;
}

void write_header_sr(int *packet, NBSRPacketID packet_id, unsigned int net_id)
{
	assert(sizeof(unsigned int) == 4);
	assert((net_id & 0xFFFFFFF) == net_id);

	packet[0] = (net_id << 4) | static_cast<unsigned int>(packet_id);
}

int get_free_buffer(mpi_context_t *mpi, int size)
{
	int data_index;

	if (!mpi->free_send_data_index.empty()) {
		data_index = mpi->free_send_data_index.front();

		mpi->free_send_data_index.pop();

		assert(mpi->pending_send_data_ref_count[data_index] == 0);

		if (mpi->pending_send_data_nbc[data_index]->size() < size) {
			mpi->pending_send_data_nbc[data_index]->resize(2*size);
		}

		zlog_level(delta_log, ROUTER_V3, "Existing data index %d\n", data_index);
	} else {
		data_index = mpi->pending_send_data_nbc.size();

		mpi->pending_send_data_nbc.emplace_back(new vector<int>(size));
		mpi->pending_send_data_ref_count.emplace_back(0);

		assert(mpi->pending_send_data_ref_count.size() == mpi->pending_send_data_nbc.size());

		zlog_level(delta_log, ROUTER_V3, "New data index %d\n", data_index);
	}

	return data_index;
}

static int get_free_req(mpi_context_t *mpi, int data_index)
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
		mpi->pending_send_req_dst.emplace_back(-1);

		assert(mpi->pending_send_req.size() == mpi->pending_send_req_data_ref.size());
		assert(mpi->pending_send_req.size() == mpi->pending_send_req_dst.size());

		zlog_level(delta_log, ROUTER_V3, "New req, req_index %d data_index %d\n", req_index, data_index);
	}

	return req_index;
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

void bootstrap_irecv_combined(mpi_context_t *mpi)
{
	for (int pid = 0; pid < mpi->comm_size; ++pid) {
		if (pid != mpi->rank) {
			assert(!mpi->received_last_update[pid]);
			assert(mpi->pending_send_req[pid] == MPI_REQUEST_NULL);

#ifdef USE_FIXED_TAG
			int error = MPI_Irecv(mpi->pending_send_data_nbc[pid]->data(), mpi->pending_send_data_nbc[pid]->size(), MPI_INT, pid, FIXED_TAG, mpi->comm, &mpi->pending_send_req[pid]);
#else
			int error = MPI_Irecv(mpi->pending_send_data_nbc[pid]->data(), mpi->pending_send_data_nbc[pid]->size(), MPI_INT, pid, MPI_ANY_TAG, mpi->comm, &mpi->pending_send_req[pid]);
#endif

			++mpi->num_pending_reqs;

			assert(error == MPI_SUCCESS);
		}
	}
}

void unicast_nonblocking(int data_index, int size, int dst, int tag, mpi_context_t *mpi)
{
	int req_index = get_free_req(mpi, data_index);
	int error;

	assert(req_index >= mpi->comm_size);
	assert(mpi->pending_send_req[req_index] == MPI_REQUEST_NULL);

	mpi->pending_send_req_dst[req_index] = dst;

	if (data_index < 0) {
		error = MPI_Isend(nullptr, 0, MPI_INT, dst, tag, mpi->comm, &mpi->pending_send_req[req_index]);
	} else {
		int *data = mpi->pending_send_data_nbc[data_index]->data();

		error = MPI_Isend(data, size, MPI_INT, dst, tag, mpi->comm, &mpi->pending_send_req[req_index]);

		assert(mpi->pending_send_data_ref_count[data_index] >= 0);
		++mpi->pending_send_data_ref_count[data_index];

		mpi->total_send_data_size += size;
		++mpi->total_send_count;
	}

	assert(error == MPI_SUCCESS);

	++mpi->num_pending_reqs;
	++mpi->num_pending_reqs_by_rank[dst];
}

void broadcast_nonblocking(int data_index, int size, int tag, mpi_context_t *mpi)
{
	assert(size <= mpi->buffer_size);

	mpi->max_send_data_size = std::max(mpi->max_send_data_size, size);

	for (int i = 0; i < mpi->comm_size; ++i) {
		if (i != mpi->rank) {
			zlog_level(delta_log, ROUTER_V2, "Sent packet from %d to %d\n", mpi->rank, i);
			
			unicast_nonblocking(data_index, size, i, tag, mpi);
		}
	}
}

static bool any_done(mpi_context_t *mpi)
{
	bool any_received_last_update = false;

	for (int pid = 0; pid < mpi->comm_size && !any_received_last_update; ++pid) {
		if (pid != mpi->rank) {
			any_received_last_update = mpi->received_last_update[pid];
		}
	}

	return any_received_last_update;
}

void selective_broadcast(int data_index, int packet_size, int tag, mpi_context_t *mpi)
{
	assert(packet_size <= mpi->buffer_size);

	if (any_done(mpi)) {
		zlog_level(delta_log, ROUTER_V3, "Selective broadcast\n");

		int num_small_unicasts = 0;

		for (int i = 0; i < mpi->comm_size; ++i) {
			if (i == mpi->rank) {
				continue;
			}

			if (!mpi->received_last_update[i]) {
				zlog_level(delta_log, ROUTER_V3, "Small unicast, index %d size %d rank %d\n", data_index, packet_size, i);

				unicast_nonblocking(data_index, packet_size, i, tag, mpi);
				++num_small_unicasts;
			} else {
				int *small_data = mpi->pending_send_data_nbc[data_index]->data();

				large_t *large = &mpi->large_packets[i];

				if (large->data_index == -1) {
					large->data_index = get_free_buffer(mpi, mpi->buffer_size);
					large->data_offset = 0;

					zlog_level(delta_log, ROUTER_V3, "First large packet, index %d rank %d\n", large->data_index, i);
				}
				
				if (large->data_offset + packet_size < mpi->buffer_size) {
					zlog_level(delta_log, ROUTER_V3, "Buffering data, index %d offset %d size %d rank %d\n", large->data_index, large->data_offset, packet_size, i);

					int *large_data = mpi->pending_send_data_nbc[large->data_index]->data();

					memcpy(&large_data[large->data_offset], small_data, packet_size*sizeof(int));

					large->data_offset += packet_size;
				} else if (large->data_offset + packet_size == mpi->buffer_size) {
					zlog_level(delta_log, ROUTER_V3, "Flushing full data, index %d size %d rank %d\n", large->data_index, mpi->buffer_size, i);

					int *large_data = mpi->pending_send_data_nbc[large->data_index]->data();

					memcpy(&large_data[large->data_offset], small_data, packet_size*sizeof(int));

					unicast_nonblocking(large->data_index, mpi->buffer_size, i, tag, mpi);

					large->data_index = -1;
					large->data_offset = -1;
				} else {
					zlog_level(delta_log, ROUTER_V3, "Flushing data, index %d size %d rank %d\n", large->data_index, large->data_offset, i);

					unicast_nonblocking(large->data_index, large->data_offset, i, tag, mpi);

					large->data_index = get_free_buffer(mpi, mpi->buffer_size);

					int *large_data = mpi->pending_send_data_nbc[large->data_index]->data();

					memcpy(large_data, small_data, packet_size*sizeof(int));

					large->data_offset = packet_size;

					zlog_level(delta_log, ROUTER_V3, "Rebuffering data, index %d size %d rank %d\n", large->data_index, packet_size, i);
				}
			}
		}

		if (num_small_unicasts == 0) {
			assert(mpi->pending_send_data_ref_count[data_index] == 0);
			mpi->free_send_data_index.push(data_index);
		}
	} else {
		zlog_level(delta_log, ROUTER_V3, "Normal broadcast\n");
		broadcast_nonblocking(data_index, packet_size, tag, mpi);
	}
}

void broadcast_trailer_nonblocking(mpi_context_t *mpi)
{
	zlog_level(delta_log, ROUTER_V2, "Sending trailer packet\n");

#ifdef USE_FIXED_TAG
	int data_index = get_free_buffer(mpi, 1);
	int *data = mpi->pending_send_data_nbc[data_index]->data();

	write_header_sr(data, NBSRPacketID::TRAILER, 0);

	selective_broadcast(data_index, 1, FIXED_TAG, mpi);

	/* flush */
	for (int i = 0; i < mpi->comm_size; ++i) {
		if (i == mpi->rank) {
			continue;
		}

		large_t *large = &mpi->large_packets[i];

		if (large->data_index != -1) {
			zlog_level(delta_log, ROUTER_V3, "Flushing last data, index %d size %d rank %d\n", large->data_index, large->data_offset, i);

			int *large_data = mpi->pending_send_data_nbc[large->data_index]->data();

			unicast_nonblocking(large->data_index, large->data_offset, i, FIXED_TAG, mpi);

			large->data_index = -1;
			large->data_offset = -1;
		}
	}
#else
	selective_broadcast(-1, 0, LAST_TAG, mpi);
#endif
}

void broadcast_rip_up_nonblocking(int net_id, mpi_context_t *mpi)
{
	zlog_level(delta_log, ROUTER_V3, "MPI sent rip up");

#ifdef USE_FIXED_TAG
	int data_index = get_free_buffer(mpi, 1);
	int *data = mpi->pending_send_data_nbc[data_index]->data();

	write_header_sr(data, NBSRPacketID::RIP_UP, net_id);

	selective_broadcast(data_index, 1, FIXED_TAG, mpi);
#else
	selective_broadcast(-1, 0, RIP_UP_TAG + net_id, mpi);
#endif
}

void complete_send_request(mpi_context_t *mpi, int completed_index)
{
	assert(mpi->pending_send_req[completed_index] == MPI_REQUEST_NULL);

	int data_index = mpi->pending_send_req_data_ref[completed_index];

	if (data_index >= 0) {
		zlog_level(delta_log, ROUTER_V3, "Reducing ref count, req_index %d data_index %d new_ref_count %d\n", completed_index, data_index, mpi->pending_send_data_ref_count[data_index]);

		--mpi->pending_send_data_ref_count[data_index];

		assert(mpi->pending_send_data_ref_count[data_index] >= 0);

		if (mpi->pending_send_data_ref_count[data_index] == 0) {
			zlog_level(delta_log, ROUTER_V3, "Freeing data, req_index %d data_index %d\n", completed_index, data_index);

			mpi->free_send_data_index.push(data_index);
		}
	}

	zlog_level(delta_log, ROUTER_V3, "Freeing req, req_index %d\n", completed_index);

	mpi->free_send_req_index.push(completed_index);
}

void progress_sends_waitall(mpi_context_t *mpi)
{
	assert(!mpi->pending_send_req.empty());

	mpi->max_req_buffer_size = std::max((unsigned long)mpi->max_req_buffer_size, mpi->pending_send_req.size());

	vector<MPI_Status> statuses(mpi->pending_send_req.size());

	int error = MPI_Waitall(mpi->pending_send_req.size(), mpi->pending_send_req.data(), statuses.data());

	assert(error == MPI_SUCCESS);

	for (int i = 0; i < statuses.size(); ++i) {
		if (statuses[i].MPI_TAG == MPI_ANY_TAG && statuses[i].MPI_SOURCE == MPI_ANY_SOURCE && statuses[i].MPI_ERROR == MPI_SUCCESS) {
			/* empty status */
		} else {
			complete_send_request(mpi, i);
		}
	}
}

void progress_sends_testsome(mpi_context_t *mpi)
{
	assert(!mpi->pending_send_req.empty());

	mpi->max_req_buffer_size = std::max((unsigned long)mpi->max_req_buffer_size, mpi->pending_send_req.size());

	if (mpi->completed_send_indices.size() < mpi->pending_send_req.size()) {
		mpi->completed_send_indices.resize(mpi->pending_send_req.size());
	}

	int num_completed;

	int error = MPI_Testsome(mpi->pending_send_req.size(), mpi->pending_send_req.data(),
					&num_completed, mpi->completed_send_indices.data(), MPI_STATUSES_IGNORE);

	assert(error == MPI_SUCCESS);
	assert(num_completed != MPI_UNDEFINED);

	for (int i = 0; i < num_completed; ++i) {
		complete_send_request(mpi, mpi->completed_send_indices[i]);
	}
}

static void handle_packet_fixed_tag(const int *data, const MPI_Status &status, const vector<net_t> &nets, congestion_t *congestion, const RRGraph &g, float pres_fac, vector<vector<RRNode>> &net_route_trees, mpi_context_t *mpi, mpi_perf_t *mpi_perf)
{
	//auto update_start = clock::now();

	int num_recvd;
	int error = MPI_Get_count(&status, MPI_INT, &num_recvd);

	assert(error == MPI_SUCCESS);
	assert(status.MPI_TAG == FIXED_TAG);

	//if (mpi_perf) {
		//mpi_perf->total_update_time += clock::now()-update_start;
	//}
	
	int offset = 0;

	while (offset < num_recvd) {
		const int *current_data = &data[offset];

		NBSRPacketID packet_id;
		int net_id;
		read_header_sr((unsigned int *)current_data, packet_id, net_id);

		switch (packet_id) {
			case NBSRPacketID::TRAILER:
				assert(net_id == 0);

				zlog_level(delta_log, ROUTER_V3, "Trailer packet from %d\n", status.MPI_SOURCE);

				assert(!mpi->received_last_update[status.MPI_SOURCE]);
				mpi->received_last_update[status.MPI_SOURCE] = true;

				++offset;

				//if (mpi_perf) {
				//mpi_perf->total_update_time += clock::now()-update_start;
				//}
				break;

			case NBSRPacketID::ROUTE:
				//#define RIP_UP_TAG LAST_TAG+1
				//#define COST_UPDATE_TAG RIP_UP_TAG+0x1000000 //support for up to 16M nets
				//
				{
					assert(net_id >= 0 && net_id < net_route_trees.size());
					assert(nets[net_id].local_id == net_id);

					zlog_level(delta_log, ROUTER_V3, "Cost update of net %d from %d\n", nets[net_id].vpr_id, status.MPI_SOURCE);

					//auto update_start = clock::now();

					//update_start = clock::now();

					int rr_node = current_data[1];
					unsigned int val;
					path_decoder_t dec;
					decoder_init(dec, (unsigned int *)&current_data[2], mpi->bit_width);
					do {
						update_one_cost_internal(rr_node, g, congestion, 1, pres_fac);
						//zlog_level(delta_log, ROUTER_V3, "MPI update, source %d node %d delta %d\n", status.MPI_SOURCE, rr_node, delta);

						//if (delta > 0) {
						//assert(find(begin(net_route_trees[net_id]), end(net_route_trees[net_id]), rr_node) == end(net_route_trees[net_id]));

						net_route_trees[net_id].push_back(rr_node);
						//} else {
						//assert(delta < 0);
						//assert(find(begin(net_route_trees[net_id]), end(net_route_trees[net_id]), rr_node) != end(net_route_trees[net_id]));
						//net_route_trees[net_id].erase(std::remove(begin(net_route_trees[net_id]), end(net_route_trees[net_id]), rr_node), end(net_route_trees[net_id]));
						//}

						val = decoder_read(dec);
						zlog_level(delta_log, ROUTER_V3, "switch id: %u\n", val);
						if (val != dec.bit_mask) {
							rr_node = get_target(g, get_edge_by_index(g, rr_node, val));
						}
					} while (val != dec.bit_mask);

					zlog_level(delta_log, ROUTER_V3, "address %X num recv %d word offset: %d bit offset %d\n", current_data, num_recvd, dec.word_offset, dec.bit_offset);

					offset += 1+1+decoder_get_num_words(dec);

					//if (mpi_perf) {
					//mpi_perf->total_update_time += clock::now()-update_start;
					//}
					break;
				} 

			case NBSRPacketID::RIP_UP:
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
				
				++offset;

				//if (mpi_perf) {
				//mpi_perf->total_update_time += clock::now()-update_start;
				//}
				break;

			default:
				assert(false);
				break;
		}
	}

	assert(offset == num_recvd);
}

static void handle_packet(const int *data, const MPI_Status &status, const vector<net_t> &nets, congestion_t *congestion, const RRGraph &g, float pres_fac, vector<vector<RRNode>> &net_route_trees, mpi_context_t *mpi, mpi_perf_t *mpi_perf)
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

		zlog_level(delta_log, ROUTER_V3, "Trailer packet from %d\n", status.MPI_SOURCE);

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
			int net_id = status.MPI_TAG - COST_UPDATE_TAG;

			assert(net_id >= 0 && net_id < net_route_trees.size());
			assert(nets[net_id].local_id == net_id);

			zlog_level(delta_log, ROUTER_V3, "Cost update of net %d from %d\n", nets[net_id].vpr_id, status.MPI_SOURCE);

			//auto update_start = clock::now();

			//update_start = clock::now();

			int rr_node = data[0];
			unsigned int val;
			path_decoder_t dec;
			decoder_init(dec, (unsigned int *)&data[1], mpi->bit_width);
			do {
				update_one_cost_internal(rr_node, g, congestion, 1, pres_fac);
				//zlog_level(delta_log, ROUTER_V3, "MPI update, source %d node %d delta %d\n", status.MPI_SOURCE, rr_node, delta);

				//if (delta > 0) {
					//assert(find(begin(net_route_trees[net_id]), end(net_route_trees[net_id]), rr_node) == end(net_route_trees[net_id]));
				
					net_route_trees[net_id].push_back(rr_node);
				//} else {
					//assert(delta < 0);
					//assert(find(begin(net_route_trees[net_id]), end(net_route_trees[net_id]), rr_node) != end(net_route_trees[net_id]));
					//net_route_trees[net_id].erase(std::remove(begin(net_route_trees[net_id]), end(net_route_trees[net_id]), rr_node), end(net_route_trees[net_id]));
				//}

				val = decoder_read(dec);
				zlog_level(delta_log, ROUTER_V3, "switch id: %u\n", val);
				if (val != dec.bit_mask) {
					rr_node = get_target(g, get_edge_by_index(g, rr_node, val));
				}
			} while (val != dec.bit_mask);

			zlog_level(delta_log, ROUTER_V3, "address %X num recv %d word offset: %d bit offset %d\n", data, num_recvd, dec.word_offset, dec.bit_offset);

			assert(num_recvd == 1+decoder_get_num_words(dec));

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

void handle_completed_request(int completed_index, const MPI_Status &status, const vector<net_t> &nets, congestion_t *congestion, const RRGraph &g, float pres_fac, vector<vector<RRNode>> &net_route_trees, mpi_context_t *mpi, mpi_perf_t *mpi_perf)
{
	if (completed_index < mpi->comm_size) {
		assert(completed_index != mpi->rank);
		assert(mpi->pending_send_req[completed_index] == MPI_REQUEST_NULL);
		assert(status.MPI_ERROR == MPI_SUCCESS);
		assert(status.MPI_SOURCE == completed_index);

		zlog_level(delta_log, ROUTER_V3, "Handling received packet, completed_index %d source %d\n", completed_index, status.MPI_SOURCE);

#ifdef USE_FIXED_TAG
		handle_packet_fixed_tag(mpi->pending_send_data_nbc[completed_index]->data(), status, nets, congestion, g, pres_fac, net_route_trees, mpi, mpi_perf);
#else 
		handle_packet(mpi->pending_send_data_nbc[completed_index]->data(), status, nets, congestion, g, pres_fac, net_route_trees, mpi, mpi_perf);
#endif

		if (!mpi->received_last_update[completed_index]) {
#ifdef USE_FIXED_TAG
			int error = MPI_Irecv(mpi->pending_send_data_nbc[completed_index]->data(), mpi->pending_send_data_nbc[completed_index]->size(), MPI_INT, completed_index, FIXED_TAG, mpi->comm, &mpi->pending_send_req[completed_index]);
#else 
			int error = MPI_Irecv(mpi->pending_send_data_nbc[completed_index]->data(), mpi->pending_send_data_nbc[completed_index]->size(), MPI_INT, completed_index, MPI_ANY_TAG, mpi->comm, &mpi->pending_send_req[completed_index]);
#endif

			assert(error == MPI_SUCCESS);

			++mpi->num_pending_reqs;
		}
	} else {
		complete_send_request(mpi, completed_index);

		--mpi->num_pending_reqs_by_rank[mpi->pending_send_req_dst[completed_index]];
		assert(mpi->num_pending_reqs_by_rank[mpi->pending_send_req_dst[completed_index]] >= 0);
	}

	--mpi->num_pending_reqs;
}

void sync_combined_wait(const vector<net_t> &nets, congestion_t *congestion, const RRGraph &g, float pres_fac, vector<vector<RRNode>> &net_route_trees, mpi_context_t *mpi, mpi_perf_t *mpi_perf)
{
	assert(!mpi->pending_send_req.empty());

	mpi->max_req_buffer_size = std::max((unsigned long)mpi->max_req_buffer_size, mpi->pending_send_req.size());
	mpi->max_active_reqs = std::max(mpi->max_active_reqs, mpi->num_pending_reqs);
	for (int i = 0; i < mpi->comm_size; ++i) {
		mpi->max_pending_send_reqs_by_rank[i] = std::max(mpi->max_pending_send_reqs_by_rank[i], mpi->num_pending_reqs_by_rank[i]);
	}
	//for (int i = 0; i < mpi->comm_size; ++i) {
		//mpi->num_pending_reqs_by_time[i].push_back(mpi->num_pending_reqs_by_rank[i]);
	//}

	vector<MPI_Status> statuses(mpi->pending_send_req.size());

	int error = MPI_Waitall(mpi->pending_send_req.size(), mpi->pending_send_req.data(), statuses.data());

	assert(error == MPI_SUCCESS);

	for (int i = 0; i < mpi->pending_send_req.size(); ++i) {
		if (statuses[i].MPI_TAG == MPI_ANY_TAG && statuses[i].MPI_SOURCE == MPI_ANY_SOURCE && statuses[i].MPI_ERROR == MPI_SUCCESS) {
			/*empty status*/
			//int num_recvd;
			//int error = MPI_Get_count(&statuses[i], MPI_INT, &num_recvd);
			//assert(error == MPI_SUCCESS);
			//assert(num_recvd == 0);
		} else {
			handle_completed_request(i, statuses[i], nets, congestion, g, pres_fac, net_route_trees, mpi, mpi_perf);
		}
	}
}

void sync_combined(const vector<net_t> &nets, congestion_t *congestion, const RRGraph &g, float pres_fac, vector<vector<RRNode>> &net_route_trees, mpi_context_t *mpi, mpi_perf_t *mpi_perf)
{
	assert(!mpi->pending_send_req.empty());

	mpi->max_req_buffer_size = std::max((unsigned long)mpi->max_req_buffer_size, mpi->pending_send_req.size());
	mpi->max_active_reqs = std::max(mpi->max_active_reqs, mpi->num_pending_reqs);
	for (int i = 0; i < mpi->comm_size; ++i) {
		mpi->max_pending_send_reqs_by_rank[i] = std::max(mpi->max_pending_send_reqs_by_rank[i], mpi->num_pending_reqs_by_rank[i]);
	}
	//for (int i = 0; i < mpi->comm_size; ++i) {
		//mpi->num_pending_reqs_by_time[i].push_back(mpi->num_pending_reqs_by_rank[i]);
	//}

	if (mpi->completed_indices.size() < mpi->pending_send_req.size()) {
		mpi->completed_indices.resize(2*mpi->pending_send_req.size());
	}

	if (mpi->completed_recv_statuses.size() < mpi->pending_send_req.size()) {
		mpi->completed_recv_statuses.resize(2*mpi->pending_send_req.size());
	}

	int num_completed;
	int error = MPI_Testsome(mpi->pending_send_req.size(), mpi->pending_send_req.data(),
			&num_completed, mpi->completed_indices.data(), mpi->completed_recv_statuses.data());
	assert(error == MPI_SUCCESS);

	while (num_completed > 0) {
		for (int i = 0; i < num_completed; ++i) {
			handle_completed_request(mpi->completed_indices[i], mpi->completed_recv_statuses[i], nets, congestion, g, pres_fac, net_route_trees, mpi, mpi_perf);
		}

		error = MPI_Testsome(mpi->pending_send_req.size(), mpi->pending_send_req.data(),
				&num_completed, mpi->completed_indices.data(), mpi->completed_recv_statuses.data());
		assert(error == MPI_SUCCESS);
	}
}

void sync_irecv(const vector<net_t> &nets, congestion_t *congestion, const RRGraph &g, float pres_fac, vector<vector<RRNode>> &net_route_trees, mpi_context_t *mpi, mpi_perf_t *mpi_perf)
{
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

			handle_packet(mpi->pending_recv_data_flat[completed_rank].data(), mpi->completed_recv_statuses[i], nets, congestion, g, pres_fac, net_route_trees, mpi, mpi_perf);

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

void encode_path(int *data, int size, int bit_width, int net_id, int first_rr_node, const vector<int> &edges)
{
#ifdef USE_FIXED_TAG
	write_header_sr(&data[0], NBSRPacketID::ROUTE, net_id);

	//zlog_level(delta_log, ROUTER_V3, "Route packet header %08X,%08X pid %d net_id %d\n", data[0], data[1], static_cast<int>(IbcastPacketID::ROUTE), net_id);
	data[1] = first_rr_node;
#else
	data[0] = first_rr_node;
#endif

	path_encoder_t enc;
#ifdef USE_FIXED_TAG
	encoder_init(enc, (unsigned int *)&data[2], bit_width);
#else
	encoder_init(enc, (unsigned int *)&data[1], bit_width);
#endif
	for (auto i = edges.crbegin(); i != edges.crend(); ++i) {
		encoder_write(enc, *i);
	}
	encoder_write_trailer(enc);

#ifdef USE_FIXED_TAG
	assert(encoder_get_num_words(enc) == size-1-1); /* 1 word for header, 1 word for first_rr_node */
#else
	assert(encoder_get_num_words(enc) == size-1); /* 1 word for first_rr_node */
#endif
}

void route_net_mpi_nonblocking_send_recv_encoded(const RRGraph &g, int vpr_id, int net_id, const source_t *source, const vector<sink_t *> &sinks, const t_router_opts *params, float pres_fac, const vector<net_t> &nets, const vector<vector<net_t *>> &partition_nets, int current_net_index, route_state_t *state, congestion_t *congestion, route_tree_t &rt, t_net_timing &net_timing, vector<vector<RRNode>> &net_route_trees, mpi_context_t *mpi, perf_t *perf, mpi_perf_t *mpi_perf, bool delayed_progress)
{
    using clock = std::chrono::high_resolution_clock;

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
	RRNode existing_opin = RRGraph::null_vertex();
	for (const auto &sink : sorted_sinks) {
		//if (delayed_progress && isink > 0 && mpi->comm_size > 1) {
			//auto progress_start = clock::now();
			////if (q_size(&mpi->pending_bcast_data_q) > 16) { 
				////progress_ibcast_blocking(nets, congestion, g, pres_fac, net_route_trees, mpi);
			////} else {
				//progress_sends_testsome(mpi);
			////}
			//if (mpi_perf) {
				//mpi_perf->total_send_testsome_time += clock::now()-progress_start;
			//}
		//}

		sprintf_rr_node(sink->rr_node, buffer);
		zlog_level(delta_log, ROUTER_V1, "Net %d current sink %d: %s criticality: %g BB: %d-%d %d-%d\n", vpr_id, sink->id, buffer, sink->criticality_fac, sink->current_bounding_box.xmin, sink->current_bounding_box.xmax, sink->current_bounding_box.ymin, sink->current_bounding_box.ymax);

		const auto &sink_rr_node = get_vertex_props(g, sink->rr_node);

		route_tree_multi_root_add_to_heap(rt, g, sink->rr_node, sink->criticality_fac, params->astar_fac, sink->current_bounding_box, heap, perf);

		int count = 0;

		unsigned long old_num_heap_pops = perf ? perf->num_heap_pops : 0;

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

				expand_neighbors(g, item.rr_node, state, congestion, sink_rr_node, sink->criticality_fac, params->astar_fac, [&g, &sink, &sink_rr_node, &v, &item, &existing_opin] (const RRNode &n) -> bool {

					/*if (trace_has_node(prev_trace, id(n))) {*/
						/*zlog_level(delta_log, ROUTER_V3, " existing node route tree ");*/
					/*}*/
					const auto &prop = get_vertex_props(g, n);

					if (existing_opin != RRGraph::null_vertex() && prop.type == OPIN && n != existing_opin) {
					return false;
					}

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

				if (count > 0 && count % params->progress_freq == 0 && mpi->comm_size > 1) {
					auto progress_start = clock::now();
					sync_combined(nets, congestion, g, pres_fac, net_route_trees, mpi, mpi_perf);
					//progress_sends_testsome(mpi);
					if (mpi_perf) {
						mpi_perf->total_sync_time += clock::now()-progress_start;
						++mpi_perf->num_syncs_while_expanding;
					}
				}

				++count;
			} else {
				zlog_level(delta_log, ROUTER_V3, "\tNot expanding neighbor because known_cost %g > %g or cost %g > %g\n", item.known_cost, state[item.rr_node].known_cost, item.cost, state[item.rr_node].cost);
			}
		}

		unsigned long num_heap_pops_per_sink = perf ? (perf->num_heap_pops - old_num_heap_pops) : 0;

		if (perf) {
			perf->max_num_heap_pops_per_sink = std::max(perf->max_num_heap_pops_per_sink, num_heap_pops_per_sink);
			perf->min_num_heap_pops_per_sink = std::min(perf->min_num_heap_pops_per_sink, num_heap_pops_per_sink);
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

			vector<int> edges;

			RRNode child_rr_node = sink->rr_node;
			RRNode first_rr_node = RRGraph::null_vertex();

			get_path(sink->rr_node, state, rt, g, [&] (const RREdge &edge) -> void
					{
					RRNode parent_rr_node = get_source(g, edge);

					const auto &parent_rr_node_p = get_vertex_props(g, parent_rr_node);

					bool update_cost = false;

					sprintf_rr_node(parent_rr_node, buffer);
					zlog_level(delta_log, ROUTER_V3, "Adding %s to route tree\n", buffer);

					RouteTreeNode parent_rt_node = route_tree_add_rr_node(rt, parent_rr_node, g);
					if (parent_rt_node != RouteTree::null_vertex()) {
						route_tree_set_node_properties(rt, parent_rt_node, parent_rr_node_p.type != IPIN && parent_rr_node_p.type != SINK, state[parent_rr_node].upstream_R, state[parent_rr_node].delay);
					} 

					assert(child_rr_node == get_target(g, edge));

					char buffer2[256];
					sprintf_rr_node(child_rr_node, buffer2);

					zlog_level(delta_log, ROUTER_V3, "Adding edge %s -> %s\n", buffer, buffer2);

					const RouteTreeEdge &rt_edge = route_tree_add_edge_between_rr_node(rt, parent_rr_node, child_rr_node); 
					auto &rt_edge_props = get_edge_props(rt.graph, rt_edge);
					rt_edge_props.rr_edge = edge;

					if (parent_rr_node_p.type == OPIN && existing_opin == RRGraph::null_vertex()) {
						existing_opin = parent_rr_node;
					}

					child_rr_node = parent_rr_node;

					if ((parent_rt_node != RouteTree::null_vertex()) || (get_vertex_props(g, parent_rr_node).type == SOURCE)) {
						update_one_cost_internal(parent_rr_node, g, congestion, 1, pres_fac);

						int id = get_edge_props(g, edge).id;

						zlog_level(delta_log, ROUTER_V3, "Adding edge id %d\n", id);

						edges.push_back(id);

						first_rr_node = parent_rr_node;
					}

					zlog_level(delta_log, ROUTER_V3, "\n");
					});

			assert(first_rr_node != RRGraph::null_vertex());

			if (mpi->comm_size > 1) {
				auto broadcast_start = clock::now();

				int payload_size = 1 + (int)ceil((float)((edges.size()+1)*mpi->bit_width) / 32);
#ifdef USE_FIXED_TAG
				int header_size = 1;
				int tag = FIXED_TAG;
#else
				int header_size = 0;
				int tag = COST_UPDATE_TAG+net_id;
#endif
				int packet_size = payload_size+header_size;
				int data_index = get_free_buffer(mpi, packet_size);

				assert(data_index >= mpi->comm_size);

				int *data = mpi->pending_send_data_nbc[data_index]->data();

				encode_path(data, packet_size, mpi->bit_width, net_id, first_rr_node, edges);

				zlog_level(delta_log, ROUTER_V3, "Broadcasting cost update, net %d sink %d\n", net_id, sink->rr_node);

				selective_broadcast(data_index, packet_size, tag, mpi);

				if (mpi_perf) {
					mpi_perf->total_broadcast_time += clock::now()-broadcast_start;
				}
			} 

			if (sink->id != -1) {
				net_timing.delay[sink->id+1] = get_vertex_props(rt.graph, route_tree_get_rt_node(rt, sink->rr_node)).delay;
			}

			zlog_level(delta_log, ROUTER_V2, "Routed sink %d\n", sink->rr_node);
		}

		if (mpi->comm_size > 1) {
#if 0
			auto sync_start = clock::now();

			sync_irecv(nets, congestion, g, pres_fac, net_route_trees, mpi, mpi_perf);

			if (mpi_perf) {
				mpi_perf->total_sync_time += clock::now()-sync_start;
			}

			auto progress_start = clock::now();

			progress_sends_testsome(mpi);

			if (mpi_perf) {
				mpi_perf->total_send_testsome_time += clock::now()-progress_start;
			}
#else
			auto sync_start = clock::now();
			sync_combined(nets, congestion, g, pres_fac, net_route_trees, mpi, mpi_perf);
			if (mpi_perf) {
				mpi_perf->total_sync_time += clock::now()-sync_start;
			}
#endif
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

