#include <mpi.h>
#include "route.h"
#include "graph.h"
#include "route_tree.h"
#include "util.h"
#include "log.h"

static MPI_Datatype occ_dt;
static MPI_Datatype acc_dt; 

MPI_Datatype get_real_occ_dt()
{
	return occ_dt;
}

MPI_Datatype get_occ_dt()
{
	return MPI_INT;
}

int get_occ_disp(RRNode rr_node)
{
	return sizeof(congestion_mpi_t)*rr_node + offsetof(congestion_mpi_t, occ);
}

void init_congestion_mpi_datatype()
{
	int bl[] = { 1, 1 };
	MPI_Aint disp[] = { offsetof(congestion_mpi_t, occ), sizeof(congestion_mpi_t) };
	MPI_Datatype dts[] = { MPI_INT, MPI_UB };

	assert(MPI_Type_create_struct(1, bl, disp, dts, &occ_dt) == MPI_SUCCESS);
	assert(MPI_Type_commit(&occ_dt) == MPI_SUCCESS);


	disp[0] = offsetof(congestion_mpi_t, acc_cost);
	dts[0] = MPI_FLOAT;
	assert(MPI_Type_create_struct(1, bl, disp, dts, &acc_dt) == MPI_SUCCESS);
	assert(MPI_Type_commit(&acc_dt) == MPI_SUCCESS);
}

void update_costs_mpi(const RRGraph &g, vector<int> &pid, int this_pid, congestion_mpi_t *congestion, MPI_Win win, float pres_fac, float acc_fac)
{
	for (const auto &rr_node : get_vertices(g)) {
		int rr_node_pid = pid[rr_node];
		int from_pid = rr_node_pid == -1 ? 0 : rr_node_pid;
		assert(MPI_Win_lock(MPI_LOCK_SHARED, from_pid, 0, win) == MPI_SUCCESS);

		if (from_pid != this_pid) {
			congestion[rr_node].occ = std::numeric_limits<int>::min();
			assert(MPI_Get(&congestion[rr_node].occ, 1, get_occ_dt(),
						from_pid,
						get_occ_disp(rr_node), 1, get_occ_dt(),
						win) == MPI_SUCCESS);
			assert(MPI_Win_flush(from_pid, win) == MPI_SUCCESS);
			assert(congestion[rr_node].occ != std::numeric_limits<int>::min());
		}

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

		assert(MPI_Win_unlock(from_pid, win) == MPI_SUCCESS);
	}
}

void update_costs_mpi_old(const RRGraph &g, vector<int> &pid, int this_pid, congestion_mpi_t *congestion, MPI_Win win, float pres_fac, float acc_fac)
{
	for (const auto &rr_node : get_vertices(g)) {
		int rr_node_pid = pid[rr_node];
		int from_pid = rr_node_pid == -1 ? 0 : rr_node_pid;
		assert(MPI_Win_lock(MPI_LOCK_EXCLUSIVE, from_pid, 0, win) == MPI_SUCCESS);

		if (from_pid != this_pid) {
			congestion[rr_node].occ = std::numeric_limits<int>::min();
			assert(MPI_Get(&congestion[rr_node], 1, occ_dt,
						from_pid,
						rr_node, 1, occ_dt,
						win) == MPI_SUCCESS);
			assert(MPI_Win_flush(from_pid, win) == MPI_SUCCESS);
			assert(congestion[rr_node].occ != std::numeric_limits<int>::min());
		}

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

		if (from_pid != this_pid) {
			assert(MPI_Put(&congestion[rr_node], 1, occ_dt,
						from_pid,
						rr_node, 1, occ_dt,
						win) == MPI_SUCCESS);
		} 

		assert(MPI_Win_unlock(from_pid, win) == MPI_SUCCESS);
	}
}

void update_one_cost_internal_mpi_new(RRNode rr_node, const RRGraph &g, const vector<int> &pid, int this_pid, int num_procs, congestion_mpi_t *congestion, int delta, float pres_fac)
{
	int rr_node_pid = pid[rr_node];

	int from_pid = rr_node_pid == -1 ? 0 : rr_node_pid;

	congestion_mpi_t acc;
	acc.occ = delta;

	congestion[rr_node].occ += delta;

	assert(congestion[rr_node].occ >= 0);

	const auto &rr_node_p = get_vertex_props(g, rr_node);

	if (congestion[rr_node].occ < rr_node_p.capacity) {
		congestion[rr_node].pres_cost = 1;
	} else {
		congestion[rr_node].pres_cost = 1 + (congestion[rr_node].occ + 1 - rr_node_p.capacity) * pres_fac;
	}

	for (int i = 0; i < num_procs; ++i) {
		if (i != this_pid) {
			MPI_Request req;
			assert(MPI_Isend(&delta, 1, MPI_INT, i, rr_node, MPI_COMM_WORLD, &req)
					== MPI_SUCCESS);
		}
	}

	char buffer[256];
	sprintf_rr_node(rr_node, buffer);
	zlog_error(delta_log, "Update cost of %s delta: %d new_occ: %d pres_fac: %g\n", buffer, delta, congestion[rr_node].occ, pres_fac);
	//zlog_level(delta_log, ROUTER_V2, "Update cost of %s delta: %d new_occ: %d pres_fac: %g\n", buffer, delta, congestion[rr_node].occ, pres_fac);
}

void update_one_cost_internal_mpi(RRNode rr_node, const RRGraph &g, const vector<int> &pid, int this_pid, congestion_mpi_t *congestion, MPI_Win win, /*int net_id, */int delta, float pres_fac)
{
	int rr_node_pid = pid[rr_node];

	int from_pid = rr_node_pid == -1 ? 0 : rr_node_pid;

	assert(MPI_Win_lock(MPI_LOCK_EXCLUSIVE, from_pid, 0, win) == MPI_SUCCESS);

	congestion_mpi_t acc;
	acc.occ = delta;

	if (from_pid != this_pid) {
		congestion[rr_node].occ = std::numeric_limits<int>::min();
		assert(MPI_Fetch_and_op(&acc.occ, &congestion[rr_node].occ, get_occ_dt(), from_pid, get_occ_disp(rr_node), MPI_SUM, win) == MPI_SUCCESS);
		assert(MPI_Win_flush(from_pid, win) == MPI_SUCCESS);
		assert(congestion[rr_node].occ != std::numeric_limits<int>::min());
	} 

	congestion[rr_node].occ += delta;

	assert(congestion[rr_node].occ >= 0);

	const auto &rr_node_p = get_vertex_props(g, rr_node);

	if (congestion[rr_node].occ < rr_node_p.capacity) {
		congestion[rr_node].pres_cost = 1;
	} else {
		congestion[rr_node].pres_cost = 1 + (congestion[rr_node].occ + 1 - rr_node_p.capacity) * pres_fac;
	}

	char buffer[256];
	sprintf_rr_node(rr_node, buffer);
	zlog_error(delta_log, "Update cost of %s delta: %d new_occ: %d pres_fac: %g\n", buffer, delta, congestion[rr_node].occ, pres_fac);
	//zlog_level(delta_log, ROUTER_V2, "Update cost of %s delta: %d new_occ: %d pres_fac: %g\n", buffer, delta, congestion[rr_node].occ, pres_fac);

	assert(MPI_Win_unlock(from_pid, win) == MPI_SUCCESS);
}

void update_one_cost_internal_mpi_old(RRNode rr_node, const RRGraph &g, const vector<int> &pid, int this_pid, congestion_mpi_t *congestion, MPI_Win win, /*int net_id, */int delta, float pres_fac)
{
	int rr_node_pid = pid[rr_node];

	int from_pid = rr_node_pid == -1 ? 0 : rr_node_pid;

	assert(MPI_Win_lock(MPI_LOCK_EXCLUSIVE, from_pid, 0, win) == MPI_SUCCESS);

	if (from_pid != this_pid) {
		congestion[rr_node].occ = std::numeric_limits<int>::min();
		assert(MPI_Get(&congestion[rr_node], 1, occ_dt,
				from_pid,
				rr_node, 1, occ_dt,
				win) == MPI_SUCCESS);
		assert(MPI_Win_flush(from_pid, win) == MPI_SUCCESS);
		assert(congestion[rr_node].occ != std::numeric_limits<int>::min());
	} 

	congestion[rr_node].occ += delta;

	assert(congestion[rr_node].occ >= 0);

	const auto &rr_node_p = get_vertex_props(g, rr_node);

	if (congestion[rr_node].occ < rr_node_p.capacity) {
		congestion[rr_node].pres_cost = 1;
	} else {
		congestion[rr_node].pres_cost = 1 + (congestion[rr_node].occ + 1 - rr_node_p.capacity) * pres_fac;
	}

	char buffer[256];
	sprintf_rr_node(rr_node, buffer);
	zlog_level(delta_log, ROUTER_V2, "Update cost of %s delta: %d new_occ: %d pres_fac: %g\n", buffer, delta, congestion[rr_node].occ, pres_fac);

	if (from_pid != this_pid) {
		assert(MPI_Put(&congestion[rr_node], 1, occ_dt,
					from_pid,
					rr_node, 1, occ_dt,
					win) == MPI_SUCCESS);
	} 

	assert(MPI_Win_unlock(from_pid, win) == MPI_SUCCESS);
}

void update_one_cost_mpi(const vector<RRNode>::const_iterator &rr_nodes_begin, const vector<RRNode>::const_iterator &rr_nodes_end, const RRGraph &g, const vector<int> &pid, int this_pid, congestion_mpi_t *congestion, MPI_Win win, int delta, float pres_fac)
{
	/*const RouteTreeNode *last = nullptr;*/
	for (auto iter = rr_nodes_begin; iter != rr_nodes_end; ++iter) {
		/* we don't update get_source because that's the link to existing route tree and the cost is handled by update_one_cost_internal */
		/*const RouteTreeNode &rt_node = get_target(rt.graph, get_edge(rt.graph, rt_edge_id));*/
		/*last = &get_target(rt.graph, get_edge(rt.graph, rt_edge_id));*/
		/*RRNode &rr_node = get_vertex(g, rt_node.rr_node);*/
		update_one_cost_internal_mpi(*iter, g, pid, this_pid, congestion, win, delta, pres_fac);
	}
	/*RRNode &rr_node = get_vertex(g, last->rr_node);*/
	/*update_one_cost_internal(rr_node, delta, pres_fac);*/
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

void update_one_cost_internal(RRNode rr_node, const RRGraph &g, congestion_mpi_t *congestion, /*int net_id, */int delta, float pres_fac)
{
	congestion[rr_node].occ += delta;

	assert(congestion[rr_node].occ >= 0);

	const auto &rr_node_p = get_vertex_props(g, rr_node);

	if (congestion[rr_node].occ < rr_node_p.capacity) {
		congestion[rr_node].pres_cost = 1;
	} else {
		congestion[rr_node].pres_cost = 1 + (congestion[rr_node].occ + 1 - rr_node_p.capacity) * pres_fac;
	}

	char buffer[256];
	sprintf_rr_node(rr_node, buffer);
	zlog_level(delta_log, ROUTER_V2, "Update cost of %s delta: %d new_occ: %d pres_fac: %g\n", buffer, delta, congestion[rr_node].occ, pres_fac);
}

void update_one_cost(const RRGraph &g, congestion_mpi_t *congestion, const vector<RRNode>::const_iterator &rr_nodes_begin, const vector<RRNode>::const_iterator &rr_nodes_end, int delta, float pres_fac)
{
	/*const RouteTreeNode *last = nullptr;*/
	for (auto iter = rr_nodes_begin; iter != rr_nodes_end; ++iter) {
		/* we don't update get_source because that's the link to existing route tree and the cost is handled by update_one_cost_internal */
		/*const RouteTreeNode &rt_node = get_target(rt.graph, get_edge(rt.graph, rt_edge_id));*/
		/*last = &get_target(rt.graph, get_edge(rt.graph, rt_edge_id));*/
		/*RRNode &rr_node = get_vertex(g, rt_node.rr_node);*/
		update_one_cost_internal(*iter, g, congestion, /*net_id,*/ delta, pres_fac);
	}
	/*RRNode &rr_node = get_vertex(g, last->rr_node);*/
	/*update_one_cost_internal(rr_node, delta, pres_fac);*/
}

//void update_one_cost(const RRGraph &g, congestion_t *congestion, route_tree_t &rt, RouteTreeNode rt_node, int delta, float pres_fac, bool lock)
//{
	//const auto &rt_node_p = get_vertex_props(rt.graph, rt_node);

	//assert(rt_node_p.valid);

	//int rr_node = rt_node_p.rr_node;

	//const auto &rr_node_p = get_vertex_props(g, rr_node);

	//update_one_cost_internal(rr_node, rr_node_p, congestion[rr_node], [>-1,<] delta, pres_fac, lock, nullptr);

	//for (const auto &e : get_out_edges(rt.graph, rt_node)) {
		//update_one_cost(g, congestion, rt, get_target(rt.graph, e), delta, pres_fac, lock);
	//}
//}

