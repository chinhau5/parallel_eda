#ifndef CONGESTION_H
#define CONGESTION_H

#include <mpi.h>

inline int get_occ(const congestion_local_t *local, int node)
{
	return local->global[node].occ + local->local[node].occ_delta;
}

inline void add_occ(congestion_local_t *local, int node, int delta)
{
	local->local[node].occ_delta += delta;
	local->dirty_nodes.insert(node);
}

inline float get_pres_cost(const congestion_local_t *local, int node)
{
	return local->local[node].pres_cost;
}

inline void set_pres_cost(congestion_local_t *local, int node, float pres_cost)
{
	local->local[node].pres_cost = pres_cost;
	assert(local->dirty_nodes.find(node) != end(local->dirty_nodes));
}

inline float get_acc_cost(const congestion_local_t *local, int node)
{
	return local->global[node].pres_cost;
}

inline void commit(congestion_local_t *local, set<int> &all_dirty_nodes)
{
	for (const auto &dirty_node : local->dirty_nodes) {
		local->global[dirty_node].occ += local->local[dirty_node].occ_delta;
		local->local[dirty_node].occ_delta = 0;
		all_dirty_nodes.insert(dirty_node);
	}
	local->dirty_nodes.clear();
}

//inline int get_occ(const congestion_global_t *global, int node)
//{
	//return global->current[node].occ;
//}

//inline int get_recalc_occ(const congestion_global_t *global, int node)
//{
	//return global->current[node].recalc_occ;
//}

//inline void add_recalc_occ(congestion_global_t *global, int node, int delta)
//{
	//global->current[node].recalc_occ += delta;
//}

//------

inline int get_occ(const congestion_locked_t *congestion, int node)
{
	return congestion[node].cong.occ;
}

inline void add_occ(congestion_locked_t *congestion, int node, int delta)
{
	congestion[node].cong.occ += delta;
}

inline int get_occ(const congestion_t *congestion, int node)
{
	return congestion[node].occ;
}

inline void add_occ(congestion_t *congestion, int node, int delta)
{
	congestion[node].occ += delta;
}

inline float get_pres_cost(const congestion_locked_t *congestion, int node)
{
	return congestion[node].cong.pres_cost;
}

inline void set_pres_cost(congestion_locked_t *congestion, int node, float pres_cost)
{
	congestion[node].cong.pres_cost = pres_cost;
}

inline float get_pres_cost(const congestion_t *congestion, int node)
{
	return congestion[node].pres_cost;
}

inline void set_pres_cost(congestion_t *congestion, int node, float pres_cost)
{
	congestion[node].pres_cost = pres_cost;
}

inline float get_acc_cost(const congestion_locked_t *congestion, int node)
{
	return congestion[node].cong.acc_cost;
}

inline void add_acc_cost(congestion_locked_t *congestion, int node, float delta)
{
	congestion[node].cong.acc_cost += delta;
}

inline float get_acc_cost(const congestion_t *congestion, int node)
{
	return congestion[node].acc_cost;
}

inline void add_acc_cost(congestion_t *congestion, int node, float delta)
{
	congestion[node].acc_cost += delta;
}

inline int get_recalc_occ(const congestion_locked_t *congestion, int node)
{
	return congestion[node].cong.recalc_occ;
}

inline void add_recalc_occ(congestion_locked_t *congestion, int node, int delta)
{
	congestion[node].cong.recalc_occ += delta;
}

inline int get_recalc_occ(const congestion_t *congestion, int node)
{
	return congestion[node].recalc_occ;
}

inline void add_recalc_occ(congestion_t *congestion, int node, int delta)
{
	congestion[node].recalc_occ += delta;
}

void update_one_cost(const RRGraph &g, congestion_locked_t *congestion, const vector<RRNode>::const_iterator &rr_nodes_begin, const vector<RRNode>::const_iterator &rr_nodes_end, /*int net_id,*/ int delta, float pres_fac, bool lock, lock_perf_t *lock_perf);

void update_one_cost(const RRGraph &g, congestion_t *congestion, const vector<RRNode>::const_iterator &rr_nodes_begin, const vector<RRNode>::const_iterator &rr_nodes_end, int delta, float pres_fac);

void update_one_cost_mpi_rma(const vector<RRNode>::const_iterator &rr_nodes_begin, const vector<RRNode>::const_iterator &rr_nodes_end, const RRGraph &g, const vector<int> &pid, int this_pid, congestion_t *congestion, MPI_Win win, int delta, float pres_fac);

void update_one_cost_mpi_send(const vector<RRNode>::const_iterator &rr_nodes_begin, const vector<RRNode>::const_iterator &rr_nodes_end, const RRGraph &g, congestion_t *congestion, int delta, float pres_fac, int this_pid, int num_procs, MPI_Comm comm, vector<ongoing_transaction_t> &transactions);

void update_one_cost_internal(RRNode rr_node, const RRGraph &g, congestion_locked_t *congestion, /*int net_id, */int delta, float pres_fac, bool lock, lock_perf_t *lock_perf);

void update_one_cost_internal(RRNode rr_node, const RRGraph &g, congestion_t *congestion, /*int net_id, */int delta, float pres_fac);

void update_one_cost_internal(RRNode rr_node, const RRGraph &g, congestion_local_t *congestion, int delta, float pres_fac);

//void update_one_cost(const RRGraph &g, congestion_t *congestion, route_tree_t &rt, const RouteTreeNode &node, int delta, float pres_fac, bool lock);

void update_one_cost_internal_mpi_rma(RRNode rr_node, const RRGraph &g, const vector<int> &pid, int this_pid, congestion_t *congestion, MPI_Win win, /*int net_id, */int delta, float pres_fac);

void update_one_cost_internal_mpi_send(RRNode rr_node, const RRGraph &g, congestion_t *congestion, int delta, float pres_fac, int this_pid, int num_procs, MPI_Comm comm);

template<typename Congestion>
void update_costs(const RRGraph &g, Congestion *congestion, float pres_fac, float acc_fac)
{
	for (const auto &rr_node : get_vertices(g)) {
		int occ = get_occ(congestion, rr_node);
		int capacity = get_vertex_props(g, rr_node).capacity;
		if (occ > capacity) {
			add_acc_cost(congestion, rr_node, (occ - capacity) * acc_fac);
			set_pres_cost(congestion, rr_node, 1. + (occ + 1 - capacity) * pres_fac);
		} else if (occ == capacity) {
			/* If occ == capacity, we don't need to increase acc_cost, but a change    *
			 * in pres_fac could have made it necessary to recompute the cost anyway.  */
			set_pres_cost(congestion, rr_node, 1. + pres_fac);
		}
	}
}


void update_costs_mpi(const RRGraph &g, vector<int> &pid, int this_pid, congestion_t *congestion, MPI_Win win, float pres_fac, float acc_fac);

MPI_Datatype get_occ_dt();
int get_occ_disp(RRNode rr_node);

void init_congestion_mpi_datatype();

#endif
