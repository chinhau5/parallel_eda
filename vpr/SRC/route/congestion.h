#ifndef CONGESTION_H
#define CONGESTION_H

#include <mpi.h>

void update_one_cost(const RRGraph &g, congestion_t *congestion, const vector<RRNode>::const_iterator &rr_nodes_begin, const vector<RRNode>::const_iterator &rr_nodes_end, /*int net_id,*/ int delta, float pres_fac, bool lock, lock_perf_t *lock_perf);

//void update_one_cost(const RRGraph &g, congestion_t *congestion, route_tree_t &rt, const RouteTreeNode &node, int delta, float pres_fac, bool lock);

void update_one_cost_mpi(const vector<RRNode>::const_iterator &rr_nodes_begin, const vector<RRNode>::const_iterator &rr_nodes_end, const RRGraph &g, const vector<int> &pid, int this_pid, congestion_mpi_t *congestion, MPI_Win win, int delta, float pres_fac);

void update_one_cost_internal(RRNode rr_node, const rr_node_property_t &rr_node_p, congestion_t &congestion, /*int net_id, */int delta, float pres_fac, bool lock, lock_perf_t *lock_perf);

void update_one_cost_internal_mpi(RRNode rr_node, const RRGraph &g, const vector<int> &pid, int this_pid, congestion_mpi_t *congestion, MPI_Win win, /*int net_id, */int delta, float pres_fac);

void update_costs(const RRGraph &g, congestion_t *congestion, float pres_fac, float acc_fac);

void update_costs_mpi(const RRGraph &g, vector<int> &pid, int this_pid, congestion_mpi_t *congestion, MPI_Win win, float pres_fac, float acc_fac);

MPI_Datatype get_occ_dt();

void init_congestion_mpi_datatype();

#endif
