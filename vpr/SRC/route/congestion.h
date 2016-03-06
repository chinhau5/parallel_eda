#ifndef CONGESTION_H
#define CONGESTION_H

void update_one_cost(const RRGraph &g, congestion_t *congestion, const vector<int>::const_iterator &rr_nodes_begin, const vector<int>::const_iterator &rr_nodes_end, /*int net_id,*/ int delta, float pres_fac, bool lock, lock_perf_t *lock_perf);

void update_one_cost(const RRGraph &g, congestion_t *congestion, route_tree_t &rt, const RouteTreeNode &node, int delta, float pres_fac, bool lock);

void update_one_cost_internal(int rr_node, const rr_node_property_t &rr_node_p, congestion_t &congestion, /*int net_id, */int delta, float pres_fac, bool lock, lock_perf_t *lock_perf);

void update_costs(const RRGraph &g, congestion_t *congestion, float pres_fac, float acc_fac);

#endif
