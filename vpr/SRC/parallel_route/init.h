#ifndef INIT_H
#define INIT_H

#include "new_rr_graph.h"

void init_graph(RRGraph &g);
void init_graph(RRGraph &g, std::vector<std::vector<RRNode>> &sink_in_nodes, std::vector<std::vector<RRNode>> &ipin_in_nodes);
void delete_graph(RRGraph &g);

void init_nets(std::vector<net_t> &nets, std::vector<net_t> &global_nets, int bb_factor, bool large_bb);
void init_nets(std::vector<net_t> &nets, std::vector<net_t> &global_nets, int bb_factor, const std::vector<std::vector<RRNode>> &sink_in_nodes);

void init_net_timing(const std::vector<net_t> &nets, const std::vector<net_t> &global_nets, t_net_timing *net_timing);

void delete_net_timing(const std::vector<net_t> &nets, const std::vector<net_t> &global_nets, t_net_timing *net_timing);

void init_partitioned_graph(RRGraph &g, int num_partitions);

#endif
