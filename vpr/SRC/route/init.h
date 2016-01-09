#ifndef INIT_H
#define INIT_H

#include "new_rr_graph.h"

void init_logging();

void init_graph(RRGraph &g);

void init_nets(std::vector<net_t> &nets, std::vector<net_t> &global_nets, int bb_factor);

void init_net_timing(const std::vector<net_t> &nets, const std::vector<net_t> &global_nets, t_net_timing *net_timing);

#endif
