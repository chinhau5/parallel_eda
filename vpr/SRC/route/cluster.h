#ifndef CLUSTER_H
#define CLUSTER_H

#include "route.h"

void create_clustered_virtual_nets(std::vector<net_t> &nets, int num_nodes_per_cluster, int sink_bb_area_threshold, std::vector<std::vector<virtual_net_t>> &virtual_nets);

#endif
