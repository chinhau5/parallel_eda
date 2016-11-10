#ifndef UTILITY_H
#define UTILITY_H

#include "new_rr_graph.h"

#define PRINT_RR_NODE

void init_sprintf_rr_node(const RRGraph *_rr_graph);

void sprintf_rr_node_impl(int inode, char *buffer);

#ifdef PRINT_RR_NODE

#define sprintf_rr_node(inode, buffer) sprintf_rr_node_impl((inode), (buffer))

#else

#define sprintf_rr_node(inode, buffer) 

#endif

#endif
