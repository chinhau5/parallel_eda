#ifndef UTILITY_H
#define UTILITY_H

#include "new_rr_graph.h"

//#define PRINT_RR_NODE

#ifdef PRINT_RR_NODE

void init_sprintf_rr_node(const RRGraph *_rr_graph);

void sprintf_rr_node(int inode, char *buffer);

#else

#define init_sprintf_rr_node(...)

#define sprintf_rr_node(...) 

#endif

#endif
