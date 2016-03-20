#ifndef UTILITY_H
#define UTILITY_H

#define PRINT_RR_NODE

#ifdef PRINT_RR_NODE

void sprintf_rr_node(int inode, char *buffer);

#else

#define sprintf_rr_node(...)

#endif

#endif
