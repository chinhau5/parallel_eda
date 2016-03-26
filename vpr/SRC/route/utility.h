#ifndef UTILITY_H
#define UTILITY_H

void sprintf_rr_node(int inode, char *buffer);

//#define PRINT_RR_NODE

#ifndef PRINT_RR_NODE

#define sprintf_rr_node(...)

#else

#endif

#endif
