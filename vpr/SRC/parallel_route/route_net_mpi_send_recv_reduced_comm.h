#ifndef ROUTE_NET_MPI_SEND_RECV_REDUCED_COMM_H

#define ROUTE_NET_MPI_SEND_RECV_REDUCED_COMM_H

#define LAST_TAG (0x1000)
#define RIP_UP_TAG (LAST_TAG+1)
#define COST_UPDATE_TAG (RIP_UP_TAG+0x1000000) //support for up to 16M nets

#endif
