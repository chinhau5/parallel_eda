#ifndef ROUTE_NONBLOCKING_SEND_RECV_ENCODED
#define ROUTE_NONBLOCKING_SEND_RECV_ENCODED

#define LAST_TAG (0x1000)
#define RIP_UP_TAG (LAST_TAG+1)
#define COST_UPDATE_TAG (RIP_UP_TAG+0x1000000) //support for up to 16M nets
#define FIXED_TAG (64)

enum class NBSRPacketID : unsigned int {
	RIP_UP,
	ROUTE,
	TRAILER,
	NUM_PACKET_IDS
};

#define USE_FIXED_TAG

#endif
