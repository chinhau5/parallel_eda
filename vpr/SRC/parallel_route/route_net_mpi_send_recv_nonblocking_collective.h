#ifndef ROUTE_NONBLOCKING_COLLECTIVE
#define ROUTE_NONBLOCKING_COLLECTIVE

enum class PacketID : unsigned int {
	RIP_UP_ALL,
	RIP_UP,
	ROUTE,
	NO_OP,
	NUM_PACKET_IDS
};

#endif
