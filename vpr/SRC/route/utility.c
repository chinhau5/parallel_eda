#include <assert.h>
#include "vpr_types.h"

static const char *rr_types[] =  {
	"SOURCE", "SINK", "IPIN", "OPIN", "CHANX", "CHANY", "INTRA_CLUSTER_EDGE"
};

/*#define PRINT_RR_NODE*/

/*#ifdef PRINT_RR_NODE*/

void sprintf_rr_node(int inode, char *buffer)
{
	extern t_rr_node *rr_node;
	if (rr_node[inode].direction == INC_DIRECTION) {
		sprintf(buffer, "%d %s (%d,%d)(%d,%d) ptc=%d dir=%d", inode, rr_types[rr_node[inode].type], rr_node[inode].xlow, rr_node[inode].ylow, rr_node[inode].xhigh, rr_node[inode].yhigh, rr_node[inode].ptc_num, rr_node[inode].direction);
	} else if (rr_node[inode].direction == DEC_DIRECTION) {
		sprintf(buffer, "%d %s (%d,%d)(%d,%d) ptc=%d dir=%d", inode, rr_types[rr_node[inode].type], rr_node[inode].xhigh, rr_node[inode].yhigh, rr_node[inode].xlow, rr_node[inode].ylow, rr_node[inode].ptc_num, rr_node[inode].direction);
	} else {
		sprintf(buffer, "%d %s (%d,%d)(%d,%d) ptc=%d dir=%d", inode, rr_types[rr_node[inode].type], rr_node[inode].xlow, rr_node[inode].ylow, rr_node[inode].xhigh, rr_node[inode].yhigh, rr_node[inode].ptc_num, rr_node[inode].direction);
	}

	/*if (rr_node[inode].direction == INC_DIRECTION) {*/
		/*if (rr_node[inode].type == CHANX) {*/
			/*sprintf(buffer, "%d %s (%d->%d,%d) %d %d", inode, rr_types[rr_node[inode].type], rr_node[inode].xlow, rr_node[inode].xhigh, rr_node[inode].ylow, rr_node[inode].ptc_num, rr_node[inode].direction);*/
		/*} else {*/
			/*assert(rr_node[inode].type == CHANY);*/
			/*sprintf(buffer, "%d %s (%d,%d->%d) %d %d", inode, rr_types[rr_node[inode].type], rr_node[inode].xlow, rr_node[inode].ylow, rr_node[inode].yhigh, rr_node[inode].ptc_num, rr_node[inode].direction); */
		/*} */
	/*} else if (rr_node[inode].direction == DEC_DIRECTION) {*/
		/*if (rr_node[inode].type == CHANX) {*/
			/*sprintf(buffer, "%d %s (%d->%d,%d) %d %d", inode, rr_types[rr_node[inode].type], rr_node[inode].xhigh, rr_node[inode].xlow, rr_node[inode].ylow, rr_node[inode].ptc_num, rr_node[inode].direction);*/
		/*} else {*/
			/*assert(rr_node[inode].type == CHANY);*/
			/*sprintf(buffer, "%d %s (%d,%d->%d) %d %d", inode, rr_types[rr_node[inode].type], rr_node[inode].xlow, rr_node[inode].yhigh, rr_node[inode].ylow, rr_node[inode].ptc_num, rr_node[inode].direction);*/
		/*}*/
	/*} else {*/
		/*sprintf(buffer, "%d %s (%d,%d)(%d,%d) %d %d", inode, rr_types[rr_node[inode].type], rr_node[inode].xlow, rr_node[inode].ylow, rr_node[inode].xhigh, rr_node[inode].yhigh, rr_node[inode].ptc_num, rr_node[inode].direction);*/
	/*}*/
}

/*#endif*/
