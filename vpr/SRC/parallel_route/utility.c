#include <assert.h>
#include "vpr_types.h"
#include "utility.h"

static const char *rr_types[] =  {
	"SOURCE", "SINK", "IPIN", "OPIN", "CHANX", "CHANY", "INTRA_CLUSTER_EDGE"
};

/*#define PRINT_RR_NODE*/

#ifdef PRINT_RR_NODE

static const RRGraph *rr_graph;

void init_sprintf_rr_node(const RRGraph *_rr_graph)
{
	rr_graph = _rr_graph;
}

void sprintf_rr_node(int inode, char *buffer)
{
	/*assert(false);*/
	/*extern t_rr_node *rr_node;*/
	/*if (rr_node[inode].direction == INC_DIRECTION) {*/
		/*sprintf(buffer, "%d %s (%d,%d)(%d,%d) ptc=%d dir=%d", inode, rr_types[rr_node[inode].type], rr_node[inode].xlow, rr_node[inode].ylow, rr_node[inode].xhigh, rr_node[inode].yhigh, rr_node[inode].ptc_num, rr_node[inode].direction);*/
	/*} else if (rr_node[inode].direction == DEC_DIRECTION) {*/
		/*sprintf(buffer, "%d %s (%d,%d)(%d,%d) ptc=%d dir=%d", inode, rr_types[rr_node[inode].type], rr_node[inode].xhigh, rr_node[inode].yhigh, rr_node[inode].xlow, rr_node[inode].ylow, rr_node[inode].ptc_num, rr_node[inode].direction);*/
	/*} else {*/
		/*sprintf(buffer, "%d %s (%d,%d)(%d,%d) ptc=%d dir=%d", inode, rr_types[rr_node[inode].type], rr_node[inode].xlow, rr_node[inode].ylow, rr_node[inode].xhigh, rr_node[inode].yhigh, rr_node[inode].ptc_num, rr_node[inode].direction);*/
	/*}*/

	const auto &rr_node_p = get_vertex_props(*rr_graph, inode);

	if (!rr_node_p.inc_direction) {
		sprintf(buffer, "%d %s (%d,%d)(%d,%d) dir=%d", inode, rr_types[rr_node_p.type], rr_node_p.xhigh, rr_node_p.yhigh, rr_node_p.xlow, rr_node_p.ylow, rr_node_p.inc_direction ? 1 : 0);
	} else {
		sprintf(buffer, "%d %s (%d,%d)(%d,%d) dir=%d", inode, rr_types[rr_node_p.type], rr_node_p.xlow, rr_node_p.ylow, rr_node_p.xhigh, rr_node_p.yhigh, rr_node_p.inc_direction ? 1 : 0);
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

#endif
