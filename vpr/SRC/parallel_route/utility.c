#include <assert.h>
#include "vpr_types.h"
#include "utility.h"

static const char *rr_types[] =  {
	"SOURCE", "SINK", "IPIN", "OPIN", "CHANX", "CHANY", "INTRA_CLUSTER_EDGE"
};

/*#define PRINT_RR_NODE*/

static const RRGraph *rr_graph;

void init_sprintf_rr_node(const RRGraph *_rr_graph)
{
	rr_graph = _rr_graph;
}

void sprintf_rr_node_impl(int inode, char *buffer)
{
	/*assert(false);*/
	extern t_rr_node *rr_node;
	extern struct s_grid_tile **grid;

	char s_pin_side[16];
	s_pin_side[0] = '\0';
	if (rr_node[inode].type == IPIN || rr_node[inode].type == OPIN) {
		auto tile = &grid[rr_node[inode].xlow][rr_node[inode].ylow];

		/*int num_pins = 0;*/
		/*for (int height = 0; height < tile->type->height; ++height) {*/
			/*for (int side = 0; side < 4; ++side) {*/
				/*for (int ipin = 0; ipin < tile->type->num_pins; ++ipin) {*/
					/*if (tile->type->pinloc[height][side][ipin] == 1[> && !tile->type->is_global_pin[ipin]<]) {*/
						/*[>printf("is_global_pin: %d\n", tile->type->is_global_pin[ipin]);<]*/
						/*++num_pins;*/
					/*}*/
				/*}*/
			/*}*/
		/*}*/
		/*printf("height: %d\n", tile->type->height);*/
		/*printf("%d %d\n", num_pins, tile->type->num_pins);*/
		/*assert(num_pins == tile->type->num_pins);*/
		/*if (num_pins == tile->type->num_pins) {*/
			/*printf("same recalc num_pins %s\n", tile->type->name);*/
		/*} else {*/
			/*printf("diff recalc num_pins %s\n", tile->type->name);*/
		/*}*/
		int pin_height = tile->type->pin_height[rr_node[inode].ptc_num];
		int num_matches = 0;
		int pin_side = -1;
		for (int side = 0; side < 4; ++side) {
			if (tile->type->pinloc[pin_height][side][rr_node[inode].ptc_num] == 1) {
				pin_side = side;
				++num_matches;
			}
		}
		/*assert(num_matches == 1);*/

		switch (pin_side) {
			case 0:
				sprintf(s_pin_side, "TOP");
				break;
			case 1:
				sprintf(s_pin_side, "RIGHT");
				break;
			case 2:
				sprintf(s_pin_side, "BOTTOM");
				break;
			case 3:
				sprintf(s_pin_side, "LEFT");
				break;
			default:
				break;
		}
	}

	if (rr_node[inode].direction == INC_DIRECTION) {
		sprintf(buffer, "%d %s (%d,%d)(%d,%d) ptc=%d dir=%d pin_side: %s", inode, rr_types[rr_node[inode].type], rr_node[inode].xlow, rr_node[inode].ylow, rr_node[inode].xhigh, rr_node[inode].yhigh, rr_node[inode].ptc_num, rr_node[inode].direction, s_pin_side);
	} else if (rr_node[inode].direction == DEC_DIRECTION) {
		sprintf(buffer, "%d %s (%d,%d)(%d,%d) ptc=%d dir=%d pin_side: %s", inode, rr_types[rr_node[inode].type], rr_node[inode].xhigh, rr_node[inode].yhigh, rr_node[inode].xlow, rr_node[inode].ylow, rr_node[inode].ptc_num, rr_node[inode].direction, s_pin_side);
	} else {
		sprintf(buffer, "%d %s (%d,%d)(%d,%d) ptc=%d dir=%d pin_side: %s", inode, rr_types[rr_node[inode].type], rr_node[inode].xlow, rr_node[inode].ylow, rr_node[inode].xhigh, rr_node[inode].yhigh, rr_node[inode].ptc_num, rr_node[inode].direction, s_pin_side);
	}

	/*const auto &rr_node_p = get_vertex_props(*rr_graph, inode);*/

	/*if (!rr_node_p.inc_direction) {*/
		/*sprintf(buffer, "%d %s (%d,%d)(%d,%d) dir=%d", inode, rr_types[rr_node_p.type], rr_node_p.xhigh, rr_node_p.yhigh, rr_node_p.xlow, rr_node_p.ylow, rr_node_p.inc_direction ? 1 : 0);*/
	/*} else {*/
		/*sprintf(buffer, "%d %s (%d,%d)(%d,%d) dir=%d", inode, rr_types[rr_node_p.type], rr_node_p.xlow, rr_node_p.ylow, rr_node_p.xhigh, rr_node_p.yhigh, rr_node_p.inc_direction ? 1 : 0);*/
	/*}*/

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
