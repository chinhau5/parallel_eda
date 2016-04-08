#include <math.h>
#include <stdio.h>
#include <assert.h>
#include <time.h>
#include <queue>
#include <zlog.h>
#include "util.h"
#include "vpr_types.h"
#include "vpr_utils.h"
/*#include "globals.h"*/
#include "route_export.h"
#include "route_common.h"
#include "route_tree_timing.h"
#include "route_timing.h"
#include "route_breadth_first.h"
#include "place_and_route.h"
#include "rr_graph.h"
#include "rr_graph_util.h"
#include "rr_graph2.h"
#include "read_xml_arch_file.h"
#include "ReadOptions.h"
#include "parallel_route_timing.h"
#include "utility.h"

extern int nx, ny;
extern int num_rr_nodes;
extern t_rr_node *rr_node; /* [0..num_rr_nodes-1]          */
extern t_rr_indexed_data *rr_indexed_data; /* [0 .. num_rr_indexed_data-1] */
extern int num_nets;
extern struct s_net *clb_net;
extern int **net_rr_terminals; /* [0..num_nets-1][0..num_pins-1] */
extern int num_blocks;
extern struct s_block *block;
extern t_type_ptr IO_TYPE;
extern struct s_grid_tile **grid; /* FPGA complex blocks grid [0..nx+1][0..ny+1] */
extern int **rr_blk_source; /* [0..num_blocks-1][0..num_class-1] */

static boolean is_parallel_route;

void print_one_route(char *route_file, int inet, t_net_route *net_route, t_rt_node **l_rr_node_to_rt_node)
{
	/* Prints out the routing to file route_file.  */

	int inode, ipin, bnum, ilow, jlow, node_block_pin, iclass;
	t_rr_type rr_type;
	struct s_trace *tptr;
	const char *name_type[] = { "SOURCE", "SINK", "IPIN", "OPIN", "CHANX", "CHANY",
			"INTRA_CLUSTER_EDGE" };
	FILE *fp;

	fp = fopen(route_file, "w");

	//fprintf(fp, "Array size: %d x %d logic blocks.\n", nx, ny);
	//fprintf(fp, "\nRouting:");
	//fprintf(fp, "\n\nNet %d (%s)\n\n", inet, clb_net[inet].name);
	tptr = net_route->l_trace_head;

	int source_node = tptr->index;
	assert(rr_node[source_node].type == SOURCE);

	int sink = 1;
	while (tptr != NULL) {
		inode = tptr->index;
		rr_type = rr_node[inode].type;

		if (rr_type == CHANX || rr_type == CHANY) {
			if (rr_node[inode].direction == INC_DIRECTION) {
				//fprintf(fp, "Node:\t%d\t%6s (%d,%d) ", inode, name_type[rr_type], rr_node[inode].xlow, rr_node[inode].ylow);
				//fprintf(fp, "to (%d,%d) ", rr_node[inode].xhigh, rr_node[inode].yhigh);
			} else {
				//fprintf(fp, "Node:\t%d\t%6s (%d,%d) ", inode, name_type[rr_type], rr_node[inode].xhigh, rr_node[inode].yhigh);
				//fprintf(fp, "to (%d,%d) ", rr_node[inode].xlow, rr_node[inode].ylow);
			}
		} else {
			ilow = rr_node[inode].xlow;
			jlow = rr_node[inode].ylow;

			//fprintf(fp, "Node:\t%d\t%6s (%d,%d) ", inode, name_type[rr_type], ilow, jlow);

			if ((ilow != rr_node[inode].xhigh)
					|| (jlow != rr_node[inode].yhigh)) {
				//fprintf(fp, "to (%d,%d) ", rr_node[inode].xhigh, rr_node[inode].yhigh);
			}
		}

		switch (rr_type) {

			case IPIN:
			case OPIN:
				if (grid[ilow][jlow].type == IO_TYPE) {
					//fprintf(fp, " Pad: ");
				} else { /* IO Pad. */
					//fprintf(fp, " Pin: ");
				}
				break;

			case CHANX:
			case CHANY:
				//fprintf(fp, " Track: ");
				break;

			case SOURCE:
			case SINK:
				if (grid[ilow][jlow].type == IO_TYPE) {
					//fprintf(fp, " Pad: ");
				} else { /* IO Pad. */
					//fprintf(fp, " Class: ");
				}
				if (rr_type == SINK) {
					fprintf(fp, " Delay: %g ", l_rr_node_to_rt_node[inode]->Tdel);
					++sink;
				}
				break;

			default:
				vpr_printf(TIO_MESSAGE_ERROR, "in print_route: Unexpected traceback element type: %d (%s).\n", 
						rr_type, name_type[rr_type]);
				exit(1);
				break;
		}

		//fprintf(fp, "%d  ", rr_node[inode].ptc_num);

		/* Uncomment line below if you're debugging and want to see the switch types *
		 * used in the routing.                                                      */
		/*          fprintf (fp, "Switch: %d", tptr->iswitch);    */

		//fprintf(fp, "\n");

		if (rr_type == SINK) {
			fprintf(fp, "\n");
		}

		tptr = tptr->next;
	}

	fclose(fp);
}

void print_route(char *route_file, t_net_route *net_route) {

	/* Prints out the routing to file route_file.  */

	int inet, inode, ipin, bnum, ilow, jlow, node_block_pin, iclass;
	t_rr_type rr_type;
	struct s_trace *tptr;
	const char *name_type[] = { "SOURCE", "SINK", "IPIN", "OPIN", "CHANX", "CHANY",
			"INTRA_CLUSTER_EDGE" };
	FILE *fp;

	fp = fopen(route_file, "w");

	fprintf(fp, "Array size: %d x %d logic blocks.\n", nx, ny);
	fprintf(fp, "\nRouting:");
	for (inet = 0; inet < num_nets; inet++) {
		if (clb_net[inet].is_global == FALSE) {
			if (clb_net[inet].num_sinks == FALSE) {
				fprintf(fp, "\n\nNet %d (%s)\n\n", inet, clb_net[inet].name);
				fprintf(fp, "\n\nUsed in local cluster only, reserved one CLB pin\n\n");
			} else {
				fprintf(fp, "\n\nNet %d (%s)\n\n", inet, clb_net[inet].name);
				tptr = net_route[inet].l_trace_head;

				int source_node = tptr->index;
				assert(rr_node[source_node].type == SOURCE);

				int sink = 1;
				while (tptr != NULL) {
					inode = tptr->index;
					rr_type = rr_node[inode].type;

					if (rr_type == CHANX || rr_type == CHANY) {
						if (rr_node[inode].direction == INC_DIRECTION) {
							fprintf(fp, "Node:\t%d\t%6s (%d,%d) ", inode, name_type[rr_type], rr_node[inode].xlow, rr_node[inode].ylow);
							fprintf(fp, "to (%d,%d) ", rr_node[inode].xhigh, rr_node[inode].yhigh);
						} else {
							fprintf(fp, "Node:\t%d\t%6s (%d,%d) ", inode, name_type[rr_type], rr_node[inode].xhigh, rr_node[inode].yhigh);
							fprintf(fp, "to (%d,%d) ", rr_node[inode].xlow, rr_node[inode].ylow);
						}
					} else {
						ilow = rr_node[inode].xlow;
						jlow = rr_node[inode].ylow;

						fprintf(fp, "Node:\t%d\t%6s (%d,%d) ", inode, name_type[rr_type], ilow, jlow);

						if ((ilow != rr_node[inode].xhigh)
								|| (jlow != rr_node[inode].yhigh))
							fprintf(fp, "to (%d,%d) ", rr_node[inode].xhigh,
									rr_node[inode].yhigh);
					}

					int current_sink_node = net_rr_terminals[inet][net_route[inet].sink_order[sink]];

					std::pair<int, int> source_pos = get_node_start(source_node); 
					std::pair<int, int> sink_pos = get_node_start(current_sink_node); 
					std::pair<int, int> current_pos = get_node_start(inode); 
					std::pair<int, int> neighbor_pos = tptr->next ? get_node_start(tptr->next->index) : current_pos;

					if (rr_type == OPIN) {
						current_pos = neighbor_pos;
					} else if (rr_type == SINK) {
						neighbor_pos = current_pos;
					}

					/*Interval<int> hor(source_pos.first, sink_pos.first);*/
					/*Interval<int> vert(source_pos.second, sink_pos.second);*/

					switch (rr_type) {

					case IPIN:
					case OPIN:
						if (grid[ilow][jlow].type == IO_TYPE) {
							fprintf(fp, " Pad: ");
						} else { /* IO Pad. */
							fprintf(fp, " Pin: ");
						}
						break;

					case CHANX:
					case CHANY:
						fprintf(fp, " Track: ");
						break;

					case SOURCE:
					case SINK:
						if (grid[ilow][jlow].type == IO_TYPE) {
							fprintf(fp, " Pad: ");
						} else { /* IO Pad. */
							fprintf(fp, " Class: ");
						}
						if (rr_type == SINK) {
							assert(current_sink_node == inode);
							++sink;
						}
						break;

					default:
						vpr_printf(TIO_MESSAGE_ERROR, "in print_route: Unexpected traceback element type: %d (%s).\n", 
								rr_type, name_type[rr_type]);
						exit(1);
						break;
					}

					fprintf(fp, "%d  ", rr_node[inode].ptc_num);

					if (rr_node[inode].occ > rr_node[inode].capacity) {
						fprintf(fp, "congested ");
					}

					int current_to_sink_manhattan_distance = abs(current_pos.first - sink_pos.first) + abs(current_pos.second - sink_pos.second);
					int neighbor_to_sink_manhattan_distance = abs(neighbor_pos.first - sink_pos.first) + abs(neighbor_pos.second - sink_pos.second);

					if (neighbor_to_sink_manhattan_distance > current_to_sink_manhattan_distance) {
						fprintf(fp, "Neighbor is non-monotonic [%d > %d]", neighbor_to_sink_manhattan_distance, current_to_sink_manhattan_distance);
					}
					if (rr_type == SINK) {
						fprintf(fp, "\n");
					}

					/* Uncomment line below if you're debugging and want to see the switch types *
					 * used in the routing.                                                      */
					/*          fprintf (fp, "Switch: %d", tptr->iswitch);    */

					fprintf(fp, "\n");

					tptr = tptr->next;
				}
			}
		}

		else { /* Global net.  Never routed. */
			fprintf(fp, "\n\nNet %d (%s): global net connecting:\n\n", inet,
					clb_net[inet].name);

			for (ipin = 0; ipin <= clb_net[inet].num_sinks; ipin++) {
				bnum = clb_net[inet].node_block[ipin];

				node_block_pin = clb_net[inet].node_block_pin[ipin];
				iclass = block[bnum].type->pin_class[node_block_pin];

				fprintf(fp, "Block %s (#%d) at (%d, %d), Pin class %d.\n",
						block[bnum].name, bnum, block[bnum].x, block[bnum].y,
						iclass);
			}
		}
	}

	fclose(fp);

	/*if (getEchoEnabled() && isEchoFileEnabled(E_ECHO_MEM)) {*/
		/*fp = my_fopen(getEchoFileName(E_ECHO_MEM), "w", 0);*/
		/*fprintf(fp, "\nNum_heap_allocated: %d   Num_trace_allocated: %d\n",*/
				/*num_heap_allocated, num_trace_allocated);*/
		/*fprintf(fp, "Num_linked_f_pointer_allocated: %d\n",*/
				/*num_linked_f_pointer_allocated);*/
		/*fclose(fp);*/
	/*}*/
}

void set_parallel_route(boolean _is_parallel_route)
{
	is_parallel_route = _is_parallel_route;
}

void thread_safe_pathfinder_update_one_cost(int inet, struct s_trace *route_segment_start,
		int add_or_sub, float pres_fac, int thread_id, int sub_iter, bool set_occ_by_thread) {

	/* This routine updates the occupancy and pres_cost of the rr_nodes that are *
	 * affected by the portion of the routing of one net that starts at          *
	 * route_segment_start.  If route_segment_start is trace_head[inet], the     *
	 * cost of all the nodes in the routing of net inet are updated.  If         *
	 * add_or_sub is -1 the net (or net portion) is ripped up, if it is 1 the    *
	 * net is added to the routing.  The size of pres_fac determines how severly *
	 * oversubscribed rr_nodes are penalized.                                    */

	struct s_trace *tptr;
	int inode, occ, capacity;
	extern zlog_category_t *route_inner_log;

	tptr = route_segment_start;
	if (tptr == NULL) /* No routing yet. */
		return;

/*	assert(!pthread_mutex_lock(&clb_net[inet].lock));*/
    
	for (;;) {
		inode = tptr->index;

		if (set_occ_by_thread) {
			rr_node[inode].occ_by_thread[sub_iter][thread_id] = add_or_sub;
		}

		assert(!pthread_mutex_lock(&rr_node[inode].lock));

		if (add_or_sub > 0) {
			if (rr_node[inode].occupant_net_id.find(inet) != rr_node[inode].occupant_net_id.end()) {
				zlog_warn(route_inner_log, "Trying to insert existing net %d into node %d\n", inet, inode);
			}
			rr_node[inode].occupant_net_id.insert(inet);
		} else {
			assert(add_or_sub < 0);
			if (rr_node[inode].occupant_net_id.find(inet) == rr_node[inode].occupant_net_id.end()) {
				zlog_warn(route_inner_log, "Trying to erase non-existing net %d from node %d\n", inet, inode);
			}
			rr_node[inode].occupant_net_id.erase(inet);
		}
		occ = rr_node[inode].occ + add_or_sub;
		if (occ < 0) {
			zlog_fatal(route_inner_log, "occ is less than 0\n");
			assert(false);
		}
		capacity = rr_node[inode].capacity;

		rr_node[inode].occ = occ;

		/* pres_cost is Pn in the Pathfinder paper. I set my pres_cost according to *
		 * the overuse that would result from having ONE MORE net use this routing  *
		 * node.                                                                    */

		if (occ < capacity) {
			rr_node[inode].pres_cost = 1.;
		} else {
			rr_node[inode].pres_cost = 1.
					+ (occ + 1 - capacity) * pres_fac;
		}

		zlog_debug(route_inner_log, "Update cost for %d\n", inode);
    
		assert(!pthread_mutex_unlock(&rr_node[inode].lock));

		if (rr_node[inode].type == SINK) {
			tptr = tptr->next; /* Skip next segment. */
			if (tptr == NULL)
				break;
		}

		tptr = tptr->next;

	} /* End while loop -- did an entire traceback. */
	
/*	assert(!pthread_mutex_unlock(&clb_net[inet].lock));*/
}

void thread_safe_pathfinder_update_cost(float pres_fac, float acc_fac) {

	/* This routine recomputes the pres_cost and acc_cost of each routing        *
	 * resource for the pathfinder algorithm after all nets have been routed.    *
	 * It updates the accumulated cost to by adding in the number of extra       *
	 * signals sharing a resource right now (i.e. after each complete iteration) *
	 * times acc_fac.  It also updates pres_cost, since pres_fac may have        *
	 * changed.  THIS ROUTINE ASSUMES THE OCCUPANCY VALUES IN RR_NODE ARE UP TO  *
	 * DATE.                                                                     */

	int inode, occ, capacity;

	for (inode = 0; inode < num_rr_nodes; inode++) {
		occ = rr_node[inode].occ;
		capacity = rr_node[inode].capacity;

		if (occ > capacity) {
			rr_node[inode].acc_cost += (occ - capacity) * acc_fac;
			rr_node[inode].pres_cost = 1.
					+ (occ + 1 - capacity) * pres_fac;
		}

		/* If occ == capacity, we don't need to increase acc_cost, but a change    *
		 * in pres_fac could have made it necessary to recompute the cost anyway.  */

		else if (occ == capacity) {
			rr_node[inode].pres_cost = 1. + pres_fac;
		}
	}
}

struct s_trace *
thread_safe_update_traceback(t_trace **l_trace_head, t_trace **l_trace_tail, const struct s_heap *hptr, t_rr_node_route_inf *l_rr_node_route_inf) {

	/* This routine adds the most recently finished wire segment to the         *
	 * traceback linked list.  The first connection starts with the net SOURCE  *
	 * and begins at the structure pointed to by l_trace_head[inet]. Each         *
	 * connection ends with a SINK.  After each SINK, the next connection       *
	 * begins (if the net has more than 2 pins).  The first element after the   *
	 * SINK gives the routing node on a previous piece of the routing, which is *
	 * the link from the existing net to this new piece of the net.             *
	 * In each traceback I start at the end of a path and trace back through    *
	 * its predecessors to the beginning.  I have stored information on the     *
	 * predecesser of each node to make traceback easy -- this sacrificies some *
	 * memory for easier code maintenance.  This routine returns a pointer to   *
	 * the first "new" node in the traceback (node not previously in trace).    */

	struct s_trace *tptr, *prevptr, *temptail, *ret_ptr;
	int inode;
	short iedge;

#ifdef DEBUG
	t_rr_type rr_type;
#endif

	inode = hptr->index;

#ifdef DEBUG
	rr_type = rr_node[inode].type;
	if (rr_type != SINK) {
		vpr_printf(TIO_MESSAGE_ERROR, "in update_traceback. Expected type = SINK (%d).\n", SINK);
		vpr_printf(TIO_MESSAGE_ERROR, "\tGot type = %d while tracing back net %d.\n", rr_type, -1);
		exit(1);
	}
#endif

	tptr = new t_trace; //alloc_trace_data(); /* SINK on the end of the connection */
	tptr->index = inode;
	tptr->iswitch = OPEN;
	tptr->next = NULL;
	temptail = tptr; /* This will become the new tail at the end */
	/* of the routine.                          */

	extern zlog_category_t *route_inner_log;
	char buffer[256];
	sprintf_rr_node(inode, buffer);
	zlog_debug(route_inner_log, "Traceback: %s\n", buffer);

	/* Now do it's predecessor. */

	inode = hptr->u.prev_node;
	iedge = hptr->prev_edge;

	/* the asserts below are not true because l_rr_node_route_inf is not updated when target_node is reached in the while (inode != target_node) loop */
	/*assert(inode == l_rr_node_route_inf[hptr->index].prev_node);*/
	/*assert(iedge == l_rr_node_route_inf[hptr->index].prev_edge);*/

	int prev_inode = tptr->index;

	while (inode != NO_PREVIOUS) {
		sprintf_rr_node(inode, buffer);
		zlog_debug(route_inner_log, "Traceback: %s\n", buffer);

		prevptr = new t_trace;//alloc_trace_data();
		prevptr->index = inode;
		prevptr->iswitch = rr_node[inode].switches[iedge];
		prevptr->next = tptr;

		tptr = prevptr;

		prev_inode = inode;

		iedge = l_rr_node_route_inf[inode].prev_edge;
		inode = l_rr_node_route_inf[inode].prev_node;
	}

	/*assert(rr_node[prev_inode].type == SOURCE);*/

	if (*l_trace_tail != NULL) {
		(*l_trace_tail)->next = tptr; /* Traceback ends with tptr */
		ret_ptr = tptr->next; /* First new segment.       */
	} else { /* This was the first "chunk" of the net's routing */
		*l_trace_head = tptr;
		ret_ptr = tptr; /* Whole traceback is new. */
	}

	*l_trace_tail = temptail;
	return (ret_ptr);
}

float get_weighted_rr_cong_cost(int inode, int num_sinks) {

	/* Returns the *congestion* cost of using this rr_node. */

	short cost_index;
	float cost;

	cost_index = rr_node[inode].cost_index;
/*	cost = rr_indexed_data[cost_index].saved_base_cost * sqrt((float)num_sinks)*/
	assert(!pthread_mutex_lock(&rr_node[inode].lock));
	cost = rr_indexed_data[cost_index].saved_base_cost
			* rr_node[inode].acc_cost
			* rr_node[inode].weighted_pres_cost;
	assert(!pthread_mutex_unlock(&rr_node[inode].lock));
	return (cost);
}

float get_rr_cong_cost(int inode, int num_sinks) {

	/* Returns the *congestion* cost of using this rr_node. */

	short cost_index;
	float cost;

	cost_index = rr_node[inode].cost_index;
/*	cost = rr_indexed_data[cost_index].saved_base_cost * sqrt((float)num_sinks)*/
	assert(!pthread_mutex_lock(&rr_node[inode].lock));
	cost = rr_indexed_data[cost_index].saved_base_cost
			* rr_node[inode].acc_cost
			* rr_node[inode].pres_cost;
	assert(!pthread_mutex_unlock(&rr_node[inode].lock));
	return (cost);
}

void node_to_heap(int thread_id,
		int inode, int prev_node, int prev_edge,
		float cost, float backward_path_cost, float R_upstream,
		const t_rr_node_route_inf *l_rr_node_route_inf,
		std::priority_queue<struct s_heap> &heap)
{

	/* Puts an rr_node on the heap, if the new cost given is lower than the     *
	 * current path_cost to this channel segment.  The index of its predecessor *
	 * is stored to make traceback easy.  The index of the edge used to get     *
	 * from its predecessor to it is also stored to make timing analysis, etc.  *
	 * easy.  The backward_path_cost and R_upstream values are used only by the *
	 * timing-driven router -- the breadth-first router ignores them.           */

	struct s_heap item;

	char buffer[256];
	sprintf_rr_node(inode, buffer);

	/*if (cost >= l_rr_node_route_inf[inode].path_cost) {*/
		/*extern zlog_category_t *route_inner_log;*/
		/*zlog_debug(route_inner_log, "\t\t[%d] Not adding %s to heap because cost %g is >= %g\n", thread_id, buffer, cost, l_rr_node_route_inf[inode].path_cost);*/
		/*return;*/
	/*}*/

	item.index = inode;
	item.u.prev_node = prev_node;
	item.prev_edge = prev_edge;
	item.cost = cost;
	item.backward_path_cost = backward_path_cost;
	item.R_upstream = R_upstream;

	extern zlog_category_t *route_inner_log;
	zlog_debug(route_inner_log, "\t\t[%d] Adding node: %s cost: %g backward_path_cost: %g prev_node: %d to heap\n", thread_id, buffer, cost, backward_path_cost, prev_node);

	heap.push(item);
}

void node_to_heap(
		int inode, int prev_node, int prev_edge,
		float cost, float backward_path_cost, float R_upstream,
		const t_rr_node_route_inf *l_rr_node_route_inf,
		std::priority_queue<struct s_heap> &heap)
{

	/* Puts an rr_node on the heap, if the new cost given is lower than the     *
	 * current path_cost to this channel segment.  The index of its predecessor *
	 * is stored to make traceback easy.  The index of the edge used to get     *
	 * from its predecessor to it is also stored to make timing analysis, etc.  *
	 * easy.  The backward_path_cost and R_upstream values are used only by the *
	 * timing-driven router -- the breadth-first router ignores them.           */

	struct s_heap item;

	if (cost >= l_rr_node_route_inf[inode].path_cost) {
		extern zlog_category_t *route_inner_log;
		zlog_debug(route_inner_log, "Not adding %d to heap because cost %g is >= %g\n", inode, cost, l_rr_node_route_inf[inode].path_cost);
		return;
	}

	item.index = inode;
	item.u.prev_node = prev_node;
	item.prev_edge = prev_edge;
	item.cost = cost;
	item.backward_path_cost = backward_path_cost;
	item.R_upstream = R_upstream;

	extern zlog_category_t *route_inner_log;
	char buffer[256];
	sprintf_rr_node(inode, buffer);
	zlog_debug(route_inner_log, "Adding node: %s cost: %g backward_path_cost: %g prev_node: %d to heap\n", buffer, cost, backward_path_cost, prev_node);

	heap.push(item);
}

void thread_safe_free_traceback(t_trace **l_trace_head, t_trace **l_trace_tail) {

	/* Puts the entire traceback (old routing) for this net on the free list *
	 * and sets the trace_head pointers etc. for the net to NULL.            */

	struct s_trace *tptr, *tempptr;

	tptr = *l_trace_head;

	while (tptr != NULL) {
		tempptr = tptr->next;
		delete tptr;
		tptr = tempptr;
	}

	*l_trace_head = NULL;
	*l_trace_tail = NULL;
}

void thread_safe_mark_ends(int inet, t_rr_node_route_inf *l_rr_node_route_inf) {

	/* Mark all the SINKs of this net as targets by setting their target flags  *
	 * to the number of times the net must connect to each SINK.  Note that     *
	 * this number can occassionally be greater than 1 -- think of connecting   *
	 * the same net to two inputs of an and-gate (and-gate inputs are logically *
	 * equivalent, so both will connect to the same SINK).                      */

	int ipin, inode;

	for (ipin = 1; ipin <= clb_net[inet].num_sinks; ipin++) {
		inode = net_rr_terminals[inet][ipin];
		l_rr_node_route_inf[inode].target_flag++;
	}
}

static void adjust_one_rr_occ_and_pcost(int inet, int inode, int add_or_sub,
		float pres_fac) {

	/* Increments or decrements (depending on add_or_sub) the occupancy of    *
	 * one rr_node, and adjusts the present cost of that node appropriately.  */

	int occ, capacity;

	occ = rr_node[inode].occ + add_or_sub;
	capacity = rr_node[inode].capacity;
	rr_node[inode].occ = occ;
	rr_node[inode].num_reservation += add_or_sub;

	if (occ < capacity) {
		rr_node[inode].pres_cost = 1.;
	} else {
		rr_node[inode].pres_cost = 1.
				+ (occ + 1 - capacity) * pres_fac;
	}
}

/* TODO: check if this is still necessary for speed */
void thread_safe_reserve_locally_used_opins(float pres_fac, boolean rip_up_local_opins,
		t_ivec ** clb_opins_used_locally, int inet) {

	/* In the past, this function implicitly allowed LUT duplication when there are free LUTs.
	 This was especially important for logical equivalence; however, now that we have a very general
	 logic cluster, it does not make sense to allow LUT duplication implicitly. we'll need to look into how we want to handle this case

	 */

	int iblk, num_local_opin, inode, from_node, iconn, num_edges, to_node;
	int iclass, ipin;
	float cost;
	const struct s_heap *heap_head_ptr;
	t_type_ptr type;

	if (rip_up_local_opins) {
		for (iblk = 0; iblk < num_blocks; iblk++) {
			type = block[iblk].type;
			for (iclass = 0; iclass < type->num_class; iclass++) {
				num_local_opin = clb_opins_used_locally[iblk][iclass].nelem;
				/* Always 0 for pads and for RECEIVER (IPIN) classes */
				for (ipin = 0; ipin < num_local_opin; ipin++) {
					inode = clb_opins_used_locally[iblk][iclass].list[ipin];
					adjust_one_rr_occ_and_pcost(inet, inode, -1, pres_fac);
				}
			}
		}
	}

	std::priority_queue<struct s_heap> heap;

	for (iblk = 0; iblk < num_blocks; iblk++) {
		type = block[iblk].type;
		for (iclass = 0; iclass < type->num_class; iclass++) {
			num_local_opin = clb_opins_used_locally[iblk][iclass].nelem;
			/* Always 0 for pads and for RECEIVER (IPIN) classes */

			if (num_local_opin != 0) { /* Have to reserve (use) some OPINs */
				from_node = rr_blk_source[iblk][iclass];
				num_edges = rr_node[from_node].num_edges;
				for (iconn = 0; iconn < num_edges; iconn++) {
					to_node = rr_node[from_node].edges[iconn];
					cost = get_rr_cong_cost(to_node, 1);

					struct s_heap item;
					item.index = to_node;
					item.cost = cost;
					item.u.prev_node = OPEN;
					item.prev_edge = OPEN;
					item.backward_path_cost = 0;
					item.R_upstream = 0;
					dzlog_debug("Adding node: %d cost: %g backward_path_cost: %g prev_node: %d to heap\n", to_node, cost, 0.0, OPEN);
					heap.push(item);
				}

				for (ipin = 0; ipin < num_local_opin; ipin++) {
					assert(!heap.empty());
					heap_head_ptr = &heap.top();
					inode = heap_head_ptr->index;
					heap.pop();
					adjust_one_rr_occ_and_pcost(inet, inode, 1, pres_fac);
					clb_opins_used_locally[iblk][iclass].list[ipin] = inode;
				}

				heap = std::priority_queue<struct s_heap>();
			}
		}
	}
}


