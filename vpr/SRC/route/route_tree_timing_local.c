#include <stdio.h>
#include <zlog.h>
#include <vector>
#include "util.h"
#include "vpr_types.h"
#include "route_common.h"
#include "route_tree_timing.h"
#include "parallel_route_timing.h"
#include "utility.h"

extern t_rr_node *rr_node; /* [0..num_rr_nodes-1]          */
extern struct s_switch_inf *switch_inf; /* [0..det_routing_arch.num_switch-1] */

static t_rt_node *
add_path_to_route_tree(const struct s_heap *sink, t_rt_node ** sink_rt_node_ptr, const t_rr_node_route_inf *l_rr_node_route_inf, t_rt_node **l_rr_node_to_rt_node, std::vector<int> &modified) {

	/* Adds the most recent wire segment, ending at the SINK indicated by sink, *
	 * to the routing tree.  It returns the first (most upstream) new rt_node,  *
	 * and (via a pointer) the rt_node of the new SINK.                         */

	int inode, remaining_connections_to_sink, no_route_throughs;
	short iedge, iswitch;
	float C_downstream;
	t_rt_node *rt_node, *downstream_rt_node, *sink_rt_node;
	t_linked_rt_edge *linked_rt_edge;

	inode = sink->index;

#ifdef DEBUG
	if (rr_node[inode].type != SINK) {
		vpr_printf(TIO_MESSAGE_ERROR, "in add_path_to_route_tree. Expected type = SINK (%d).\n", SINK);
		vpr_printf(TIO_MESSAGE_INFO, "Got type = %d.", rr_node[inode].type);
		exit(1);
	}
#endif

	remaining_connections_to_sink = l_rr_node_route_inf[inode].target_flag;
	sink_rt_node = new t_rt_node; //alloc_rt_node();
	sink_rt_node->u.child_list = NULL;
	sink_rt_node->inode = inode;
	C_downstream = rr_node[inode].C;
	sink_rt_node->C_downstream = C_downstream;
	l_rr_node_to_rt_node[inode] = sink_rt_node;
	modified.push_back(inode);

	char buffer[256];
	extern zlog_category_t *route_inner_log;
	sprintf_rr_node(inode, buffer);
	zlog_debug(route_inner_log, "Adding %s to route tree\n", buffer);

	/* In the code below I'm marking SINKs and IPINs as not to be re-expanded.  *
	 * Undefine NO_ROUTE_THROUGHS if you want route-throughs or ipin doglegs.   *
	 * It makes the code more efficient (though not vastly) to prune this way   *
	 * when there aren't route-throughs or ipin doglegs.                        */

#define NO_ROUTE_THROUGHS 1	/* Can't route through unused CLB outputs */
	no_route_throughs = 1;
	if (no_route_throughs == 1)
		sink_rt_node->re_expand = FALSE;
	else {
		if (remaining_connections_to_sink == 0) { /* Usual case */
			sink_rt_node->re_expand = TRUE;
		}

		/* Weird case.  This net connects several times to the same SINK.  Thus I   *
		 * can't re_expand this node as part of the partial routing for subsequent  *
		 * connections, since I need to reach it again via another path.            */

		else {
			sink_rt_node->re_expand = FALSE;
		}
	}

	/* Now do it's predecessor. */

	downstream_rt_node = sink_rt_node;
	inode = sink->u.prev_node;
	iedge = sink->prev_edge;
	iswitch = rr_node[inode].switches[iedge];

	/* For all "new" nodes in the path */

	while (l_rr_node_route_inf[inode].prev_node != NO_PREVIOUS) {
		linked_rt_edge = new t_linked_rt_edge; //alloc_linked_rt_edge();
		linked_rt_edge->child = downstream_rt_node;
		linked_rt_edge->iswitch = iswitch;
		linked_rt_edge->next = NULL;

		rt_node = new t_rt_node; //alloc_rt_node();
		downstream_rt_node->parent_node = rt_node;
		downstream_rt_node->parent_switch = iswitch;

		rt_node->u.child_list = linked_rt_edge;
		rt_node->inode = inode;

		if (switch_inf[iswitch].buffered == FALSE)
			C_downstream += rr_node[inode].C;
		else
			C_downstream = rr_node[inode].C;

		rt_node->C_downstream = C_downstream;
		l_rr_node_to_rt_node[inode] = rt_node;
		modified.push_back(inode);

		sprintf_rr_node(inode, buffer);
		zlog_debug(route_inner_log, "Adding %s to route tree\n", buffer);

		if (no_route_throughs == 1)
			if (rr_node[inode].type == IPIN)
				rt_node->re_expand = FALSE;
			else
				rt_node->re_expand = TRUE;

		else {
			if (remaining_connections_to_sink == 0) { /* Normal case */
				rt_node->re_expand = TRUE;
			} else { /* This is the IPIN before a multiply-connected SINK */
				rt_node->re_expand = FALSE;

				/* Reset flag so wire segments get reused */

				remaining_connections_to_sink = 0;
			}
		}

		downstream_rt_node = rt_node;
		iedge = l_rr_node_route_inf[inode].prev_edge;
		inode = l_rr_node_route_inf[inode].prev_node;
		iswitch = rr_node[inode].switches[iedge];
	}

	/* Inode is the join point to the old routing */

	rt_node = l_rr_node_to_rt_node[inode];

	linked_rt_edge = new t_linked_rt_edge; //alloc_linked_rt_edge();
	linked_rt_edge->child = downstream_rt_node;
	linked_rt_edge->iswitch = iswitch;
	linked_rt_edge->next = rt_node->u.child_list;
	rt_node->u.child_list = linked_rt_edge;

	downstream_rt_node->parent_node = rt_node;
	downstream_rt_node->parent_switch = iswitch;

	*sink_rt_node_ptr = sink_rt_node;
	return (downstream_rt_node);
}

static void load_new_path_R_upstream(t_rt_node * start_of_new_path_rt_node) {

	/* Sets the R_upstream values of all the nodes in the new path to the       *
	 * correct value.                                                           */

	float R_upstream;
	int inode;
	short iswitch;
	t_rt_node *rt_node, *parent_rt_node;
	t_linked_rt_edge *linked_rt_edge;

	rt_node = start_of_new_path_rt_node;
	iswitch = rt_node->parent_switch;
	inode = rt_node->inode;
	parent_rt_node = rt_node->parent_node;

	R_upstream = switch_inf[iswitch].R + rr_node[inode].R;

	if (switch_inf[iswitch].buffered == FALSE)
		R_upstream += parent_rt_node->R_upstream;

	rt_node->R_upstream = R_upstream;

	/* Note:  the traversal below makes use of the fact that this new path      *
	 * really is a path (not a tree with branches) to do a traversal without    *
	 * recursion, etc.                                                          */

	linked_rt_edge = rt_node->u.child_list;

	while (linked_rt_edge != NULL) { /* While SINK not reached. */

#ifdef DEBUG
		if (linked_rt_edge->next != NULL) {
			vpr_printf(TIO_MESSAGE_ERROR, "in load_new_path_R_upstream: new routing addition is a tree (not a path).\n");
			exit(1);
		}
#endif

		rt_node = linked_rt_edge->child;
		iswitch = linked_rt_edge->iswitch;
		inode = rt_node->inode;

		if (switch_inf[iswitch].buffered)
			R_upstream = switch_inf[iswitch].R + rr_node[inode].R;
		else
			R_upstream += switch_inf[iswitch].R + rr_node[inode].R;

		rt_node->R_upstream = R_upstream;
		linked_rt_edge = rt_node->u.child_list;
	}
}

static void load_rt_subtree_Tdel(t_rt_node * subtree_rt_root, float Tarrival) {

	/* Updates the Tdel values of the subtree rooted at subtree_rt_root by      *
	 * by calling itself recursively.  The C_downstream values of all the nodes *
	 * must be correct before this routine is called.  Tarrival is the time at  *
	 * at which the signal arrives at this node's *input*.                      */

	int inode;
	short iswitch;
	t_rt_node *child_node;
	t_linked_rt_edge *linked_rt_edge;
	float Tdel, Tchild;

	inode = subtree_rt_root->inode;

	/* Assuming the downstream connections are, on average, connected halfway    *
	 * along a wire segment's length.  See discussion in net_delay.c if you want *
	 * to change this.                                                           */

	Tdel = Tarrival + 0.5 * subtree_rt_root->C_downstream * rr_node[inode].R;
	subtree_rt_root->Tdel = Tdel;

	/* Now expand the children of this node to load their Tdel values (depth-   *
	 * first pre-order traversal).                                              */

	linked_rt_edge = subtree_rt_root->u.child_list;

	while (linked_rt_edge != NULL) {
		iswitch = linked_rt_edge->iswitch;
		child_node = linked_rt_edge->child;

		Tchild = Tdel + switch_inf[iswitch].R * child_node->C_downstream;
		Tchild += switch_inf[iswitch].Tdel; /* Intrinsic switch delay. */
		load_rt_subtree_Tdel(child_node, Tchild);

		linked_rt_edge = linked_rt_edge->next;
	}
}

static t_rt_node *
update_unbuffered_ancestors_C_downstream(t_rt_node * start_of_new_path_rt_node) {

	/* Updates the C_downstream values for the ancestors of the new path.  Once *
	 * a buffered switch is found amongst the ancestors, no more ancestors are  *
	 * affected.  Returns the root of the "unbuffered subtree" whose Tdel       *
	 * values are affected by the new path's addition.                          */

	t_rt_node *rt_node, *parent_rt_node;
	short iswitch;
	float C_downstream_addition;

	rt_node = start_of_new_path_rt_node;
	C_downstream_addition = rt_node->C_downstream;
	parent_rt_node = rt_node->parent_node;
	iswitch = rt_node->parent_switch;

	while (parent_rt_node != NULL && switch_inf[iswitch].buffered == FALSE) {
		rt_node = parent_rt_node;
		rt_node->C_downstream += C_downstream_addition;
		parent_rt_node = rt_node->parent_node;
		iswitch = rt_node->parent_switch;
	}

	return (rt_node);
}

t_rt_node *
thread_safe_update_route_tree(const struct s_heap * sink, const t_rr_node_route_inf *l_rr_node_route_inf, t_rt_node **l_rr_node_to_rt_node, std::vector<int> &modified)
{
	/* Adds the most recently finished wire segment to the routing tree, and    *
	 * updates the Tdel, etc. numbers for the rest of the routing tree.  sink   *
	 * is the heap pointer of the SINK that was reached.  This routine returns  *
	 * a pointer to the rt_node of the SINK that it adds to the routing.        */

	t_rt_node *start_of_new_path_rt_node, *sink_rt_node;
	t_rt_node *unbuffered_subtree_rt_root, *subtree_parent_rt_node;
	float Tdel_start;
	short iswitch;

	start_of_new_path_rt_node = add_path_to_route_tree(sink, &sink_rt_node, l_rr_node_route_inf, l_rr_node_to_rt_node, modified);
	load_new_path_R_upstream(start_of_new_path_rt_node);
	unbuffered_subtree_rt_root = update_unbuffered_ancestors_C_downstream(
			start_of_new_path_rt_node);

	subtree_parent_rt_node = unbuffered_subtree_rt_root->parent_node;

	if (subtree_parent_rt_node != NULL) { /* Parent exists. */
		Tdel_start = subtree_parent_rt_node->Tdel;
		iswitch = unbuffered_subtree_rt_root->parent_switch;
		Tdel_start += switch_inf[iswitch].R
				* unbuffered_subtree_rt_root->C_downstream;
		Tdel_start += switch_inf[iswitch].Tdel;
	} else { /* Subtree starts at SOURCE */
		Tdel_start = 0.;
	}

	load_rt_subtree_Tdel(unbuffered_subtree_rt_root, Tdel_start);

	return (sink_rt_node);
}

t_rt_node *
init_route_tree_to_source(int source_inode, t_rt_node **l_rr_node_to_rt_node, std::vector<int> &modified) {

	/* Initializes the routing tree to just the net source, and returns the root *
	 * node of the rt_tree (which is just the net source).                       */

	t_rt_node *rt_root;

	rt_root = new t_rt_node;
	rt_root->u.child_list = NULL;
	rt_root->parent_node = NULL;
	rt_root->parent_switch = OPEN;
	rt_root->re_expand = TRUE;
	rt_root->inode = source_inode;
	rt_root->C_downstream = rr_node[source_inode].C;
	rt_root->R_upstream = rr_node[source_inode].R;
	rt_root->Tdel = 0.5 * rr_node[source_inode].R * rr_node[source_inode].C;

	l_rr_node_to_rt_node[source_inode] = rt_root;
	modified.push_back(source_inode);

	return (rt_root);
}

void alloc_per_thread_route_tree_timing_structs(void) {

	/* Allocates any structures needed to build the routing trees. */
/**/
/*	if (rr_node_to_rt_node != NULL || rt_node_free_list != NULL*/
/*			|| rt_node_free_list != NULL) {*/
/*		vpr_printf(TIO_MESSAGE_ERROR, "in alloc_route_tree_timing_structs: old structures already exist.\n");*/
/*		exit(1);*/
/*	}*/
/**/
/*	rr_node_to_rt_node = (t_rt_node ***) my_malloc(*/
/*			get_num_threads() * sizeof(t_rt_node **));*/
/**/
/*	for (int i = 0; i < get_num_threads(); ++i) {*/
/*		rr_node_to_rt_node[i] = (t_rt_node **) my_malloc(*/
/*				num_rr_nodes * sizeof(t_rt_node *));*/
/*	}*/
}
