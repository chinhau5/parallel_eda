#ifndef ROUTE_COMMON_TYPES_H
#define ROUTE_COMMON_TYPES_H

/************ Defines and types shared by all route files ********************/

struct s_heap {
	int index;
	int level;
	float cost;
	union {
		int prev_node;
		struct s_heap *next;
	} u;
	int prev_edge;
	float backward_path_cost;
	float R_upstream;

	bool operator<(const s_heap &other) const {
		return this->cost > other.cost;
	}
};


typedef struct {
	float old_tcost;
	float old_bcost;
	float new_tcost;
	float new_bcost;
	int prev_node;
	int prev_edge;
	int target_node;
	bool written;
} transaction_t;

/* Used by the heap as its fundamental data structure.                      * 
 * index:   Index (ID) of this routing resource node.                       * 
 * cost:    Cost up to and including this node.                             * 
 * u.prev_node:  Index (ID) of the predecessor to this node for             * 
 *          use in traceback.  NO_PREVIOUS if none.                         * 
 * u.next:  pointer to the next s_heap structure in the free                * 
 *          linked list.  Not used when on the heap.                        * 
 * prev_edge:  Index of the edge (between 0 and num_edges-1) used to        *
 *             connect the previous node to this one.  NO_PREVIOUS if       *
 *             there is no previous node.                                   *
 * backward_path_cost:  Used only by the timing-driven router.  The "known" *
 *                      cost of the path up to and including this node.     *
 *                      In this case, the .cost member contains not only    *
 *                      the known backward cost but also an expected cost   *
 *                      to the target.                                      *
 * R_upstream: Used only by the timing-driven router.  Stores the upstream  *
 *             resistance to ground from this node, including the           *
 *             resistance of the node itself (rr_node[index].R).            */

typedef struct {
	int prev_node;
	float path_cost;
	float backward_path_cost;
	short prev_edge;
	short target_flag;
	/* the members below should be moved to rr_node and protected by the same lock
	 * that protects rr_node.occ */
	float pres_cost;
	float acc_cost;
	bool modified;
	std::vector<transaction_t> transactions;
} t_rr_node_route_inf;

/* Extra information about each rr_node needed only during routing (i.e.    *
 * during the maze expansion).                                              *
 *                                                                          *
 * prev_node:  Index of the previous node used to reach this one;           *
 *             used to generate the traceback.  If there is no              *
 *             predecessor, prev_node = NO_PREVIOUS.                        *
 * pres_cost:  Present congestion cost term for this node.                  *
 * acc_cost:   Accumulated cost term from previous Pathfinder iterations.   *
 * path_cost:  Total cost of the path up to and including this node +       *
 *             the expected cost to the target if the timing_driven router  *
 *             is being used.                                               *
 * backward_path_cost:  Total cost of the path up to and including this     *
 *                      node.  Not used by breadth-first router.            *
 * prev_edge:  Index of the edge (from 0 to num_edges-1) that was used      *
 *             to reach this node from the previous node.  If there is      *
 *             no predecessor, prev_edge = NO_PREVIOUS.                     *
 * target_flag:  Is this node a target (sink) for the current routing?      *
 *               Number of times this node must be reached to fully route.  */
#endif
