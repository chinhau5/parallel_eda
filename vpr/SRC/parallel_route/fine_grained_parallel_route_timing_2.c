#include <stdio.h>
#include <math.h>
#include <errno.h>
#include <time.h>
#include <sys/time.h>
/*#include <mach/mach_time.h>*/
#include <string.h>
#include <assert.h>
#include <pthread.h>
#include <fcntl.h>
#include <unistd.h>
#include <semaphore.h>
#include <queue>
#include <mutex>
#include <random>
#include <algorithm>
#include <functional>
#include <boost/timer/timer.hpp>
#include <zlog.h>
#include "tbb/tbb.h"
#include "tbb/task_scheduler_init.h"
#include "util.h"
#include "vpr_types.h"
/*#include "globals.h"*/
//#include "route_export.h"
#include "route_common_local.h"
#include "route_tree_timing_local.h"
/*#include "route_timing.h"*/
#include "heapsort.h"
#include "path_delay.h"
#include "net_delay.h"
#include "stats.h"
#include "ReadOptions.h"
#include "rr_graph_util.h"
#include "rr_graph2.h"
#include "rr_graph.h"
#include "parallel_route_timing.h"
#include "barrier.h"
#include "bounding_box.h"
#include "utility.h"
using namespace std;

extern int num_nets;
extern struct s_net *clb_net;
extern int nx, ny;
extern int num_rr_nodes;
extern t_rr_node *rr_node; /* [0..num_rr_nodes-1]          */
extern t_rr_indexed_data *rr_indexed_data; /* [0 .. num_rr_indexed_data-1] */
extern struct s_trace **trace_head, **trace_tail;
extern int **net_rr_terminals; /* [0..num_nets-1][0..num_pins-1] */
extern struct s_switch_inf *switch_inf; /* [0..det_routing_arch.num_switch-1] */
extern struct s_bb *route_bb;
extern struct s_block *block;

#define PARALLEL_ROUTE

/* Macro used below to ensure that fractions are rounded up, but floating   *
 * point values very close to an integer are rounded to that integer.       */

#define ROUND_UP(x) (ceil (x - 0.001))

#define ERROR_TOL 0.0001


static const char *rr_types[] =  {
	"SOURCE", "SINK", "IPIN", "OPIN", "CHANX", "CHANY", "INTRA_CLUSTER_EDGE"
};
/*truct mach_timebase_info {*/
		/*uint32_t	numer;*/
			/*uint32_t	denom;*/
/*};*/

/*typedef struct mach_timebase_info	*mach_timebase_info_t;*/
/*typedef struct mach_timebase_info	mach_timebase_info_data_t;*/


/*extern "C" {*/
/*kern_return_t		mach_timebase_info(*/
								/*mach_timebase_info_t	info);*/

/*kern_return_t		mach_wait_until(*/
								/*uint64_t		deadline);*/


/*uint64_t			mach_absolute_time(void);*/
/*uint64_t			mach_approximate_time(void);*/
/*}*/

/*void reserve_locally_used_opins(float pres_fac, boolean rip_up_local_opins,*/
		/*t_ivec ** clb_opins_used_locally);*/

boolean feasible_routing(void);

void print_route(char *name);

void update_net_delays_from_route_tree_new(t_net_timing *net_timing,
		t_rt_node ** rt_node_of_sink, int inet);

/******************** Subroutines local to route_timing.c ********************/

static int get_max_pins_per_net(void);

static float get_timing_driven_expected_cost(int inode, int target_node,
		float criticality_fac, float R_upstream, float base_cost_scale_factor);

static int get_expected_segs_to_target(int inode, int target_node,
		int *num_segs_ortho_dir_ptr);

static void timing_driven_check_net_delays(t_net_timing *net_timing);

static int mark_node_expansion_by_bin(int inet, int target_node,
		t_rt_node * rt_node);

/************************ Subroutine definitions *****************************/

extern zlog_category_t *next_net_log;
extern zlog_category_t *route_inner_log;
extern zlog_category_t *route_outer_log;
extern zlog_category_t *check_route_log;
extern zlog_category_t *congestion_log;
extern zlog_category_t *net_sort_log;
extern zlog_category_t *net_route_log;

typedef struct s_next_net {
	int current_net;
	int num_local_nets;
	int num_local_nets_routed;
	pthread_mutex_t lock;
	bool reversed;
	int num_threads;
} t_next_net;

static void init_next_net(t_next_net *next_net, int num_threads, bool reversed_order)
{
	next_net->current_net = 0;
	next_net->num_local_nets = 0;
	for (int i = 0; i < num_nets; ++i) {
		if (!clb_net[i].is_global) {
			++next_net->num_local_nets;
		}
	}
	zlog_info(next_net_log, "Num local nets: %d\n", next_net->num_local_nets);
	next_net->num_local_nets_routed = 0;
	if (pthread_mutex_init(&next_net->lock, NULL)) {
		printf("Failed to initialize next net mutex\n");
		exit(-1);
	}
	next_net->reversed = reversed_order;
	next_net->num_threads = num_threads;
}

static void reset_next_net(t_next_net *next_net)
{
	if (next_net->reversed) {
		next_net->current_net = num_nets-1;
	} else {
		next_net->current_net = 0;
	}
	next_net->num_local_nets_routed = 0;
}

static int get_next_net(t_next_net *next_net, int *net_index, int *real_index, bool *final)
{
	int net;
	pthread_mutex_lock(&next_net->lock);
	bool valid_net = true;
	bool is_global = true;
	while (is_global && valid_net) {
		if ((!next_net->reversed && next_net->current_net >= num_nets) || (next_net->reversed && next_net->current_net < 0)){
			valid_net = false;
		} else if (clb_net[net_index[next_net->current_net]].is_global) {
			if (next_net->reversed) {
				--next_net->current_net;
			} else {
				++next_net->current_net;
			}
		} else {
			is_global = false;
		}
	}
	if (valid_net && !is_global) {
		net = net_index[next_net->current_net];
		*real_index = next_net->num_local_nets_routed;
		if (next_net->reversed) {
			--next_net->current_net;
		} else {
			++next_net->current_net;
		}
		if ((next_net->num_local_nets - 1 - next_net->num_local_nets_routed) < (next_net->num_local_nets % next_net->num_threads)) {
			*final = true;
		} else {
			*final = false;
		}
		++next_net->num_local_nets_routed;
	} else {
		net = -1;
		*real_index = -1;
		*final = true;
	}
	pthread_mutex_unlock(&next_net->lock);
	return net;
}

typedef struct thread_info {
	int thread_index;

	sem_t *start_route_sem;
	sem_t *complete_route_sem;
	my_pthread_barrier_t *barrier;

	int *iter;
	int *net_index;
	t_router_opts *router_opts;
	float *pres_fac;
	t_net_timing *net_timing;
	t_net_route *net_route;
	int **sink_order; //[0..inet-1][1..num_sinks]
	timeval *net_route_time;

	bool successful;

	t_next_net *next_net;

	double total_stall_time;
	
} thread_info;

static void alloc_per_thread_timing_driven_route_structs(float **pin_criticality_ptr,
		int **sink_order_ptr) {

	/* Allocates all the structures needed only by the timing-driven router.   */

	int max_pins_per_net;

	max_pins_per_net = get_max_pins_per_net();
	*pin_criticality_ptr = new float[max_pins_per_net];
	*sink_order_ptr = new int[max_pins_per_net];
}

static bool compare_sink(const t_sink &a, const t_sink &b)
{
	return make_tuple(a.x, a.y, a.inode) < make_tuple(b.x, b.y, b.inode);
}

static double get_estimated_workload(int num_sinks, int bb_area)
{
	return -16.652675 + bb_area*0.241073 + num_sinks*27.626815;
}

static int get_max_pins_per_net(void) {

	/* Returns the largest number of pins on any non-global net.    */

	int inet, max_pins_per_net;

	max_pins_per_net = 0;
	for (inet = 0; inet < num_nets; inet++) {
		if (clb_net[inet].is_global == FALSE) {
			max_pins_per_net = std::max(max_pins_per_net,
					(clb_net[inet].num_sinks + 1));
		}
	}

	return (max_pins_per_net);
}

static void check_net_route(const t_trace *l_trace_head, int inet) {
	const t_trace *tptr = l_trace_head;
	int inode, prev_inode = -1;
	vector<int> visited_nodes;
	char buffer[256];

	zlog_debug(check_route_log, "Net %d visited nodes:\n", inet);
	while (tptr) {
		inode = tptr->index;
		visited_nodes.push_back(inode);

		sprintf_rr_node(inode, buffer);
		zlog_debug(check_route_log, "%s\n", buffer);

		if (prev_inode != -1) {
			t_rr_node *prev_node = &rr_node[prev_inode];
			bool found = false;
			for (int i = 0; i < prev_node->num_edges && !found; ++i) {
				found = prev_node->edges[i] == inode;
			}
			assert(found);
		} else {
			assert(rr_node[inode].type == SOURCE);
		}

		prev_inode = inode;
		tptr = tptr->next;

		if (rr_node[inode].type == SINK) {
			if (tptr != NULL) {
				prev_inode = tptr->index;
				tptr = tptr->next;
			}
		}
	}

	sort(visited_nodes.begin(), visited_nodes.end());
	for (int i = 0; i < visited_nodes.size()-1; ++i) {
		if (visited_nodes[i] == visited_nodes[i+1]) {
			zlog_fatal(check_route_log, "Node %d repeated\n", visited_nodes[i]);
		}
	}
	set<int> visited_nodes_set;
	for (const auto &n : visited_nodes) {
		visited_nodes_set.insert(n);
	}

	zlog_debug(check_route_log, "%d vs %d\n", visited_nodes_set.size(), visited_nodes.size());
	/*assert(visited_nodes_set.size() == visited_nodes.size());*/

	assert(visited_nodes_set.find(net_rr_terminals[inet][0]) != visited_nodes_set.end());
	for (int i = 0; i < clb_net[inet].num_sinks; ++i) {
		assert(visited_nodes_set.find(net_rr_terminals[inet][i+1]) != visited_nodes_set.end());
	}
}

typedef struct s_simple_net {
	int source_inode;
	vector<int> sink_inodes;

	int current_sink;
	struct s_bb current_bounding_box;

	t_trace *trace_head;
	t_trace *trace_tail;

	t_rt_node **rr_node_to_rt_node;

	t_net_timing timing;
} t_simple_net;

static void add_route_tree_to_heap(
		const t_rt_node * rt_node,
		int target_node, float target_criticality, float astar_fac, const t_rr_node_route_inf *l_rr_node_route_inf, int inet,
		std::priority_queue<struct s_heap> &heap, int &num_heap_pushes) {

	/* Puts the entire partial routing below and including rt_node onto the heap *
	 * (except for those parts marked as not to be expanded) by calling itself   *
	 * recursively.                                                              */

	int inode;
	t_rt_node *child_node;
	t_linked_rt_edge *linked_rt_edge;
	float tot_cost, backward_path_cost, R_upstream;
	char buffer[256];

	/* Pre-order depth-first traversal */

	if (rt_node->re_expand) {
		inode = rt_node->inode;
		backward_path_cost = target_criticality * rt_node->Tdel;
		R_upstream = rt_node->R_upstream;
		tot_cost = backward_path_cost
				+ astar_fac
						* get_timing_driven_expected_cost(inode, target_node,
								target_criticality, R_upstream, 1);

		node_to_heap(inode, NO_PREVIOUS, NO_PREVIOUS,
				tot_cost, backward_path_cost, R_upstream, l_rr_node_route_inf, heap);
		++num_heap_pushes;
	} else {
		sprintf_rr_node(rt_node->inode, buffer);
		zlog_debug(route_inner_log, "Not expanding %s\n", buffer);
	}

	linked_rt_edge = rt_node->u.child_list;

	while (linked_rt_edge != NULL) {
		child_node = linked_rt_edge->child;
		add_route_tree_to_heap(child_node, target_node, target_criticality,
				astar_fac, l_rr_node_route_inf, inet, heap, num_heap_pushes);
		linked_rt_edge = linked_rt_edge->next;
	}
}

/*std::ostream &operator<<(std::ostream &out, const t_rr_node &node)*/
/*{*/
	/*out */
		/*<< rr_types[node.type] << " ";*/
	/*if (node.direction == INC_DIRECTION) {*/
		/*if (node.type == CHANX) {*/
			/*out */
				/*<< "(" << node.xlow << "->" << node.xhigh << ", " << node.ylow << ")";*/
		/*} else if (node.type == CHANY) {*/
			/*out */
				/*<< "(" << node.xlow << ", " << node.ylow << "->" << node.yhigh << ")";*/
		/*} else {*/
			/*out */
				/*<< "(" << node.xlow << ", " << node.ylow << ")"*/
				/*<< "(" << node.xhigh << ", " << node.yhigh << ")";*/
		/*}*/
	/*} else {*/
		/*if (node.type == CHANX) {*/
			/*out */
				/*<< "(" << node.xhigh << "->" << node.xlow << ", " << node.ylow << ")";*/
		/*} else if (node.type == CHANY) {*/
			/*out */
				/*<< "(" << node.xlow << ", " << node.yhigh << "->" << node.ylow << ")";*/
		/*} else {*/
			/*out */
				/*<< "(" << node.xhigh << ", " << node.yhigh << ")"*/
				/*<< "(" << node.xlow << ", " << node.ylow << ")";*/
		/*}*/
	/*}*/
	/*return out;*/
/*}*/

static int timing_driven_expand_neighbours(
		int thread_id,
		const struct s_heap *current, int target_node,
		int inet,
		float bend_cost, float astar_fac, float criticality_fac, int highfanout_rlim, bool monotonic,
		const t_rr_node_route_inf *l_rr_node_route_inf,
		std::priority_queue<struct s_heap> &heap) {

	/* Puts all the rr_nodes adjacent to current on the heap.  rr_nodes outside *
	 * the expanded bounding box specified in route_bb are not added to the     *
	 * heap.                                                                    */

	int iconn, to_node, num_edges, inode, iswitch, target_x, target_y;
	t_rr_type from_type, to_type;
	float new_tot_cost, old_back_pcost, new_back_pcost, R_upstream;
	float new_R_upstream, Tdel;
	char buffer[256];

	inode = current->index;
	old_back_pcost = current->backward_path_cost;
	R_upstream = current->R_upstream;
	num_edges = rr_node[inode].num_edges;

	target_x = rr_node[target_node].xhigh;
	target_y = rr_node[target_node].yhigh;

	assert(rr_node[target_node].xlow == rr_node[target_node].xhigh &&
			rr_node[target_node].ylow == rr_node[target_node].yhigh);

	int num_pushes = 0;

	for (iconn = 0; iconn < num_edges; iconn++) {
		to_node = rr_node[inode].edges[iconn];

		sprintf_rr_node(to_node, buffer);
		zlog_debug(route_inner_log, "\t[%d] Neighbor: %s \n", thread_id, buffer);

		std::pair<int, int> current_pos = get_node_start(inode); 
		std::pair<int, int> sink_pos = get_node_end(target_node); 

		/*Interval<int> hor(current_pos.first, sink_pos.first);*/
		/*Interval<int> vert(current_pos.second, sink_pos.second);*/

		std::pair<int, int> neighbor_pos = get_node_start(to_node); 

		if (rr_node[inode].type == OPIN) {
			assert(rr_node[inode].xlow == rr_node[inode].xhigh &&
					rr_node[inode].ylow == rr_node[inode].yhigh);

			current_pos = neighbor_pos;
		}
		
		int current_to_sink_manhattan_distance = abs(current_pos.first - sink_pos.first) + abs(current_pos.second - sink_pos.second);
		int neighbor_to_sink_manhattan_distance = abs(neighbor_pos.first - sink_pos.first) + abs(neighbor_pos.second - sink_pos.second);

		if (neighbor_to_sink_manhattan_distance > current_to_sink_manhattan_distance) {
			/*zlog_debug(route_inner_log, "Neighbor is non-monotonic [%d > %d] ", neighbor_to_sink_manhattan_distance, current_to_sink_manhattan_distance);*/
			/*continue;*/
		}

		/*if ((!hor.contains(neighbor_pos.first) || !vert.contains(neighbor_pos.second)) && monotonic) {*/
			/*DVLOG(2) << rr_node[inode] << " -> " << rr_node[to_node] << " is outside of " << rr_node[target_node];*/
			/*continue;*/
		/*}*/

		if (rr_node[to_node].xhigh < route_bb[inet].xmin
				|| rr_node[to_node].xlow > route_bb[inet].xmax
				|| rr_node[to_node].yhigh < route_bb[inet].ymin
				|| rr_node[to_node].ylow > route_bb[inet].ymax) {
			zlog_debug(route_inner_log, "\t\t[%d] Outside of bounding box\n", thread_id);
			continue; /* Node is outside (expanded) bounding box. */
		}

		if (clb_net[inet].num_sinks >= HIGH_FANOUT_NET_LIM) {
			if (rr_node[to_node].xhigh < target_x - highfanout_rlim
					|| rr_node[to_node].xlow > target_x + highfanout_rlim
					|| rr_node[to_node].yhigh < target_y - highfanout_rlim
					|| rr_node[to_node].ylow > target_y + highfanout_rlim) {
				zlog_debug(route_inner_log, "\t\t[%d] Outside of high fanout bin\n", thread_id);
				continue; /* Node is outside high fanout bin. */
			}
		}

		/* Prune away IPINs that lead to blocks other than the target one.  Avoids  *
		 * the issue of how to cost them properly so they don't get expanded before *
		 * more promising routes, but makes route-throughs (via CLBs) impossible.   *
		 * Change this if you want to investigate route-throughs.                   */

		to_type = rr_node[to_node].type;
		/*assert(to_type != IPIN || (rr_node[to_node].xlow == rr_node[to_node].xhigh &&*/
					/*rr_node[to_node].ylow == rr_node[to_node].yhigh));*/
		if (to_type == IPIN
				&& (rr_node[to_node].xhigh != target_x
						|| rr_node[to_node].yhigh != target_y)) {
			zlog_debug(route_inner_log, "\t\t[%d] Not the target IPIN\n", thread_id);
			continue;
		}

		/* new_back_pcost stores the "known" part of the cost to this node -- the   *
		 * congestion cost of all the routing resources back to the existing route  *
		 * plus the known delay of the total path back to the source.  new_tot_cost *
		 * is this "known" backward cost + an expected cost to get to the target.   */

		float cong_cost = get_rr_cong_cost(to_node, clb_net[inet].num_sinks);
		new_back_pcost = old_back_pcost
				+ (1. - criticality_fac) * cong_cost;

		iswitch = rr_node[inode].switches[iconn];
		if (switch_inf[iswitch].buffered) {
			new_R_upstream = switch_inf[iswitch].R;
		} else {
			new_R_upstream = R_upstream + switch_inf[iswitch].R;
		}

		Tdel = rr_node[to_node].C * (new_R_upstream + 0.5 * rr_node[to_node].R);
		Tdel += switch_inf[iswitch].Tdel;
		new_R_upstream += rr_node[to_node].R;
		new_back_pcost += criticality_fac * Tdel;

		if (bend_cost != 0.) {
			from_type = rr_node[inode].type;
			to_type = rr_node[to_node].type;
			if ((from_type == CHANX && to_type == CHANY)
					|| (from_type == CHANY && to_type == CHANX))
				new_back_pcost += bend_cost;
		}

		zlog_debug(route_inner_log, "\t\t[%d] old_b: %g crit: %g cong: %g Tdel: %g bend: %g\n", thread_id,
				old_back_pcost, criticality_fac, cong_cost, Tdel, bend_cost);

		new_tot_cost = new_back_pcost
				+ astar_fac
						* get_timing_driven_expected_cost(to_node, target_node,
								criticality_fac, new_R_upstream, 1);

		node_to_heap(thread_id, to_node, inode, iconn, new_tot_cost, new_back_pcost, new_R_upstream, l_rr_node_route_inf, heap);
		++num_pushes;

	} /* End for all neighbours */
	return num_pushes;
}

static float get_timing_driven_expected_cost(int inode, int target_node,
		float criticality_fac, float R_upstream, float base_cost_scale_factor) {

	/* Determines the expected cost (due to both delay and resouce cost) to reach *
	 * the target node from inode.  It doesn't include the cost of inode --       *
	 * that's already in the "known" path_cost.                                   */

	t_rr_type rr_type;
	int cost_index, ortho_cost_index, num_segs_same_dir, num_segs_ortho_dir;
	float expected_cost, cong_cost, Tdel;

	rr_type = rr_node[inode].type;

	if (rr_type == CHANX || rr_type == CHANY) {
		num_segs_same_dir = get_expected_segs_to_target(inode, target_node,
				&num_segs_ortho_dir);
		cost_index = rr_node[inode].cost_index;
		ortho_cost_index = rr_indexed_data[cost_index].ortho_cost_index;

		cong_cost = num_segs_same_dir * rr_indexed_data[cost_index].saved_base_cost
				+ num_segs_ortho_dir
						* rr_indexed_data[ortho_cost_index].saved_base_cost;
		cong_cost += rr_indexed_data[IPIN_COST_INDEX].base_cost
				+ rr_indexed_data[SINK_COST_INDEX].base_cost;

		Tdel =
				num_segs_same_dir * rr_indexed_data[cost_index].T_linear
						+ num_segs_ortho_dir
								* rr_indexed_data[ortho_cost_index].T_linear
						+ num_segs_same_dir * num_segs_same_dir
								* rr_indexed_data[cost_index].T_quadratic
						+ num_segs_ortho_dir * num_segs_ortho_dir
								* rr_indexed_data[ortho_cost_index].T_quadratic
						+ R_upstream
								* (num_segs_same_dir
										* rr_indexed_data[cost_index].C_load
										+ num_segs_ortho_dir
												* rr_indexed_data[ortho_cost_index].C_load);

		Tdel += rr_indexed_data[IPIN_COST_INDEX].T_linear;

		expected_cost = criticality_fac * Tdel
				+ (1. - criticality_fac) * cong_cost;
		return (expected_cost);
	}

	else if (rr_type == IPIN) { /* Change if you're allowing route-throughs */
		return (rr_indexed_data[SINK_COST_INDEX].base_cost);
	}

	else { /* Change this if you want to investigate route-throughs */
		return (0.);
	}
}

static int get_expected_segs_to_target(int inode, int target_node,
		int *num_segs_ortho_dir_ptr) {

	/* Returns the number of segments the same type as inode that will be needed *
	 * to reach target_node (not including inode) in each direction (the same    *
	 * direction (horizontal or vertical) as inode and the orthogonal direction).*/

	t_rr_type rr_type;
	int target_x, target_y, num_segs_same_dir, cost_index, ortho_cost_index;
	int no_need_to_pass_by_clb;
	float inv_length, ortho_inv_length, ylow, yhigh, xlow, xhigh;

	target_x = rr_node[target_node].xlow;
	target_y = rr_node[target_node].ylow;
	cost_index = rr_node[inode].cost_index;
	inv_length = rr_indexed_data[cost_index].inv_length;
	ortho_cost_index = rr_indexed_data[cost_index].ortho_cost_index;
	ortho_inv_length = rr_indexed_data[ortho_cost_index].inv_length;
	rr_type = rr_node[inode].type;

	if (rr_type == CHANX) {
		ylow = rr_node[inode].ylow;
		xhigh = rr_node[inode].xhigh;
		xlow = rr_node[inode].xlow;

		/* Count vertical (orthogonal to inode) segs first. */

		if (ylow > target_y) { /* Coming from a row above target? */
			*num_segs_ortho_dir_ptr =
					(int)(ROUND_UP((ylow - target_y + 1.) * ortho_inv_length));
			no_need_to_pass_by_clb = 1;
		} else if (ylow < target_y - 1) { /* Below the CLB bottom? */
			*num_segs_ortho_dir_ptr = (int)(ROUND_UP((target_y - ylow) *
					ortho_inv_length));
			no_need_to_pass_by_clb = 1;
		} else { /* In a row that passes by target CLB */
			*num_segs_ortho_dir_ptr = 0;
			no_need_to_pass_by_clb = 0;
		}

		/* Now count horizontal (same dir. as inode) segs. */

		if (xlow > target_x + no_need_to_pass_by_clb) {
			num_segs_same_dir = (int)(ROUND_UP((xlow - no_need_to_pass_by_clb -
							target_x) * inv_length));
		} else if (xhigh < target_x - no_need_to_pass_by_clb) {
			num_segs_same_dir = (int)(ROUND_UP((target_x - no_need_to_pass_by_clb -
							xhigh) * inv_length));
		} else {
			num_segs_same_dir = 0;
		}
	}

	else { /* inode is a CHANY */
		ylow = rr_node[inode].ylow;
		yhigh = rr_node[inode].yhigh;
		xlow = rr_node[inode].xlow;

		/* Count horizontal (orthogonal to inode) segs first. */

		if (xlow > target_x) { /* Coming from a column right of target? */
			*num_segs_ortho_dir_ptr = (int)(
					ROUND_UP((xlow - target_x + 1.) * ortho_inv_length));
			no_need_to_pass_by_clb = 1;
		} else if (xlow < target_x - 1) { /* Left of and not adjacent to the CLB? */
			*num_segs_ortho_dir_ptr = (int)(ROUND_UP((target_x - xlow) *
					ortho_inv_length));
			no_need_to_pass_by_clb = 1;
		} else { /* In a column that passes by target CLB */
			*num_segs_ortho_dir_ptr = 0;
			no_need_to_pass_by_clb = 0;
		}

		/* Now count vertical (same dir. as inode) segs. */

		if (ylow > target_y + no_need_to_pass_by_clb) {
			num_segs_same_dir = (int)(ROUND_UP((ylow - no_need_to_pass_by_clb -
							target_y) * inv_length));
		} else if (yhigh < target_y - no_need_to_pass_by_clb) {
			num_segs_same_dir = (int)(ROUND_UP((target_y - no_need_to_pass_by_clb -
							yhigh) * inv_length));
		} else {
			num_segs_same_dir = 0;
		}
	}

	return (num_segs_same_dir);
}

/* Nets that have high fanout can take a very long time to route.  Each sink should be routed contained within a bin instead of the entire bounding box to speed things up */
static int mark_node_expansion_by_bin(int inet, int target_node,
		t_rt_node * rt_node) {
	int target_x, target_y;
	int rlim = 1;
	int inode;
	float area;
	boolean success;
	t_linked_rt_edge *linked_rt_edge;
	t_rt_node * child_node;

	target_x = rr_node[target_node].xlow;
	target_y = rr_node[target_node].ylow;

	if (clb_net[inet].num_sinks < HIGH_FANOUT_NET_LIM) {
		/* This algorithm only applies to high fanout nets */
		return 1;
	}

	area = (route_bb[inet].xmax - route_bb[inet].xmin)
			* (route_bb[inet].ymax - route_bb[inet].ymin);
	if (area <= 0) {
		area = 1;
	}

	rlim = (int)(ceil(sqrt((float) area / (float) clb_net[inet].num_sinks)));
	if (rt_node == NULL || rt_node->u.child_list == NULL) {
		/* If unknown traceback, set radius of bin to be size of chip */
		rlim = std::max(nx + 2, ny + 2);
		return rlim;
	}

	success = FALSE;
	/* determine quickly a feasible bin radius to route sink for high fanout nets 
	 this is necessary to prevent super long runtimes for high fanout nets; in best case, a reduction in complexity from O(N^2logN) to O(NlogN) (Swartz fast router)
	 */
	linked_rt_edge = rt_node->u.child_list;
	while (success == FALSE && linked_rt_edge != NULL) {
		while (linked_rt_edge != NULL && success == FALSE) {
			child_node = linked_rt_edge->child;
			inode = child_node->inode;
			if (!(rr_node[inode].type == IPIN || rr_node[inode].type == SINK)) {
				if (rr_node[inode].xlow <= target_x + rlim
						&& rr_node[inode].xhigh >= target_x - rlim
						&& rr_node[inode].ylow <= target_y + rlim
						&& rr_node[inode].yhigh >= target_y - rlim) {
					success = TRUE;
				}
			}
			linked_rt_edge = linked_rt_edge->next;
		}

		if (success == FALSE) {
			if (rlim > std::max(nx + 2, ny + 2)) {
				vpr_printf(TIO_MESSAGE_ERROR, "VPR internal error, net %s has paths that are not found in traceback.\n",
						clb_net[inet].name);
				exit(1);
			}
			/* if sink not in bin, increase bin size until fit */
			rlim *= 2;
		} else {
			/* Sometimes might just catch a wire in the end segment, need to give it some channel space to explore */
			rlim += 4;
		}
		linked_rt_edge = rt_node->u.child_list;
	}

	/* redetermine expansion based on rlim */
	linked_rt_edge = rt_node->u.child_list;
	/*printf("Starting to print\n");*/
	/*printf("%d\n", rt_node->inode);*/
	while (linked_rt_edge != NULL) {
		child_node = linked_rt_edge->child;
		inode = child_node->inode;
		/*printf("%d\n", inode);*/
		if (!(rr_node[inode].type == IPIN || rr_node[inode].type == SINK)) {
			if (rr_node[inode].xlow <= target_x + rlim
					&& rr_node[inode].xhigh >= target_x - rlim
					&& rr_node[inode].ylow <= target_y + rlim
					&& rr_node[inode].yhigh >= target_y - rlim) {
				child_node->re_expand = TRUE;
			} else {
				child_node->re_expand = FALSE;
			}
		}
		linked_rt_edge = linked_rt_edge->next;
	}
	return rlim;
}

void distribute_work(int inode, int num_threads, vector<int> &works)
{
	std::deque<int> inodes;
	inodes.push_back(inode);
	bool done = false;
	while (inodes.size() < num_threads) {
		int current = inodes.front();
		inodes.pop_front();
		for (int i = 0; i <  rr_node[current].num_edges; ++i) {
			int neighbor = rr_node[current].edges[i];
			inodes.push_back(neighbor);
		}
	}
	for (const auto &i : inodes) {
		works.push_back(i);
	}
}

void test_distributed_work()
{
	vector<int> works;
	distribute_work(7652, 8, works);
	for (const auto &i : works) {
		printf("work: %d\n", i);
	}
}

typedef struct {
	std::priority_queue<struct s_heap> heap;
	int worker_id;
} RouteWorkerFeederItem;

class GlobalRouteWorker {
	public:
		/*std::priority_queue<struct s_heap> *heaps;*/
		std::mutex &global_heap_lock;
		std::priority_queue<struct s_heap> &global_heap;

		std::mutex &l_rr_node_route_inf_lock;

		tbb::atomic<bool> &found_sink;
		tbb::atomic<int> &num_found_sink;
		struct s_heap &sink_node;

		std::mutex &spawn_lock;
		int &heap_index; /* access must be protected by mutex */
		int &num_spawned_workers;

		std::vector<int> &modified_inodes;
		std::mutex &modified_inodes_lock;

		t_rr_node_route_inf *l_rr_node_route_inf;
		t_rt_node **l_rr_node_to_rt_node;
		int target_node;
		int inet;
		float bend_cost;
		float target_criticality;
		float astar_fac;
		int highfanout_rlim;
		int num_threads;

		typedef RouteWorkerFeederItem argument_type;

		GlobalRouteWorker(std::mutex &global_heap_lock, std::priority_queue<struct s_heap> &global_heap,
				std::mutex &l_rr_node_route_inf_lock,
				tbb::atomic<bool> &found_sink, struct s_heap &sink_node, tbb::atomic<int> &num_found_sink,
				std::mutex &spawn_lock, int &heap_index, int &num_spawned_workers,
				std::vector<int> &modified_inodes, std::mutex &modified_inodes_lock
				)
			: global_heap(global_heap), global_heap_lock(global_heap_lock),
				l_rr_node_route_inf_lock(l_rr_node_route_inf_lock),
			found_sink(found_sink), sink_node(sink_node), num_found_sink(num_found_sink),
			spawn_lock(spawn_lock), heap_index(heap_index), num_spawned_workers(num_spawned_workers),
			modified_inodes(modified_inodes), modified_inodes_lock(modified_inodes_lock)
		{
		}

		void operator()(argument_type item, tbb::parallel_do_feeder<argument_type> &feeder) const
		{
			/*char buffer[256];*/
			/*struct s_heap current_heap_item;// = item.heap_item;*/

			/*std::priority_queue<struct s_heap> local_heap;*/

			/*global_heap_lock.lock();*/
			/*current_heap_item = global_heap.top();*/
			/*global_heap.pop();*/
			/*global_heap_lock.unlock();*/

			/*int inode = current_heap_item.index;*/
			
			/*while (!found_sink && inode != target_node) {*/

				/*sprintf_rr_node(inode, buffer);*/

				/*l_rr_node_route_inf_lock.lock();*/

				/*float old_tcost = l_rr_node_route_inf[inode].path_cost;*/
				/*float old_back_cost;*/

				/*if (old_tcost > 0.99 * HUGE_POSITIVE_FLOAT) [> First time touched. <]*/
					/*old_back_cost = HUGE_POSITIVE_FLOAT;*/
				/*else*/
					/*old_back_cost = l_rr_node_route_inf[inode].backward_path_cost;*/

				/*float new_tcost = current_heap_item.cost;*/
				/*float new_back_cost = current_heap_item.backward_path_cost;*/

				/*bool expand_neighbours = false;*/
				/*bool in_route_tree = l_rr_node_to_rt_node[inode] != NULL;*/
				/*bool added_from_route_tree = current_heap_item.u.prev_node == NO_PREVIOUS;*/

				/*zlog_debug(route_inner_log, "Current: %s old_t: %g new_t: %g route_inf_b: %g old_b: %g new_b: %g old_prev_node: %d new_prev_node: %d in_route_tree: %d added_from_route_tree %d\n", buffer, old_tcost, new_tcost, l_rr_node_route_inf[inode].backward_path_cost, old_back_cost, new_back_cost, l_rr_node_route_inf[inode].prev_node, current_heap_item.u.prev_node, in_route_tree ? 1 : 0, added_from_route_tree ? 1 : 0);*/

				/*if (in_route_tree && !added_from_route_tree) {*/
					/*expand_neighbours = false;*/
				/*} else if (in_route_tree && added_from_route_tree) {*/
					/*assert(current_heap_item.u.prev_node == NO_PREVIOUS);*/
					/*assert(current_heap_item.prev_edge == NO_PREVIOUS);*/

					/*l_rr_node_route_inf[inode].prev_node = current_heap_item.u.prev_node;*/
					/*l_rr_node_route_inf[inode].prev_edge = current_heap_item.prev_edge;*/
					/*l_rr_node_route_inf[inode].path_cost = current_heap_item.cost; [> astar cost changes, so we can't check for equality here <]*/
					/*l_rr_node_route_inf[inode].backward_path_cost = current_heap_item.backward_path_cost; [> when after adding to heap, the backward path cost does not include the congestion cost <]*/
					/*zlog_debug(route_inner_log, "\tForcing expansion of existing route tree node\n");*/
					/*expand_neighbours = true;*/
				/*} else if (old_tcost > new_tcost && old_back_cost > new_back_cost) {*/
					/*if (in_route_tree && !added_from_route_tree && current_heap_item.u.prev_node != l_rr_node_route_inf[inode].prev_node) {*/
						/*zlog_warn(route_inner_log, "Trying to drive an existing route tree node %s with a new driver\n", buffer);*/
					/*} else {*/
						/*l_rr_node_route_inf[inode].prev_node = current_heap_item.u.prev_node;*/
						/*l_rr_node_route_inf[inode].prev_edge = current_heap_item.prev_edge;*/
						/*l_rr_node_route_inf[inode].path_cost = new_tcost;*/
						/*l_rr_node_route_inf[inode].backward_path_cost = new_back_cost;*/

						/*if (old_tcost > 0.99 * HUGE_POSITIVE_FLOAT) [> First time touched. <] {*/
							/*modified_inodes_lock.lock();*/
							/*modified_inodes.push_back(inode);*/
							/*modified_inodes_lock.unlock();*/
						/*}*/

						/*zlog_debug(route_inner_log, "\tNormal expansion\n");*/
						/*expand_neighbours = true;*/
					/*}*/
				/*}*/

				/*l_rr_node_route_inf_lock.unlock();*/

				/*int num_heap_pushes = 0;*/

				/*[>std::priority_queue<struct s_heap> local_heap = heaps[local_heap_index];<]*/

				/*if (expand_neighbours) {*/
					/*timing_driven_expand_neighbours(*/
							/*&current_heap_item, target_node,*/
							/*inet,*/
							/*bend_cost, astar_fac, target_criticality, highfanout_rlim, false,*/
							/*l_rr_node_route_inf,*/
							/*local_heap);*/
				/*}*/

				/*[>if (heap.empty()) {<]*/
				/*[>zlog_fatal(route_inner_log, "Failed to find a path to %d\n", target_node);<]*/
				/*[>return FALSE;<]*/
				/*[>}<]*/

				/*[> prioritize best item for current task <]*/
				/*if (!local_heap.empty()) {*/
					/*current_heap_item = local_heap.top();*/
					/*local_heap.pop();*/

					/*inode = current_heap_item.index;*/
				/*} else {*/
					/*spawn_lock.lock();*/
					/*--num_spawned_workers;*/
					/*spawn_lock.unlock();*/
					/*break;*/
				/*}*/

				/*spawn_lock.lock();*/
				/*while (num_spawned_workers < num_threads && !local_heap.empty()) {*/
					/*current_heap_item = local_heap.top();*/
					/*local_heap.pop();*/

					/*feeder.add({ local_heap, 0 });*/
					/*++num_spawned_workers;*/
				/*}*/
				/*spawn_lock.unlock();*/
			/*}*/
			/*if (inode == target_node) {*/
				/*found_sink = true;*/
				/*++num_found_sink;*/
				/*sink_node = current_heap_item;*/
			/*}*/
			/*zlog_debug(route_inner_log, "RouteWorker done\n");*/
		}
};

class RouteWorker2 {
	public:
		/*std::priority_queue<struct s_heap> *heaps;*/

		std::mutex &l_rr_node_route_inf_lock;

		std::mutex &found_sinks_lock;
		vector<struct s_heap> &found_sinks;

		std::mutex &spawn_lock;
		int &heap_index; /* access must be protected by mutex */
		int &num_spawned_workers;
		int &worker_id;

		tbb::atomic<int> **visit_count;

		std::vector<int> &modified_inodes;
		std::mutex &modified_inodes_lock;

		std::mutex &push_lock;
		int &num_heap_pushes;

		t_rr_node_route_inf *l_rr_node_route_inf;
		t_rt_node **l_rr_node_to_rt_node;
		int target_node;
		int inet;
		float bend_cost;
		float target_criticality;
		float astar_fac;
		int highfanout_rlim;
		int num_threads;

		typedef RouteWorkerFeederItem argument_type;

		RouteWorker2(std::mutex &l_rr_node_route_inf_lock,
				std::mutex &found_sinks_lock, vector<struct s_heap> &found_sinks,
				std::mutex &spawn_lock, int &heap_index, int &num_spawned_workers, int &worker_id,
				std::vector<int> &modified_inodes, std::mutex &modified_inodes_lock,
				std::mutex &push_lock, int &num_heap_pushes
				)
			: l_rr_node_route_inf_lock(l_rr_node_route_inf_lock),
			found_sinks_lock(found_sinks_lock), found_sinks(found_sinks),
			spawn_lock(spawn_lock), heap_index(heap_index), num_spawned_workers(num_spawned_workers), worker_id(worker_id),
			modified_inodes(modified_inodes), modified_inodes_lock(modified_inodes_lock),
			push_lock(push_lock), num_heap_pushes(num_heap_pushes)
		{
		}

		void operator()(argument_type item, tbb::parallel_do_feeder<argument_type> &feeder) const
		{
			char buffer[256];

			std::priority_queue<struct s_heap> &local_heap = item.heap;

			assert(!local_heap.empty());
			struct s_heap current_heap_item = local_heap.top();
			local_heap.pop();
			int inode = current_heap_item.index;
			auto start = get_node_start(inode);
			++visit_count[start.first][start.second];

			sprintf_rr_node(inode, buffer);
			zlog_debug(route_inner_log, "[%d] Start item %s\n", item.worker_id, buffer);
			
			while (inode != target_node) {
				l_rr_node_route_inf_lock.lock();

				float old_tcost = l_rr_node_route_inf[inode].path_cost;
				float old_back_cost;

				if (!l_rr_node_route_inf[inode].modified) /* First time touched. */
					old_back_cost = HUGE_POSITIVE_FLOAT;
				else
					old_back_cost = l_rr_node_route_inf[inode].backward_path_cost;

				float new_tcost = current_heap_item.cost;
				float new_back_cost = current_heap_item.backward_path_cost;

				bool expand_neighbours = false;
				bool in_route_tree = l_rr_node_to_rt_node[inode] != NULL;
				bool added_from_route_tree = current_heap_item.u.prev_node == NO_PREVIOUS;

				sprintf_rr_node(inode, buffer);
				zlog_debug(route_inner_log, "[%d] Current: %s old_t: %g new_t: %g route_inf_b: %g old_b: %g new_b: %g old_prev_node: %d new_prev_node: %d in_route_tree: %d added_from_route_tree %d\n", item.worker_id, buffer, old_tcost, new_tcost, l_rr_node_route_inf[inode].backward_path_cost, old_back_cost, new_back_cost, l_rr_node_route_inf[inode].prev_node, current_heap_item.u.prev_node, in_route_tree ? 1 : 0, added_from_route_tree ? 1 : 0);

				/*if (in_route_tree && !added_from_route_tree) {*/
					/*expand_neighbours = false;*/
				/*} else if (in_route_tree && added_from_route_tree) {*/
					/*assert(current_heap_item.u.prev_node == NO_PREVIOUS);*/
					/*assert(current_heap_item.prev_edge == NO_PREVIOUS);*/

					/*l_rr_node_route_inf[inode].prev_node = current_heap_item.u.prev_node;*/
					/*l_rr_node_route_inf[inode].prev_edge = current_heap_item.prev_edge;*/
					/*l_rr_node_route_inf[inode].path_cost = current_heap_item.cost; [> astar cost changes, so we can't check for equality here <]*/
					/*l_rr_node_route_inf[inode].backward_path_cost = current_heap_item.backward_path_cost; [> when after adding to heap, the backward path cost does not include the congestion cost <]*/
					/*zlog_debug(route_inner_log, "\tForcing expansion of existing route tree node\n");*/
					/*expand_neighbours = true;*/
				/*} else*/
				bool written = false;

				if (old_back_cost > new_back_cost && old_tcost > new_tcost) {
					if (in_route_tree && !added_from_route_tree && current_heap_item.u.prev_node != l_rr_node_route_inf[inode].prev_node) {
						/*zlog_warn(route_inner_log, "\t[%d] Trying to drive an existing route tree node %s with a new driver\n", item.worker_id, buffer);*/
					} else {
						l_rr_node_route_inf[inode].prev_node = current_heap_item.u.prev_node;
						l_rr_node_route_inf[inode].prev_edge = current_heap_item.prev_edge;
						l_rr_node_route_inf[inode].path_cost = new_tcost;
						l_rr_node_route_inf[inode].backward_path_cost = new_back_cost;

						written = true;

						if (!l_rr_node_route_inf[inode].modified) /* First time touched. */ {
							modified_inodes_lock.lock();
							modified_inodes.push_back(inode);
							modified_inodes_lock.unlock();
							l_rr_node_route_inf[inode].modified = true;
						}

						zlog_debug(route_inner_log, "\t[%d] Normal expansion\n", item.worker_id);
						expand_neighbours = true;
					}
				}

				transaction_t trans;
				trans.target_node = target_node;
				trans.old_tcost = old_tcost;
				trans.old_bcost = old_back_cost;
				trans.new_bcost = new_back_cost;
				trans.new_tcost = new_tcost;
				trans.prev_node = current_heap_item.u.prev_node;
				trans.prev_edge = current_heap_item.prev_edge;
				trans.written = written;

				l_rr_node_route_inf[inode].transactions.push_back(trans);

				l_rr_node_route_inf_lock.unlock();

				/*std::priority_queue<struct s_heap> local_heap = heaps[local_heap_index];*/
				int local_num_heap_pushes = 0;

				if (expand_neighbours) {
					local_num_heap_pushes = timing_driven_expand_neighbours(
							item.worker_id,
							&current_heap_item, target_node,
							inet,
							bend_cost, astar_fac, target_criticality, highfanout_rlim, false,
							l_rr_node_route_inf,
							local_heap);
				}

				push_lock.lock();
				num_heap_pushes += local_num_heap_pushes;
				push_lock.unlock();

				/*if (heap.empty()) {*/
				/*zlog_fatal(route_inner_log, "Failed to find a path to %d\n", target_node);*/
				/*return FALSE;*/
				/*}*/

				/* prioritize best item for current task */
				if (!local_heap.empty()) {
					current_heap_item = local_heap.top();
					local_heap.pop();

					inode = current_heap_item.index;
					start = get_node_start(inode);
					++visit_count[start.first][start.second];
				} else {
					zlog_debug(route_inner_log, "[%d] Breaking because of empty heap\n", item.worker_id);
					spawn_lock.lock();
					--num_spawned_workers;
					spawn_lock.unlock();
					break;
				}

				spawn_lock.lock();
				while (num_spawned_workers < num_threads && !local_heap.empty()) {
					std::priority_queue<struct s_heap> new_heap;
					new_heap.push(local_heap.top());
					local_heap.pop();

					feeder.add({ new_heap, worker_id++ });

					sprintf_rr_node(new_heap.top().index, buffer);
					zlog_debug(route_inner_log, "[%d] Adding to feeder inode %s\n", item.worker_id, buffer);

					++num_spawned_workers;
				}
				spawn_lock.unlock();
			}

			if (inode == target_node) {
				assert(current_heap_item.index == inode);

				sprintf_rr_node(inode, buffer);
				zlog_debug(route_inner_log, "[%d] Found %s\n", item.worker_id, buffer);

				found_sinks_lock.lock();
				found_sinks.push_back(current_heap_item);
				found_sinks_lock.unlock();

				/*l_rr_node_route_inf_lock.lock();*/
				/*int sink_inode = current_heap_item.index;*/
				/*if (current_heap_item.backward_path_cost < l_rr_node_route_inf[sink_inode].backward_path_cost*/
						/*&& current_heap_item.cost < l_rr_node_route_inf[sink_inode].path_cost) {*/
					/*l_rr_node_route_inf[sink_inode].backward_path_cost = current_heap_item.backward_path_cost;*/
					/*l_rr_node_route_inf[sink_inode].path_cost = current_heap_item.cost;*/
					/*l_rr_node_route_inf[sink_inode].prev_node = current_heap_item.u.prev_node;*/
					/*l_rr_node_route_inf[sink_inode].prev_edge = current_heap_item.prev_edge;*/

					/*modified_inodes_lock.lock();*/
					/*modified_inodes.push_back(sink_inode);*/
					/*modified_inodes_lock.unlock();*/

					/*sprintf_rr_node(inode, buffer);*/
					/*zlog_info(route_inner_log, "[%d] Found %s\n", item.worker_id, buffer);*/

					/*sink_lock.lock();*/
					/*found_sink = true;*/
					/*++num_found_sink;*/
					/*sink_node = current_heap_item;*/
					/*sink_lock.unlock();*/
				/*} else {*/
					/*zlog_debug(route_inner_log, "[%d] Not setting found sink because old_b %g new_b %g old_t %g new_t %g prev_node: %d prev_edge: %d\n", item.worker_id, l_rr_node_route_inf[sink_inode].backward_path_cost, current_heap_item.backward_path_cost, l_rr_node_route_inf[sink_inode].path_cost, current_heap_item.cost, l_rr_node_route_inf[sink_inode].prev_node, l_rr_node_route_inf[sink_inode].prev_edge);*/
				/*}*/
				/*l_rr_node_route_inf_lock.unlock();*/

			}
			zlog_debug(route_inner_log, "[%d] Worker done\n", item.worker_id);
		}
};

void check_route_tree_internal(t_rt_node *rt_node, set<int> &visited_sinks)
{
	if (rr_node[rt_node->inode].type == SINK) {
		assert(visited_sinks.find(rt_node->inode) == visited_sinks.end());
		visited_sinks.insert(rt_node->inode);
	}
	t_linked_rt_edge *edge = rt_node->u.child_list;
	while (edge) {
		check_route_tree_internal(edge->child, visited_sinks);
		edge = edge->next;
	}
}

void check_route_tree(t_rt_node *root, int inet)
{
	set<int> visited_sinks;
	check_route_tree_internal(root, visited_sinks);
	assert(visited_sinks.size() == clb_net[inet].num_sinks);
	for (int i = 1; i <= clb_net[inet].num_sinks; ++i) {
		assert(visited_sinks.find(net_rr_terminals[inet][i]) != visited_sinks.end());
	}
}

static boolean tbb_parallel_timing_driven_route_net(
		int inet,

		/* route parameters */
		int iter,
		const t_router_opts *opts,
		float pres_fac,

		/* preallocated for efficiency */
		float *pin_criticality,
		int *sink_order,
		t_rr_node_route_inf *l_rr_node_route_inf,
		struct s_rt_node **l_rr_node_to_rt_node,

		/* outputs */
		t_net_route *net_route,
		t_net_timing *net_timing,
		timeval *time)
{

	/* Returns TRUE as long is found some way to hook up this net, even if that *
	 * way resulted in overuse of resources (congestion).  If there is no way   *
	 * to route this net, even ignoring congestion, it returns FALSE.  In this  *
	 * case the rr_graph is disconnected and you can give up. If slacks = NULL, *
	 * give each net a dummy criticality of 0.									*/

	int ipin, num_sinks, itarget, target_pin, target_node, inode;
	float target_criticality, old_tcost, new_tcost, largest_criticality,
		old_back_cost, new_back_cost;
	t_rt_node *rt_root;
	struct s_heap current;
	struct s_trace *new_route_start_tptr;
	int highfanout_rlim;

	/* Rip-up any old routing. */
	thread_safe_pathfinder_update_one_cost(inet, net_route->l_trace_head, -1, pres_fac, 0, 0, false);
	thread_safe_free_traceback(&net_route->l_trace_head, &net_route->l_trace_tail);

	for (int i = 0; i < num_rr_nodes; ++i) {
		l_rr_node_to_rt_node[i] = NULL;
	}

	for (ipin = 1; ipin <= clb_net[inet].num_sinks; ipin++) { 
		if (!net_timing) {
			/* Use criticality of 1. This makes all nets critical.  Note: There is a big difference between setting pin criticality to 0
			compared to 1.  If pin criticality is set to 0, then the current path delay is completely ignored during routing.  By setting
			pin criticality to 1, the current path delay to the pin will always be considered and optimized for */
			pin_criticality[ipin] = 1.0;
		} else { 
#ifdef PATH_COUNTING
			/* Pin criticality is based on a weighted sum of timing and path criticalities. */	
			pin_criticality[ipin] =		 ROUTE_PATH_WEIGHT	* net_timing->path_criticality[ipin]
								  + (1 - ROUTE_PATH_WEIGHT) * net_timing->timing_criticality[ipin]; 
#else
			/* Pin criticality is based on only timing criticality. */
			pin_criticality[ipin] = net_timing->timing_criticality[ipin];
#endif
			/* Currently, pin criticality is between 0 and 1. Now shift it downwards 
			by 1 - max_criticality (max_criticality is 0.99 by default, so shift down
			by 0.01) and cut off at 0.  This means that all pins with small criticalities 
			(<0.01) get criticality 0 and are ignored entirely, and everything
			else becomes a bit less critical. This effect becomes more pronounced if
			max_criticality is set lower. */
			assert(pin_criticality[ipin] > -0.01 && pin_criticality[ipin] < 1.01);
			pin_criticality[ipin] = std::max(pin_criticality[ipin] - (1.0 - opts->max_criticality), 0.0);

			/* Take pin criticality to some power (1 by default). */
			pin_criticality[ipin] = pow(pin_criticality[ipin], opts->criticality_exp);
			
			/* Cut off pin criticality at max_criticality. */
			pin_criticality[ipin] = std::min(pin_criticality[ipin], opts->max_criticality);
		}
	}

	num_sinks = clb_net[inet].num_sinks;
	heapsort(sink_order, pin_criticality, num_sinks, 0);

	/* Update base costs according to fanout and criticality rules */

	/* we don't need this for parallel route */
	/* base costs are calculated dynamically based on num_sinks */
/*	largest_criticality = pin_criticality[sink_order[1]];*/
	//update_rr_base_costs(inet, largest_criticality);

	thread_safe_mark_ends(inet, l_rr_node_route_inf); /* Only needed to check for multiply-connected SINKs */

	vector<int> modified_route_tree;

	rt_root = init_route_tree_to_source(net_rr_terminals[inet][0], l_rr_node_to_rt_node, modified_route_tree);

	std::vector<int> modified_inodes;
	std::mutex modified_inodes_lock;
	char buffer[256];

	net_route->num_heap_pushes = 0;
	net_route->num_heap_pops = 0;

	timeval t1, t2;
	assert(!gettimeofday(&t1, NULL));

	for (int i = 0; i < num_rr_nodes; ++i) {
		assert(l_rr_node_route_inf[i].transactions.size() == 0);
	}

	tbb::atomic<int> **visit_count = new tbb::atomic<int>*[nx+2];
	for (int x = 0; x < nx+2; ++x) {
		visit_count[x] = new tbb::atomic<int>[ny+2];
	}

	sprintf_rr_node(net_rr_terminals[inet][0], buffer);
	zlog_debug(route_inner_log, "Source: %s\n", buffer);

	for (itarget = 1; itarget <= num_sinks; itarget++) {
		for (int i = 0; i < num_rr_nodes; ++i) {
			/*if (l_rr_node_to_rt_node[i]) {*/
			/*l_rr_node_route_inf[i].path_cost = HUGE_POSITIVE_FLOAT;*/
			/*} else {*/
			/*assert(l_rr_node_route_inf[i].path_cost > 0.999*HUGE_POSITIVE_FLOAT);*/
			/*}*/
			/*assert(l_rr_node_route_inf[i].path_cost > 0.99*HUGE_POSITIVE_FLOAT);*/
			if (!(l_rr_node_route_inf[i].path_cost > 0.99*HUGE_POSITIVE_FLOAT)) {
				sprintf_rr_node(i, buffer);
				zlog_fatal(route_inner_log, "%s doesn't have HUGE_POSITIVE_FLOAT\n", buffer);
			}
			assert(l_rr_node_route_inf[i].path_cost > 0.99*HUGE_POSITIVE_FLOAT);
			assert(l_rr_node_route_inf[i].modified == false);
		}

		target_pin = sink_order[itarget];
		target_node = net_rr_terminals[inet][target_pin];

		sprintf_rr_node(target_node, buffer);
		zlog_debug(route_inner_log, "Sink %d: %s\n", itarget, buffer);

		target_criticality = pin_criticality[target_pin];

		highfanout_rlim = mark_node_expansion_by_bin(inet, target_node,
				rt_root);

		int num_threads = opts->num_threads;

		/*std::priority_queue<struct s_heap> *heaps = new std::priority_queue<struct s_heap>[num_threads];*/

		std::priority_queue<struct s_heap> temp_heap;

		zlog_debug(route_inner_log, "Adding route tree to heap\n");
		add_route_tree_to_heap(rt_root, target_node, target_criticality,
				opts->astar_fac, l_rr_node_route_inf, inet, temp_heap, net_route->num_heap_pushes);

		if (temp_heap.empty()) { /* Infeasible routing.  No possible path for net. */
			zlog_fatal(route_inner_log, "Cannot route net #%d (%s) to sink #%d -- no possible path.\n",
					   inet, clb_net[inet].name, itarget);
			for (const auto &modified : modified_inodes) {
				l_rr_node_route_inf[modified].path_cost = HUGE_POSITIVE_FLOAT;
				/*l_rr_node_route_inf[modified].prev_node = NO_PREVIOUS;*/
			}
/*			free_route_tree(rt_root);*/
			return (FALSE);
		}

		std::mutex spawn_lock;
		std::mutex push_lock;
		vector<struct s_heap> found_sinks;
		std::mutex found_sinks_lock;
		std::mutex l_rr_node_route_inf_lock;
		int worker_id = 0;

		for (int x = 0; x < nx+2; ++x) {
			for (int y = 0; y < ny+2; ++y) {
				visit_count[x][y] = 0;
			}
		}

		vector<RouteWorker2::argument_type> items;
		/*int num_spawned_workers = 0;*/
		int heap_index = 0;
		/*while (!temp_heap.empty() && num_spawned_workers < num_threads) {*/
		/*items.push_back({ temp_heap.top(), num_spawned_workers });*/
		/*temp_heap.pop();*/
		/*++num_spawned_workers;*/
		/*}*/
		items.push_back({ temp_heap, worker_id++ });
		int num_spawned_workers = items.size();

		RouteWorker2 worker(
				l_rr_node_route_inf_lock,
				found_sinks_lock, found_sinks, 
				spawn_lock, heap_index, num_spawned_workers, worker_id,
				modified_inodes, modified_inodes_lock,
				push_lock, net_route->num_heap_pushes
				);
		/*worker.heaps = heaps;*/
		worker.visit_count = visit_count;
		worker.l_rr_node_route_inf = l_rr_node_route_inf;
		worker.l_rr_node_to_rt_node = l_rr_node_to_rt_node;
		worker.target_node = target_node;
		worker.inet = inet;
		worker.bend_cost = opts->bend_cost;
		worker.target_criticality = target_criticality; 
		worker.astar_fac = opts->astar_fac;
		worker.highfanout_rlim = highfanout_rlim;
		worker.num_threads = num_threads;

		tbb::parallel_do(items.begin(), items.end(), worker);

		/* NB:  In the code below I keep two records of the partial routing:  the   *
		 * traceback and the route_tree.  The route_tree enables fast recomputation *
		 * of the Elmore delay to each node in the partial routing.  The traceback  *
		 * lets me reuse all the routines written for breadth-first routing, which  *
		 * all take a traceback structure as input.  Before this routine exits the  *
		 * route_tree structure is destroyed; only the traceback is needed at that  *
		 * point.                                                                   */

		sprintf(buffer,"heatmap_net_%d_sink_%d.txt", inet, itarget);
		FILE *heat_map = fopen(buffer, "w");
		for (int x = 0; x < nx+2; ++x) {
			for (int y = 0; y < ny+2; ++y) {
				for (int i = 0; i < visit_count[x][y]; ++i) {
					fprintf(heat_map, "%d %d\n", x, y);
				}
			}
		}
		fclose(heat_map);

		assert(found_sinks.size() > 0);
		/*assert(num_found_sink == 1);*/

		float min_b_cost = HUGE_POSITIVE_FLOAT;
		for (const auto &sink : found_sinks) {
			if (sink.backward_path_cost < min_b_cost) {
				min_b_cost = sink.backward_path_cost;
				current = sink;
			}
		}

		inode = current.index;
		if (inode != target_node) {
			char buffer2[256];
			sprintf_rr_node(inode, buffer);
			sprintf_rr_node(target_node, buffer2);
			zlog_fatal(route_inner_log, "INODE %s != TARGET %s\n", buffer, buffer2);
		}
		assert(inode == target_node);
		l_rr_node_route_inf[inode].target_flag--; /* Connected to this SINK. */

		new_route_start_tptr = thread_safe_update_traceback(&net_route->l_trace_head, &net_route->l_trace_tail, &current, l_rr_node_route_inf);
		/*rt_node_of_sink[target_pin] = thread_safe_update_route_tree(&current, l_rr_node_route_inf, l_rr_node_to_rt_node);*/
		t_rt_node *sink_rt_node = NULL;//thread_safe_update_route_tree(&current, l_rr_node_route_inf, l_rr_node_to_rt_node);
		assert(sink_rt_node == l_rr_node_to_rt_node[current.index]);
/*		free_heap_data(current);*/
		thread_safe_pathfinder_update_one_cost(inet, new_route_start_tptr, 1, pres_fac, 0, 0, false);

		/* clear the heap */
		/*for (int i = 0; i < num_threads; ++i) {*/
			/*heaps[i] = std::priority_queue<struct s_heap>();*/
		/*}*/
		/*for (int modified = 0; modified < num_rr_nodes; ++modified) {*/
		for (const auto &modified : modified_inodes) {
			/*if (!l_rr_node_to_rt_node[modified]) {*/
				l_rr_node_route_inf[modified].path_cost = HUGE_POSITIVE_FLOAT;
				l_rr_node_route_inf[modified].modified = false;
			/*}*/
			/*l_rr_node_route_inf[modified].prev_node = NO_PREVIOUS;*/
		}
		modified_inodes.clear();
	}

	assert(!gettimeofday(&t2, NULL));
	time->tv_sec = t2.tv_sec - t1.tv_sec;
	time->tv_usec = t2.tv_usec - t1.tv_usec;

	/*sprintf(buffer, "transactions/trans_iter_%d_net_%d.txt", iter, inet);*/
	/*FILE *file = fopen(buffer, "w");*/
	/*assert(file);*/
	for (int i = 0; i < num_rr_nodes; ++i) {
		if (l_rr_node_route_inf[i].transactions.size() > 0) {
			/*fprintf(file, "Transactions for node %d:\n", i);*/
			/*for (const auto &trans : l_rr_node_route_inf[i].transactions) {*/
				/*fprintf(file, "sink: %6d old_b: %10.4g new_b: %10.4g old_t: %10.4g new_t: %10.4g prev_node: %6d prev_edge: %2d written: %d\n", trans.target_node, trans.old_bcost, trans.new_bcost, trans.old_tcost, trans.new_tcost, trans.prev_node, trans.prev_edge, trans.written ? 1 : 0);*/
			/*}*/
			/*fprintf(file, "\n");*/

			l_rr_node_route_inf[i].transactions.clear();
		}
	}
	/*fclose(file);*/

	/*FILE *prev_index_file = fopen("prev_index.txt", "r");*/
	/*assert(prev_index_file);*/
	/*int prev_index;*/
	/*fscanf(prev_index_file, "%d", &prev_index);*/
	/*fclose(prev_index_file);*/

	/*prev_index_file = fopen("prev_index.txt", "w");*/
	/*fprintf(prev_index_file, "%d", prev_index+1);*/
	
	/*sprintf(buffer, "lol.%d.txt", prev_index+1);*/
	/*print_one_route(buffer, inet, net_route, l_rr_node_to_rt_node);*/
	print_one_route("lol.txt", inet, net_route, l_rr_node_to_rt_node);

	for (int x = 0; x < nx+2; ++x) {
		delete [] visit_count[x];
	}
	delete [] visit_count;

	/*exit(0);*/

	/* For later timing analysis. */

	update_net_delays_from_route_tree_new(net_timing, l_rr_node_to_rt_node, inet);
	check_net_route(net_route->l_trace_head, inet);
	check_route_tree(rt_root, inet);
	/*free_route_tree_new(rt_root);*/
	return (TRUE);
}

static void timing_driven_check_net_delays(t_net_timing *net_timing) {

	/* Checks that the net delays computed incrementally during timing driven    *
	 * routing match those computed from scratch by the net_delay.c module.      */

	int inet, ipin;
	float **net_delay_check;

	t_chunk list_head_net_delay_check_ch = {NULL, 0, NULL};

	/*struct s_linked_vptr *ch_list_head_net_delay_check;*/

	net_delay_check = alloc_net_delay(&list_head_net_delay_check_ch, clb_net,
			num_nets);
	load_net_delay_from_routing(net_delay_check, clb_net, num_nets);

	for (inet = 0; inet < num_nets; inet++) {
		for (ipin = 1; ipin <= clb_net[inet].num_sinks; ipin++) {
			if (net_delay_check[inet][ipin] == 0.) { /* Should be only GLOBAL nets */
				if (fabs(net_timing[inet].delay[ipin]) > ERROR_TOL) {
					vpr_printf(TIO_MESSAGE_ERROR, "in timing_driven_check_net_delays: net %d pin %d.\n",
							inet, ipin);
					vpr_printf(TIO_MESSAGE_ERROR, "\tIncremental calc. net_delay is %g, but from scratch net delay is %g.\n",
							net_timing[inet].delay[ipin], net_delay_check[inet][ipin]);
					exit(1);
				}
			} else {
				if (fabs(1.0 - net_timing[inet].delay[ipin] / net_delay_check[inet][ipin]) > ERROR_TOL) {
					vpr_printf(TIO_MESSAGE_ERROR, "in timing_driven_check_net_delays: net %d pin %d.\n",
							inet, ipin);
					vpr_printf(TIO_MESSAGE_ERROR, "\tIncremental calc. net_delay is %g, but from scratch net delay is %g.\n",
							net_timing[inet].delay[ipin], net_delay_check[inet][ipin]);
					exit(1);
				}
			}
		}
	}

	free_net_delay(net_delay_check, &list_head_net_delay_check_ch);
	vpr_printf(TIO_MESSAGE_INFO, "Completed net delay value cross check successfully.\n");
}

static boolean try_fine_grained_parallel_timing_driven_route(struct s_router_opts router_opts,
		t_net_timing *net_timing, t_ivec ** clb_opins_used_locally, boolean timing_analysis_enabled, int round, std::mt19937 &mt) {

	/* Timing-driven routing algorithm.  The timing graph (includes slack)   *
	 * must have already been allocated, and net_delay must have been allocated. *
	 * Returns TRUE if the routing succeeds, FALSE otherwise.                    */

	int itry, inet, ipin, i, bends, wirelength, total_wirelength, available_wirelength, 
		segments, *net_index, *sink_order /* [1..max_pins_per_net-1] */;
	boolean success, is_routable, rip_up_local_opins;
	float *pin_criticality /* [1..max_pins_per_net-1] */, pres_fac, *sinks, 
		critical_path_delay, init_timing_criticality_val;
	t_rt_node **rt_node_of_sink; /* [1..max_pins_per_net-1] */
	/*clock_t begin,end;*/
	timeval t1, t2;
	sinks = (float*)my_malloc(sizeof(float) * num_nets);
	net_index = (int*)my_malloc(sizeof(int) * num_nets);

	for (i = 0; i < num_nets; i++) {
		sinks[i] = clb_net[i].num_sinks;
		net_index[i] = i;
	}
	heapsort(net_index, sinks, num_nets, 1);

	for (int inet = 0; inet < num_nets; ++inet) {
		for (int ipin = 0; ipin <= clb_net[inet].num_sinks; ++ipin) {
			int x = block[clb_net[inet].node_block[ipin]].x;
			int y = block[clb_net[inet].node_block[ipin]].y
				+ block[clb_net[inet].node_block[ipin]].type->pin_height[clb_net[inet].node_block_pin[ipin]];
			int inode = net_rr_terminals[inet][ipin];


			clb_net[inet].sorted_sinks.emplace_back(x, y, inode);
		}
		std::sort(clb_net[inet].sorted_sinks.begin(), clb_net[inet].sorted_sinks.end(), compare_sink); 
		zlog_debug(net_sort_log, "Net %d has sorted sinks:\n");
		for (const auto &s : clb_net[inet].sorted_sinks) {
			zlog_debug(net_sort_log, "%d %d %d\n", s.x, s.y, s.inode);
		}
	}
	/*TODO: need to alloc route tree timing struct */	

	/* these structures are allocated per thread, so no need to do it here */
/*	alloc_timing_driven_route_structs(&pin_criticality, &sink_order,*/
/*			&rt_node_of_sink);*/

	/* First do one routing iteration ignoring congestion to	
	get reasonable net delay estimates. Set criticalities to 1 
	when timing analysis is on to optimize timing, and to 0 
	when timing analysis is off to optimize routability. */

	if (timing_analysis_enabled) {
		init_timing_criticality_val = 1.;
	} else {
		init_timing_criticality_val = 0.;
	}

	for (inet = 0; inet < num_nets; inet++) {
		if (clb_net[inet].is_global == FALSE) {
			for (ipin = 1; ipin <= clb_net[inet].num_sinks; ipin++)
				net_timing[inet].timing_criticality[ipin] = init_timing_criticality_val;
#ifdef PATH_COUNTING
				net_timing[inet].path_criticality[ipin] = init_timing_criticality_val;
#endif		
		} else { 
			/* Set delay of global signals to zero. Non-global net 
			delays are set by update_net_delays_from_route_tree() 
			inside timing_driven_route_net(), which is only called
			for non-global nets. */
			for (ipin = 1; ipin <= clb_net[inet].num_sinks; ipin++) {
				net_timing[inet].delay[ipin] = 0.;
			}
		}
	}

	pres_fac = router_opts.first_iter_pres_fac; /* Typically 0 -> ignore cong. */


    int num_threads = router_opts.num_threads;

	int *prev_occ = new int[num_rr_nodes];

	for (int i = 0; i < num_rr_nodes; ++i) {
		prev_occ[i] = 0;
	}

	t_net_route *prev_net_route = new t_net_route[num_nets];
	for (int i = 0; i < num_nets; ++i) {
		prev_net_route[i].l_trace_head = NULL;
		prev_net_route[i].l_trace_tail = NULL;
	}

	for (int i = 0; i < num_rr_nodes; ++i) {
		/*assert(rr_node[i].occ == 0);*/
		/*assert(rr_node[i].occupant_net_id.empty());*/
		/*assert(rr_node[i].driver_nets.empty());*/
		/*assert(rr_node[i].num_reservation == 0);*/
		rr_node[i].occ = 0;
		rr_node[i].occupant_net_id.clear();
		rr_node[i].driver_nets.clear();
		rr_node[i].num_reservation = 0;
		rr_node[i].acc_cost = 1;
		rr_node[i].pres_cost = 1;
		rr_node[i].weighted_pres_cost = 1;
	}

	alloc_per_thread_timing_driven_route_structs(&pin_criticality, &sink_order);

	t_rr_node_route_inf *l_rr_node_route_inf = new t_rr_node_route_inf[num_rr_nodes];

	for (int inode = 0; inode < num_rr_nodes; inode++) {
		l_rr_node_route_inf[inode].prev_node = NO_PREVIOUS;
		l_rr_node_route_inf[inode].prev_edge = NO_PREVIOUS;
		l_rr_node_route_inf[inode].path_cost = HUGE_POSITIVE_FLOAT;
		l_rr_node_route_inf[inode].backward_path_cost = HUGE_POSITIVE_FLOAT;
		l_rr_node_route_inf[inode].target_flag = 0;
		l_rr_node_route_inf[inode].modified = false;
	}

	t_rt_node **l_rr_node_to_rt_node = new t_rt_node*[num_rr_nodes];

	t_net_route *net_route = new t_net_route[num_nets];

	for (int i = 0; i < num_nets; ++i) {
		net_route[i].l_trace_head = NULL;
		net_route[i].l_trace_tail = NULL;
	}

	timeval *net_route_time = new timeval[num_nets];

	char buffer[256];
	sprintf(buffer, "%d", -1);
	assert(!zlog_put_mdc("tid", buffer));

	for (itry = 1; itry <= router_opts.max_router_iterations; itry++) {
		/*begin = clock();*/

		sprintf(buffer, "%d", itry);
		assert(!zlog_put_mdc("iter", buffer));

		zlog_info(route_outer_log, "\n");
		zlog_info(route_outer_log, "Routing iteration: %d\n", itry);
		
		/*struct mach_timebase_info tbi;*/
		/*mach_timebase_info(&tbi);*/
		
		/*sprintf(buffer, "route_inner_%d_%d_test.txt", info->thread_index, *info->iter);*/
		/*FILE *file = fopen(buffer, "w");*/
		/*assert(file);*/

		assert(!gettimeofday(&t1, NULL));

		for (i = 0; i < num_nets; i++) {
			inet = net_index[i];
			if (clb_net[inet].is_global == FALSE) { /* Skip global nets. */
				timeval route_time_start, route_time_end;
				zlog_info(route_inner_log, "Routing net %d\n", inet);

				/*fprintf(file, "%d Routing net %d real %d final %d final_no_route: %d\n", info->thread_index, inet, real_net_index, final ? 1 : 0, final_no_route ? 1 : 0);*/
				/*fflush(file);*/

				is_routable = tbb_parallel_timing_driven_route_net(
						inet,

						itry,
						&router_opts,
						pres_fac,

						pin_criticality,
						sink_order,
						l_rr_node_route_inf,
						l_rr_node_to_rt_node,

						&net_route[inet],
						&net_timing[inet],
						&net_route_time[inet]
						);

				if (!is_routable) {
					zlog_fatal(route_outer_log, "Routing failed.\n");
					/*			free_timing_driven_route_structs(pin_criticality,*/
					/*					sink_order, rt_node_of_sink);*/
					/*			free(net_index);*/
					/*			free(sinks);*/
					exit(-1);
				}

				time_t rawtime;
				struct tm * timeinfo;
				char buffer [80];

				time (&rawtime);
				timeinfo = localtime (&rawtime);

				strftime(buffer, 256, "route_of_net_%d.%H:%M:%S.txt", timeinfo);
				print_one_route(buffer, inet, &net_route[inet], l_rr_node_to_rt_node);
			}
		}

		assert(!gettimeofday(&t2, NULL));

		double route_time = t2.tv_sec - t1.tv_sec;
		route_time += (double)(t2.tv_usec - t1.tv_usec)/1e6;
		zlog_info(route_outer_log, "Route time: %g\n", route_time);

		for (int i = 0; i < num_nets; ++i) {
			trace_head[i] = net_route[i].l_trace_head;
			trace_tail[i] = net_route[i].l_trace_tail;
		}

		if (itry == 1) {
			/* Early exit code for cases where it is obvious that a successful route will not be found 
			 Heuristic: If total wirelength used in first routing iteration is X% of total available wirelength, exit
			 */
			total_wirelength = 0;
			available_wirelength = 0;

			for (i = 0; i < num_rr_nodes; i++) {
				if (rr_node[i].type == CHANX || rr_node[i].type == CHANY) {
					available_wirelength += 1 + rr_node[i].xhigh
							- rr_node[i].xlow + rr_node[i].yhigh
							- rr_node[i].ylow;
				}
			}

			for (inet = 0; inet < num_nets; inet++) {
				if (clb_net[inet].is_global == FALSE
						&& clb_net[inet].num_sinks != 0) { /* Globals don't count. */
					get_num_bends_and_length(inet, &bends, &wirelength,
							&segments);

					total_wirelength += wirelength;
				}
			}
			zlog_info(route_outer_log, "Wire length after first iteration %d, total available wire length %d, ratio %g\n",
					total_wirelength, available_wirelength,
					(float) (total_wirelength) / (float) (available_wirelength));
			if ((float) (total_wirelength) / (float) (available_wirelength)> FIRST_ITER_WIRELENTH_LIMIT) {
				zlog_info(route_outer_log, "Wire length usage ratio exceeds limit of %g, fail routing.\n",
						FIRST_ITER_WIRELENTH_LIMIT);
/*				free_timing_driven_route_structs(pin_criticality, sink_order,*/
/*						rt_node_of_sink);*/
				free(net_index);
				free(sinks);
				return FALSE;
			}
		}

		/* Make sure any CLB OPINs used up by subblocks being hooked directly     *
		 * to them are reserved for that purpose.                                 */

		if (itry == 1)
			rip_up_local_opins = FALSE;
		else
			rip_up_local_opins = TRUE;

		thread_safe_reserve_locally_used_opins(pres_fac, rip_up_local_opins,
				clb_opins_used_locally, -1);

		// ------------------------ //
		
		const float prev_acc_cost_weight = 0.9;
		for (int i = 0; i < num_rr_nodes; ++i) {
			rr_node[i].weighted_pres_cost +=
				prev_acc_cost_weight * rr_node[i].weighted_pres_cost
				+ pres_fac * (rr_node[i].occ + 1 - rr_node[i].capacity);
		}	

		// ------------------------ //

		char buffer[256];

		sprintf(buffer, "heap_pushes_%d_%d.txt", round, itry);
		FILE *heap_pushes = fopen(buffer, "w");
		assert(heap_pushes);
		vector<pair<int, int>> sorted_num_heap_pushes;
		for (int i = 0; i < num_nets; ++i) {
			/*zlog_info(route_outer_log, "*/
			if (!clb_net[i].is_global) {
				sorted_num_heap_pushes.push_back(make_pair(net_route[i].num_heap_pushes, i));
			}
		}
		std::sort(sorted_num_heap_pushes.begin(), sorted_num_heap_pushes.end());

		zlog_debug(route_outer_log, "Sorted num_heap_pushes:\n");
		/*for (const auto &heap_push : sorted_num_heap_pushes) {*/
		for (auto iter = sorted_num_heap_pushes.rbegin(); iter != sorted_num_heap_pushes.rend(); ++iter) {
			int inet = iter->second;
			int area = get_bounding_box_area(&route_bb[inet]);

			double elapsed = net_route_time[inet].tv_sec;
			elapsed += (double)net_route_time[inet].tv_usec/1e6;

			zlog_debug(route_outer_log, "Net: %d BB Area: %d Num sinks: %d Num pushes: %d Num pops: %d Pop/Push: %g Route time: %g\n", inet, area, clb_net[inet].num_sinks, iter->first, net_route[inet].num_heap_pops, (double)net_route[inet].num_heap_pops*100/iter->first, elapsed);
			fprintf(heap_pushes, "%d %d %d\n", area, clb_net[inet].num_sinks, iter->first);
		}
		fclose(heap_pushes);

		// ------------------------ //

		int num_different_occs = 0;
		for (int i = 0; i < num_rr_nodes; ++i) {
			if (rr_node[i].occ != prev_occ[i]) {
				++num_different_occs;
			}
			prev_occ[i] = rr_node[i].occ;
		}
		zlog_info(route_outer_log, "Number of nodes with changed occ: %d\n", num_different_occs);

		// ------------------------ //

		sprintf(buffer, "routes_%d_%d.txt", round, itry);
		/*print_route(buffer, net_route, sink_order);*/

		int num_congested_nodes_by_type[NUM_RR_TYPES] = { 0 };

		std::set<int> congested_nets;
		for (int inode = 0; inode < num_rr_nodes; ++inode) {
			if (rr_node[inode].occ > rr_node[inode].capacity) {
				sprintf_rr_node(inode, buffer);
				zlog_debug(congestion_log, "Node %s is congested by %d nets:\n", buffer, rr_node[inode].occupant_net_id.size());
				for (const auto &n : rr_node[inode].occupant_net_id) {
					zlog_debug(congestion_log, "%d\n", n);
					congested_nets.insert(n);
				}
				++num_congested_nodes_by_type[rr_node[inode].type];
			}
			assert(rr_node[inode].occ - rr_node[inode].num_reservation == rr_node[inode].occupant_net_id.size());
		}
		for (int i = 0; i < NUM_RR_TYPES; ++i) {
			if (num_congested_nodes_by_type[i] > 0) {
				zlog_info(route_outer_log, "%d %s(s) congested\n", num_congested_nodes_by_type[i], rr_types[i]);
				zlog_info(congestion_log, "%d %s(s) congested\n", num_congested_nodes_by_type[i], rr_types[i]);
			}
		}
		zlog_info(route_outer_log, "%d / %d (%.2f) congested nets\n", congested_nets.size(), num_nets, (float)congested_nets.size()*100/num_nets);
		zlog_info(congestion_log, "%d / %d (%.2f) congested nets\n", congested_nets.size(), num_nets, (float)congested_nets.size()*100/num_nets);

		// ------------------------ //

		struct net_stats {
			int area;
			int num_segments;
			int num_congested_segments;
			int index;
		};

		for (int i = 0; i < num_rr_nodes; ++i) {
			rr_node[i].driver_nets.clear();
		}

		std::vector<net_stats> congested_net_stats;
		int congested_nets_total_area = 0;
		for (const auto &net : congested_nets) {
			int area = (route_bb[net].xmax - route_bb[net].xmin) * (route_bb[net].ymax - route_bb[net].ymin);
			congested_nets_total_area += area;

			net_stats stats;
			stats.index = net;
			stats.area = area;

			struct s_trace *tptr = net_route[net].l_trace_head;
			struct s_trace *prev_tptr = NULL;
			assert(tptr != NULL);

			stats.num_segments = 0;
			stats.num_congested_segments = 0;
			while (tptr != NULL) {
				int inode = tptr->index;

				if (rr_node[inode].occ > rr_node[inode].capacity) {
					++stats.num_congested_segments;
					if (prev_tptr) {
						rr_node[inode].driver_nets.insert(make_pair(net, prev_tptr->index));
					}
				}

				prev_tptr = tptr;
				tptr = tptr->next;

				if (rr_node[inode].type == SINK && tptr != NULL) {
					prev_tptr = tptr;
					tptr = tptr->next; /* skip the joining segment */
				} 

				++stats.num_segments;
			}

			congested_net_stats.push_back(stats);
		}

		for (int i = 0; i < num_rr_nodes; ++i) {
			if (rr_node[i].occ > rr_node[i].capacity) {
				assert(rr_node[i].driver_nets.size() == rr_node[i].occupant_net_id.size());
				for (const auto &driver : rr_node[i].driver_nets) {
					assert(rr_node[i].occupant_net_id.find(driver.first) !=
							rr_node[i].occupant_net_id.end());
				}
			}
		}

		std::sort(congested_net_stats.begin(), congested_net_stats.end(),
				[](const net_stats &a, const net_stats &b) -> bool
				{
					return (float)a.num_congested_segments/a.num_segments < (float)b.num_congested_segments/b.num_segments;
				});

		zlog_debug(congestion_log, "Sorted congested nets based on bb area:\n");
		for (const auto &c : congested_net_stats) {
			zlog_debug(congestion_log, "Net %6d Area: %6d (%5.1f) Congested segments: %6d /%6d (%5.1f)\n", c.index, c.area, (float)c.area*100/(nx*ny), c.num_congested_segments, c.num_segments, (float)c.num_congested_segments*100/c.num_segments);
		}

		double congestion_nets_average_area = (double)congested_nets_total_area/congested_nets.size();
		zlog_info(congestion_log, "Congestion nets average area: %g/%d (%g)\n", congestion_nets_average_area, nx*ny, congestion_nets_average_area*100/(nx*ny));

		sprintf(buffer, "rr_occ_%d.txt", itry);
		FILE *rr_occ = fopen(buffer, "w");
		for (int i = 0; i < num_rr_nodes; ++i) {
			fprintf(rr_occ, "%d %d\n", i, rr_node[i].occ);
		}
		fclose(rr_occ);

		if (itry > 1) {
			/*printf("inet: %d\n", inet);*/
			/*prev_net_route[inet].l_trace_head = tinfo[i].net_route[inet].l_trace_head;*/
			/*prev_net_route[inet].l_trace_tail = tinfo[i].net_route[inet].l_trace_tail;*/
		}

		/* Pathfinder guys quit after finding a feasible route. I may want to keep *
		 * going longer, trying to improve timing.  Think about this some.         */

		success = feasible_routing();
		if (success) {
			zlog_info(route_outer_log, "Successfully routed after %d routing iterations.\n", itry);
/*			free_timing_driven_route_structs(pin_criticality, sink_order, rt_node_of_sink);*/
#ifdef DEBUG
			timing_driven_check_net_delays(net_timing);
#endif
			free(net_index);
			free(sinks);
			return (TRUE);
		}

		if (itry == 1) {
			pres_fac = router_opts.initial_pres_fac;
			thread_safe_pathfinder_update_cost(pres_fac, 0.); /* Acc_fac=0 for first iter. */
		} else {
			pres_fac *= router_opts.pres_fac_mult;

			/* Avoid overflow for high iteration counts, even if acc_cost is big */
			pres_fac = std::min(pres_fac, static_cast<float>(HUGE_POSITIVE_FLOAT / 1e5));

			thread_safe_pathfinder_update_cost(pres_fac, router_opts.acc_fac);
		}

		if (timing_analysis_enabled) {		
			/* Update slack values by doing another timing analysis.                 *
			 * Timing_driven_route_net updated the net delay values.                 */

			load_timing_graph_net_delays_new(net_timing);

	#ifdef HACK_LUT_PIN_SWAPPING
			do_timing_analysis_new(net_timing, FALSE, TRUE, FALSE);
	#else
			do_timing_analysis_new(net_timing, FALSE, FALSE, FALSE);
	#endif

			/* Print critical path delay - convert to nanoseconds. */
			critical_path_delay = get_critical_path_delay();
			zlog_info(route_outer_log, "Critical path: %g ns\n", critical_path_delay);
		} else {
			/* If timing analysis is not enabled, make sure that the criticalities and the 	*
			 * net_delays stay as 0 so that wirelength can be optimized. 			*/
			
			for (inet = 0; inet < num_nets; inet++) {
				for (ipin = 1; ipin <= clb_net[inet].num_sinks; ipin++) {
					net_timing[inet].timing_criticality[ipin] = 0.;
#ifdef PATH_COUNTING 		
					net_timing[inet].path_criticality[ipin] = 0.; 		
#endif
					net_timing[inet].delay[ipin] = 0.;
				}
			}
		}
		
		/*end = clock();*/
		/*#ifdef CLOCKS_PER_SEC*/
			/*vpr_printf(TIO_MESSAGE_INFO, "Routing iteration took %g seconds.\n", (float)(end - begin) / CLOCKS_PER_SEC);*/
		/*#else*/
			/*vpr_printf(TIO_MESSAGE_INFO, "Routing iteration took %g seconds.\n", (float)(end - begin) / CLK_PER_SEC);*/
		/*#endif*/
		
		fflush(stdout);
	}

	zlog_fatal(route_outer_log, "Routing failed.\n");
/*	free_timing_driven_route_structs(pin_criticality, sink_order,*/
/*			rt_node_of_sink);*/
	free(net_index);
	free(sinks);
	return (FALSE);
}

boolean try_fine_grained_parallel_timing_driven_route_top_2(struct s_router_opts router_opts,
		t_net_timing *net_timing, t_ivec ** clb_opins_used_locally, boolean timing_analysis_enabled)
{
	/*test_distributed_work();*/
	/*exit(0);*/
	int num_repetitions = 1;
	std::mt19937 mt(1039);
	tbb::task_scheduler_init sch_init(router_opts.num_threads);

	for (int i = 0; i < num_repetitions; ++i) {
		try_fine_grained_parallel_timing_driven_route(router_opts, net_timing, clb_opins_used_locally, timing_analysis_enabled, i, mt);
	}
}
