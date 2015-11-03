#include <stdio.h>
#include <math.h>
#include <errno.h>
#include <time.h>
#include <sys/time.h>
#include <string.h>
#include <assert.h>
#include <pthread.h>
#include <fcntl.h>
#include <unistd.h>
#include <semaphore.h>
#include <queue>
#include <random>
#include <algorithm>
#include <functional>
#include <zlog.h>
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

static const char *rr_types[] =  {
	"SOURCE", "SINK", "IPIN", "OPIN", "CHANX", "CHANY", "INTRA_CLUSTER_EDGE"
};

/*void reserve_locally_used_opins(float pres_fac, boolean rip_up_local_opins,*/
		/*t_ivec ** clb_opins_used_locally);*/

boolean feasible_routing(void);

void print_route(char *name);

void update_net_delays_from_route_tree_new(t_net_timing *net_timing,
		t_rt_node ** rt_node_of_sink, int inet);

/******************** Subroutines local to route_timing.c ********************/

typedef struct s_inode_with_pos {
	s_inode_with_pos()
		: original_clb_net_index(-1), inode(-1), x(-1), y(-1)
	{
	}

	s_inode_with_pos(int original_clb_net_index, int inode, int x, int y) :
		original_clb_net_index(original_clb_net_index), inode(inode), x(x), y(y)
	{
	}

	int original_clb_net_index;
	int inode;
	int x;
	int y;
} t_inode_with_pos;

typedef struct s_simple_net {
	int index;
	t_inode_with_pos source;

	/*t_sink_direction current_sink_direction;*/
	vector<t_inode_with_pos> sorted_sinks_to_left_of_source_inode;
	vector<t_inode_with_pos> sorted_sinks_to_right_of_source_inode;
	vector<t_inode_with_pos> sorted_sinks;

	vector<t_inode_with_pos> sinks_in_direction[4];

	t_inode_with_pos current_source;
	int current_sink_index;
	t_inode_with_pos *current_sink;

	struct s_bb current_bounding_box;

	vector<struct s_simple_net *> non_overlapping_nets;
	vector<struct s_simple_net *> overlapping_nets;
	int current_non_overlapping_net_index;

	t_trace *trace_head;
	t_trace *trace_tail;

	t_rt_node **rr_node_to_rt_node;
} t_simple_net;

typedef struct s_sink_direction {
	bool horizontal_increasing;
	bool vertical_increasing;
	s_sink_direction(const t_simple_net &net) {
		horizontal_increasing = net.current_sink->x >= net.current_source.x;
		vertical_increasing = net.current_sink->y >= net.current_source.y;
	}
	bool operator==(const s_sink_direction &other) const {
		return other.horizontal_increasing == horizontal_increasing && other.vertical_increasing == vertical_increasing;
	}
	/*bool operator<(const s_sink_direction &other) const {*/
		/*if (*/
	/*}*/
} t_sink_direction;

typedef struct s_bb_overlap_cost {
	t_sink_direction direction;
	int area;
	/*bool operator<(const*/
} t_bb_overlap_cost;

static boolean simple_route_net(t_simple_net *net,
		t_router_opts *opts, float pres_fac,
		t_rr_node_route_inf *l_rr_node_route_inf, 
		t_net_timing *net_timing,
		int *num_heap_pushes);

static int get_max_pins_per_net(void);

static void add_route_tree_to_heap(t_rt_node * rt_node, int target_node,
		float target_criticality, float astar_fac, t_rr_node_route_inf *l_rr_node_route_inf,  std::priority_queue<struct s_heap> &heap, int inet);

static void timing_driven_expand_neighbours(const struct s_heap *current, int inet,
		float bend_cost, float astar_fac, float criticality_fac, int target_node,
		int highfanout_rlim, bool monotonic,
		int *num_heap_pushes,
		t_rr_node_route_inf *l_rr_node_route_inf,
		std::priority_queue<struct s_heap> &heap);

static float get_timing_driven_expected_cost(int inode, int target_node,
		float criticality_fac, float R_upstream, float base_cost_scale_factor);

static int get_expected_segs_to_target(int inode, int target_node,
		int *num_segs_ortho_dir_ptr);

static void timing_driven_check_net_delays(t_net_timing *net_timing);

static int mark_node_expansion_by_bin(int inet, int target_node,
		t_rt_node * rt_node);

static boolean parallel_timing_driven_route_net(int inet, t_router_opts *opts,
		float pres_fac,
		float *pin_criticality, int *sink_order,
		t_rt_node **l_rr_node_to_rt_node,  
		t_rr_node_route_inf *l_rr_node_route_inf, t_net_timing *net_timing,
		int *num_heap_pushes);

/************************ Subroutine definitions *****************************/

extern zlog_category_t *route_inner_log;
static zlog_category_t *next_net_log;
static zlog_category_t *route_outer_log;
static zlog_category_t *check_route_log;
static zlog_category_t *congestion_log;
static zlog_category_t *net_sort_log;
static zlog_category_t *net_partition_log;
static zlog_category_t *init_net_log;

void init_advanced_parallel_route_logging()
{
	next_net_log = zlog_get_category("next_net");
	assert(next_net_log);
	route_outer_log = zlog_get_category("route_outer");
	assert(route_outer_log);
	check_route_log = zlog_get_category("check_route");
	assert(check_route_log);
	congestion_log = zlog_get_category("congestion");
	assert(congestion_log);
	net_sort_log = zlog_get_category("net_sort");
	assert(net_sort_log);
	net_partition_log = zlog_get_category("net_partition");
	assert(net_partition_log);
	init_net_log = zlog_get_category("init_net");
	assert(init_net_log);
}

typedef struct s_next_net {
	int current_net;
	int num_local_nets;
	int num_local_nets_routed;
	pthread_mutex_t lock;
	bool reversed;
	int num_threads;
} t_next_net;

typedef struct thread_info {
	int thread_index;

	sem_t *start_route_sem;
	sem_t *complete_route_sem;
	my_pthread_barrier_t *barrier;

	int *net_index;
	t_router_opts *router_opts;
	float *pres_fac;
	t_net_timing *net_timing;
	int **sink_order; //[0..inet-1][1..num_sinks]

	bool successful;

	t_next_net *next_net;

	double total_stall_time;
	
} thread_info;

struct compare_to_source {
	const t_inode_with_pos &source;
	compare_to_source(const t_inode_with_pos &source) :
		source(source)
	{
	}

	bool operator()(const struct s_inode_with_pos &sink_a, const struct s_inode_with_pos &sink_b) const
	{
		return abs(source.x - sink_a.x) + abs(source.y - sink_a.y)
			< abs(source.x - sink_b.x) + abs(source.y - sink_b.y);
	}
};

static void write_hmetis_graph_file(const char *filename, t_simple_net *nets, int l_num_nets)
{
	FILE *file = fopen(filename, "w");
	assert(file);

	int num_edges = 0;
	for (int i = 0; i < l_num_nets; ++i) {
		num_edges += nets[i].overlapping_nets.size() > 0 ? 1 : 0;
	}

	/*num_edges, num_vertices*/
	/*fprintf(file, "%d %d 10\n", num_edges, l_num_nets);*/
	fprintf(file, "%d %d 10\n", num_edges, l_num_nets);
	int inet = 0;
	for (int inet = 0; inet < l_num_nets; ++inet) {
		fprintf(file, "%% Net %d overlaps with %d net(s)\n", inet, nets[inet].overlapping_nets.size());

		if (nets[inet].overlapping_nets.size() > 0) {
			fprintf(file, "%d ", inet+1);
			for (const auto &overlap : nets[inet].overlapping_nets) {
				fprintf(file, "%d ", overlap->index+1);
			}
			fprintf(file, "\n");
		}
	}

	for (int i = 0; i < l_num_nets; ++i) {
		int area = get_bounding_box_area(&nets[i].current_bounding_box);
		fprintf(file, "%d\n", area);
	}

	fclose(file);
}

static void run_hmetis(int num_partitions, const char *graph_filename)
{
	char buffer[256];
	sprintf(buffer, "./shmetis %s %d 1", graph_filename, num_partitions);
	FILE *pipe = popen(buffer, "r");
	while (!feof(pipe)) {
		fgets(buffer, 256, pipe);
		printf("%s", buffer);
	}
	pclose(pipe);
}

void dispatch_work(t_simple_net *l_nets, int l_num_nets, thread_info *info, int num_threads)
{
	std::sort(l_nets, l_nets + l_num_nets,
			[](const t_simple_net &net_a, const t_simple_net &net_b) -> bool {
				return net_a.non_overlapping_nets.size() < net_b.non_overlapping_nets.size();
			});
	bool *done = new bool[l_num_nets];
	int to_thread = 0;
	bool has_work = true;
	/*while (has_work) {*/
		/*for (int i = 0; i < l_num_nets; ++i) {*/
			/*if (l_nets[i].non_overlapping_nets[i]*/
		/*}*/
		/*to_thread = (to_thread+1) % num_threads;*/
	/*}*/
	delete [] done;
}

void init_nets_sink_in_same_direction(t_simple_net *nets)
{
	for (int inet = 0; inet < num_nets; ++inet) {
		t_simple_net *net = &nets[inet];

		net->index = inet;

		zlog_debug(init_net_log, "Net %d\n", inet);

		int source_x = block[clb_net[inet].node_block[0]].x;
		int source_y = block[clb_net[inet].node_block[0]].y
			+ block[clb_net[inet].node_block[0]].type->pin_height[clb_net[inet].node_block_pin[0]];

		net->source.original_clb_net_index = 0;
		net->source.inode = net_rr_terminals[inet][0];
		net->source.x = source_x;
		net->source.y = source_y;

		zlog_debug(init_net_log, "Original source (%d,%d)\n", source_x, source_y);

		net->current_source = net->source;

		/*for (int isink = 1; isink <= clb_net[inet].num_sinks; ++isink) {*/
			/*int sink_x = block[clb_net[inet].node_block[isink]].x;*/
			/*int sink_y = block[clb_net[inet].node_block[isink]].y*/
				/*+ block[clb_net[inet].node_block[isink]].type->pin_height[clb_net[inet].node_block_pin[isink]];*/

			/*if (sink_x > source_x) {*/
				/*net->sorted_sinks_to_right_of_source_inode.emplace_back(*/
						/*isink, net_rr_terminals[inet][isink], sink_x, sink_y);*/
			/*} else {*/
				/*net->sorted_sinks_to_left_of_source_inode.emplace_back(*/
						/*isink, net_rr_terminals[inet][isink], sink_x, sink_y);*/
			/*}*/

			/*zlog_debug(init_net_log, "Original sink (%d,%d)\n", sink_x, sink_y);*/
		/*}*/


			for (int isink = 1; isink <= clb_net[inet].num_sinks; ++isink) {
				int sink_x = block[clb_net[inet].node_block[isink]].x;
				int sink_y = block[clb_net[inet].node_block[isink]].y
					+ block[clb_net[inet].node_block[isink]].type->pin_height[clb_net[inet].node_block_pin[isink]];

				bool hi = sink_x >= source_x;
				bool vi = sink_y >= source_y;
				int direction = 0;
				direction |= hi ? 1 : 0;
				direction |= vi ? 2 : 0;

				net->sinks_in_direction[direction].emplace_back(
						isink, net_rr_terminals[inet][isink], sink_x, sink_y);

				zlog_debug(init_net_log, "Original sink (%d,%d)\n", sink_x, sink_y);
			}

		vector<pair<int, int>> num_sinks_in_direction;
		for (int dir = 0; dir < 4; ++dir) {
			num_sinks_in_direction.emplace_back(net->sinks_in_direction[dir].size(), dir);
		}
		std::sort(num_sinks_in_direction.begin(), num_sinks_in_direction.end());
		zlog_debug(init_net_log, "Sorted percentage number of sinks in direction: ");
		for (const auto &n : num_sinks_in_direction) {
			zlog_debug(init_net_log, "%.2f ", (float)n.first/clb_net[inet].num_sinks);
		}
		zlog_debug(init_net_log, "\n");

		/*for (const auto &s: net->sorted_sinks_to_left_of_source_inode) {*/
			/*zlog_debug(init_net_log, "Sorted sinks to left of source (%d,%d)\n", s.x, s.y);*/
		/*}*/
		/*for (const auto &s: net->sorted_sinks_to_right_of_source_inode) {*/
			/*zlog_debug(init_net_log, "Sorted sinks to right of source (%d,%d)\n", s.x, s.y);*/
		/*}*/

		net->current_sink_index = 0;
		net->current_sink = &net->sorted_sinks[net->current_sink_index];

		/*if (net->sorted_sinks_to_left_of_source_inode.size() == 0) {*/
			/*net->done_sinks_to_left_of_source_inode = true;*/
			/*net->current_sink = &net->sorted_sinks_to_right_of_source_inode[net->current_sink_index];*/
		/*} else {*/
			/*net->current_sink = &net->sorted_sinks_to_left_of_source_inode[net->current_sink_index];*/
		/*}*/

		const int bb_factor = 1;

		net->current_bounding_box.xmin = std::min(net->current_source.x, net->current_sink->x);
		net->current_bounding_box.xmax = std::max(net->current_source.x, net->current_sink->x);
		net->current_bounding_box.ymin = std::min(net->current_source.y, net->current_sink->y);
		net->current_bounding_box.ymax = std::max(net->current_source.y, net->current_sink->y);
		--net->current_bounding_box.xmin;
		--net->current_bounding_box.ymin;

		net->current_bounding_box.xmin = std::max(net->current_bounding_box.xmin - bb_factor, 0);
		net->current_bounding_box.xmax = std::min(net->current_bounding_box.xmax + bb_factor, nx + 1);
		net->current_bounding_box.ymin = std::max(net->current_bounding_box.ymin - bb_factor, 0);
		net->current_bounding_box.ymax = std::min(net->current_bounding_box.ymax + bb_factor, ny + 1);
	}

	/*for (int net_a = 0; net_a < num_nets; ++net_a) {*/
		/*for (int net_b = net_a + 1; net_b < num_nets; ++net_b) {*/
			/*struct s_bb *net_a_bb = &nets[net_a].current_bounding_box;*/
			/*struct s_bb *net_b_bb = &nets[net_b].current_bounding_box;*/

			/*t_sink_direction net_a_sink_direction(nets[net_a]);*/
			/*t_sink_direction net_b_sink_direction(nets[net_b]);*/

			/*zlog_debug(init_net_log, "Net %d src %d,%d sink %d,%d lrbt %d %d %d %d (hi %d vi %d) ", net_a,*/
					/*nets[net_a].current_source.x,*/
					/*nets[net_a].current_source.y,*/
					/*nets[net_a].current_sink->x,*/
					/*nets[net_a].current_sink->y,*/
					/*net_a_bb->xmin,*/
					/*net_a_bb->xmax,*/
					/*net_a_bb->ymin,*/
					/*net_a_bb->ymax,*/
					/*net_a_sink_direction.horizontal_increasing, net_a_sink_direction.vertical_increasing);*/

			/*zlog_debug(init_net_log, "Net %d src %d,%d sink %d,%d lrbt %d %d %d %d (hi %d vi %d) ", net_b,*/
					/*nets[net_b].current_source.x,*/
					/*nets[net_b].current_source.y,*/
					/*nets[net_b].current_sink->x,*/
					/*nets[net_b].current_sink->y,*/
					/*net_b_bb->xmin,*/
					/*net_b_bb->xmax,*/
					/*net_b_bb->ymin,*/
					/*net_b_bb->ymax, */
					/*net_b_sink_direction.horizontal_increasing, net_b_sink_direction.vertical_increasing);*/

			/*int net_a_bb_area = get_bounding_box_area(&nets[net_a].current_bounding_box);*/
			/*int net_b_bb_area = get_bounding_box_area(&nets[net_b].current_bounding_box);*/

			/*[> relax <]*/
			/*bool same_sink_direction = net_a_sink_direction.horizontal_increasing == net_b_sink_direction.horizontal_increasing*/
					/*&& net_a_sink_direction.vertical_increasing == net_b_sink_direction.vertical_increasing;*/
			/*bool overlap = bounding_box_overlap(nets[net_a].current_bounding_box, nets[net_b].current_bounding_box);*/
			/*bool overlap_area_more_than_threshold = (float)get_bounding_box_overlap_area(nets[net_a].current_bounding_box, nets[net_b].current_bounding_box)*/
				/*> 0.5f * (net_a_bb_area + net_b_bb_area) / 2;*/


			/*if (overlap && same_sink_direction && overlap_area_more_than_threshold) {*/
				/*nets[net_a].overlapping_nets.push_back(&nets[net_b]);*/
				/*if (overlap) {*/
					/*zlog_debug(init_net_log, "overlap ");*/
				/*}*/
				/*if (same_sink_direction) {*/
					/*zlog_debug(init_net_log, "same_sink_direction ");*/
				/*}*/
			/*}*/
			/*zlog_debug(init_net_log, "\n");*/

			/*[>if (!bounding_box_overlap(nets[net_a].current_bounding_box, nets[net_b].current_bounding_box)<]*/
					/*[>|| nets[net_a].done_sinks_to_left_of_source_inode != nets[net_b].done_sinks_to_left_of_source_inode) {<]*/
				/*[>nets[net_a].non_overlapping_nets.push_back(&nets[net_b]);<]*/
				/*[>zlog_debug(init_net_log, "no overlap\n");<]*/
			/*[>} else {<]*/
				/*[>nets[net_a].overlapping_nets.push_back(&nets[net_b]);<]*/
				/*[>zlog_debug(init_net_log, "overlap\n");<]*/
			/*[>}<]*/
		/*}*/
	/*}*/
}

void init_nets(t_simple_net *nets)
{
	for (int inet = 0; inet < num_nets; ++inet) {
		t_simple_net *net = &nets[inet];

		net->index = inet;

		zlog_debug(init_net_log, "Net %d\n", inet);

		int source_x = block[clb_net[inet].node_block[0]].x;
		int source_y = block[clb_net[inet].node_block[0]].y
			+ block[clb_net[inet].node_block[0]].type->pin_height[clb_net[inet].node_block_pin[0]];

		net->source.original_clb_net_index = 0;
		net->source.inode = net_rr_terminals[inet][0];
		net->source.x = source_x;
		net->source.y = source_y;

		zlog_debug(init_net_log, "Original source (%d,%d)\n", source_x, source_y);

		net->current_source = net->source;

		/*for (int isink = 1; isink <= clb_net[inet].num_sinks; ++isink) {*/
			/*int sink_x = block[clb_net[inet].node_block[isink]].x;*/
			/*int sink_y = block[clb_net[inet].node_block[isink]].y*/
				/*+ block[clb_net[inet].node_block[isink]].type->pin_height[clb_net[inet].node_block_pin[isink]];*/

			/*if (sink_x > source_x) {*/
				/*net->sorted_sinks_to_right_of_source_inode.emplace_back(*/
						/*isink, net_rr_terminals[inet][isink], sink_x, sink_y);*/
			/*} else {*/
				/*net->sorted_sinks_to_left_of_source_inode.emplace_back(*/
						/*isink, net_rr_terminals[inet][isink], sink_x, sink_y);*/
			/*}*/

			/*zlog_debug(init_net_log, "Original sink (%d,%d)\n", sink_x, sink_y);*/
		/*}*/

		for (int isink = 1; isink <= clb_net[inet].num_sinks; ++isink) {
			int sink_x = block[clb_net[inet].node_block[isink]].x;
			int sink_y = block[clb_net[inet].node_block[isink]].y
				+ block[clb_net[inet].node_block[isink]].type->pin_height[clb_net[inet].node_block_pin[isink]];

			net->sorted_sinks.emplace_back(
					isink, net_rr_terminals[inet][isink], sink_x, sink_y);

			zlog_debug(init_net_log, "Original sink (%d,%d)\n", sink_x, sink_y);
		}

		/*std::sort(net->sorted_sinks_to_right_of_source_inode.begin(),*/
				/*net->sorted_sinks_to_right_of_source_inode.end(),*/
				/*compare_to_source(net->source)); */
		/*std::sort(net->sorted_sinks_to_left_of_source_inode.begin(),*/
				/*net->sorted_sinks_to_left_of_source_inode.end(),*/
				/*compare_to_source(net->source)); */

		std::sort(net->sorted_sinks.begin(),
				net->sorted_sinks.end(),
				compare_to_source(net->source)); 


		zlog_debug(init_net_log, "Sorted sinks\n");
		for (const auto &s: net->sorted_sinks) {
			zlog_debug(init_net_log, "(%d,%d) dist_to_source: %d\n", s.x, s.y, abs(net->current_source.x - s.x) + abs(net->current_source.y - s.y));
		}

		/*for (const auto &s: net->sorted_sinks_to_left_of_source_inode) {*/
			/*zlog_debug(init_net_log, "Sorted sinks to left of source (%d,%d)\n", s.x, s.y);*/
		/*}*/
		/*for (const auto &s: net->sorted_sinks_to_right_of_source_inode) {*/
			/*zlog_debug(init_net_log, "Sorted sinks to right of source (%d,%d)\n", s.x, s.y);*/
		/*}*/

		net->current_sink_index = 0;
		net->current_sink = &net->sorted_sinks[net->current_sink_index];

		/*if (net->sorted_sinks_to_left_of_source_inode.size() == 0) {*/
			/*net->done_sinks_to_left_of_source_inode = true;*/
			/*net->current_sink = &net->sorted_sinks_to_right_of_source_inode[net->current_sink_index];*/
		/*} else {*/
			/*net->current_sink = &net->sorted_sinks_to_left_of_source_inode[net->current_sink_index];*/
		/*}*/

		const int bb_factor = 1;

		net->current_bounding_box.xmin = std::min(net->current_source.x, net->current_sink->x);
		net->current_bounding_box.xmax = std::max(net->current_source.x, net->current_sink->x);
		net->current_bounding_box.ymin = std::min(net->current_source.y, net->current_sink->y);
		net->current_bounding_box.ymax = std::max(net->current_source.y, net->current_sink->y);
		--net->current_bounding_box.xmin;
		--net->current_bounding_box.ymin;

		net->current_bounding_box.xmin = std::max(net->current_bounding_box.xmin - bb_factor, 0);
		net->current_bounding_box.xmax = std::min(net->current_bounding_box.xmax + bb_factor, nx + 1);
		net->current_bounding_box.ymin = std::max(net->current_bounding_box.ymin - bb_factor, 0);
		net->current_bounding_box.ymax = std::min(net->current_bounding_box.ymax + bb_factor, ny + 1);
	}

	for (int net_a = 0; net_a < num_nets; ++net_a) {
		for (int net_b = net_a + 1; net_b < num_nets; ++net_b) {
			struct s_bb *net_a_bb = &nets[net_a].current_bounding_box;
			struct s_bb *net_b_bb = &nets[net_b].current_bounding_box;

			t_sink_direction net_a_sink_direction(nets[net_a]);
			t_sink_direction net_b_sink_direction(nets[net_b]);

			zlog_debug(init_net_log, "Net %d src %d,%d sink %d,%d lrbt %d %d %d %d (hi %d vi %d) ", net_a,
					nets[net_a].current_source.x,
					nets[net_a].current_source.y,
					nets[net_a].current_sink->x,
					nets[net_a].current_sink->y,
					net_a_bb->xmin,
					net_a_bb->xmax,
					net_a_bb->ymin,
					net_a_bb->ymax,
					net_a_sink_direction.horizontal_increasing, net_a_sink_direction.vertical_increasing);

			zlog_debug(init_net_log, "Net %d src %d,%d sink %d,%d lrbt %d %d %d %d (hi %d vi %d) ", net_b,
					nets[net_b].current_source.x,
					nets[net_b].current_source.y,
					nets[net_b].current_sink->x,
					nets[net_b].current_sink->y,
					net_b_bb->xmin,
					net_b_bb->xmax,
					net_b_bb->ymin,
					net_b_bb->ymax, 
					net_b_sink_direction.horizontal_increasing, net_b_sink_direction.vertical_increasing);

			int net_a_bb_area = get_bounding_box_area(&nets[net_a].current_bounding_box);
			int net_b_bb_area = get_bounding_box_area(&nets[net_b].current_bounding_box);

			/* relax */
			bool same_sink_direction = net_a_sink_direction.horizontal_increasing == net_b_sink_direction.horizontal_increasing
					&& net_a_sink_direction.vertical_increasing == net_b_sink_direction.vertical_increasing;
			bool overlap = bounding_box_overlap(nets[net_a].current_bounding_box, nets[net_b].current_bounding_box);
			bool overlap_area_more_than_threshold = (float)get_bounding_box_overlap_area(nets[net_a].current_bounding_box, nets[net_b].current_bounding_box)
				> 0.5f * (net_a_bb_area + net_b_bb_area) / 2;


			if (overlap && same_sink_direction && overlap_area_more_than_threshold) {
				nets[net_a].overlapping_nets.push_back(&nets[net_b]);
				if (overlap) {
					zlog_debug(init_net_log, "overlap ");
				}
				if (same_sink_direction) {
					zlog_debug(init_net_log, "same_sink_direction ");
				}
			}
			zlog_debug(init_net_log, "\n");

			/*if (!bounding_box_overlap(nets[net_a].current_bounding_box, nets[net_b].current_bounding_box)*/
					/*|| nets[net_a].done_sinks_to_left_of_source_inode != nets[net_b].done_sinks_to_left_of_source_inode) {*/
				/*nets[net_a].non_overlapping_nets.push_back(&nets[net_b]);*/
				/*zlog_debug(init_net_log, "no overlap\n");*/
			/*} else {*/
				/*nets[net_a].overlapping_nets.push_back(&nets[net_b]);*/
				/*zlog_debug(init_net_log, "overlap\n");*/
			/*}*/
		}
	}
}

static void alloc_per_thread_timing_driven_route_structs(float **pin_criticality_ptr,
		int **sink_order_ptr, t_rt_node *** rt_node_of_sink_ptr) {

	/* Allocates all the structures needed only by the timing-driven router.   */

	int max_pins_per_net;
	float *pin_criticality;
	int *sink_order;
	t_rt_node **rt_node_of_sink;

	max_pins_per_net = get_max_pins_per_net();

	pin_criticality = (float *) my_malloc(
			(max_pins_per_net - 1) * sizeof(float));
	*pin_criticality_ptr = pin_criticality - 1; /* First sink is pin #1. */

	sink_order = (int *) my_malloc((max_pins_per_net - 1) * sizeof(int));
	*sink_order_ptr = sink_order - 1;

	rt_node_of_sink = (t_rt_node **) my_malloc(
			(max_pins_per_net - 1) * sizeof(t_rt_node *));
	*rt_node_of_sink_ptr = rt_node_of_sink - 1;
}

static void *worker_thread(void *arg)
{
	/*thread_info *info = (thread_info *)arg;*/

	/*float *pin_criticality;*/
	/*int *sink_order;*/
	/*t_rt_node ** rt_node_of_sink;*/

	/*alloc_per_thread_timing_driven_route_structs(&pin_criticality, &sink_order, &rt_node_of_sink);*/

	/*t_rr_node_route_inf *l_rr_node_route_inf = (t_rr_node_route_inf *)malloc(num_rr_nodes * sizeof(t_rr_node_route_inf));*/

	/*for (int inode = 0; inode < num_rr_nodes; inode++) {*/
		/*l_rr_node_route_inf[inode].prev_node = NO_PREVIOUS;*/
		/*l_rr_node_route_inf[inode].prev_edge = NO_PREVIOUS;*/
		/*l_rr_node_route_inf[inode].path_cost = HUGE_POSITIVE_FLOAT;*/
		/*l_rr_node_route_inf[inode].backward_path_cost = HUGE_POSITIVE_FLOAT;*/
		/*l_rr_node_route_inf[inode].target_flag = 0;*/
	/*}*/

	/*t_rt_node **l_rr_node_to_rt_node = (t_rt_node **) malloc(*/
			/*num_rr_nodes * sizeof(t_rt_node *));*/

	/*boolean is_routable;*/
	/*int num_heap_pushes;*/


/*#ifdef PARALLEL_ROUTE*/
	/*while (true) {*/
		/*char error_str[256];*/
		/*zlog_debug(route_outer_log, "[%d waiting for signal]\n", info->thread_index);*/
		/*[> wait for main thread to start this worker thread <]*/
		/*while (sem_wait(info->start_route_sem) == -1 && errno == EINTR) {*/
			/*printf("Thread %d wait for start route interrupted\n", info->thread_index);*/
		/*}*/
		/*[>if (sem_wait(info->start_route_sem)) {	<]*/
			/*[>strerror_r(errno, error_str, sizeof(error_str));<]*/
			/*[>printf("Error waiting for main thread's signal: %s\n", error_str);<]*/
			/*[>exit(-1);<]*/
		/*[>}<]*/

		/*if (info->successful) {*/
			/*break;			*/
		/*}*/
/*#endif*/

		/*zlog_debug(route_outer_log, "[%d start route]\n", info->thread_index);*/

		/*double total_stall_time = 0;*/
		/*bool final;*/
		/*bool final_no_route = false;*/
		/*bool routed_last = false;*/
		/*int real_net_index;*/
		/*t_simple_net *current_net;*/
		/*while (current_net) {*/
			/*zlog_debug(route_outer_log, "%d Routing net %d real %d final %d final_no_route: %d\n");*/


			/*thread_safe_pathfinder_update_one_cost(inet, trace_head[inet], -1, *info->pres_fac);*/
			/*thread_safe_free_traceback(inet);*/

			/*pthread_barrier_wait(info->barrier, info->thread_index);*/

			/*num_heap_pushes = 0;*/
			/*is_routable = simple_route_net(current_net,*/
					/*info->router_opts, *info->pres_fac,*/
					/*l_rr_node_route_inf,*/
					/*info->net_timing*/
					/*&num_heap_pushes);*/

			/*pthread_barrier_wait(info->barrier, info->thread_index);*/

			/*thread_safe_pathfinder_update_one_cost(inet, trace_head[inet], 1, *info->pres_fac);*/

			/*[>if (!final) {<]*/
				/*[>timeval t1, t2;<]*/
				/*[>assert(!gettimeofday(&t1, NULL));<]*/

				/*[>if (pthread_barrier_wait(info->barrier, info->thread_index)) {<]*/
					/*[>assert(!gettimeofday(&t2, NULL));<]*/
					/*[>double stall_time = t2.tv_sec - t1.tv_sec;<]*/
					/*[>stall_time += (double)(t2.tv_usec - t1.tv_usec)/1e6;<]*/
					/*[>info->total_stall_time += stall_time;<]*/
				/*[>}<]*/
				/*[>[>printf("Stall time of thread %d: %g\n", info->thread_index, elapsed);<]<]*/
			/*[>} <]*/


			/*[> Impossible to route? (disconnected rr_graph) <]*/

			/*if (!is_routable) {*/
				/*zlog_fatal(route_outer_log, "Routing failed.\n");*/
				/*[>			free_timing_driven_route_structs(pin_criticality,<]*/
				/*[>					sink_order, rt_node_of_sink);<]*/
				/*[>			free(net_index);<]*/
				/*[>			free(sinks);<]*/
				/*exit(-1);*/
			/*}*/

			/*for (int i = 0; i < clb_net[inet].num_sinks; ++i) {*/
				/*info->sink_order[inet][i] = sink_order[i+1]; //sink order is offset by 1 due to allocation*/
			/*}*/

			/*zlog_debug(route_outer_log, "%d Done routing routing net %d\n", info->thread_index, inet);*/

			/*[>if (inet < 0 && final) {<]*/
				/*[>final_no_route = true;<]*/
			/*[>}<]*/
		/*}*/

		/*if (inet < 0) {*/
			/*pthread_barrier_wait(info->barrier, info->thread_index);*/
			/*pthread_barrier_wait(info->barrier, info->thread_index);*/
			/*zlog_debug(route_outer_log, "%d Final no route\n", info->thread_index);*/
		/*}*/

/*#ifdef PARALLEL_ROUTE*/
		/*zlog_debug(route_outer_log, "[%d complete route]\n", info->thread_index);*/
		/*assert(!sem_post(info->complete_route_sem));*/
	/*}*/
/*#endif*/

	/*return NULL;*/
}

static void init_thread_info(thread_info *tinfo, int num_threads,
		my_pthread_barrier_t *barrier,
		int *net_index,
		t_router_opts *router_opts,
		float *pres_fac,
		t_net_timing *net_timing,
		t_next_net *next_net)
{
	char str[256];
	sem_t *start_route_sem = sem_open("start_route_sem", O_CREAT | O_EXCL, S_IWUSR | S_IRUSR, 0);
	if (start_route_sem == SEM_FAILED) {
		if (errno == EEXIST) {
			printf("Existing start route sem sem, recreating\n");
			assert(!sem_unlink("start_route_sem"));
			start_route_sem = sem_open("start_route_sem", O_CREAT | O_EXCL, S_IWUSR | S_IRUSR, 0);
			assert(start_route_sem != SEM_FAILED);
		} else {
			strerror_r(errno, str, sizeof(str));
			printf("failed to create start route sem: %s\n", str);
			exit(-1);
		}
	}

	int **sink_order = (int **)malloc(sizeof(int *) * num_nets);
	int max_num_pins = get_max_pins_per_net();
	for (int i = 0; i < num_nets; ++i) {
		sink_order[i] = (int *)malloc(sizeof(int) * max_num_pins);
	}

	for (int i = 0; i < num_threads; ++i) {
		tinfo[i].thread_index = i;

		tinfo[i].start_route_sem = start_route_sem;
		sprintf(str, "complete_route_sem_%d", i);
		tinfo[i].complete_route_sem = sem_open(str, O_CREAT | O_EXCL, S_IWUSR | S_IRUSR, 0);
		if (tinfo[i].complete_route_sem == SEM_FAILED) {
			if (errno == EEXIST) {
				printf("Existing complete route sem, recreating\n");
				assert(!sem_unlink(str));
				tinfo[i].complete_route_sem = sem_open(str, O_CREAT | O_EXCL, S_IWUSR | S_IRUSR, 0);
				assert(tinfo[i].complete_route_sem != SEM_FAILED);
			} else {
				strerror_r(errno, str, sizeof(str));
				printf("failed to create start route sem: %s\n", str);
				exit(-1);
			}
		}
		tinfo[i].barrier = barrier;

		tinfo[i].net_index = net_index;
		tinfo[i].router_opts = router_opts;
		tinfo[i].pres_fac = pres_fac;
		tinfo[i].net_timing = net_timing;
		tinfo[i].sink_order = sink_order;

		tinfo[i].successful = false;

		tinfo[i].next_net = next_net;
	}
}

static void create_worker_threads(int num_threads, thread_info *tinfo, pthread_t **tids)
{
	*tids = (pthread_t *)malloc(num_threads * sizeof(pthread_t));
	
	for (int i = 0; i < num_threads; ++i) {
		if (pthread_create(&(*tids)[i], NULL, worker_thread, &tinfo[i])) {
			printf("Failed to create worker thread %d\n", i);
		}
	}
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

static void check_net_route(int inet) {
	struct s_trace *tptr = trace_head[inet];
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
			zlog_info(check_route_log, "Node %d repeated\n", visited_nodes[i]);
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


static void simple_expand_neighbours(const struct s_heap *current, 
		float bend_cost, float astar_fac, float criticality_fac, int target_node,
		int *num_heap_pushes, struct s_bb *l_route_bb,
		t_rr_node_route_inf *l_rr_node_route_inf,
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

	for (iconn = 0; iconn < num_edges; iconn++) {
		to_node = rr_node[inode].edges[iconn];

		sprintf_rr_node(to_node, buffer);
		zlog_debug(route_inner_log, "\tNeighbor: %s ", buffer);

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
			zlog_debug(route_inner_log, "Neighbor is non-monotonic [%d > %d] ", neighbor_to_sink_manhattan_distance, current_to_sink_manhattan_distance);
			/*continue;*/
		}

		/*if ((!hor.contains(neighbor_pos.first) || !vert.contains(neighbor_pos.second)) && monotonic) {*/
			/*DVLOG(2) << rr_node[inode] << " -> " << rr_node[to_node] << " is outside of " << rr_node[target_node];*/
			/*continue;*/
		/*}*/

		if (rr_node[to_node].xhigh < l_route_bb->xmin
				|| rr_node[to_node].xlow > l_route_bb->xmax
				|| rr_node[to_node].yhigh < l_route_bb->ymin
				|| rr_node[to_node].ylow > l_route_bb->ymax) {
			zlog_debug(route_inner_log, "Outside of bounding box\n");
			continue; /* Node is outside (expanded) bounding box. */
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
			zlog_debug(route_inner_log, "Not the target IPIN\n");
			continue;
		}

		/* new_back_pcost stores the "known" part of the cost to this node -- the   *
		 * congestion cost of all the routing resources back to the existing route  *
		 * plus the known delay of the total path back to the source.  new_tot_cost *
		 * is this "known" backward cost + an expected cost to get to the target.   */

		new_back_pcost = old_back_pcost
				+ (1. - criticality_fac) * get_rr_cong_cost(to_node, 1);

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

		new_tot_cost = new_back_pcost
				+ astar_fac
						* get_timing_driven_expected_cost(to_node, target_node,
								criticality_fac, new_R_upstream, 1);

		node_to_heap(to_node, inode, iconn, new_tot_cost, new_back_pcost,
				new_R_upstream, l_rr_node_route_inf, heap);
		++(*num_heap_pushes);

	} /* End for all neighbours */
}

const t_inode_with_pos *get_current_sink(const t_simple_net *net)
{
	const t_inode_with_pos *sink;
	/*if (!net->done_sinks_to_left_of_source_inode) {*/
		/*sink = &net->sorted_sinks_to_left_of_source_inode[net->current_sink_index];*/
	/*} else {*/
		/*sink = &net->sorted_sinks_to_right_of_source_inode[net->current_sink_index];*/
	/*}*/
	return sink;
}

boolean simple_route_net(t_simple_net *net,
		t_router_opts *opts, float pres_fac,
		t_rr_node_route_inf *l_rr_node_route_inf, 
		t_net_timing *net_timing,
		int *num_heap_pushes) {

	/*Returns TRUE as long is found some way to hook up this net, even if that **/
	/*way resulted in overuse of resources (congestion).  If there is no way   **/
	/*to route this net, even ignoring congestion, it returns FALSE.  In this  **/
	/*case the rr_graph is disconnected and you can give up. If slacks = NULL, **/
	/*give each net a dummy criticality of 0.									*/

	int ipin, num_sinks, itarget, target_pin, target_node, inode;
	float target_criticality, old_tcost, new_tcost, largest_criticality,
		  old_back_cost, new_back_cost;
	struct s_heap current;
	struct s_trace *new_route_start_tptr;
	int highfanout_rlim;
	std::priority_queue<struct s_heap> heap;

	/* Rip-up any old routing. */
	/*thread_safe_pathfinder_update_one_cost(inet, trace_head[inet], -1, pres_fac);*/
	/*thread_safe_free_traceback(inet);*/

	t_rt_node *rt_root = (t_rt_node *)malloc(sizeof(t_rt_node));
	rt_root->u.child_list = NULL;
	rt_root->parent_node = NULL;
	rt_root->parent_switch = OPEN;
	rt_root->re_expand = TRUE;
	rt_root->inode = net->current_source.inode;
	rt_root->C_downstream = rr_node[rt_root->inode].C;
	rt_root->R_upstream = rr_node[rt_root->inode].R;
	rt_root->Tdel = 0.5 * rr_node[rt_root->inode].R * rr_node[rt_root->inode].C;

	net->rr_node_to_rt_node[rt_root->inode] = rt_root;

	std::vector<int> modified_inodes;
	char buffer[256];

	const t_inode_with_pos *sink = get_current_sink(net);
	target_node = sink->inode;

	sprintf_rr_node(target_node, buffer);
	zlog_debug(route_inner_log, "Sink: %s\n", buffer);

	target_criticality = 1;

	float backward_path_cost = target_criticality * rt_root->Tdel;
	float tot_cost = backward_path_cost
		+ opts->astar_fac
		* get_timing_driven_expected_cost(rt_root->inode, target_node,
				target_criticality, rt_root->R_upstream, 1);

	node_to_heap(rt_root->inode, NO_PREVIOUS, NO_PREVIOUS, tot_cost, backward_path_cost, rt_root->R_upstream, l_rr_node_route_inf, heap);

	current = heap.top();
	heap.pop();

	inode = current.index;
	assert(inode >= 0 && inode < num_rr_nodes);

	while (inode != target_node) {
		old_tcost = l_rr_node_route_inf[inode].path_cost;
		/*old_back_cost = l_rr_node_route_inf[inode].backward_path_cost;*/

		if (old_tcost > 0.99 * HUGE_POSITIVE_FLOAT) /* First time touched. */
			old_back_cost = HUGE_POSITIVE_FLOAT;
		else
			old_back_cost = l_rr_node_route_inf[inode].backward_path_cost;

		new_tcost = current.cost;
		new_back_cost = current.backward_path_cost;

		sprintf_rr_node(inode, buffer);
		zlog_debug(route_inner_log, "Current: %s old_t: %g new_t: %g old_b: %g new_b: %g old_prev_node: %d new_prev_node: %d\n", buffer, old_tcost, new_tcost, old_back_cost, new_back_cost, l_rr_node_route_inf[inode].prev_node, current.u.prev_node);

		if (old_tcost > new_tcost && old_back_cost > new_back_cost) {
			bool in_route_tree = net->rr_node_to_rt_node[inode] != NULL;
			bool not_in_previous_route_tree = current.u.prev_node != NO_PREVIOUS;

			if (in_route_tree && not_in_previous_route_tree && current.u.prev_node != l_rr_node_route_inf[inode].prev_node) {
				zlog_warn(route_inner_log, "Trying to drive an existing route tree node %s with a new driver\n", buffer);
			} else {
				l_rr_node_route_inf[inode].prev_node = current.u.prev_node;
				l_rr_node_route_inf[inode].prev_edge = current.prev_edge;
				l_rr_node_route_inf[inode].path_cost = new_tcost;
				l_rr_node_route_inf[inode].backward_path_cost = new_back_cost;

				if (old_tcost > 0.99 * HUGE_POSITIVE_FLOAT) /* First time touched. */
					modified_inodes.push_back(inode);

				struct s_bb l_route_bb;
				simple_expand_neighbours(&current, opts->bend_cost, opts->astar_fac,
						target_criticality, target_node, 
						num_heap_pushes,
						&l_route_bb,
						l_rr_node_route_inf,
						heap);
			}
		}

		if (heap.empty()) {
			zlog_fatal(route_inner_log, "Failed to find a path to %d\n", target_node);
			return FALSE;
		}

		current = heap.top();
		heap.pop();

		inode = current.index;
		assert(inode >= 0 && inode < num_rr_nodes);
	}

	new_route_start_tptr = thread_safe_update_traceback(&net->trace_head, &net->trace_tail, &current, l_rr_node_route_inf);
	/*rt_node_of_sink[target_pin] = thread_safe_update_route_tree(&current, l_rr_node_route_inf, l_rr_node_to_rt_node);*/
	t_rt_node *sink_rt_node = NULL;//thread_safe_update_route_tree(&current, l_rr_node_route_inf, net->rr_node_to_rt_node);
	assert(sink_rt_node == net->rr_node_to_rt_node[current.index]);
	/*		free_heap_data(current);*/
	/*thread_safe_pathfinder_update_one_cost(inet, new_route_start_tptr, 1, pres_fac);*/

	for (const auto &modified : modified_inodes) {
		l_rr_node_route_inf[modified].path_cost = HUGE_POSITIVE_FLOAT;
		/*l_rr_node_route_inf[modified].prev_node = NO_PREVIOUS;*/
	}

	/* For later timing analysis. */
	net_timing->delay[sink->original_clb_net_index] = net->rr_node_to_rt_node[target_node]->Tdel;
	/*check_net_route(inet);*/
	/*	free_route_tree(rt_root);*/
	return (TRUE);
}

static boolean parallel_timing_driven_route_net(int inet, t_router_opts *opts,
		float pres_fac,
		float *pin_criticality, int *sink_order,
		t_rt_node **l_rr_node_to_rt_node,  
		t_rr_node_route_inf *l_rr_node_route_inf, t_net_timing *net_timing,
		int *num_heap_pushes) {

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
	std::priority_queue<struct s_heap> heap;

	/* Rip-up any old routing. */
	/*thread_safe_pathfinder_update_one_cost(inet, trace_head[inet], -1, pres_fac);*/
	/*thread_safe_free_traceback(inet);*/

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
			pin_criticality[ipin] =		 ROUTE_PATH_WEIGHT	* net_timing[inet].path_criticality[ipin]
								  + (1 - ROUTE_PATH_WEIGHT) * net_timing[inet].timing_criticality[ipin]; 
#else
			/* Pin criticality is based on only timing criticality. */
			pin_criticality[ipin] = net_timing[inet].timing_criticality[ipin];
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
	char buffer[256];

	for (itarget = 1; itarget <= num_sinks; itarget++) {
		target_pin = sink_order[itarget];
		target_node = net_rr_terminals[inet][target_pin];

		sprintf_rr_node(target_node, buffer);
		zlog_debug(route_inner_log, "Sink: %s\n", buffer);

		target_criticality = pin_criticality[target_pin];

		highfanout_rlim = mark_node_expansion_by_bin(inet, target_node,
				rt_root);

		zlog_debug(route_inner_log, "Adding route tree to heap\n");
		add_route_tree_to_heap(rt_root, target_node, target_criticality,
				opts->astar_fac, l_rr_node_route_inf, heap, inet);

		if (heap.empty()) { /* Infeasible routing.  No possible path for net. */
			vpr_printf(TIO_MESSAGE_INFO, "Cannot route net #%d (%s) to sink #%d -- no possible path.\n",
					   inet, clb_net[inet].name, itarget);
			for (const auto &modified : modified_inodes) {
				l_rr_node_route_inf[modified].path_cost = HUGE_POSITIVE_FLOAT;
				/*l_rr_node_route_inf[modified].prev_node = NO_PREVIOUS;*/
			}
/*			free_route_tree(rt_root);*/
			return (FALSE);
		}

		current = heap.top();
		heap.pop();

		inode = current.index;
		assert(inode >= 0 && inode < num_rr_nodes);

		while (inode != target_node) {
			old_tcost = l_rr_node_route_inf[inode].path_cost;
			/*old_back_cost = l_rr_node_route_inf[inode].backward_path_cost;*/

			if (old_tcost > 0.99 * HUGE_POSITIVE_FLOAT) /* First time touched. */
				old_back_cost = HUGE_POSITIVE_FLOAT;
			else
				old_back_cost = l_rr_node_route_inf[inode].backward_path_cost;

			new_tcost = current.cost;
			new_back_cost = current.backward_path_cost;

			sprintf_rr_node(inode, buffer);
			zlog_debug(route_inner_log, "Current: %s old_t: %g new_t: %g old_b: %g new_b: %g old_prev_node: %d new_prev_node: %d\n", buffer, old_tcost, new_tcost, old_back_cost, new_back_cost, l_rr_node_route_inf[inode].prev_node, current.u.prev_node);

			/* I only re-expand a node if both the "known" backward cost is lower  *
			 * in the new expansion (this is necessary to prevent loops from       *
			 * forming in the routing and causing havoc) *and* the expected total  *
			 * cost to the sink is lower than the old value.  Different R_upstream *
			 * values could make a path with lower back_path_cost less desirable   *
			 * than one with higher cost.  Test whether or not I should disallow   *
			 * re-expansion based on a higher total cost.                          */


			/*assert(current.u.prev_node != NO_PREVIOUS || (current.u.prev_node == NO_PREVIOUS && l_rr_node_to_rt_node[inode] != NULL));*/

			if (old_tcost > new_tcost && old_back_cost > new_back_cost) {
				bool in_route_tree = l_rr_node_to_rt_node[inode] != NULL;
				bool not_in_previous_route_tree = current.u.prev_node != NO_PREVIOUS;

				if (in_route_tree && not_in_previous_route_tree && current.u.prev_node != l_rr_node_route_inf[inode].prev_node) {
					zlog_warn(route_inner_log, "Trying to drive an existing route tree node %s with a new driver\n", buffer);
				} else {
					l_rr_node_route_inf[inode].prev_node = current.u.prev_node;
					l_rr_node_route_inf[inode].prev_edge = current.prev_edge;
					l_rr_node_route_inf[inode].path_cost = new_tcost;
					l_rr_node_route_inf[inode].backward_path_cost = new_back_cost;

					if (old_tcost > 0.99 * HUGE_POSITIVE_FLOAT) /* First time touched. */
						modified_inodes.push_back(inode);

					timing_driven_expand_neighbours(&current, inet, opts->bend_cost, opts->astar_fac,
							target_criticality, target_node, 
							highfanout_rlim, false,
							num_heap_pushes,
							l_rr_node_route_inf,
							heap);
				}
			}

			if (heap.empty()) {
				zlog_fatal(route_inner_log, "Failed to find a path to %d\n", target_node);
				return FALSE;
			}

			current = heap.top();
			heap.pop();

			inode = current.index;
			assert(inode >= 0 && inode < num_rr_nodes);
		}

		/* NB:  In the code below I keep two records of the partial routing:  the   *
		 * traceback and the route_tree.  The route_tree enables fast recomputation *
		 * of the Elmore delay to each node in the partial routing.  The traceback  *
		 * lets me reuse all the routines written for breadth-first routing, which  *
		 * all take a traceback structure as input.  Before this routine exits the  *
		 * route_tree structure is destroyed; only the traceback is needed at that  *
		 * point.                                                                   */

		l_rr_node_route_inf[inode].target_flag--; /* Connected to this SINK. */

		new_route_start_tptr = thread_safe_update_traceback(&trace_head[inet], &trace_tail[inet], &current, l_rr_node_route_inf);
		/*rt_node_of_sink[target_pin] = thread_safe_update_route_tree(&current, l_rr_node_route_inf, l_rr_node_to_rt_node);*/
		t_rt_node *sink_rt_node = NULL;//thread_safe_update_route_tree(&current, l_rr_node_route_inf, l_rr_node_to_rt_node);
		assert(sink_rt_node == l_rr_node_to_rt_node[current.index]);
/*		free_heap_data(current);*/
		/*thread_safe_pathfinder_update_one_cost(inet, new_route_start_tptr, 1, pres_fac);*/

		/* clear the heap */
		heap = std::priority_queue<struct s_heap>();
		for (const auto &modified : modified_inodes) {
			l_rr_node_route_inf[modified].path_cost = HUGE_POSITIVE_FLOAT;
			/*l_rr_node_route_inf[modified].prev_node = NO_PREVIOUS;*/
		}
		modified_inodes.clear();
	}

	/* For later timing analysis. */

	update_net_delays_from_route_tree_new(net_timing, l_rr_node_to_rt_node, inet);
	check_net_route(inet);
/*	free_route_tree(rt_root);*/
	return (TRUE);
}

static void add_route_tree_to_heap(t_rt_node * rt_node, int target_node,
		float target_criticality, float astar_fac, t_rr_node_route_inf *l_rr_node_route_inf, std::priority_queue<struct s_heap> &heap, int inet) {

	/* Puts the entire partial routing below and including rt_node onto the heap *
	 * (except for those parts marked as not to be expanded) by calling itself   *
	 * recursively.                                                              */

	int inode;
	t_rt_node *child_node;
	t_linked_rt_edge *linked_rt_edge;
	float tot_cost, backward_path_cost, R_upstream;

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
	}

	linked_rt_edge = rt_node->u.child_list;

	while (linked_rt_edge != NULL) {
		child_node = linked_rt_edge->child;
		add_route_tree_to_heap(child_node, target_node, target_criticality,
				astar_fac, l_rr_node_route_inf, heap, inet);
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

static void timing_driven_expand_neighbours(const struct s_heap *current, int inet,
		float bend_cost, float astar_fac, float criticality_fac, int target_node,
		int highfanout_rlim, bool monotonic,
		int *num_heap_pushes,
		t_rr_node_route_inf *l_rr_node_route_inf,
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

	for (iconn = 0; iconn < num_edges; iconn++) {
		to_node = rr_node[inode].edges[iconn];

		sprintf_rr_node(to_node, buffer);
		zlog_debug(route_inner_log, "\tNeighbor: %s ", buffer);

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
			zlog_debug(route_inner_log, "Neighbor is non-monotonic [%d > %d] ", neighbor_to_sink_manhattan_distance, current_to_sink_manhattan_distance);
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
			zlog_debug(route_inner_log, "Outside of bounding box\n");
			continue; /* Node is outside (expanded) bounding box. */
		}

		if (clb_net[inet].num_sinks >= HIGH_FANOUT_NET_LIM) {
			if (rr_node[to_node].xhigh < target_x - highfanout_rlim
					|| rr_node[to_node].xlow > target_x + highfanout_rlim
					|| rr_node[to_node].yhigh < target_y - highfanout_rlim
					|| rr_node[to_node].ylow > target_y + highfanout_rlim) {
				zlog_debug(route_inner_log, "Outside of high fanout bin\n");
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
			zlog_debug(route_inner_log, "Not the target IPIN\n");
			continue;
		}

		/* new_back_pcost stores the "known" part of the cost to this node -- the   *
		 * congestion cost of all the routing resources back to the existing route  *
		 * plus the known delay of the total path back to the source.  new_tot_cost *
		 * is this "known" backward cost + an expected cost to get to the target.   */

		new_back_pcost = old_back_pcost
				+ (1. - criticality_fac) * get_rr_cong_cost(to_node, clb_net[inet].num_sinks);

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

		new_tot_cost = new_back_pcost
				+ astar_fac
						* get_timing_driven_expected_cost(to_node, target_node,
								criticality_fac, new_R_upstream, 1);

		node_to_heap(to_node, inode, iconn, new_tot_cost, new_back_pcost, new_R_upstream, l_rr_node_route_inf, heap);
		++(*num_heap_pushes);

	} /* End for all neighbours */
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

/* Macro used below to ensure that fractions are rounded up, but floating   *
 * point values very close to an integer are rounded to that integer.       */

#define ROUND_UP(x) (ceil (x - 0.001))

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
	while (linked_rt_edge != NULL) {
		child_node = linked_rt_edge->child;
		inode = child_node->inode;
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

#define ERROR_TOL 0.0001

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

boolean try_advanced_parallel_timing_driven_route(struct s_router_opts router_opts,
		t_net_timing *net_timing, t_ivec ** clb_opins_used_locally, boolean timing_analysis_enabled) {

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

	t_simple_net *nets = new t_simple_net[num_nets];

	init_nets(nets);
	write_hmetis_graph_file("test.graph", nets, num_nets);
	run_hmetis(8, "test.graph");
	exit(0);

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

	my_pthread_barrier_t barrier;
	assert(!my_pthread_barrier_init(&barrier, NULL, num_threads));

	thread_info *tinfo = (thread_info *)malloc(sizeof(thread_info) * num_threads);
	int *ordered_index = new int[num_nets];
	for (int i = 0; i < num_nets; ++i) {
		ordered_index[i] = i;
	}
	t_next_net next_net;
	init_thread_info(tinfo, num_threads, &barrier, net_index, &router_opts, &pres_fac, net_timing, &next_net);

#ifdef PARALLEL_ROUTE
	pthread_t *tids = NULL;
	create_worker_threads(num_threads, tinfo, &tids);
#endif
	
	int num_heap_pushes = 0;

	for (itry = 1; itry <= router_opts.max_router_iterations; itry++) {
		/*begin = clock();*/
		assert(!gettimeofday(&t1, NULL));
		vpr_printf(TIO_MESSAGE_INFO, "\n");
		vpr_printf(TIO_MESSAGE_INFO, "Routing iteration: %d\n", itry);

#ifdef PARALLEL_ROUTE
		/*dispatch_work_to_threads();*/
		/*start_worker_threads();*/
		/*wait_for_worker_threads();*/
		/* start the worker threads */
		for (i = 0; i < num_threads; ++i) {
			assert(!sem_post(tinfo[i].start_route_sem));
		}
		
/*		printf("Signaled threads\n");*/

		/* wait for the worker threads to complete */
		for (i = 0; i < num_threads; ++i) {
			/* to allow for debugger to work */
			while (sem_wait(tinfo[i].complete_route_sem) == -1 && errno == EINTR) {
				printf("Thread wait for complete route interrupted\n");
			}
		}

/*		printf("Threads completed\n");*/
#else
/*		printf("Sequential route\n");*/
		worker_thread(tinfo);
#endif


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
			vpr_printf(TIO_MESSAGE_INFO, "Wire length after first iteration %d, total available wire length %d, ratio %g\n",
					total_wirelength, available_wirelength,
					(float) (total_wirelength) / (float) (available_wirelength));
			if ((float) (total_wirelength) / (float) (available_wirelength)> FIRST_ITER_WIRELENTH_LIMIT) {
				vpr_printf(TIO_MESSAGE_INFO, "Wire length usage ratio exceeds limit of %g, fail routing.\n",
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

		char buffer[256];
		sprintf(buffer, "routes_%d.txt", itry);
		//print_route(buffer, tinfo[0].sink_order);

		/* Pathfinder guys quit after finding a feasible route. I may want to keep *
		 * going longer, trying to improve timing.  Think about this some.         */

		success = feasible_routing();
		if (success) {
			vpr_printf(TIO_MESSAGE_INFO, "Successfully routed after %d routing iterations.\n", itry);
/*			free_timing_driven_route_structs(pin_criticality, sink_order, rt_node_of_sink);*/
#ifdef DEBUG
			timing_driven_check_net_delays(net_timing);
#endif
			free(net_index);
			free(sinks);
			return (TRUE);
		}

		set<int> congested_nets;
		for (int inode = 0; inode < num_rr_nodes; ++inode) {
			if (rr_node[inode].occ > rr_node[inode].capacity) {
				zlog_debug(congestion_log, "Node %d is congested by the following nets:\n", inode);
				for (const auto &n : rr_node[inode].occupant_net_id) {
					zlog_debug(congestion_log, "%d\n", n);
					congested_nets.insert(n);
				}
			}
			assert(rr_node[inode].occ - rr_node[inode].num_reservation == rr_node[inode].occupant_net_id.size());
		}
		zlog_info(congestion_log, "There are %d congested nets:\n", congested_nets.size());
		int congested_nets_total_area = 0;
		for (const auto &c : congested_nets) {
			int area = get_bounding_box_area(&route_bb[c]);
			congested_nets_total_area += area;
			zlog_info(congestion_log, "Net %d Area: %d\n", c, area);
		}
		double congestion_nets_average_area = (double)congested_nets_total_area/congested_nets.size();
		zlog_info(congestion_log, "Congestion nets average area: %g/%d [%g]\n", congestion_nets_average_area, nx*ny, congestion_nets_average_area*100/(nx*ny));

		for (int i = 0; i < num_threads; ++i) {
			tinfo[i].successful = false;
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
			vpr_printf(TIO_MESSAGE_INFO, "Critical path: %g ns\n", critical_path_delay);
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
		assert(!gettimeofday(&t2, NULL));
		double elapsed = t2.tv_sec - t1.tv_sec;
		elapsed += (double)(t2.tv_usec - t1.tv_usec)/1e6;
		vpr_printf(TIO_MESSAGE_INFO, "Routing iteration took %g seconds.\n", elapsed);
		
		fflush(stdout);
	}

	vpr_printf(TIO_MESSAGE_INFO, "Routing failed.\n");
/*	free_timing_driven_route_structs(pin_criticality, sink_order,*/
/*			rt_node_of_sink);*/
	free(net_index);
	free(sinks);
	return (FALSE);
}
