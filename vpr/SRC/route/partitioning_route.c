#include <assert.h>
#include <fcntl.h>
#include <unistd.h>
#include <semaphore.h>
#include <sstream>
#include <random>
#include <memory>
#include <ctime>
#include <chrono>
#include <mutex>
/*#include <boost/numeric/interval.hpp>*/
/*#define TBB_USE_DEBUG 1*/
#include <tbb/task_scheduler_init.h>
#include <tbb/concurrent_queue.h>
#include <tbb/tbb.h>
#ifdef __linux__
#include <sys/syscall.h>
#include <ittnotify.h>
#endif
#include <boost/timer/timer.hpp>

#include "vpr_types.h"
#include "path_delay.h"

#include "log.h"
#include "barrier.h"
#include "graph.h"
#include "route.h"
#include "route_tree.h"
#include "trace.h"
#include "scheduler.h"
#include "geometry.h"
#include "quadtree.h"
#include "utility.h"
#include "args.h"
#include "init.h"

#ifdef __linux__
__itt_domain* pD = __itt_domain_create("Domain");
__itt_domain* dispatch_domain = __itt_domain_create("Dispatch");
__itt_domain* update_domain = __itt_domain_create("Update");
__itt_string_handle *shMyTask = __itt_string_handle_create("My Task");
__itt_string_handle *shMainTask = __itt_string_handle_create("Main Task");
#endif

/* TODO: check whether nets are global before routing */
using namespace boost::timer;

/*#define zlog_level(cat, level, ...)*/
/*#define zlog_debug(cat, level, ...)*/

bool operator<(const sink_t &a, const sink_t &b)
{
	return a.criticality_fac > b.criticality_fac;
}

void load_overlapping_nets_rtree(vector<net_t *> &nets);

/*typedef struct heap_item_t {*/
	/*int node;*/
	/*float cost;*/
	/*float known_cost;*/
	/*float upstream_R;*/
/*} heap_item_t;*/

zlog_category_t *sort_log;
zlog_category_t *delta_log;
zlog_category_t *rr_log;
zlog_category_t *net_log;
zlog_category_t *schedule_log;
zlog_category_t *scheduler_log;
zlog_category_t *independent_log;
zlog_category_t *static_log;
zlog_category_t *dynamic_log;

/*int get_num_nets()*/
/*{*/
	/*return num_nets;*/
/*}*/

/*int get_net_rr_terminals(int inet, int term)*/
/*{*/
	/*assert(inet < num_nets && term <= clb_net[inet].num_sinks);*/
	/*return net_rr_terminals[inet][term];*/
/*}*/

/*struct s_net *get_net(int inet)*/
/*{*/
	/*return &clb_net[inet];*/
/*}*/

bool operator<(const route_state_t &a, const route_state_t &b)
{
	return a.cost > b.cost;
}

float get_delay(const RREdge &e, const RRNode &v, float unbuffered_upstream_R)
{
	float upstream_R = e.properties.R;
	if (!e.properties.buffered) {
		upstream_R += unbuffered_upstream_R;
	}

	float delay = e.properties.switch_delay;
	delay += v.properties.C * (upstream_R + 0.5 * v.properties.R);

	/*zlog_level(delta_log, ROUTER_V3, " [edge_delay: %g edge_R: %g node_R: %g node_C: %g] ", e.properties.switch_delay, e.properties.R, v.properties.R, v.properties.C);*/

	return delay;
}

float get_congestion_cost(const RRNode &v)
{
	extern t_rr_indexed_data *rr_indexed_data;
	/*zlog_level(delta_log, ROUTER_V3, " [pres: %g acc: %g] ", v.properties.pres_cost, v.properties.acc_cost);*/
	return rr_indexed_data[v.properties.cost_index].base_cost * v.properties.acc_cost * v.properties.pres_cost;
}

float get_known_cost(const RRGraph &g, const RREdge &e, float criticality_fac, float unbuffered_upstream_R)
{
	const auto &target = get_target(g, e);

	float delay = get_delay(e, target, unbuffered_upstream_R);
	float congestion = get_congestion_cost(target);

	zlog_level(delta_log, ROUTER_V3, " [delay: %g congestion %g crit_fac: %g] ", delay, congestion, criticality_fac);

	return criticality_fac * delay + (1 - criticality_fac) * congestion;
}

/* Macro used below to ensure that fractions are rounded up, but floating   *
 * point values very close to an integer are rounded to that integer.       */

#define ROUND_UP(x) (ceil (x - 0.001))

static int get_expected_segs_to_target(const RRNode &current, const RRNode &target,
		int *num_segs_ortho_dir_ptr) {

	/* Returns the number of segments the same type as inode that will be needed *
	 * to reach target_node (not including inode) in each direction (the same    *
	 * direction (horizontal or vertical) as inode and the orthogonal direction).*/

	t_rr_type rr_type;
	int target_x, target_y, num_segs_same_dir, cost_index, ortho_cost_index;
	int no_need_to_pass_by_clb;
	float inv_length, ortho_inv_length, ylow, yhigh, xlow, xhigh;

	extern t_rr_indexed_data *rr_indexed_data;

	target_x = target.properties.xlow;
	target_y = target.properties.ylow;
	cost_index = current.properties.cost_index;
	inv_length = rr_indexed_data[cost_index].inv_length;
	ortho_cost_index = rr_indexed_data[cost_index].ortho_cost_index;
	ortho_inv_length = rr_indexed_data[ortho_cost_index].inv_length;
	rr_type = current.properties.type;

	if (rr_type == CHANX) {
		ylow = current.properties.ylow;
		xhigh = current.properties.xhigh;
		xlow = current.properties.xlow;

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
		ylow = current.properties.ylow;
		yhigh = current.properties.yhigh;
		xlow = current.properties.xlow;

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

float get_timing_driven_expected_cost(const RRNode &current, const RRNode &target,
		float criticality_fac, float R_upstream) {

	/* Determines the expected cost (due to both delay and resouce cost) to reach *
	 * the target node from inode.  It doesn't include the cost of inode --       *
	 * that's already in the "known" path_cost.                                   */

	t_rr_type rr_type;
	int cost_index, ortho_cost_index, num_segs_same_dir, num_segs_ortho_dir;
	float expected_cost, cong_cost, Tdel;

	extern t_rr_indexed_data *rr_indexed_data;

	rr_type = current.properties.type;

	if (rr_type == CHANX || rr_type == CHANY) {
		num_segs_same_dir = get_expected_segs_to_target(current, target,
				&num_segs_ortho_dir);
		cost_index = current.properties.cost_index;
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

template<typename ShouldExpandFunc>
void expand_neighbors_fast(const RRGraph &g, const route_state_t &current, const RRNode &target, float criticality_fac, float astar_fac, std::priority_queue<route_state_t> &heap, const ShouldExpandFunc &should_expand, perf_t *perf)
{
	for_all_out_edges(g, get_vertex(g, current.rr_node), [&heap, &g, &current, &target, &criticality_fac, &astar_fac, &should_expand, &perf] (const RREdge &e) -> void {
			auto &neighbor = get_target(g, e);

			char buffer[256];
			sprintf_rr_node(id(neighbor), buffer);
			zlog_level(delta_log, ROUTER_V3, "\tNeighbor: %s ", buffer);
			
			if (!should_expand(neighbor)) {
				return;
			}

			route_state_t item;

			item.rr_node = id(neighbor);
			item.prev_edge = &e;

			float unbuffered_upstream_R = current.upstream_R;
			float upstream_R = e.properties.R + neighbor.properties.R;
			if (!e.properties.buffered) {
				upstream_R += unbuffered_upstream_R;
			}
			item.upstream_R = upstream_R;
			
			item.delay = current.delay + get_delay(e, neighbor, unbuffered_upstream_R);

			float known_cost = current.known_cost + get_known_cost(g, e, criticality_fac, unbuffered_upstream_R);
			item.known_cost = known_cost;

			float expected_cost = get_timing_driven_expected_cost(get_vertex(g, item.rr_node), target, criticality_fac, upstream_R);
			item.cost = known_cost + astar_fac * expected_cost;

			heap.push(item);

			if (perf) {
				++perf->num_heap_pushes;
			}
			
			zlog_level(delta_log, ROUTER_V3, "added with prev_edge: %X upstream_R: %g delay: %g known_cost: %g expected_cost: %g cost: %g\n", item.prev_edge, item.upstream_R, item.delay, item.known_cost, expected_cost, item.cost);
	});
}

template<typename ShouldExpandFunc>
void expand_neighbors(const RRGraph &g, const route_state_t &current, const RRNode &target, float criticality_fac, float astar_fac, std::priority_queue<route_state_t> &heap, const ShouldExpandFunc &should_expand, perf_t *perf, bool lock)
{
	for_all_out_edges(g, get_vertex(g, current.rr_node), [&heap, &g, &current, &target, &criticality_fac, &astar_fac, &should_expand, &perf, &lock] (const RREdge &e) -> void {
			auto &neighbor = get_target(g, e);

			char buffer[256];
			sprintf_rr_node(id(neighbor), buffer);
			zlog_level(delta_log, ROUTER_V3, "\tNeighbor: %s ", buffer);

			if (perf) {
				++perf->num_neighbor_visits;
			}
			
			if (!should_expand(neighbor)) {
				return;
			}

			route_state_t item;

			item.rr_node = id(neighbor);
			item.prev_edge = &e;

			float unbuffered_upstream_R = current.upstream_R;
			float upstream_R = e.properties.R + neighbor.properties.R;
			if (!e.properties.buffered) {
				upstream_R += unbuffered_upstream_R;
			}
			item.upstream_R = upstream_R;
			
			if (lock) {
				neighbor.properties.lock->lock();
			}
			float congestion = get_congestion_cost(neighbor);
			if (lock) {
				neighbor.properties.lock->unlock();
			}
			float delay = get_delay(e, neighbor, unbuffered_upstream_R);

			item.delay = current.delay + delay;

			float known_cost = criticality_fac * delay + (1 - criticality_fac) * congestion;
			item.known_cost = current.known_cost + known_cost;

			float expected_cost = get_timing_driven_expected_cost(get_vertex(g, item.rr_node), target, criticality_fac, upstream_R);
			item.cost = item.known_cost + astar_fac * expected_cost;

			heap.push(item);

			if (perf) {
				++perf->num_heap_pushes;
			}

			zlog_level(delta_log, ROUTER_V3, " [cost: %g known_cost: %g][occ/cap: %d/%d pres: %g acc: %g][edge_delay: %g edge_R: %g node_R: %g node_C: %g] \n",
					item.cost, item.known_cost, 
					neighbor.properties.occ, neighbor.properties.capacity, neighbor.properties.pres_cost, neighbor.properties.acc_cost,
					e.properties.switch_delay, e.properties.R, neighbor.properties.R, neighbor.properties.C);
	});
}

void node_to_heap(int inode, float cost, int prev_node, int prev_edge,
		float backward_path_cost, float R_upstream);

struct s_heap *get_heap_head();

void test_heap()
{
	std::priority_queue<route_state_t> heap;
	std::priority_queue<std::unique_ptr<route_state_t>> ptr_heap;
	vector<int> random_costs(1e6);

	std::mt19937 mt(time(NULL));
	std::uniform_real_distribution<float> uni(0, 10);
	for (int i = 0; i < 1e6; ++i) {
		random_costs[i] = uni(mt);
		/*printf("%d\n", i);*/
	}

	cpu_timer timer;
	timer.start();
	for (const auto &cost : random_costs) {
		route_state_t item;
		item.cost = cost;
		heap.push(item);
	}
	timer.stop();
	cpu_times elapsed = timer.elapsed();
	printf("Pushing 1 million non-ptr items took %g\n", elapsed.wall / 1e9);

	timer.start();
	for (const auto &cost : random_costs) {
		std::unique_ptr<route_state_t> item(new route_state_t);
		item->cost = cost;
		ptr_heap.push(std::move(item));
	}
	timer.stop();
	elapsed = timer.elapsed();
	printf("Pushing 1 million ptr items took %g\n", elapsed.wall / 1e9);

	timer.start();
	for (int i = 0; i < random_costs.size(); ++i) {
		heap.pop();
	}
	timer.stop();
	elapsed = timer.elapsed();
	printf("Popping 1 million non-ptr items took %g\n", elapsed.wall / 1e9);

	timer.start();
	for (int i = 0; i < random_costs.size(); ++i) {
		ptr_heap.pop();
	}
	timer.stop();
	elapsed = timer.elapsed();
	printf("Popping 1 million ptr items took %g\n", elapsed.wall / 1e9);

	timer.start();
	for (const auto &cost : random_costs) {
		node_to_heap(0, cost, 0, 0, 0, 0);
	}
	timer.stop();
	elapsed = timer.elapsed();
	printf("Pushing 1 million vpr items took %g\n", elapsed.wall / 1e9);

	timer.start();
	for (const auto &cost : random_costs) {
		get_heap_head();
	}
	timer.stop();
	elapsed = timer.elapsed();
	printf("Popping 1 million vpr items took %g\n", elapsed.wall / 1e9);

	timer.start();
	for (int i = 0; i < 50000; ++i) {
		heap = std::priority_queue<route_state_t>();
	}
	timer.stop();
	elapsed = timer.elapsed();
	printf("Clearing heap 50000 times took %g\n", elapsed.wall / 1e9);
}

void overused_stats(const trace_t &trace, const route_tree_t &rt, const net_t &net, const RRGraph &g)
{
	for (const auto &item : trace.segments) {
		const Segment &segment = item.second;
		int sink_rr_node = item.first;

		int num_overused_nodes = 0;
		for (const auto &rr_node : segment) {
			auto &v = get_vertex(g, rr_node);
			if (v.properties.occ > v.properties.capacity) {
				++num_overused_nodes;
			}
		}

		int isink = -1;
		int num_matches = 0;
		for (int i = 0; i < net.sinks.size(); ++i) {
			if (net.sinks[i].rr_node == sink_rr_node) {
				isink = i;
				++num_matches;
			}
		}
		assert(num_matches == 1 && isink != -1);
		auto &sink = net.sinks[isink];

		char buffer[256];
		sprintf_rr_node(sink_rr_node, buffer);

		if (num_overused_nodes > 0) {
			const bounding_box_t &bb = sink.current_bounding_box;

			zlog_level(delta_log, ROUTER_V2, "Net %d sink %d: %s (distance rank %d/%lu BB %d-%d %d-%d) segment has %d/%lu (%g) overused nodes\n", net.vpr_id, sink.id, buffer, sink.distance_to_source_rank, net.sinks.size(), bb.xmin, bb.xmax, bb.ymin, bb.ymax, num_overused_nodes, segment.size(), (float)num_overused_nodes/segment.size()*100);

			bool is_only_path = route_tree_is_only_path(rt, sink_rr_node, g);
			if (is_only_path) {
				zlog_level(delta_log, ROUTER_V2, "IS_ONLY_PATH\n");
			}
		} 		
	}
}

void adjust_bb_factor(const trace_t &trace, net_t &net, const RRGraph &g, int threshold)
{
	assert(net.sinks.size() == trace.segments.size());

	for (const auto &item : trace.segments) {
		const Segment &segment = item.second;
		int sink_rr_node = item.first;

		int num_overused_nodes = 0;
		for (const auto &rr_node : segment) {
			auto &v = get_vertex(g, rr_node);
			if (v.properties.occ > v.properties.capacity) {
				++num_overused_nodes;
			}
		}

		int isink = -1;
		int num_matches = 0;
		for (int i = 0; i < net.sinks.size(); ++i) {
			if (net.sinks[i].rr_node == sink_rr_node) {
				isink = i;
				++num_matches;
			}
		}
		assert(num_matches == 1 && isink != -1);
		auto &sink = net.sinks[isink];

		char buffer[256];
		sprintf_rr_node(sink_rr_node, buffer);

		if (num_overused_nodes > 0) {
			const bounding_box_t &bb = sink.current_bounding_box;

			++sink.congested_iterations;
			if (sink.congested_iterations >= threshold) {
				sink.congested_iterations = 0;
				zlog_level(delta_log, ROUTER_V2, "Increasing net %d %s bb_factor by 1\n", net.vpr_id, buffer);
				++sink.bb_factor;
			}
		} else {
			/*if (sink.congested_iterations > 0) {*/
				/*zlog_level(delta_log, ROUTER_V2, "Decreasing net %d %s bb_factor by 1\n", net.vpr_id, buffer);*/
				/*--sink.congested_iterations;*/
			/*}*/
		}

		assert(sink.congested_iterations >= 0);
	}
}

void update_costs(RRGraph &g, float pres_fac, float acc_fac)
{
	for_all_vertices(g, [&pres_fac, &acc_fac] (RRNode &v) -> void {
			int occ = v.properties.occ;
			int capacity = v.properties.capacity;
			if (occ > capacity) {
				v.properties.acc_cost += (occ - capacity) * acc_fac;
				v.properties.pres_cost = 1. + (occ + 1 - capacity) * pres_fac;
			} else if (occ == capacity) {
				/* If occ == capacity, we don't need to increase acc_cost, but a change    *
				 * in pres_fac could have made it necessary to recompute the cost anyway.  */
				v.properties.pres_cost = 1. + pres_fac;
			}
			});
}

void update_one_cost_internal(RRNode &rr_node, int delta, float pres_fac, bool lock)
{
	if (lock) {
		rr_node.properties.lock->lock();
	}
	
	rr_node.properties.occ += delta;

	assert(rr_node.properties.occ >= 0);

	if (rr_node.properties.occ < rr_node.properties.capacity) {
		rr_node.properties.pres_cost = 1;
	} else {
		rr_node.properties.pres_cost = 1 + (rr_node.properties.occ + 1 - rr_node.properties.capacity) * pres_fac;
	}

	if (lock) {
		rr_node.properties.lock->unlock();
	}
		
	char buffer[256];
	sprintf_rr_node(id(rr_node), buffer);
	zlog_level(delta_log, ROUTER_V2, "Update cost of %s delta: %d new_occ: %d pres_fac: %g\n", buffer, delta, rr_node.properties.occ, pres_fac);
}

void update_one_cost(RRGraph &g, const vector<int>::const_iterator &rr_nodes_begin, const vector<int>::const_iterator &rr_nodes_end, int delta, float pres_fac, bool lock)
{
	/*const RouteTreeNode *last = nullptr;*/
	for (auto iter = rr_nodes_begin; iter != rr_nodes_end; ++iter) {
		/* we don't update get_source because that's the link to existing route tree and the cost is handled by update_one_cost_internal */
		/*const RouteTreeNode &rt_node = get_target(rt.graph, get_edge(rt.graph, rt_edge_id));*/
		/*last = &get_target(rt.graph, get_edge(rt.graph, rt_edge_id));*/
		/*RRNode &rr_node = get_vertex(g, rt_node.properties.rr_node);*/
		RRNode &rr_node = get_vertex(g, *iter);
		update_one_cost_internal(rr_node, delta, pres_fac, lock);
	}
	/*RRNode &rr_node = get_vertex(g, last->properties.rr_node);*/
	/*update_one_cost_internal(rr_node, delta, pres_fac);*/
}

void update_one_cost(RRGraph &g, route_tree_t &rt, const RouteTreeNode &node, int delta, float pres_fac, bool lock)
{
	assert(node.properties.valid);

	RRNode &rr_node = get_vertex(g, node.properties.rr_node);

	update_one_cost_internal(rr_node, delta, pres_fac, lock);

	for_all_out_edges(rt.graph, node, [&g, &rt, &delta, &pres_fac, &lock] (const RouteTreeEdge &e) -> void {
			update_one_cost(g, rt, get_target(rt.graph, e), delta, pres_fac, lock);
			});
}

void update_R()
{
}

/*void route_tree_update_timing(route_tree_t &rt, const RouteTreeNode &node)*/
/*{*/
	/*for_all_out_edges(rt.graph, node, [] (const RouteTreeEdge &e) -> void {*/
			/*e.properties.rr_edge;*/
	/*});*/
/*}*/

void update_sink_criticalities(net_t &net, const t_net_timing &net_timing, const route_parameters_t &params)
{
	for (int ipin = 1; ipin <= net.sinks.size(); ipin++) { 
		float pin_criticality;
		/*if (!net_timing) {*/
			/* Use criticality of 1. This makes all nets critical.  Note: There is a big difference between setting pin criticality to 0*/
			/*compared to 1.  If pin criticality is set to 0, then the current path delay is completely ignored during routing.  By setting*/
			/*pin criticality to 1, the current path delay to the pin will always be considered and optimized for */
			/*pin_criticality = 1.0;*/
		/*} else { */
#ifdef PATH_COUNTING
			/* Pin criticality is based on a weighted sum of timing and path criticalities. */	
			pin_criticality =		 ROUTE_PATH_WEIGHT	* net_timing.path_criticality[ipin]
								  + (1 - ROUTE_PATH_WEIGHT) * net_timing.timing_criticality[ipin]; 
#else
			/* Pin criticality is based on only timing criticality. */
			pin_criticality = net_timing.timing_criticality[ipin];
#endif
			/* Currently, pin criticality is between 0 and 1. Now shift it downwards 
			by 1 - max_criticality (max_criticality is 0.99 by default, so shift down
			by 0.01) and cut off at 0.  This means that all pins with small criticalities 
			(<0.01) get criticality 0 and are ignored entirely, and everything
			else becomes a bit less critical. This effect becomes more pronounced if
			max_criticality is set lower. */
			assert(pin_criticality > -0.01 && pin_criticality < 1.01);
			pin_criticality = std::max(pin_criticality - (1.0 - params.max_criticality), 0.0);

			/* Take pin criticality to some power (1 by default). */
			pin_criticality = pow(pin_criticality, params.criticality_exp);
			
			/* Cut off pin criticality at max_criticality. */
			pin_criticality = std::min(pin_criticality, params.max_criticality);
		/*}*/

		net.sinks[ipin-1].criticality_fac = pin_criticality;
	}
}

void get_overused_nodes(const route_tree_t &rt, const RouteTreeNode &node, RRGraph &g, vector<int> &overused_rr_node)
{
	assert(node.properties.valid);

	const auto &rr_node = get_vertex(g, node.properties.rr_node);

	if (rr_node.properties.occ > rr_node.properties.capacity) {
		overused_rr_node.push_back(id(rr_node));
	}
	
	for_all_out_edges(rt.graph, node, [&rt, &g, &overused_rr_node] (const RouteTreeEdge &e) -> void {
			const auto &neighbor = get_target(rt.graph, e);
			get_overused_nodes(rt, neighbor, g, overused_rr_node);
			});
}

void check_route_tree_internal(const route_tree_t &rt, const RouteTreeNode &node, RRGraph &g, vector<int> &visited_sinks, vector<int> &visited_nodes)
{
	assert(node.properties.valid);

	auto &rr_node = get_vertex(g, node.properties.rr_node);
	if (rr_node.properties.type == SINK) {
		visited_sinks.push_back(id(rr_node));
		char buffer[256];
		sprintf_rr_node(id(rr_node), buffer);
		/*zlog_level(delta_log, ROUTER_V2, "route_tree_check: %s\n", buffer);*/
	}

	visited_nodes.push_back(node.properties.rr_node);

	if (rr_node.properties.type == SOURCE) {
		rr_node.properties.recalc_occ += num_out_edges(rt.graph, node);
	} else {
		++rr_node.properties.recalc_occ;
	}

	for (const auto &branch : route_tree_get_branches(rt, node)) {
		/*for_all_out_edges(rt.graph, node, [&rt, &g, &visited_sinks, &visited_nodes] (const RouteTreeEdge &e) -> void {*/
		const auto &child = get_target(rt.graph, branch);
		check_route_tree_internal(rt, child, g, visited_sinks, visited_nodes);
	}
}

void check_route_tree(const route_tree_t &rt, const net_t &net, RRGraph &g)
{
	const auto &rt_root = *route_tree_get_rt_node(rt, net.source.rr_node);
	vector<int> sinks;
	for (const auto &s : net.sinks) {
		sinks.push_back(s.rr_node);
	}
	vector<int> visited_sinks;

	vector<int> visited_nodes;
	check_route_tree_internal(rt, rt_root, g, visited_sinks, visited_nodes);
	vector<int> duplicated_nodes;
	sort(begin(visited_nodes), end(visited_nodes));
	int current = visited_nodes[0];
	bool added_duplicated = false;
	for (int i = 1; i < visited_nodes.size(); ++i) {
		/*zlog_info(delta_log, "visited nodes %d\n", visited_nodes[i]);*/
		if (visited_nodes[i] == current) {
			if (!added_duplicated) {
				duplicated_nodes.push_back(visited_nodes[i]);
				added_duplicated = true;
			}
		} else {
			current = visited_nodes[i];
			added_duplicated = false;
		}
	}
	char buffer[256];
	if (!duplicated_nodes.empty()) {
		zlog_error(delta_log, "Error: net %d visited_nodes has duplicates: \n", net.vpr_id);
		for (const auto &d : duplicated_nodes) {
			sprintf_rr_node(d, buffer);
			zlog_error(delta_log, "%s\n", buffer);
		}
		write_graph(rt.graph, "duplicate.dot",
				[&duplicated_nodes] (const RouteTreeNode &rt_node) -> string {
				char s_rr_node[256];
				char buffer[256];
				if (find(begin(duplicated_nodes), end(duplicated_nodes), rt_node.properties.rr_node) !=
						end(duplicated_nodes)) {
				sprintf_rr_node(rt_node.properties.rr_node, s_rr_node);
				sprintf(buffer, "label=\"%s\" style=filled fillcolor=red", s_rr_node);
				} else {
				sprintf_rr_node(rt_node.properties.rr_node, s_rr_node);
				sprintf(buffer, "label=\"%s\"", s_rr_node);
				}
				return buffer;
				},
				[] (const RouteTreeEdge &rt_edge) -> string {
				return "";
				},
				[] (const RouteTreeNode &rt_node) -> bool {
				return false;
				});
		assert(false);
	}

	sort(visited_sinks.begin(), visited_sinks.end());
	sort(sinks.begin(), sinks.end());
	
	if (visited_sinks != sinks) {
		zlog_error(delta_log, "Error: Visited %lu sinks out of %lu sinks of net %d\n", visited_sinks.size(), sinks.size(), net.vpr_id);
		vector<int> only_in_required;
		vector<int> only_in_visited;
		vector<int> sym_difference;
		set_difference(sinks.begin(), sinks.end(), visited_sinks.begin(), visited_sinks.end(), back_inserter(only_in_required));
		set_difference(visited_sinks.begin(), visited_sinks.end(), sinks.begin(), sinks.end(), back_inserter(only_in_visited));
		/*set_symmetric_difference(sinks.begin(), sinks.end(), visited_sinks.begin(), visited_sinks.end(), back_inserter(sym_difference));*/
		/*assert(difference == sym_difference);*/

		write_graph(rt.graph, "error.dot",
				[] (const RouteTreeNode &rt_node) -> string {
					char buffer[256];
					sprintf_rr_node(rt_node.properties.rr_node, buffer);
					return buffer;
				},
				[] (const RouteTreeEdge &rt_edge) -> string {
					return "";
				},
				[] (const RouteTreeNode &rt_node) -> bool {
					return !rt_node.properties.valid;
				});


		zlog_error(delta_log, "Only in required: ");
		for (const int d : only_in_required) {
			zlog_error(delta_log, "%d ", d);
		}
		zlog_error(delta_log, "\n");

		zlog_error(delta_log, "Only in visited: ");
		for (const int d : only_in_visited) {
			zlog_error(delta_log, "%d ", d);
		}
		zlog_error(delta_log, "\n");

		zlog_error(delta_log, "Visited sinks: ");
		for (const int d : visited_sinks) {
			zlog_error(delta_log, "%d ", d);
		}
		zlog_error(delta_log, "\n");

		assert(false);
	}
}

void route_net_fast(RRGraph &g, const net_t &net, const route_parameters_t &params, route_state_t *state, route_tree_t &rt, trace_t &trace, t_net_timing &net_timing, perf_t *perf)
{
	std::priority_queue<route_state_t> heap;

	vector<int> modified;

	zlog_level(delta_log, ROUTER_V1, "Routing VPR net id %d\n", net.vpr_id);

	if (route_tree_empty(rt)) {
		/* special case */
		/*update_one_cost_internal(get_vertex(g, net.current_source.rr_node), 1, params.pres_fac);*/

		/*route_tree_set_source(rt, get_vertex(g, net.source.rr_node));*/
	} else {
		RouteTreeNode rt_root = get_vertex(rt.graph, rt.root_rt_node_id);
		if (rt_root.properties.rr_node != net.source.rr_node) {
			char root[256];
			char source[256];
			sprintf_rr_node(rt_root.properties.rr_node, root);
			sprintf_rr_node(net.source.rr_node, source);
			zlog_error(delta_log, "Error: Non empty route tree node has root %s that is different from net source %s\n",
					root, source);
			assert(false);
		}
	}

	vector<sink_t> sorted_sinks = net.sinks;
	std::sort(sorted_sinks.begin(), sorted_sinks.end());

	char buffer[256];

	sprintf_rr_node(net.source.rr_node, buffer);
	zlog_level(delta_log, ROUTER_V1, "Source: %s\n", buffer);

	int isink = 0;
	for (const auto &sink : sorted_sinks) {
		sprintf_rr_node(sink.rr_node, buffer);
		zlog_level(delta_log, ROUTER_V1, "Sink %d: %s criticality: %g BB: %d-%d %d-%d\n", sink.id, buffer, sink.criticality_fac, sink.current_bounding_box.xmin, sink.current_bounding_box.xmax, sink.current_bounding_box.ymin, sink.current_bounding_box.ymax);
		route_tree_add_to_heap(rt, g, get_vertex(g, sink.rr_node), sink.criticality_fac, params.astar_fac, sink.current_bounding_box, heap, perf);

		bool found_sink = false;
		while (!heap.empty() && !found_sink) {
			auto item = heap.top();
			heap.pop();

			sprintf_rr_node(item.rr_node, buffer);
			zlog_level(delta_log, ROUTER_V2, "Current: %s prev: %d old_delay: %g new_delay: %g old_known: %g new_known: %g old_cost: %g new_cost: %g\n", buffer, item.prev_edge ? id(get_source(g, *item.prev_edge)) : -1, state[item.rr_node].delay, item.delay, state[item.rr_node].known_cost, item.known_cost, state[item.rr_node].cost, item.cost);

			if (item.rr_node == sink.rr_node) {
				state[item.rr_node] = item;
				modified.push_back(item.rr_node);
				found_sink = true;
			} else if (item.known_cost < state[item.rr_node].known_cost && item.cost < state[item.rr_node].cost) {
				/*if (route_tree_get_rt_node(rt, item.rr_node)) {*/
					/*sprintf_rr_node(item.rr_node, buffer);*/
					/*zlog_warn(delta_log, "Warning: Found a lower cost path to existing route tree node %s\n", buffer);*/
					/*assert(false);*/
				/*}*/
				const auto &sink_vertex = get_vertex(g, sink.rr_node);

				state[item.rr_node] = item;
				modified.push_back(item.rr_node);

				for (const auto &e_i : get_vertex(g, item.rr_node).edges) {
					auto &e = get_edge(g, e_i);
					auto &neighbor = get_target(g, e);

					char buffer[256];
					sprintf_rr_node(id(neighbor), buffer);
					zlog_level(delta_log, ROUTER_V3, "\tNeighbor: %s ", buffer);

					const auto &prop = neighbor.properties;

					if (prop.xhigh < sink.current_bounding_box.xmin
							|| prop.xlow > sink.current_bounding_box.xmax
							|| prop.yhigh < sink.current_bounding_box.ymin
							|| prop.ylow > sink.current_bounding_box.ymax) {
						zlog_level(delta_log, ROUTER_V3, "outside of bounding box\n");
						continue;
					}

					if (prop.type == IPIN
							&& (prop.xhigh != sink_vertex.properties.xhigh 
								|| prop.yhigh != sink_vertex.properties.yhigh)) {
						zlog_level(delta_log, ROUTER_V3, "not target IPIN\n");
						continue;
					}

					route_state_t new_item;

					new_item.rr_node = id(neighbor);
					new_item.prev_edge = &e;

					float unbuffered_upstream_R = item.upstream_R;
					float upstream_R = e.properties.R + neighbor.properties.R;
					if (!e.properties.buffered) {
						upstream_R += unbuffered_upstream_R;
					}
					new_item.upstream_R = upstream_R;

					new_item.delay = item.delay + get_delay(e, neighbor, unbuffered_upstream_R);

					float known_cost = item.known_cost + get_known_cost(g, e, sink.criticality_fac, unbuffered_upstream_R);
					new_item.known_cost = known_cost;

					float expected_cost = get_timing_driven_expected_cost(get_vertex(g, new_item.rr_node), sink_vertex, sink.criticality_fac, upstream_R);
					new_item.cost = known_cost + params.astar_fac * expected_cost;

					heap.push(new_item);

					if (perf) {
						++perf->num_heap_pushes;
					}

					zlog_level(delta_log, ROUTER_V3, "added with prev_edge: %X upstream_R: %g delay: %g known_cost: %g expected_cost: %g cost: %g\n", new_item.prev_edge, new_item.upstream_R, new_item.delay, new_item.known_cost, expected_cost, new_item.cost);
				}
			} else {
				zlog_level(delta_log, ROUTER_V3, "\tNot expanding neighbor because known_cost %g > %g or cost %g > %g\n", item.known_cost, state[item.rr_node].known_cost, item.cost, state[item.rr_node].cost);
			}
		}

		assert(found_sink);

		const auto &new_path = trace_add_path(trace, g, state, sink.rr_node, net.vpr_id);
		/* new_path.second-1 because we do not update the last node because it's the node in existing route tree (which was already updated previously */
		/*int new_path_first_rr_node = *(new_path.end()-1);*/
		/* if the first node is a connection to existing route, don't update its cost again */
		/*if (route_tree_get_rt_node(rt, new_path_first_rr_node)) {*/
			/*update_one_cost(g, new_path.begin(), new_path.end()-1, 1, params.pres_fac);*/
		/*} else {*/
			/*update_one_cost(g, new_path.begin(), new_path.end(), 1, params.pres_fac);*/
		/*}*/
		update_one_cost(g, new_path.begin(), new_path.end(), 1, params.pres_fac, false);

		/* important to update route tree only after updating cost because the previous check relies on whether the
		 * first node is already in the route tree */
		route_tree_add_path(rt, g, state, new_path, net.vpr_id);

		net_timing.delay[sink.id+1] = route_tree_get_rt_node(rt, sink.rr_node)->properties.delay;

		for (const auto &m : modified)  {
			state[m].known_cost = std::numeric_limits<float>::max();
			state[m].cost = std::numeric_limits<float>::max();
		}

		heap = std::priority_queue<route_state_t>();
	}

	/*delete [] state;*/

	/*check_route_tree(rt, net, g);*/
}

void route_net_2(RRGraph &g, int vpr_id, const source_t *source, const vector<sink_t *> &sinks, const route_parameters_t &params, route_state_t *state, route_tree_t &rt, t_net_timing &net_timing, perf_t *perf, bool lock)
{
	std::priority_queue<route_state_t> heap;

	vector<int> modified;

	zlog_level(delta_log, ROUTER_V1, "Routing net %d (%lu sinks)\n", vpr_id, sinks.size());

	vector<sink_t *> sorted_sinks = sinks;
	std::sort(sorted_sinks.begin(), sorted_sinks.end());

	char buffer[256];

	sprintf_rr_node(source->rr_node, buffer);
	zlog_level(delta_log, ROUTER_V1, "Source: %s\n", buffer);

	if (route_tree_empty(rt)) {
		/* special case */
		zlog_level(delta_log, ROUTER_V2, "Empty route tree. Setting root to %s\n", buffer);

		auto &source_rr_node = get_vertex(g, source->rr_node);
		RouteTreeNode *root_rt_node = route_tree_add_rr_node(rt, source_rr_node);
		route_tree_set_node_properties(*root_rt_node, true, nullptr, source_rr_node.properties.R, 0.5 * source_rr_node.properties.R * source_rr_node.properties.C);
		route_tree_set_root(rt, source->rr_node);

		/*update_one_cost_internal(source_rr_node, 1, params.pres_fac);*/
	} else {
		const RouteTreeNode &rt_root = get_vertex(rt.graph, rt.root_rt_node_id);
		if (source && rt_root.properties.rr_node != source->rr_node) {
			char root[256];
			char source_str[256];
			sprintf_rr_node(rt_root.properties.rr_node, root);
			sprintf_rr_node(source->rr_node, source_str);
			zlog_error(delta_log, "Error: Non empty route tree node has root %s that is different from net source %s\n",
					root, source_str);
			assert(false);
		}
	}

	int isink = 0;
	for (const auto &sink : sorted_sinks) {
		sprintf_rr_node(sink->rr_node, buffer);
		zlog_level(delta_log, ROUTER_V1, "Net %d Sink %d: %s criticality: %g BB: %d-%d %d-%d Prev BB: %d-%d %d-%d\n", vpr_id, sink->id, buffer, sink->criticality_fac, sink->current_bounding_box.xmin, sink->current_bounding_box.xmax, sink->current_bounding_box.ymin, sink->current_bounding_box.ymax, sink->previous_bounding_box.xmin, sink->previous_bounding_box.xmax, sink->previous_bounding_box.ymin, sink->previous_bounding_box.ymax);

		const auto &sink_rr_node = get_vertex(g, sink->rr_node);

		route_tree_add_to_heap(rt, g, sink_rr_node, sink->criticality_fac, params.astar_fac, sink->current_bounding_box, heap, perf);

		bool found_sink = false;
		while (!heap.empty() && !found_sink) {
			auto item = heap.top();
			heap.pop();

			if (perf) {
				++perf->num_heap_pops;
			}

			sprintf_rr_node(item.rr_node, buffer);
			const auto &v = get_vertex(g, item.rr_node);
			zlog_level(delta_log, ROUTER_V3, "Current: %s occ/cap: %d/%d prev: %d new_cost: %g new_known: %g new_delay: %g old_cost: %g old_known: %g old_delay: %g\n", buffer, v.properties.occ, v.properties.capacity, item.prev_edge ? id(get_source(g, *item.prev_edge)) : -1, item.cost, item.known_cost, item.delay, state[item.rr_node].cost, state[item.rr_node].known_cost, state[item.rr_node].delay);

			if (item.rr_node == sink->rr_node) {
				state[item.rr_node] = item;
				modified.push_back(item.rr_node);
				found_sink = true;
			} else if (item.known_cost < state[item.rr_node].known_cost && item.cost < state[item.rr_node].cost) {
				/*if (route_tree_get_rt_node(rt, item.rr_node)) {*/
					/*sprintf_rr_node(item.rr_node, buffer);*/
					/*zlog_warn(delta_log, "Warning: Found a lower cost path to existing route tree node %s\n", buffer);*/
					/*assert(false);*/
				/*}*/

				state[item.rr_node] = item;
				modified.push_back(item.rr_node);

				expand_neighbors(g, item, sink_rr_node, sink->criticality_fac, params.astar_fac, heap, [&sink, &sink_rr_node] (const RRNode &v) -> bool {

					/*if (trace_has_node(prev_trace, id(v))) {*/
						/*zlog_level(delta_log, ROUTER_V3, " existing node route tree ");*/
					/*}*/
					const auto &prop = v.properties;

					if (prop.xhigh < sink->current_bounding_box.xmin
							|| prop.xlow > sink->current_bounding_box.xmax
							|| prop.yhigh < sink->current_bounding_box.ymin
							|| prop.ylow > sink->current_bounding_box.ymax) {
					zlog_level(delta_log, ROUTER_V3, "outside of bounding box\n");
					return false;
					}

					if (prop.type == IPIN
							&& (prop.xhigh != sink_rr_node.properties.xhigh 
								|| prop.yhigh != sink_rr_node.properties.yhigh)) {
					zlog_level(delta_log, ROUTER_V3, "not target IPIN\n");
					return false;
					}

					/*if (net.sinks.size() >= HIGH_FANOUT_NET_LIM) {*/
					/*if (prop.xhigh < target_x - highfanout_rlim*/
							/*|| prop.xlow > target_x + highfanout_rlim*/
							/*|| prop.yhigh < target_y - highfanout_rlim*/
							/*|| prop.ylow > target_y + highfanout_rlim) {*/
						/*return false;*/
					/*}*/
					/*}*/
					return true;
				}, perf, lock);
			} else {
				zlog_level(delta_log, ROUTER_V3, "\tNot expanding neighbor because known_cost %g > %g or cost %g > %g\n", item.known_cost, state[item.rr_node].known_cost, item.cost, state[item.rr_node].cost);
			}
		}

		if (!found_sink) {
			zlog_error(delta_log, "Error: Failed to find sink\n");
			assert(false);
		}

		vector<int> added_rr_nodes;
		route_tree_add_path(rt, g, state, sink->rr_node, vpr_id, added_rr_nodes);

		update_one_cost(g, added_rr_nodes.begin(), added_rr_nodes.end(), 1, params.pres_fac, lock);

		net_timing.delay[sink->id+1] = route_tree_get_rt_node(rt, sink->rr_node)->properties.delay;

		for (const auto &m : modified)  {
			state[m].known_cost = std::numeric_limits<float>::max();
			state[m].cost = std::numeric_limits<float>::max();
		}

		heap = std::priority_queue<route_state_t>();
	}

	zlog_level(delta_log, ROUTER_V1, "\n");

	/*delete [] state;*/

	/*check_route_tree(rt, net, g);*/
}

void route_net(RRGraph &g, int vpr_id, const source_t *source, const vector<sink_t *> &sinks, const route_parameters_t &params, route_state_t *state, route_tree_t &rt, trace_t &trace, const trace_t &prev_trace, t_net_timing &net_timing, perf_t *perf)
{
	std::priority_queue<route_state_t> heap;

	vector<int> modified;

	zlog_level(delta_log, ROUTER_V1, "Routing VPR net id %d\n", vpr_id);

	if (route_tree_empty(rt)) {
		/* special case */
		/*update_one_cost_internal(get_vertex(g, net.current_source.rr_node), 1, params.pres_fac);*/

		/*route_tree_set_source(rt, get_vertex(g, source->rr_node));*/
		auto &source_rr_node = get_vertex(g, source->rr_node);
		RouteTreeNode *root_rt_node = route_tree_add_rr_node(rt, source_rr_node);
		route_tree_set_node_properties(*root_rt_node, true, nullptr, source_rr_node.properties.R, 0.5 * source_rr_node.properties.R * source_rr_node.properties.C);
		route_tree_set_root(rt, source->rr_node);
	} else {
		RouteTreeNode rt_root = get_vertex(rt.graph, rt.root_rt_node_id);
		if (rt_root.properties.rr_node != source->rr_node) {
			char root[256];
			char source_str[256];
			sprintf_rr_node(rt_root.properties.rr_node, root);
			sprintf_rr_node(source->rr_node, source_str);
			zlog_error(delta_log, "Error: Non empty route tree node has root %s that is different from net source %s\n",
					root, source_str);
			assert(false);
		}
	}

	vector<sink_t *> sorted_sinks = sinks;
	std::sort(sorted_sinks.begin(), sorted_sinks.end());

	char buffer[256];

	sprintf_rr_node(source->rr_node, buffer);
	zlog_level(delta_log, ROUTER_V1, "Source: %s\n", buffer);

	int isink = 0;
	for (const auto &sink : sorted_sinks) {
		sprintf_rr_node(sink->rr_node, buffer);
		zlog_level(delta_log, ROUTER_V1, "Net %d Sink %d: %s criticality: %g BB: %d-%d %d-%d Prev BB: %d-%d %d-%d\n", vpr_id, sink->id, buffer, sink->criticality_fac, sink->current_bounding_box.xmin, sink->current_bounding_box.xmax, sink->current_bounding_box.ymin, sink->current_bounding_box.ymax, sink->previous_bounding_box.xmin, sink->previous_bounding_box.xmax, sink->previous_bounding_box.ymin, sink->previous_bounding_box.ymax);
		route_tree_add_to_heap(rt, g, get_vertex(g, sink->rr_node), sink->criticality_fac, params.astar_fac, sink->current_bounding_box, heap, perf);

		bool found_sink = false;
		while (!heap.empty() && !found_sink) {
			auto item = heap.top();
			heap.pop();

			sprintf_rr_node(item.rr_node, buffer);
			const auto &v = get_vertex(g, item.rr_node);
			zlog_level(delta_log, ROUTER_V2, "Current: %s occ/cap: %d/%d prev: %d old_cost: %g new_cost: %g old_delay: %g new_delay: %g old_known: %g new_known: %g \n", buffer, v.properties.occ, v.properties.capacity, item.prev_edge ? id(get_source(g, *item.prev_edge)) : -1, state[item.rr_node].cost, item.cost, state[item.rr_node].delay, item.delay, state[item.rr_node].known_cost, item.known_cost);

			if (item.rr_node == sink->rr_node) {
				state[item.rr_node] = item;
				modified.push_back(item.rr_node);
				found_sink = true;
			} else if (item.known_cost < state[item.rr_node].known_cost && item.cost < state[item.rr_node].cost) {
				/*if (route_tree_get_rt_node(rt, item.rr_node)) {*/
					/*sprintf_rr_node(item.rr_node, buffer);*/
					/*zlog_warn(delta_log, "Warning: Found a lower cost path to existing route tree node %s\n", buffer);*/
					/*assert(false);*/
				/*}*/
				const auto &sink_vertex = get_vertex(g, sink->rr_node);

				state[item.rr_node] = item;
				modified.push_back(item.rr_node);

				expand_neighbors(g, item, sink_vertex, sink->criticality_fac, params.astar_fac, heap, [&sink, &prev_trace, &sink_vertex] (const RRNode &v) -> bool {

					if (trace_has_node(prev_trace, id(v))) {
						zlog_level(delta_log, ROUTER_V3, " existing node route tree ");
					}
					const auto &prop = v.properties;

					if (prop.xhigh < sink->current_bounding_box.xmin
							|| prop.xlow > sink->current_bounding_box.xmax
							|| prop.yhigh < sink->current_bounding_box.ymin
							|| prop.ylow > sink->current_bounding_box.ymax) {
					zlog_level(delta_log, ROUTER_V3, "outside of bounding box\n");
					return false;
					}

					if (prop.type == IPIN
							&& (prop.xhigh != sink_vertex.properties.xhigh 
								|| prop.yhigh != sink_vertex.properties.yhigh)) {
					zlog_level(delta_log, ROUTER_V3, "not target IPIN\n");
					return false;
					}

					/*if (net.sinks.size() >= HIGH_FANOUT_NET_LIM) {*/
					/*if (prop.xhigh < target_x - highfanout_rlim*/
							/*|| prop.xlow > target_x + highfanout_rlim*/
							/*|| prop.yhigh < target_y - highfanout_rlim*/
							/*|| prop.ylow > target_y + highfanout_rlim) {*/
						/*return false;*/
					/*}*/
					/*}*/
					return true;
				}, perf, false);
			} else {
				zlog_level(delta_log, ROUTER_V3, "\tNot expanding neighbor because known_cost %g > %g or cost %g > %g\n", item.known_cost, state[item.rr_node].known_cost, item.cost, state[item.rr_node].cost);
			}
		}

		if (!found_sink) {
			zlog_error(delta_log, "Failed to find sink\n");
			assert(false);
		}

		const auto &new_path = trace_add_path(trace, g, state, sink->rr_node, vpr_id);
		/* new_path.second-1 because we do not update the last node because it's the node in existing route tree (which was already updated previously */
		/*int new_path_first_rr_node = *(new_path.end()-1);*/
		/* if the first node is a connection to existing route, don't update its cost again */
		/*if (route_tree_get_rt_node(rt, new_path_first_rr_node)) {*/
			/*update_one_cost(g, new_path.begin(), new_path.end()-1, 1, params.pres_fac);*/
		/*} else {*/
			/*update_one_cost(g, new_path.begin(), new_path.end(), 1, params.pres_fac);*/
		/*}*/
		update_one_cost(g, new_path.begin(), new_path.end(), 1, params.pres_fac, false);

		/* important to update route tree only after updating cost because the previous check relies on whether the
		 * first node is already in the route tree */
		route_tree_add_path(rt, g, state, new_path, vpr_id);

		net_timing.delay[sink->id+1] = route_tree_get_rt_node(rt, sink->rr_node)->properties.delay;

		for (const auto &m : modified)  {
			state[m].known_cost = std::numeric_limits<float>::max();
			state[m].cost = std::numeric_limits<float>::max();
		}

		heap = std::priority_queue<route_state_t>();
	}

	/*delete [] state;*/

	/*check_route_tree(rt, net, g);*/
}

void route_net(RRGraph &g, const net_t &net, const route_parameters_t &params, route_state_t *state, route_tree_t &rt, trace_t &trace, t_net_timing &net_timing, perf_t *perf)
{
	std::priority_queue<route_state_t> heap;

	vector<int> modified;

	zlog_level(delta_log, ROUTER_V1, "Routing VPR net id %d\n", net.vpr_id);

	if (route_tree_empty(rt)) {
		/* special case */
		/*update_one_cost_internal(get_vertex(g, net.current_source.rr_node), 1, params.pres_fac);*/

		/*route_tree_set_source(rt, get_vertex(g, net.source.rr_node));*/
	} else {
		RouteTreeNode rt_root = get_vertex(rt.graph, rt.root_rt_node_id);
		if (rt_root.properties.rr_node != net.source.rr_node) {
			char root[256];
			char source[256];
			sprintf_rr_node(rt_root.properties.rr_node, root);
			sprintf_rr_node(net.source.rr_node, source);
			zlog_error(delta_log, "Error: Non empty route tree node has root %s that is different from net source %s\n",
					root, source);
			assert(false);
		}
	}

	vector<sink_t> sorted_sinks = net.sinks;
	std::sort(sorted_sinks.begin(), sorted_sinks.end());

	char buffer[256];

	sprintf_rr_node(net.source.rr_node, buffer);
	zlog_level(delta_log, ROUTER_V1, "Source: %s\n", buffer);

	int isink = 0;
	for (const auto &sink : sorted_sinks) {
		sprintf_rr_node(sink.rr_node, buffer);
		zlog_level(delta_log, ROUTER_V1, "Sink %d: %s criticality: %g BB: %d-%d %d-%d\n", sink.id, buffer, sink.criticality_fac, sink.current_bounding_box.xmin, sink.current_bounding_box.xmax, sink.current_bounding_box.ymin, sink.current_bounding_box.ymax);
		route_tree_add_to_heap(rt, g, get_vertex(g, sink.rr_node), sink.criticality_fac, params.astar_fac, sink.current_bounding_box, heap, perf);

		bool found_sink = false;
		while (!heap.empty() && !found_sink) {
			auto item = heap.top();
			heap.pop();

			sprintf_rr_node(item.rr_node, buffer);
			zlog_level(delta_log, ROUTER_V2, "Current: %s prev: %d old_delay: %g new_delay: %g old_known: %g new_known: %g old_cost: %g new_cost: %g\n", buffer, item.prev_edge ? id(get_source(g, *item.prev_edge)) : -1, state[item.rr_node].delay, item.delay, state[item.rr_node].known_cost, item.known_cost, state[item.rr_node].cost, item.cost);

			if (item.rr_node == sink.rr_node) {
				state[item.rr_node] = item;
				modified.push_back(item.rr_node);
				found_sink = true;
			} else if (item.known_cost < state[item.rr_node].known_cost && item.cost < state[item.rr_node].cost) {
				/*if (route_tree_get_rt_node(rt, item.rr_node)) {*/
					/*sprintf_rr_node(item.rr_node, buffer);*/
					/*zlog_warn(delta_log, "Warning: Found a lower cost path to existing route tree node %s\n", buffer);*/
					/*assert(false);*/
				/*}*/
				const auto &sink_vertex = get_vertex(g, sink.rr_node);

				state[item.rr_node] = item;
				modified.push_back(item.rr_node);

				expand_neighbors(g, item, sink_vertex, sink.criticality_fac, params.astar_fac, heap, [&net, &sink, &sink_vertex] (const RRNode &v) -> bool {
					const auto &prop = v.properties;

					if (prop.xhigh < sink.current_bounding_box.xmin
							|| prop.xlow > sink.current_bounding_box.xmax
							|| prop.yhigh < sink.current_bounding_box.ymin
							|| prop.ylow > sink.current_bounding_box.ymax) {
					zlog_level(delta_log, ROUTER_V3, "outside of bounding box\n");
					return false;
					}

					if (prop.type == IPIN
							&& (prop.xhigh != sink_vertex.properties.xhigh 
								|| prop.yhigh != sink_vertex.properties.yhigh)) {
					zlog_level(delta_log, ROUTER_V3, "not target IPIN\n");
					return false;
					}

					/*if (net.sinks.size() >= HIGH_FANOUT_NET_LIM) {*/
					/*if (prop.xhigh < target_x - highfanout_rlim*/
							/*|| prop.xlow > target_x + highfanout_rlim*/
							/*|| prop.yhigh < target_y - highfanout_rlim*/
							/*|| prop.ylow > target_y + highfanout_rlim) {*/
						/*return false;*/
					/*}*/
					/*}*/
					return true;
				}, perf, false);
			} else {
				zlog_level(delta_log, ROUTER_V3, "\tNot expanding neighbor because known_cost %g > %g or cost %g > %g\n", item.known_cost, state[item.rr_node].known_cost, item.cost, state[item.rr_node].cost);
			}
		}

		assert(found_sink);

		const auto &new_path = trace_add_path(trace, g, state, sink.rr_node, net.vpr_id);
		/* new_path.second-1 because we do not update the last node because it's the node in existing route tree (which was already updated previously */
		/*int new_path_first_rr_node = *(new_path.end()-1);*/
		/* if the first node is a connection to existing route, don't update its cost again */
		/*if (route_tree_get_rt_node(rt, new_path_first_rr_node)) {*/
			/*update_one_cost(g, new_path.begin(), new_path.end()-1, 1, params.pres_fac);*/
		/*} else {*/
			/*update_one_cost(g, new_path.begin(), new_path.end(), 1, params.pres_fac);*/
		/*}*/
		update_one_cost(g, new_path.begin(), new_path.end(), 1, params.pres_fac, false);

		/* important to update route tree only after updating cost because the previous check relies on whether the
		 * first node is already in the route tree */
		route_tree_add_path(rt, g, state, new_path, net.vpr_id);

		net_timing.delay[sink.id+1] = route_tree_get_rt_node(rt, sink.rr_node)->properties.delay;

		for (const auto &m : modified)  {
			state[m].known_cost = std::numeric_limits<float>::max();
			state[m].cost = std::numeric_limits<float>::max();
		}

		heap = std::priority_queue<route_state_t>();
	}

	/*delete [] state;*/

	/*check_route_tree(rt, net, g);*/
}

void fine_grained_route_nets(RRGraph &g, vector<net_t> &nets, const route_parameters_t &params, route_tree_t *route_trees, t_net_timing *net_timing)
{
	int inet = 0;
	for (auto &net : nets) {
		/*route_net(g, net, params, route_trees[net.local_id], net_timing[net.vpr_id]);*/
		++inet;
	}
}

bool use_net_level_parallelism()
{
	return true;
}

/*template<typename BoxA, typename BoxB>*/
/*bounding_box_t get_bounding_box(const BoxA &a, const BoxB &b, int bb_factor)*/
/*{*/
	/*bounding_box_t bb;*/

	/*bb.xmin = std::min(a.xmin, b.xmin);*/
	/*bb.xmax = std::max(a.xmax, b.xmax);*/
	/*bb.ymin = std::min(a.ymin, b.ymin);*/
	/*bb.ymax = std::max(a.ymax, b.ymax);*/

	/*bb.xmin -= 1;*/
	/*bb.ymin -= 1;*/

	/*extern int nx;*/
	/*extern int ny;*/

	/*bb.xmin = std::max(bb.xmin - bb_factor, 0);*/
	/*bb.xmax = std::min(bb.xmax + bb_factor, nx + 1);*/
	/*bb.ymin = std::max(bb.ymin - bb_factor, 0);*/
	/*bb.ymax = std::min(bb.ymax + bb_factor, ny + 1);*/

	/*return bb;*/
/*}*/

void test_interval()
{
	/*using namespace boost::numeric;*/
	/*interval<int> i1(1,2);*/
	/*interval<int> i2(3,3);*/

	/*if (overlap(i1, i2)) {*/
		/*printf("1 overlap\n");*/
	/*}*/
}

void test_scheduler()
{
	/*using namespace boost::numeric;*/
	/*vector<pair<interval<int>, int>> intervals = {*/
		/*{ { 1, 5 }, 0 },*/
		/*{ { 6, 7 }, 1 },*/
		/*{ { 3, 10 }, 2 },*/
		/*{ { 8, 9 }, 3 }*/
	/*};*/
	/*vector<pair<interval<int>, int> *> chosen;*/
	/*max_independent_intervals(intervals, chosen);*/
	/*for (const auto &c : chosen) {*/
		/*printf("chosen: [%d, %d]\n", c->first.lower(), c->first.upper());*/
	/*}*/

	vector<bounding_box_t> bbs = {
		{ 1, 5, 3, 4 },
		{ 9, 10, 3, 4 },
		{ 1, 5, 3, 4 },
		{ 9, 10, 3, 4 }
	};
	vector<pair<const bounding_box_t *, int> *> bbs_ptr;

	for (const auto &b : bbs) {
		bbs_ptr.push_back(new pair<const bounding_box_t *, int>(&b, 0));
	}
	vector<pair<const bounding_box_t *, int> *> chosen_bbs;
	max_independent_rectangles(bbs_ptr, chosen_bbs);
	for (const auto &c : chosen_bbs) {
		printf("chosen: x %d-%d y %d-%d]\n", c->first->xmin, c->first->xmax, c->first->ymin, c->first->ymax);
	}
}

void dump_net_bounding_boxes(const vector<net_t> &nets)
{
	vector<net_t> sorted_nets = nets;

	sort(sorted_nets.begin(), sorted_nets.end(), [] (const net_t &a, const net_t &b) -> bool {
			return get_bounding_box_area(a.bounding_box) > get_bounding_box_area(b.bounding_box);
			});

	int rank = 0;
	for (const auto &net : sorted_nets) {
		char filename[256];
		sprintf(filename, "/Volumes/DATA/net_bb/%d.txt", rank);
		FILE *file = fopen(filename, "w");
		for (int x = net.bounding_box.xmin; x <= net.bounding_box.xmax; ++x) {
			for (int y = net.bounding_box.ymin; y <= net.bounding_box.ymax; ++y) {
				fprintf(file, "%d %d\n", x, y);
			}
		}
		fclose(file);
		++rank;
	}
}

void write_metis_file(const graph_t<int, int> &dep_g, const char *filename)
{
	FILE *file = fopen("metis.graph", "w");
	assert(file);
	assert(num_edges(dep_g) % 2 == 0);
	fprintf(file, "%d %d\n", num_vertices(dep_g), num_edges(dep_g)/2);
	for_all_vertices(dep_g, [&dep_g, &file] (const vertex_t<int, int> &v) -> void {
			for_all_out_edges(dep_g, v, [&dep_g, &file] (const edge_t<int> &e) -> void {
					fprintf(file, "%d ", id(get_target(dep_g, e))+1);
					});
			fprintf(file, "\n");
			});
	fclose(file);
}

void load_overlapping_nets_vec(vector<net_t *> &nets)
{
	for (int i = 0; i < nets.size(); ++i) {
		nets[i]->overlapping_nets_vec.clear();
		nets[i]->non_overlapping_nets_vec.clear();
	}

	for (int i = 0; i < nets.size(); ++i) {
		for (int j = i+1; j < nets.size(); ++j) {
			if (box_overlap(nets[i]->sinks[nets[i]->current_sink_index].current_bounding_box, nets[j]->sinks[nets[j]->current_sink_index].current_bounding_box)) {
				nets[i]->overlapping_nets_vec.push_back(nets[j]);
				nets[j]->overlapping_nets_vec.push_back(nets[i]);
			} else {
				nets[i]->non_overlapping_nets_vec.push_back(nets[j]);
				nets[j]->non_overlapping_nets_vec.push_back(nets[i]);
			}
		}
		/*zlog_debug(net_log, "Net %i overlaps with %d nets\n", i, num_overlaps);*/
	}
}

void load_overlapping_nets(vector<net_t *> &nets)
{
	for (auto &net : nets) {
		/*for (int i = 0; i < n.num_local_nets; ++i) {*/
			/*n.overlapping_nets[i] = false;*/
			/*n.non_overlapping_nets[i] = false;*/
			/*n.num_overlapping_nets = 0;*/
			/*n.num_non_overlapping_nets = 0;*/
		/*}*/
		net->overlapping_nets = vector<bool>(nets.size(), false);
		net->non_overlapping_nets = vector<bool>(nets.size(), false);
		net->num_overlapping_nets = 0;
		net->num_non_overlapping_nets = 0;
	}

	for (int i = 0; i < nets.size(); ++i) {
		for (int j = i+1; j < nets.size(); ++j) {
			if (box_overlap(nets[i]->sinks[nets[i]->current_sink_index].current_bounding_box, nets[j]->sinks[nets[j]->current_sink_index].current_bounding_box)) {
				assert(nets[j]->current_local_id < nets[i]->overlapping_nets.size());
				nets[i]->overlapping_nets[nets[j]->current_local_id] = true;
				assert(nets[i]->current_local_id < nets[j]->overlapping_nets.size());
				nets[j]->overlapping_nets[nets[i]->current_local_id] = true;

				++nets[i]->num_overlapping_nets;
				++nets[j]->num_overlapping_nets;
			} else {
				assert(nets[j]->current_local_id < nets[i]->non_overlapping_nets.size());
				nets[i]->non_overlapping_nets[nets[j]->current_local_id] = true;
				assert(nets[i]->current_local_id < nets[j]->non_overlapping_nets.size());
				nets[j]->non_overlapping_nets[nets[i]->current_local_id] = true;

				++nets[i]->num_non_overlapping_nets;
				++nets[j]->num_non_overlapping_nets;
			}
		}
		/*zlog_debug(net_log, "Net %i overlaps with %d nets\n", i, num_overlaps);*/
	}
}

int get_quadrant(const bounding_box_t &outer_box, const bounding_box_t &inner_box, float &distance)
{
	point_t<float> inner_centroid;
	inner_centroid.x = (float)(inner_box.xmin + inner_box.xmax)/2;
	inner_centroid.y = (float)(inner_box.ymin + inner_box.ymax)/2;

	point_t<float> outer_centroid;
	outer_centroid.x = (float)(outer_box.xmin + outer_box.xmax)/2;
	outer_centroid.y = (float)(outer_box.ymin + outer_box.ymax)/2;

	float x_offset = inner_centroid.x - outer_centroid.x;
	float y_offset = inner_centroid.y - outer_centroid.y;

	distance = sqrt(x_offset*x_offset + y_offset*y_offset);

	int quadrant;
	if (x_offset > 0) {
		if (y_offset > 0) {
			quadrant = 0;
		} else {
			quadrant = 1;
		}
	} else {
		if (y_offset > 0) {
			quadrant = 3;
		} else {
			quadrant = 2;
		}
	}
	return quadrant;
}

void test_quadtree(const vector<net_t> &nets)
{
	extern int nx;
	extern int ny;
	bounding_box_t bb;
	bb.xmin = 0;
	bb.ymin = 0;
	bb.xmax = nx+2;
	bb.ymax = ny+2;
	QuadTree<const net_t *> qt(bb, 4);
	auto start = std::chrono::high_resolution_clock::now();

	for (const auto &net : nets) {
		qt.insert({ net.sinks[0].current_bounding_box, &net });
	}

	auto elapsed = std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::high_resolution_clock::now()-start);
	printf("quadtree insertion of %lu nets took %lld us\n", nets.size(), elapsed.count());

	int num_items = 0;
	qt.verify(num_items);
	assert(num_items == nets.size());

	start = std::chrono::high_resolution_clock::now();
	for (int i = 0; i < nets.size(); ++i) {
		vector<const pair<bounding_box_t, const net_t *> *> res;
		qt.get_overlapping_boxes(nets[i].sinks[0].current_bounding_box, res);
		vector<const net_t *> v1;
		for (const auto &r : res) {
			if (r->second != &nets[i]) {
				v1.push_back(r->second);
			}
		}
		sort(v1.begin(), v1.end());

		vector<const net_t *> v2;
		for (int j = 0; j < nets.size(); ++j) {
			if (j != i && box_overlap(nets[i].sinks[0].current_bounding_box,
						nets[j].sinks[0].current_bounding_box)) {
				v2.push_back(&nets[j]);
			}
		}
		sort(v2.begin(), v2.end());
		assert(v1 == v2);
	}

	elapsed = std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::high_resolution_clock::now()-start);
	printf("quadtree queries of %lu nets took %lld us\n", nets.size(), elapsed.count());
}

void test_quadtree()
{
	bounding_box_t bb;
	bb.xmin = 0;
	bb.xmax = 20;
	bb.ymin = 0;
	bb.ymax = 20;
	bounding_box_t item1;
	item1.xmin = 0;
	item1.xmax = 20;
	item1.ymin = 0;
	item1.ymax = 20;
	QuadTree<int> qt(bb, 1);
	qt.insert({ item1, 0 });
	item1.xmin = 0;
	item1.xmax = 9;
	item1.ymin = 0;
	item1.ymax = 9;
	qt.insert({ item1, 0 });
	item1.xmin = 11;
	item1.xmax = 20;
	item1.ymin = 0;
	item1.ymax = 9;
	qt.insert({ item1, 0 });

	vector<const pair<bounding_box_t, int> *> res;
	qt.get_overlapping_boxes(bb, res);
}

void test_quadrant()
{
	bounding_box_t outer_box;
	outer_box.xmin = 0;
	outer_box.ymin = 0;
	outer_box.xmax = 20;
	outer_box.ymax = 20;

	float distance;

	assert(get_quadrant(outer_box, { 16, 17, 16, 17 }, distance) == 0);
	assert(get_quadrant(outer_box, { 16, 17, 6, 7 }, distance) == 1);
	assert(get_quadrant(outer_box, { 6, 7, 6, 7 }, distance) == 2);
	assert(get_quadrant(outer_box, { 6, 7, 16, 17 }, distance) == 3);
	assert(get_quadrant(outer_box, { 9, 11, 9, 11 }, distance) == 2);
}

void choose_sinks(vector<net_t> &nets, vector<vector<bool>> &sink_scheduled, vector<sink_t *> &scheduled_sinks)
{
	int current_quadrant = 0;

	int quadrant_utilization[4];
	for (int i = 0; i < 4; ++i) {
		quadrant_utilization[i] = 0;
	}

	bounding_box_t fpga_box;
	extern int nx;
	extern int ny;
	fpga_box.xmin = 0;
	fpga_box.ymin = 0;
	fpga_box.xmax = nx+2;
	fpga_box.ymax = ny+2;

	for (auto &net : nets) {
		zlog_debug(sort_log, "Net %d\n", net.vpr_id);
		int min_utilization_quadrant = -1;
		int min_utilization = std::numeric_limits<int>::max();
		float max_distance = std::numeric_limits<float>::lowest();
		sink_t *scheduled_sink = nullptr;
		zlog_debug(sort_log, "Max distance: %g\n", max_distance);
		for (auto sink = begin(net.sinks); sink != end(net.sinks); ++sink) {	
			if (!sink_scheduled[net.current_local_id][sink->id]) {
				float distance;
				int sink_quadrant = get_quadrant(fpga_box, sink->current_bounding_box, distance);
				zlog_debug(sort_log, "Sink %d bounding box %d-%d %d-%d quadrant %d distance %g", sink->id, sink->current_bounding_box.xmin, sink->current_bounding_box.xmax, sink->current_bounding_box.ymin, sink->current_bounding_box.ymax, sink_quadrant, distance);
				/*if (quadrant_utilization[sink_quadrant] < min_utilization) {*/
				if (distance > max_distance) {
					min_utilization_quadrant = sink_quadrant;
					/*min_utilization = quadrant_utilization[sink_quadrant];*/
					max_distance = distance;
					scheduled_sink = &(*sink);
					zlog_debug(sort_log, " UPDATING");
				}
				zlog_debug(sort_log, "\n");
			} else {
				zlog_debug(sort_log, "Sink %d scheduled\n", sink->id);
			}
		}

		assert(min_utilization_quadrant >= 0 && min_utilization_quadrant < 4);
		++quadrant_utilization[min_utilization_quadrant];
		sink_scheduled[net.current_local_id][scheduled_sink->id] = true;
		scheduled_sinks.push_back(scheduled_sink);

		zlog_debug(sort_log, "Scheduled sink %d quadrant %d\n", scheduled_sink->id, min_utilization_quadrant);
	}

	int max_imbalance = 0;
	for (int i = 0; i < 4; ++i) {
		for (int j = i+1; j < 4; ++j) {
			int imba = abs(quadrant_utilization[i]-quadrant_utilization[j]);
			if (imba > max_imbalance) {
				max_imbalance = imba;
			}
		}
	}
	for (int i = 0; i < 4; ++i) {
		zlog_info(sort_log, "Quadrant %d utilization: %d\n", i, quadrant_utilization[i]);
	}
	zlog_info(sort_log, "Max quadrant imbalance: %d\n", max_imbalance);
}

void init_overlapping_nets(vector<net_t> &nets, graph_t<int, int> &dep_g)
{
	/*add_vertex(dep_g, nets.size());*/
	/*for (auto &n : nets) {*/
		/*n.overlapping_nets = vector<bool>(nets.size(), false);*/
		/*n.non_overlapping_nets = vector<bool>(nets.size(), false);*/
	/*}*/

	/*for (int i = 0; i < nets.size(); ++i) {*/
		/*int num_overlaps = 0;*/
		/*for (int j = i+1; j < nets.size(); ++j) {*/
			/*if (box_overlap(bb_a, bb_b)) {*/
				/*++num_overlaps;*/
				/*add_edge(dep_g, i, j);*/
				/*add_edge(dep_g, j, i);*/

				/*nets[i].overlapping_nets[nets[i].id] = true;*/
				/*nets[j].overlapping_nets[nets[i].id] = true;*/
			/*} else {*/
				/*nets[i].non_overlapping_nets[nets[j].id] = true;*/
				/*nets[j].non_overlapping_nets[nets[i].id] = true;*/
			/*}*/
		/*}*/
		/*zlog_debug(net_log, "Net %i overlaps with %d nets\n", i, num_overlaps);*/
	/*}*/
}

typedef struct partition_t {
	vector<const net_t *> nets;
} partition_t;

void init_nets_partition(vector<net_t> &nets, const char *part_file, int num_parts, vector<partition_t> &partitions)
{
	FILE *file = fopen(part_file, "r");
	int i = 0;
	printf("Partitions\n");
	while (!feof(file)) {
		fscanf(file, "%d\n", &nets[i].pid);
		printf("Net %d pid %d\n", i, nets[i].pid);	
		assert(nets[i].pid < num_parts && nets[i].pid >= 0);
		partitions[nets[i].pid].nets.push_back(&nets[i]);
		++i;
	}
	
	assert(i == nets.size());
}

template<typename T>
inline bool in_set(const set<T> &s, const T &item)
{
	return s.find(item) != s.end();
}

net_t *schedule_net_greedy(const vector<net_t *> &nets, vector<bool> &scheduled_nets, vector<const net_t *> &current_time_scheduled_nets, int time)
{
	int max_num_non_overlapping_nets = std::numeric_limits<int>::min();	
	net_t *res = nullptr;

	/*zlog_debug(schedule_log, "-- Start scheduling for time %d\n", time);*/

	zlog_debug(schedule_log, "Scheduled nets: ");
	for (const auto &n : nets) {
		if (scheduled_nets[n->current_local_id]) {
			zlog_debug(schedule_log, "%d ", n->current_local_id);
		}
	}
	zlog_debug(schedule_log, "\n");

	zlog_debug(schedule_log, "Current time scheduled nets: ");
	for (const auto &n : current_time_scheduled_nets) {
		/*zlog_debug(schedule_log, "[%d BB: %d] ", n->current_local_id, get_bounding_box_area(n->current_bounding_box));*/
	}
	zlog_debug(schedule_log, "\n");

	for (const auto &net : nets) {
		bool net_scheduled = scheduled_nets[net->current_local_id];
		/*zlog_debug(schedule_log, "Checking net %d non_overlapping_nets %d net_scheduled %d ", net->id, num_non_overlapping_nets, net_scheduled ? 1 : 0);*/

		if (!net_scheduled && net->num_non_overlapping_nets > max_num_non_overlapping_nets) {
			bool overlap_scheduled_nets;

			/*if (current_time_scheduled_nets.empty()) {*/
				/*overlap_scheduled_nets = false;*/
			/*} else {*/
				/*overlap_scheduled_nets = any_of(current_time_scheduled_nets.begin(), current_time_scheduled_nets.end(), [&net] (const net_t *current_time_net) -> bool {*/
					/*return net->overlapping_nets[current_time_net->id];*/
					/*});*/
			/*}*/
			overlap_scheduled_nets = any_of(current_time_scheduled_nets.begin(), current_time_scheduled_nets.end(), [&net] (const net_t *current_time_net) -> bool {
					return net->overlapping_nets[current_time_net->current_local_id];
					});

			zlog_debug(schedule_log, "Net %d non_overlapping_nets %d not scheduled and overlap_scheduled_nets %d\n", net->current_local_id, net->num_non_overlapping_nets, overlap_scheduled_nets ? 1 : 0);

			if (!overlap_scheduled_nets) {
				max_num_non_overlapping_nets = net->num_non_overlapping_nets;
				res = net;
				zlog_debug(schedule_log, "\tRecording net %d\n", net->current_local_id);
			}
		}

		/*zlog_debug(schedule_log, "\n");*/
	}

	if (res) {
		res->schedule = time;
		scheduled_nets[res->current_local_id] = true;
		current_time_scheduled_nets.push_back(res);
		zlog_debug(schedule_log, "Scheduled net %d\n", res->current_local_id);
	}

	/*zlog_debug(schedule_log, "-- End scheduling for time %d\n", time);*/

	return res;
}

void schedule_nets_num_sinks(vector<net_t *> &nets, vector<vector<const net_t *>> &net_scheduled_at_time)
{
	sort(nets.begin(), nets.end(), [] (const net_t *a, const net_t *b) -> bool {
			return a->sinks.size() > b->sinks.size();
			});

	if (nets.size() > 0) {
		net_scheduled_at_time.push_back(vector<const net_t *>());
		for (const auto &net : nets) {
			net_scheduled_at_time[0].push_back(net);
		}
	}
}

void schedule_nets_greedy(vector<net_t *> &nets, vector<vector<const net_t *>> &net_scheduled_at_time)
{
	cpu_timer timer;
	timer.start();
	load_overlapping_nets(nets);
	timer.stop();
	cpu_times elapsed = timer.elapsed();
	/*zlog_info(delta_log, "Scheduling took %s\n", format(elapsed).c_str());*/
	zlog_info(delta_log, "Load overlapping net took %g\n", elapsed.wall / 1e9);

	/* we only schedule nets that has sinks to be routed */
	vector<bool> scheduled_nets(nets.size(), false);
	
	vector<const net_t *> current_time_scheduled_nets;

	timer.start();

	int time = 0;

	zlog_debug(schedule_log, "-- Start scheduling nets for time %d\n", time);
	net_t *scheduled_net = schedule_net_greedy(nets, scheduled_nets, current_time_scheduled_nets, time);
	if (!scheduled_net) {
		/* no more nets to schedule */
		return;
	}
	net_scheduled_at_time.push_back(vector<const net_t *>());
	net_scheduled_at_time[time].push_back(scheduled_net);
	int num_scheduled_nets = 1;

	while (num_scheduled_nets < nets.size()) {
		vector<net_t *> temp;
		for (const auto &net : nets) {
			if (scheduled_net->non_overlapping_nets[net->current_local_id]) {
				temp.push_back(net);
			}
		}
		scheduled_net = schedule_net_greedy(temp, scheduled_nets, current_time_scheduled_nets, time);
		if (!scheduled_net) {
			current_time_scheduled_nets.clear();
			++time;

			zlog_debug(schedule_log, "-- Start scheduling nets for time %d\n", time);
			scheduled_net = schedule_net_greedy(nets, scheduled_nets, current_time_scheduled_nets, time);
			assert(scheduled_net);
			net_scheduled_at_time.push_back(vector<const net_t *>());
			net_scheduled_at_time[time].push_back(scheduled_net);
		} else {
			net_scheduled_at_time[time].push_back(scheduled_net);
		}
		++num_scheduled_nets;
	}
	for (int i = 0; i < nets.size(); ++i) {
		/*zlog_info(schedule_log, "Net %d scheduled %d\n", i, scheduled_nets[i] ? 1 : 0);*/
	}
	assert(all_of(scheduled_nets.begin(), scheduled_nets.end(), [] (const bool &item) -> bool { return item; }));

	timer.stop();
	elapsed = timer.elapsed();
	/*zlog_info(delta_log, "Scheduling took %s\n", format(elapsed).c_str());*/
	zlog_info(delta_log, "Scheduling took %g\n", elapsed.wall / 1e9);

	/*delete [] scheduled_nets;*/
}

void verify_scheduling()
{
}

template<typename Nets>
net_t *schedule_net(const Nets &nets, set<int> &scheduled_net, set<int> &scheduled_partition, int time)
{
	int max_num_overlapping_nets = std::numeric_limits<int>::min();	
	net_t *res = nullptr;

	zlog_debug(schedule_log, "-- Start scheduling for time %d\n", time);

	zlog_debug(schedule_log, "Scheduled partitions: ");
	for (const auto &p : scheduled_partition) {
		zlog_debug(schedule_log, "%d ", p);
	}
	zlog_debug(schedule_log, "\n");

	zlog_debug(schedule_log, "Scheduled nets: ");
	for (const auto &n : scheduled_net) {
		zlog_debug(schedule_log, "%d ", n);
	}
	zlog_debug(schedule_log, "\n");

	for (const auto &net : nets) {
		bool net_scheduled = in_set(scheduled_net, net->local_id);
		bool part_scheduled = in_set(scheduled_partition, net->pid);

		zlog_debug(schedule_log, "Checking net %d pid %d non_overlapping_nets %d net_scheduled %d part_scheduled %d\n", net->local_id, net->pid, net->num_non_overlapping_nets, net_scheduled ? 1 : 0, part_scheduled ? 1 : 0);
		
		if (!net_scheduled && !part_scheduled
				&& net->num_non_overlapping_nets > max_num_overlapping_nets) {
			max_num_overlapping_nets = net->num_non_overlapping_nets ;
			res = net;
			zlog_debug(schedule_log, "\tRecording net %d\n", net->local_id);
		}
	}

	if (res) {
		res->schedule = time;
		scheduled_net.insert(res->local_id);
		scheduled_partition.insert(res->pid);
		zlog_debug(schedule_log, "Scheduled net %d pid %d\n", res->local_id, res->pid);
	}

	zlog_debug(schedule_log, "-- End scheduling for time %d\n", time);

	return res;
}

void schedule_nets(vector<net_t> &nets, int num_parts)
{
	for (auto &n : nets) {
		n.schedule = -1;
	}

	vector<net_t *> bootstrap;
	for (auto &n : nets) {
		bootstrap.push_back(&n);
	}
	set<int> scheduled_nets;
	set<int> scheduled_partitions;

	int time = 0;
	net_t *scheduled_net = schedule_net(bootstrap, scheduled_nets, scheduled_partitions, time);
	while (scheduled_nets.size() < nets.size()) {
		/*scheduled_net = schedule_net(scheduled_net->non_overlapping_nets, scheduled_nets, scheduled_partitions, time);*/
		if (!scheduled_net) {

			scheduled_partitions.clear();
			++time;
			scheduled_net = schedule_net(bootstrap, scheduled_nets, scheduled_partitions, time);
			assert(scheduled_net);
		} else if (scheduled_partitions.size() == num_parts) {
			scheduled_partitions.clear();
			++time;
		}
	}
}

/*void schedule_nets(vector<net_t> &nets, int num_parts)*/
/*{*/
	/*for (auto &n : nets) {*/
		/*n.schedule = -1;*/
	/*}*/

	/*int max_time;*/
	/*set<int> scheduled;*/
	/*for (int time = 0; time < max_time; ++time) {*/
		/*const net_t *cur_part_min = &nets[0];*/

		/*for (int i = 1; i < nets.size(); ++i) {*/
			/*auto iter = scheduled.find(nets[i].id);*/
			/*if (nets[i].non_overlapping_nets.size() < cur_part_min->non_overlapping_nets.size()*/
					/*&& iter == scheduled.end()) {*/
				/*cur_part_min = &nets[i];*/
			/*}*/
		/*}*/
		/*cur_part_min->schedule = time;*/
		/*scheduled.insert(cur_part_min->id);*/

		/*for (int i = 0; i < num_parts; ++i) { */
			/*const net_t *other_part_min = nullptr;*/
			/*for (int j = 0; j < cur_part_min->non_overlapping_nets.size(); ++j) {*/
				/*if (other_part_min == nullptr) {*/
					/*if (cur_part_min->non_overlapping_nets[j]->pid != cur_part_min->pid) {*/
						/*other_part_min = &cur_part_min->non_overlapping_nets[j];*/
					/*}*/
				/*} else if (cur_part_min->non_overlapping_nets[j]->pid != cur_part_min->pid && cur_part_min->non_overlapping_nets[j].size() < other_part_min->non_overlapping_nets.size()) {*/
					/*other_part_min = &cur_part_min->non_overlapping_nets[j];*/
				/*}*/
			/*}*/
		/*}*/
	/*}*/
/*}*/

void schedule_nets_faster(vector<net_t *> &nets, vector<net_t *> &output)
{
	/*load_overlapping_nets_vec(nets);*/

	/*sort(nets.begin(), nets.end(), [] (const net_t *a, const net_t *b) { return a->non_overlapping_nets_vec.size() < b->non_overlapping_nets_vec.size(); });*/
	sort(nets.begin(), nets.end(), [] (const net_t *a, const net_t *b) { return get_bounding_box_area(a->sinks[a->current_sink_index].current_bounding_box) < get_bounding_box_area(b->sinks[b->current_sink_index].current_bounding_box); });

	/*printf("Num overlapping_nets:\n");*/
	/*for (const auto &net : nets) {*/
		/*printf("Net %d: %d\n", net->current_local_id, net->overlapping_nets_vec.size());*/
	/*}*/

	for (int i = 0; i < nets.size(); ++i) {
		bool overlap_scheduled_nets = any_of(begin(output), end(output), [&i, &nets] (const net_t *scheduled_net) -> bool {
				return box_overlap(nets[i]->sinks[nets[i]->current_sink_index].current_bounding_box,
						scheduled_net->sinks[scheduled_net->current_sink_index].current_bounding_box);
				});

		/*for (const auto &scheduled_net : output) { */
			/*if (box_overlap(nets[i]->sinks[nets[i]->current_sink_index].current_bounding_box,*/
						/*scheduled_net->sinks[scheduled_net->current_sink_index].current_bounding_box)) {*/
				/*zlog_debug(schedule_log, "Net %d overlaps with net %d\n", nets[i]->vpr_id, scheduled_net->vpr_id);*/
			/*}*/
		/*}*/

		if (!overlap_scheduled_nets) {
			output.push_back(nets[i]);
		}
	}

	/*verify_ind(output);*/
}

void schedule_nets_fast(vector<net_t *> &nets, vector<vector<const net_t *>> &net_scheduled_at_time)
{
	/*load_overlapping_nets_vec(nets);*/

	/*sort(nets.begin(), nets.end(), [] (const net_t *a, const net_t *b) { return a->non_overlapping_nets_vec.size() < b->non_overlapping_nets_vec.size(); });*/
	sort(nets.begin(), nets.end(), [] (const net_t *a, const net_t *b) { return get_bounding_box_area(a->sinks[a->current_sink_index].current_bounding_box) > get_bounding_box_area(b->sinks[b->current_sink_index].current_bounding_box); });

	zlog_debug(schedule_log, "Nets sorted by current bounding box size:\n");
	for (const auto &net : nets) {
		zlog_debug(schedule_log, "Net %d Area %d\n", net->vpr_id, get_bounding_box_area(net->sinks[net->current_sink_index].current_bounding_box));
	}

	/*printf("Num overlapping_nets:\n");*/
	/*for (const auto &net : nets) {*/
		/*printf("Net %d: %d\n", net->current_local_id, net->overlapping_nets_vec.size());*/
	/*}*/

	vector<bool> scheduled_nets(nets.size(), false);
	net_scheduled_at_time.reserve(nets.size()); //worst case
	int num_scheduled_nets = 0;
	while (num_scheduled_nets < nets.size()) {
		net_scheduled_at_time.push_back(vector<const net_t *>());

		for (int i = 0; i < nets.size(); ++i) {
			if (!scheduled_nets[i]) {
				bool overlap_scheduled_nets = any_of(net_scheduled_at_time.back().begin(), net_scheduled_at_time.back().end(), [&i, &nets] (const net_t *scheduled_net) -> bool {
						return box_overlap(nets[i]->sinks[nets[i]->current_sink_index].current_bounding_box,
								scheduled_net->sinks[scheduled_net->current_sink_index].current_bounding_box);
					});

				for (const auto &scheduled_net : net_scheduled_at_time.back()) { 
					if (box_overlap(nets[i]->sinks[nets[i]->current_sink_index].current_bounding_box,
								scheduled_net->sinks[scheduled_net->current_sink_index].current_bounding_box)) {
						zlog_debug(schedule_log, "Net %d overlaps with net %d\n", nets[i]->vpr_id, scheduled_net->vpr_id);
					}
				}

				if (!overlap_scheduled_nets) {
					scheduled_nets[i] = true;
					net_scheduled_at_time.back().push_back(nets[i]);
					++num_scheduled_nets;
				}
			}
		}
	}

	assert(all_of(scheduled_nets.begin(), scheduled_nets.end(), [] (const bool &item) -> bool { return item; }));

	for (const auto &at_time : net_scheduled_at_time) {
		verify_ind(at_time);
	}
}

void analyze_timing(t_net_timing *net_timing) 
{
	load_timing_graph_net_delays_new(net_timing); 
#ifdef HACK_LUT_PIN_SWAPPING
	do_timing_analysis_new(net_timing, FALSE, TRUE, FALSE);
#else
	do_timing_analysis_new(net_timing, FALSE, FALSE, FALSE);
#endif

	float critical_path_delay = get_critical_path_delay();
	zlog_info(delta_log, "Critical path: %g ns\n", critical_path_delay);
	printf("Critical path: %g ns\n", critical_path_delay);
}

bool feasible_routing(const RRGraph &g)
{
	bool feasible = true;
	for (const auto &v : get_vertices(g)) {
		if (v.properties.occ > v.properties.capacity) {
			feasible = false;
			break;
		}
	}

	return feasible;
}

void run_metis(const char *graph_filename, int num_partitions)
{
	char buffer[256];
	sprintf(buffer, "./gpmetis %s %d", graph_filename, num_partitions);
	FILE *pipe = popen(buffer, "r");
	while (!feof(pipe)) {
		fgets(buffer, 256, pipe);
		printf("%s", buffer);
	}
	pclose(pipe);
}

void set_current_source(net_t &net, route_tree_t &rt, const RRGraph &g, const sink_t &target, float astar_fac)
{
	/*net.previous_source = net.current_source;*/
	/*net.previous_source_valid = true;*/
	/*assert(net.previous_sink_valid);*/
	/*net.previous_bounding_box = get_bounding_box(net.previous_source, net.sinks[net.previous_sink_index]);*/
	/*net.previous_bounding_box_valid = true;*/

	/*if (route_tree_empty(rt)) {*/
		/*net.current_source = net.source;*/
	/*} else {*/
		/*std::priority_queue<route_state_t> heap;*/
		/*route_tree_add_to_heap(rt, g, get_vertex(g, target.rr_node), target.criticality_fac, astar_fac, net.previous_bounding_box, heap);*/
		/*net.current_source.rr_node = heap.top().rr_node;*/
		/*extern struct s_net *clb_net;*/
		/*extern struct s_block *block;*/
		/*int b = clb_net[net.id].node_block[0];*/
		/*int p = clb_net[net.id].node_block_pin[0];*/
		/*net.current_source.x = block[b].x;*/
		/*net.current_source.y = block[b].y + block[b].type->pin_height[p];*/
	/*}*/
	/*net.current_bounding_box = get_bounding_box(net.current_source, net.sinks[net.current_sink_index]);*/
}

void adjust_bounding_box(sink_t &current_sink)
{
	if (current_sink.current_bounding_box.xmin != current_sink.previous_bounding_box.xmin ||
			current_sink.current_bounding_box.ymin != current_sink.previous_bounding_box.ymin || 
			current_sink.current_bounding_box.xmax != current_sink.previous_bounding_box.xmax ||
			current_sink.current_bounding_box.ymax != current_sink.previous_bounding_box.ymax) {
		/*char buffer[256];*/
		/*sprintf_rr_node(current_sink.rr_node, buffer);*/
		zlog_level(scheduler_log, ROUTER_V2, "Net %d sink %d bounding box differs. Prev: %d-%d %d-%d Cur: %d-%d %d-%d\n", current_sink.net->vpr_id, current_sink.id,
				current_sink.current_bounding_box.xmin, current_sink.current_bounding_box.xmax,
				current_sink.current_bounding_box.ymin, current_sink.current_bounding_box.ymax,
				current_sink.previous_bounding_box.xmin, current_sink.previous_bounding_box.xmax,
				current_sink.previous_bounding_box.ymin, current_sink.previous_bounding_box.ymax);
	}

	/*zlog_level(delta_log, ROUTER_V1, "Net %d sink %d bounding boxes. Prev: %d-%d %d-%d Cur: %d-%d %d-%d\n", current_sink.net->vpr_id, current_sink.id,*/
			/*current_sink.current_bounding_box.xmin, current_sink.current_bounding_box.ymin,*/
			/*current_sink.current_bounding_box.xmax, current_sink.current_bounding_box.ymax,*/
			/*current_sink.previous_bounding_box.xmin, current_sink.previous_bounding_box.ymin,*/
			/*current_sink.previous_bounding_box.xmax, current_sink.previous_bounding_box.ymax);*/

	current_sink.current_bounding_box.xmin = std::min(current_sink.current_bounding_box.xmin, current_sink.previous_bounding_box.xmin);
	current_sink.current_bounding_box.ymin = std::min(current_sink.current_bounding_box.ymin, current_sink.previous_bounding_box.ymin);
	current_sink.current_bounding_box.xmax = std::max(current_sink.current_bounding_box.xmax, current_sink.previous_bounding_box.xmax);
	current_sink.current_bounding_box.ymax = std::max(current_sink.current_bounding_box.ymax, current_sink.previous_bounding_box.ymax);
}

void adjust_bounding_box(net_t &net)
{
	assert(net.current_sink_index >= 0 && net.current_sink_index < net.sinks.size());
	auto &current_sink = net.sinks[net.current_sink_index];

	if (current_sink.current_bounding_box.xmin != current_sink.previous_bounding_box.xmin ||
			current_sink.current_bounding_box.ymin != current_sink.previous_bounding_box.ymin || 
			current_sink.current_bounding_box.xmax != current_sink.previous_bounding_box.xmax ||
			current_sink.current_bounding_box.ymax != current_sink.previous_bounding_box.ymax) {
		char buffer[256];
		sprintf_rr_node(current_sink.rr_node, buffer);
		zlog_level(delta_log, ROUTER_V2, "Net %d %s bounding box differs. Prev: %d-%d %d-%d Cur: %d-%d %d-%d\n", net.vpr_id, buffer,
				current_sink.current_bounding_box.xmin, current_sink.current_bounding_box.ymin,
				current_sink.current_bounding_box.xmax, current_sink.current_bounding_box.ymax,
				current_sink.previous_bounding_box.xmin, current_sink.previous_bounding_box.ymin,
				current_sink.previous_bounding_box.xmax, current_sink.previous_bounding_box.ymax);
	}

	current_sink.current_bounding_box.xmin = std::min(current_sink.current_bounding_box.xmin, current_sink.previous_bounding_box.xmin);
	current_sink.current_bounding_box.ymin = std::min(current_sink.current_bounding_box.ymin, current_sink.previous_bounding_box.ymin);
	current_sink.current_bounding_box.xmax = std::max(current_sink.current_bounding_box.xmax, current_sink.previous_bounding_box.xmax);
	current_sink.current_bounding_box.ymax = std::max(current_sink.current_bounding_box.ymax, current_sink.previous_bounding_box.ymax);
}

int adjust_bounding_box_2(net_t &net)
{
	int num_changed = 0;
	for (auto &sink : net.sinks) {
		if (!net.sink_routed[sink.id]) {
			zlog_level(delta_log, ROUTER_V2, "Adjust bounding box for net %d sink %d. Current x %d-%d y %d-%d Previous x %d-%d y %d-%d\n", net.vpr_id, sink.id, sink.current_bounding_box.xmin, sink.current_bounding_box.xmax, sink.current_bounding_box.ymin, sink.current_bounding_box.ymax, sink.previous_bounding_box.xmin, sink.previous_bounding_box.xmax, sink.previous_bounding_box.ymin, sink.previous_bounding_box.ymax);

			bool changed = false;
			if (sink.current_bounding_box.xmin != sink.previous_bounding_box.xmin) {
				changed = true;
			}
			if (sink.current_bounding_box.xmax != sink.previous_bounding_box.xmax) {
				changed = true;
			}
			if (sink.current_bounding_box.ymin != sink.previous_bounding_box.ymin) {
				changed = true;
			}
			if (sink.current_bounding_box.ymax != sink.previous_bounding_box.ymax) {
				changed = true;
			}

			if (changed) {
				++num_changed;
			}

			sink.current_bounding_box.xmin = std::min(sink.current_bounding_box.xmin, sink.previous_bounding_box.xmin);
			sink.current_bounding_box.ymin = std::min(sink.current_bounding_box.ymin, sink.previous_bounding_box.ymin);
			sink.current_bounding_box.xmax = std::max(sink.current_bounding_box.xmax, sink.previous_bounding_box.xmax);
			sink.current_bounding_box.ymax = std::max(sink.current_bounding_box.ymax, sink.previous_bounding_box.ymax);
		}
	}
	return num_changed;
}

void update_sink_bounding_boxes_3(net_t &net, route_tree_t &rt, const RRGraph &g, float astar_fac)
{
	/*net.previous_sink_index = net.current_sink_index;*/
	/*net.previous_source = net.current_source;*/
	/*net.previous_bounding_box = get_bounding_box(net.previous_source, net.sinks[net.previous_sink_index]);*/

	/*current_sink.previous_bounding_box = current_sink.current_bounding_box;*/
	assert(!route_tree_empty(rt));

	net.has_sink = false;
	for (auto &sink : net.sinks) {
		if (!net.sink_routed[sink.id]) {
			extern int nx;
			extern int ny;
			bounding_box_t fpga_bounding_box;
			fpga_bounding_box.xmin = 0;
			fpga_bounding_box.xmax = nx+2;
			fpga_bounding_box.ymin = 0;
			fpga_bounding_box.ymax = ny+2;

			auto &sink_vertex = get_vertex(g, sink.rr_node);

			pair<int, int> min_metric = make_pair(std::numeric_limits<int>::max(), std::numeric_limits<int>::max());
			const RRNode *new_from = nullptr;
			bounding_box_t new_bounding_box;
			for_all_vertices(rt.graph, [&min_metric, &new_from, &new_bounding_box, &sink, &g] (const RouteTreeNode &rt_node) -> void {
					if (!rt_node.properties.valid) {
					return;
					}

					const auto &from_node = get_vertex(g, rt_node.properties.rr_node);

					if (from_node.properties.type == IPIN || from_node.properties.type == SINK) {
					return;
					}

					point_t<int> bottom_left, top_right;

					bottom_left.x = std::min(sink.x, from_node.properties.xlow);
					bottom_left.y = std::min(sink.y, from_node.properties.ylow);
					top_right.x = std::max(sink.x, from_node.properties.xhigh);
					top_right.y = std::max(sink.y, from_node.properties.yhigh);

					new_bounding_box = get_bounding_box(bottom_left, top_right, sink.bb_factor);

					int area = get_bounding_box_area(new_bounding_box);
					int slack = from_node.properties.occ - from_node.properties.capacity;
					auto new_metric = make_pair(slack, area);
					if (new_metric < min_metric) {
						new_from = &from_node;
						min_metric = new_metric;
					}
			});

			char buffer[256];

			if (new_from == nullptr) {
				sprintf_rr_node(sink.rr_node, buffer);
				zlog_error(delta_log, "Error: Unable to find node in route tree that is closest to %s\n", buffer);
				assert(false);
			}

			sprintf_rr_node(id(*new_from), buffer);
			zlog_level(delta_log, ROUTER_V2, "New source %s\n", buffer);

			sink.current_bounding_box = new_bounding_box;

			net.has_sink = true;
		}
	}
}

void update_virtual_net_current_sinks(virtual_net_t &virtual_net, const route_tree_t &rt)
{
	virtual_net.current_sinks.clear();
	virtual_net.sink_bounding_box = bg::make_inverse<box>();

	for (const auto &sink : virtual_net.sinks) {
		char buffer[256];
		sprintf_rr_node(sink->rr_node, buffer);

		const RouteTreeNode *sink_rt_node = route_tree_get_rt_node(rt, sink->rr_node);
		if (!sink_rt_node || sink_rt_node->properties.pending_rip_up) {
			zlog_level(scheduler_log, ROUTER_V3, "Adding sink %s to net %d virtual net %d\n", buffer, virtual_net.sinks[0]->net->vpr_id, virtual_net.id);
			virtual_net.current_sinks.push_back(sink);

			bg::expand(virtual_net.sink_bounding_box, point(sink->x, sink->y));
		} else {
			zlog_level(scheduler_log, ROUTER_V3, "NOT adding sink %s to net %d virtual net %d because it's already routed\n", buffer, virtual_net.sinks[0]->net->vpr_id, virtual_net.id);
		}
	}

	bg::centroid(virtual_net.sink_bounding_box, virtual_net.centroid);
	/*std::sort(begin(virtual_net.current_sinks), end(virtual_net.current_sinks), [] (const sink_t *a, const sink_t *b) -> bool {*/
			/*return bg::distance(virtual_net.source->*/
}

void update_virtual_net_scheduler_bounding_box(virtual_net_t &virtual_net, const box &initial_scheduler_bounding_box)
{
	zlog_level(delta_log, ROUTER_V3, "Initial scheduler bounding box %s\n", 
			static_cast<const std::ostringstream &>(std::ostringstream() << bg::dsv(initial_scheduler_bounding_box)).str().c_str());

	virtual_net.scheduler_bounding_box = initial_scheduler_bounding_box;

	bg::expand(virtual_net.scheduler_bounding_box, virtual_net.current_bounding_box);

	zlog_level(delta_log, ROUTER_V3, "Scheduler bounding box after expanding with virtual net bounding box %s: %s\n",
			static_cast<const std::ostringstream &>(std::ostringstream() << bg::dsv(virtual_net.current_bounding_box)).str().c_str(),
			static_cast<const std::ostringstream &>(std::ostringstream() << bg::dsv(virtual_net.scheduler_bounding_box)).str().c_str());
}

void update_virtual_net_bounding_box(virtual_net_t &virtual_net, route_tree_t &rt, const RRGraph &g, float astar_fac, perf_t *perf)
{
	/*net.previous_sink_index = net.current_sink_index;*/
	/*net.previous_source = net.current_source;*/
	/*net.previous_bounding_box = get_bounding_box(net.previous_source, net.sinks[net.previous_sink_index]);*/

	/*current_sink.previous_bounding_box = current_sink.current_bounding_box;*/
	/*assert(!route_tree_empty(rt));*/

	int min_area = std::numeric_limits<int>::max();
	const RRNode *new_from = nullptr;

	char buffer[256];

	auto centroid_start = std::chrono::high_resolution_clock::now();

	if (perf) {
		perf->total_centroid_time += std::chrono::high_resolution_clock::now()-centroid_start;
	}

	net_t &net = *virtual_net.current_sinks[0]->net;

	/*zlog_level(delta_log, ROUTER_V3, "Getting nearest node from route tree to centroid %d,%d\n", c.get<0>(), c.get<1>());  */

	auto get_nearest_start = std::chrono::high_resolution_clock::now();

	int num_iters = 0;
	RouteTreeNode *nearest_rt_node = route_tree_get_nearest_node(rt, virtual_net.centroid, g, &num_iters);
	net.num_nearest_iters += num_iters;
	net.total_point_tree_size += rt.point_tree.size();

	if (perf) {
		perf->total_get_nearest_time += std::chrono::high_resolution_clock::now()-get_nearest_start;
	}

	auto verification_start = std::chrono::high_resolution_clock::now();

	if (nearest_rt_node == nullptr) {
		const auto &nodes = route_tree_get_nodes(rt);
		if (route_tree_empty(rt) || all_of(begin(nodes), end(nodes), [] (const RouteTreeNode &node) -> bool { return node.properties.pending_rip_up; })) {
			virtual_net.nearest_rr_node = net.source.rr_node;
			zlog_level(delta_log, ROUTER_V3, "Net %d empty route tree. Setting virtual net %d bounding box to start from source\n", net.vpr_id, virtual_net.id);
		} else {
			zlog_error(delta_log, "Error: Unable to find node in route tree that is closest to virtual net %d centroid %d,%d\n", virtual_net.id, virtual_net.centroid.get<0>(), virtual_net.centroid.get<1>());
			assert(false);
		}
	} else {
		virtual_net.nearest_rr_node = nearest_rt_node->properties.rr_node;
	}

	if (perf) {
		perf->total_verificaion_time += std::chrono::high_resolution_clock::now()-verification_start;
	}

	/*sprintf_rr_node(virtual_net.nearest_rr_node, buffer);*/
	/*zlog_level(delta_log, ROUTER_V3, "Net %d virtual net %d new source %s\n", net.vpr_id, virtual_net.id, buffer);*/

	auto expansion_start = std::chrono::high_resolution_clock::now();

	virtual_net.current_bounding_box = bg::make_inverse<box>();
	auto &rr_node = get_vertex(g, virtual_net.nearest_rr_node);
	bg::expand(virtual_net.current_bounding_box, segment(point(rr_node.properties.xlow, rr_node.properties.ylow), point(rr_node.properties.xhigh, rr_node.properties.yhigh)));
	bg::expand(virtual_net.current_bounding_box, virtual_net.sink_bounding_box);

	int bb_factor = 0;
	for (const auto &sink : virtual_net.current_sinks) {
		bb_factor += sink->bb_factor;
	}
	bb_factor /= virtual_net.current_sinks.size();

	/* expand bounding box in VPR style */
	int xmin = virtual_net.current_bounding_box.min_corner().get<0>();
	int ymin = virtual_net.current_bounding_box.min_corner().get<1>();
	int xmax = virtual_net.current_bounding_box.max_corner().get<0>();
	int ymax = virtual_net.current_bounding_box.max_corner().get<1>();

	extern int nx;
	extern int ny;

	xmin = std::max(xmin - 1 - bb_factor, 0);
	xmax = std::min(xmax + bb_factor, nx + 1);
	ymin = std::max(ymin - 1 - bb_factor, 0);
	ymax = std::min(ymax + bb_factor, ny + 1);

	virtual_net.current_bounding_box.min_corner().set<0>(xmin);
	virtual_net.current_bounding_box.min_corner().set<1>(ymin);
	virtual_net.current_bounding_box.max_corner().set<0>(xmax);
	virtual_net.current_bounding_box.max_corner().set<1>(ymax);

	/* assign all sinks to have the same bounding box as the virtual net */
	for (const auto &sink : virtual_net.current_sinks) {
		sink->current_bounding_box.xmin = virtual_net.current_bounding_box.min_corner().get<0>();
		sink->current_bounding_box.ymin = virtual_net.current_bounding_box.min_corner().get<1>();
		sink->current_bounding_box.xmax = virtual_net.current_bounding_box.max_corner().get<0>();
		sink->current_bounding_box.ymax = virtual_net.current_bounding_box.max_corner().get<1>();
	}

	if (perf) {
		perf->total_expansion_time += std::chrono::high_resolution_clock::now()-expansion_start;
	}

	zlog_level(delta_log, ROUTER_V3, "Net %d virtual net %d sink bounding box %s current bounding box %s\n", net.vpr_id, virtual_net.id, 
			static_cast<const std::ostringstream &>(std::ostringstream() << bg::dsv(virtual_net.sink_bounding_box)).str().c_str(),
			static_cast<const std::ostringstream &>(std::ostringstream() << bg::dsv(virtual_net.current_bounding_box)).str().c_str());
}

void update_sink_bounding_boxes_5(vector<virtual_net_t> &virtual_nets, route_tree_t &rt, const RRGraph &g, float astar_fac)
{
	/*net.previous_sink_index = net.current_sink_index;*/
	/*net.previous_source = net.current_source;*/
	/*net.previous_bounding_box = get_bounding_box(net.previous_source, net.sinks[net.previous_sink_index]);*/

	/*current_sink.previous_bounding_box = current_sink.current_bounding_box;*/
	/*assert(!route_tree_empty(rt));*/

	for (auto &vnet : virtual_nets) {
		if (!vnet.routed) {
			int min_area = std::numeric_limits<int>::max();
			const RRNode *new_from = nullptr;

			const net_t &net = *vnet.sinks[0]->net;

			point centroid;
			zlog_level(scheduler_log, ROUTER_V2, "Getting nearest node from route tree to centroid %d,%d\n", centroid.get<0>(), centroid.get<1>());  
			RouteTreeNode *nearest_rt_node = route_tree_get_nearest_node(rt, centroid, g, nullptr);

			char buffer[256];

			if (nearest_rt_node == nullptr) {
				const auto &nodes = route_tree_get_nodes(rt);
				if (route_tree_empty(rt) || all_of(begin(nodes), end(nodes), [] (const RouteTreeNode &node) -> bool { return node.properties.pending_rip_up; })) {
					vnet.nearest_rr_node = net.source.rr_node;
					zlog_level(scheduler_log, ROUTER_V2, "Net %d empty route tree. Setting virtual net %d bounding box to start from source\n", net.vpr_id, vnet.id);
				} else {
					zlog_error(scheduler_log, "Error: Unable to find node in route tree that is closest to virtual net %d centroid %d,%d\n", vnet.id, centroid.get<0>(), centroid.get<1>());					assert(false);
				}
			} else {
				vnet.nearest_rr_node = nearest_rt_node->properties.rr_node;
			}

			sprintf_rr_node(vnet.nearest_rr_node, buffer);
			zlog_level(scheduler_log, ROUTER_V2, "Net %d virtual net %d new source %s\n", net.vpr_id, vnet.id, buffer);
		} else {
			zlog_level(scheduler_log, ROUTER_V2, "Invalid virtual net %d of net %d\n", vnet.id, vnet.sinks[0]->net->vpr_id);
		}
	}
}

void update_sink_bounding_boxes_4(net_t &net, route_tree_t &rt, const RRGraph &g, float astar_fac)
{
	/*net.previous_sink_index = net.current_sink_index;*/
	/*net.previous_source = net.current_source;*/
	/*net.previous_bounding_box = get_bounding_box(net.previous_source, net.sinks[net.previous_sink_index]);*/

	/*current_sink.previous_bounding_box = current_sink.current_bounding_box;*/
	/*assert(!route_tree_empty(rt));*/

	net.has_sink = false;
	for (auto &sink : net.sinks) {
		if (!net.sink_routed[sink.id]) {
			int min_area = std::numeric_limits<int>::max();
			const RRNode *new_from = nullptr;

			sink.scheduler_bounding_box.xmin = std::numeric_limits<int>::max();
			sink.scheduler_bounding_box.xmax = std::numeric_limits<int>::min();
			sink.scheduler_bounding_box.ymin = std::numeric_limits<int>::max();
			sink.scheduler_bounding_box.ymax = std::numeric_limits<int>::min();

			for (const auto &rt_node : route_tree_get_nodes(rt)) {
				const auto &from_node = get_vertex(g, rt_node.properties.rr_node);

				char s_source[256];
				sprintf_rr_node(id(from_node), s_source);

				if (rt_node.properties.pending_rip_up) {
					sink.scheduler_bounding_box.xmin = std::min(sink.scheduler_bounding_box.xmin, from_node.properties.xlow);
					sink.scheduler_bounding_box.xmax = std::max(sink.scheduler_bounding_box.xmax, from_node.properties.xhigh);
					sink.scheduler_bounding_box.ymin = std::min(sink.scheduler_bounding_box.ymin, from_node.properties.ylow);
					sink.scheduler_bounding_box.ymax = std::max(sink.scheduler_bounding_box.ymax, from_node.properties.yhigh);
					zlog_level(scheduler_log, ROUTER_V3, "Route tree node %s pending rip up. Skipping.\n", s_source);
					continue;
				}

				if (from_node.properties.type == IPIN || from_node.properties.type == SINK) {
					zlog_level(scheduler_log, ROUTER_V3, "Route tree node %s is IPIN/SINK. Skipping.\n", s_source);
					continue;
				}

				point_t<int> bottom_left, top_right;

				bottom_left.x = std::min(sink.x, from_node.properties.xlow);
				bottom_left.y = std::min(sink.y, from_node.properties.ylow);
				top_right.x = std::max(sink.x, from_node.properties.xhigh);
				top_right.y = std::max(sink.y, from_node.properties.yhigh);

				const auto &bounding_box = get_bounding_box(bottom_left, top_right, sink.bb_factor);

				char s_sink[256];
				sprintf_rr_node(sink.rr_node, s_sink);
				zlog_level(scheduler_log, ROUTER_V3, "%s -> %s %d-%d %d-%d (bb_factor %d) ", s_source, s_sink, bounding_box.xmin, bounding_box.xmax, bounding_box.ymin, bounding_box.ymax, sink.bb_factor);

				int area = get_bounding_box_area(bounding_box);
				if (area < min_area) {
					new_from = &from_node;
					sink.current_bounding_box = bounding_box;
					min_area = area;
					zlog_level(scheduler_log, ROUTER_V3, "updated min_area to %d", min_area);

					sink.scheduler_bounding_box.xmin = std::min(sink.scheduler_bounding_box.xmin, bounding_box.xmin);
					sink.scheduler_bounding_box.xmax = std::max(sink.scheduler_bounding_box.xmax, bounding_box.xmax);
					sink.scheduler_bounding_box.ymin = std::min(sink.scheduler_bounding_box.ymin, bounding_box.ymin);
					sink.scheduler_bounding_box.ymax = std::max(sink.scheduler_bounding_box.ymax, bounding_box.ymax);
				}

				zlog_level(scheduler_log, ROUTER_V3, "\n");
			}

			char buffer[256];

			if (new_from == nullptr) {
				const auto &nodes = route_tree_get_nodes(rt);
				if (route_tree_empty(rt) || all_of(begin(nodes), end(nodes), [] (const RouteTreeNode &node) -> bool { return node.properties.pending_rip_up; })) {
					new_from = &get_vertex(g, net.source.rr_node);
					sink.current_bounding_box = get_bounding_box(net.source, sink, sink.bb_factor);
					zlog_level(scheduler_log, ROUTER_V2, "Net %d empty route tree. Setting sink %d bounding box to start from source\n", net.vpr_id, sink.id);
				} else {
					sprintf_rr_node(sink.rr_node, buffer);
					zlog_error(scheduler_log, "Error: Unable to find node in route tree that is closest to %s\n", buffer);
					assert(false);
				}
			}

			sprintf_rr_node(id(*new_from), buffer);
			zlog_level(scheduler_log, ROUTER_V2, "Net %d sink %d new source %s %d-%d %d-%d\n", net.vpr_id, sink.id, buffer, sink.current_bounding_box.xmin, sink.current_bounding_box.xmax, sink.current_bounding_box.ymin, sink.current_bounding_box.ymax);

			net.has_sink = true;
		}
	}
}

void update_sink_bounding_boxes_2(net_t &net, route_tree_t &rt, const RRGraph &g, float astar_fac)
{
	/*net.previous_sink_index = net.current_sink_index;*/
	/*net.previous_source = net.current_source;*/
	/*net.previous_bounding_box = get_bounding_box(net.previous_source, net.sinks[net.previous_sink_index]);*/

	/*current_sink.previous_bounding_box = current_sink.current_bounding_box;*/
	assert(!route_tree_empty(rt));

	net.has_sink = false;
	for (auto &sink : net.sinks) {
		if (!net.sink_routed[sink.id]) {
			extern int nx;
			extern int ny;
			bounding_box_t fpga_bounding_box;
			fpga_bounding_box.xmin = 0;
			fpga_bounding_box.xmax = nx+2;
			fpga_bounding_box.ymin = 0;
			fpga_bounding_box.ymax = ny+2;

			auto &sink_vertex = get_vertex(g, sink.rr_node);

			int min_area = std::numeric_limits<int>::max();
			const RRNode *new_from = nullptr;
			for_all_vertices(rt.graph, [&min_area, &new_from, &sink, &g] (const RouteTreeNode &rt_node) -> void {
					if (!rt_node.properties.valid/* || rt_node.properties.owner == sink.rr_node*/) {
					return;
					}

					const auto &from_node = get_vertex(g, rt_node.properties.rr_node);

					if (from_node.properties.type == IPIN || from_node.properties.type == SINK) {
					return;
					}

					point_t<int> bottom_left, top_right;

					bottom_left.x = std::min(sink.x, from_node.properties.xlow);
					bottom_left.y = std::min(sink.y, from_node.properties.ylow);
					top_right.x = std::max(sink.x, from_node.properties.xhigh);
					top_right.y = std::max(sink.y, from_node.properties.yhigh);

					const auto &bounding_box = get_bounding_box(bottom_left, top_right, sink.bb_factor);

					char s_sink[256];
					char s_source[256];
					sprintf_rr_node(id(from_node), s_source);
					sprintf_rr_node(sink.rr_node, s_sink);
					zlog_level(scheduler_log, ROUTER_V3, "%s -> %s %d-%d %d-%d (bb_factor %d) ", s_source, s_sink, bounding_box.xmin, bounding_box.xmax, bounding_box.ymin, bounding_box.ymax, sink.bb_factor);

					int area = get_bounding_box_area(bounding_box);
					if (area < min_area) {
						new_from = &from_node;
						sink.current_bounding_box = bounding_box;
						min_area = area;
						zlog_level(scheduler_log, ROUTER_V3, "updated min_area to %d", min_area);
					}

					zlog_level(scheduler_log, ROUTER_V3, "\n");
			});

			char buffer[256];

			if (new_from == nullptr) {
				sprintf_rr_node(sink.rr_node, buffer);
				zlog_error(delta_log, "Error: Unable to find node in route tree that is closest to %s\n", buffer);
				assert(false);
			}

			sprintf_rr_node(id(*new_from), buffer);
			zlog_level(scheduler_log, ROUTER_V2, "Net %d sink %d new source %s %d-%d %d-%d\n", net.vpr_id, sink.id, buffer, sink.current_bounding_box.xmin, sink.current_bounding_box.xmax, sink.current_bounding_box.ymin, sink.current_bounding_box.ymax);

			net.has_sink = true;
		}
	}
}

void update_sink_bounding_boxes(net_t &net, route_tree_t &rt, const RRGraph &g, float astar_fac, perf_t &perf)
{
	/*net.previous_sink_index = net.current_sink_index;*/
	/*net.previous_source = net.current_source;*/
	/*net.previous_bounding_box = get_bounding_box(net.previous_source, net.sinks[net.previous_sink_index]);*/

	/*current_sink.previous_bounding_box = current_sink.current_bounding_box;*/

	assert(!route_tree_empty(rt));

	net.has_sink = false;
	for (auto &sink : net.sinks) {
		if (!net.sink_routed[sink.id]) {
			std::priority_queue<route_state_t> heap;
			extern int nx;
			extern int ny;
			bounding_box_t fpga_bounding_box;
			fpga_bounding_box.xmin = 0;
			fpga_bounding_box.xmax = nx+2;
			fpga_bounding_box.ymin = 0;
			fpga_bounding_box.ymax = ny+2;
			route_tree_add_to_heap(rt, g, get_vertex(g, sink.rr_node), sink.criticality_fac, astar_fac, fpga_bounding_box, heap, &perf);
			const RRNode *current_rr_node = &get_vertex(g, heap.top().rr_node);
			while (current_rr_node->properties.type != CHANX && current_rr_node->properties.type != CHANY) {
				heap.pop();
				current_rr_node = &get_vertex(g, heap.top().rr_node);
			}
			sink.source.rr_node = id(*current_rr_node);
			if (current_rr_node->properties.inc_direction) {
				sink.source.x = current_rr_node->properties.xlow; 
				sink.source.y = current_rr_node->properties.ylow;  
			} else {
				sink.source.x = current_rr_node->properties.xhigh; 
				sink.source.y = current_rr_node->properties.yhigh;  
			}
			/*sink.previous_bounding_box = sink.current_bounding_box;*/
			sink.current_bounding_box = get_bounding_box(sink.source, sink, sink.bb_factor);
			char source_s[256];
			char sink_s[256];
			sprintf_rr_node(sink.rr_node, sink_s);
			sprintf_rr_node(sink.source.rr_node, source_s);
			zlog_level(delta_log, ROUTER_V2, "Updating source of sink %s to %s. New bouding box %d-%d %d-%d\n", sink_s, source_s,
					sink.current_bounding_box.xmin, sink.current_bounding_box.xmax, sink.current_bounding_box.ymin, sink.current_bounding_box.ymax);
			net.has_sink = true;
		}
	}
}

void set_to_next_sink(net_t &net, route_tree_t &rt, const RRGraph &g, float astar_fac, perf_t &perf)
{
	/*net.previous_sink_index = net.current_sink_index;*/
	/*net.previous_source = net.current_source;*/
	/*net.previous_bounding_box = get_bounding_box(net.previous_source, net.sinks[net.previous_sink_index]);*/

	assert(net.current_sink_index >= 0 && net.current_sink_index < net.sinks.size());

	auto &old_sink = net.sinks[net.current_sink_index];
	old_sink.previous_bounding_box = old_sink.current_bounding_box;

	net.current_sink_index++;
	if (net.current_sink_index >= net.sinks.size()) {
		net.has_sink = false;
	} else {
		net.has_sink = true;
	}

	if (!net.has_sink) {
		return;
	}

	zlog_level(delta_log, ROUTER_V2, "Set to next sink for net %d\n", net.vpr_id);

	assert(!route_tree_empty(rt));
	
	std::priority_queue<route_state_t> heap;
	extern int nx;
	extern int ny;
	bounding_box_t fpga_bounding_box;
	fpga_bounding_box.xmin = 0;
	fpga_bounding_box.xmax = nx+2;
	fpga_bounding_box.ymin = 0;
	fpga_bounding_box.ymax = ny+2;
	auto &new_sink = net.sinks[net.current_sink_index];
	/*route_tree_add_to_heap(rt, g, get_vertex(g, new_sink.rr_node), new_sink.criticality_fac, astar_fac, fpga_bounding_box, heap, &perf);*/
	/*const RRNode *current_rr_node = &get_vertex(g, heap.top().rr_node);*/
	/*zlog_info(delta_log, "Set to next sink for net %d\n", net.vpr_id);*/
	/*while (current_rr_node->properties.type != CHANX && current_rr_node->properties.type != CHANY) {*/
		/*heap.pop();*/
		/*current_rr_node = &get_vertex(g, heap.top().rr_node);*/
		/*zlog_info(delta_log, "Current: %d", net.vpr_id);*/
	/*}*/

	auto &sink_vertex = get_vertex(g, new_sink.rr_node);
	/*if (sink_vertex.properties.xlow != sink_vertex.properties.xhigh || sink_vertex.properties.ylow != sink_vertex.properties.yhigh) {*/

		/*assert(false);*/
	/*}*/

	int min_area = std::numeric_limits<int>::max();
	const RRNode *new_from = nullptr;
	bounding_box_t new_bounding_box;
	for_all_vertices(rt.graph, [&min_area, &new_from, &new_bounding_box, &new_sink, &g] (const RouteTreeNode &rt_node) -> void {
			if (!rt_node.properties.valid) {
				return;
			}

			const auto &from_node = get_vertex(g, rt_node.properties.rr_node);

			if (from_node.properties.type == IPIN || from_node.properties.type == SINK) {
				return;
			}

			point_t<int> bottom_left, top_right;
			
			bottom_left.x = std::min(new_sink.x, from_node.properties.xlow);
			bottom_left.y = std::min(new_sink.y, from_node.properties.ylow);
			top_right.x = std::max(new_sink.x, from_node.properties.xhigh);
			top_right.y = std::max(new_sink.y, from_node.properties.yhigh);

			new_bounding_box = get_bounding_box(bottom_left, top_right, new_sink.bb_factor);

			int area = get_bounding_box_area(new_bounding_box);
			if (area < min_area) {
				new_from = &from_node;
				min_area = area;
			}
			});

	char buffer[256];

	if (new_from == nullptr) {
		sprintf_rr_node(new_sink.rr_node, buffer);
		zlog_error(delta_log, "Error: Unable to find node in route tree that is closest to %s\n", buffer);
		assert(false);
	}

	sprintf_rr_node(id(*new_from), buffer);
	zlog_level(delta_log, ROUTER_V2, "New source %s\n", buffer);
	/*source_t new_source;*/
	/*int min_distance = std::numeric_limits<int>::max();*/
	/*for_all_vertices(rt.graph, [] (const RouteTreeNode &v) -> void {*/
			/*if (v.properties.valid) {*/
			/*v.properties.xlow*/
			/*}*/
			/*});*/
	/*const sink_t &target = net.sinks[net.current_sink_index];*/
	/*route_tree_add_to_heap(rt, g, get_vertex(g, target.rr_node), target.criticality_fac, astar_fac, fpga_bounding_box, heap, NULL);*/
	/*net.current_source.rr_node = heap.top().rr_node;*/
	/*assert(get_vertex(g, net.current_source.rr_node).properties.type == SOURCE);*/
	/*extern struct s_net *clb_net;*/
	/*extern struct s_block *block;*/
	/*int b = clb_net[net.vpr_id].node_block[0];*/
	/*int p = clb_net[net.vpr_id].node_block_pin[0];*/
	/* TODO: fix the code below, the calculation of coordinates is wrong */
	/*net.current_source.x = block[b].x;*/
	/*net.current_source.y = block[b].y + block[b].type->pin_height[p];*/

	new_sink.current_bounding_box = new_bounding_box;

	/*adjust_bounding_box(net);*/
}

void reset_current_source_sink_3(net_t &net)
{
	for (auto &sink : net.sinks) {
		sink.source = net.source;
		sink.current_bounding_box = get_bounding_box(sink.source, sink, sink.bb_factor);
	}
}

void reset_current_source_sink_2(net_t &net)
{
	net.current_sink = nullptr;

	for (auto &sink : net.sinks) {
		sink.source = net.source;
		sink.current_bounding_box = get_bounding_box(sink.source, sink, sink.bb_factor);
	}

	fill(net.sink_routed.begin(), net.sink_routed.end(), false);
	net.has_sink = true;
}

void reset_current_source_sink(net_t &net)
{
	/*net.previous_source_valid = false;*/
	/*net.previous_sink_valid = false;*/
	/*net.previous_bounding_box_valid = false;*/

	net.current_source = net.source;
	/*net.previous_source.rr_node = -1;*/
	/*net.previous_source.x = -1;*/
	/*net.previous_source.y = -1;*/

	net.current_sink_index = 0;
	/*net.previous_sink_index = -1;*/
	assert(net.sinks.size() > 0);
	net.has_sink = true;

	auto &sink = net.sinks[net.current_sink_index];
	char buffer[256];
	sprintf_rr_node(sink.rr_node, buffer);
	zlog_level(delta_log, ROUTER_V2, "Resetting %s bounding box with bb_factor %d\n", buffer, sink.bb_factor);

	sink.current_bounding_box = get_bounding_box(net.current_source, sink, sink.bb_factor);
}

/*void temp_route_nets(RRGraph &g, vector<net_t *> &nets, const route_parameters_t &params, route_tree_t *route_trees, t_net_timing *net_timing)*/
/*{*/
	/*[> rip up previous routing <]*/
	/*if (route_tree_initiated(rt)) {*/
		/*update_one_cost(g, rt, route_tree_get_rt_node(rt, net.source.rr_node), -1, params.pres_fac);*/
	/*} */

	/*route_tree_clear(rt);*/
/*}*/

/*bool route_tree_has_path(const route_tree_t &rt, int sink_rr_node)*/
/*{*/
	/*return route_tree_get_path_to_sink(rt, sink_rr_node) ? true : false;*/
/*}*/

void test_graph()
{
	graph_t<int, int> g;
	/*add_vertex(g, 5);*/
	/*add_edge(g, get_vertex(g, 0), get_vertex(g, 1));*/
	/*add_edge(g, get_vertex(g, 0), get_vertex(g, 2));*/
	/*add_edge(g, get_vertex(g, 1), get_vertex(g, 3));*/
	/*add_edge(g, get_vertex(g, 1), get_vertex(g, 4));*/
	/*remove_edge(g, get_edge(g, 0));*/

	const int num_iters = 10000;
	add_vertex(g, num_iters);
	for (int i = 0; i < num_iters; ++i) {
		get_vertex(g, i).properties = i;
	}
	/*for (int i = 0; i < num_iters; ++i) {*/
		/*add_edge(g, get_vertex(g, 0), get_vertex(g, 1)).properties = i;*/
	/*}*/

	unsigned long long sum1 = 0;
	auto time = std::chrono::high_resolution_clock::now();

	for (auto &v : get_vertices(g)) {
		/*printf("%d\n", v.properties);*/
		v.properties = 1;
		sum1 += v.properties;
	}

	auto elapsed = std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::high_resolution_clock::now()-time);
	printf("new traversing vertices took %lld\n", elapsed.count());

	unsigned long long sum2 = 0;
	time = std::chrono::high_resolution_clock::now();

	for_all_vertices(g, [&sum2] (const vertex_t<int, int> &v) -> void {
		sum2 += v.properties;
	});

	elapsed = std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::high_resolution_clock::now()-time);
	printf("old traversing vertices took %lld\n", elapsed.count());

	assert(sum1 == sum2);
}

void test_rtree();

void sort_sinks(vector<net_t> &nets)
{
	for (auto &net : nets) {
		std::sort(net.sinks.begin(), net.sinks.end(), [&net] (const sink_t &a, const sink_t &b) -> bool {
				int a_to_src_dist = abs(a.x - net.source.x) + abs(a.y - net.source.y);
				int b_to_src_dist = abs(b.x - net.source.x) + abs(b.y - net.source.y);
				return a_to_src_dist > b_to_src_dist;
				});
	}
}

bool partitioning_route_bounding_box(t_router_opts *opts)
{
	const int max_iter = 50;

	init_logging();

	/*test_scheduler();*/

	/*test_rtree();*/

	/*test_graph();*/

	RRGraph g;
	init_graph(g);

	vector<net_t> nets;
	vector<net_t> global_nets;
	init_nets(nets, global_nets, opts->bb_factor);

	vector<net_t *> nets_ptr(nets.size());
	for (int i = 0; i < nets.size(); ++i) {
		nets_ptr[i] = &nets[i];
	}
	extern s_bb *route_bb;
	sort(nets_ptr.begin(), nets_ptr.end(), [] (const net_t *a, const net_t *b) -> bool {
			return get_bounding_box_area(route_bb[a->vpr_id]) > get_bounding_box_area(route_bb[b->vpr_id]);
			});
	int rank = 0;
	for (auto &net : nets_ptr) {
		net->bb_area_rank = rank++;
	}

	extern int nx;
	extern int ny;
	zlog_info(delta_log, "Num nets: %lu Num global nets: %lu nx: %d ny: %d\n", nets.size(), global_nets.size(), nx, ny);
	
	t_net_timing *net_timing = new t_net_timing[nets.size()+global_nets.size()];
	init_net_timing(nets, global_nets, net_timing);

	vector<route_tree_t> route_trees(nets.size());
	for (int i = 0; i < nets.size(); ++i) {
		route_tree_init(route_trees[i]);
	}

	vector<trace_t> traces(nets.size());
	for (int i = 0; i < nets.size(); ++i) {
		trace_init(traces[i]);
	}

	route_state_t *state = new route_state_t[num_vertices(g)];

	for (int i = 0; i < num_vertices(g); ++i) {
		state[i].rr_node = -1;
		state[i].known_cost = std::numeric_limits<float>::max();
		state[i].cost = std::numeric_limits<float>::max();
		state[i].prev_edge = nullptr;
		state[i].upstream_R = -1;
		state[i].delay = std::numeric_limits<float>::max();
	}


	/*route_tree_t **route_trees = new route_tree_t *[nets.size()+global_nets.size()];*/
	/*for (int i = 0; i < nets.size(); ++i) {*/
		/*route_trees[nets[i].id] = new route_tree_t[nets[i].sinks.size()+1];*/
		/*for (int j = 0; j <= nets[i].sinks.size(); ++j) {*/
			/*route_tree_init(route_trees[nets[i].id][j]);*/
		/*}*/
	/*}*/
	/*for (int i = 0; i < global_nets.size(); ++i) {*/
		/*route_trees[nets[i].id] = new route_tree_t[nets[i].sinks.size()+1];*/
		/*for (int j = 0; j <= nets[i].sinks.size(); ++j) {*/
			/*route_tree_init(route_trees[nets[i].id][j]);*/
		/*}*/
	/*}*/

	route_parameters_t params;
	params.criticality_exp = opts->criticality_exp;
	params.astar_fac = opts->astar_fac;
	params.max_criticality = opts->max_criticality;
	params.bend_cost = opts->bend_cost;
	params.pres_fac = opts->first_iter_pres_fac; /* Typically 0 -> ignore cong. */

	for_all_vertices(g, [] (RRNode &v) -> void {
			v.properties.acc_cost = 1;
			v.properties.pres_cost = 1;
			v.properties.occ = 0;
			});

	bool routed = false;
	for (int iter = 0; iter < max_iter && !routed; ++iter) {
		char iter_s[64];
		sprintf(iter_s, "%d", iter);
		assert(!zlog_put_mdc("iter", iter_s));

		zlog_info(delta_log, "Routing iteration: %d\n", iter);

		cpu_timer iter_timer;
		iter_timer.start();

		for (int i = 0; i < nets.size(); ++i) {
			route_tree_clear(route_trees[i]);
		}

		for (auto &net : nets) {
			update_sink_criticalities(net, net_timing[net.vpr_id], params);
		}

		cpu_timer schedule_timer;
		schedule_timer.start();

		vector<net_t *> nets_to_route;
		int current_local_id = 0;
		for (auto &net : nets) {
			if (net.has_sink) {
				nets_to_route.push_back(&net);
				net.current_local_id = current_local_id++;
			}
		}
		if (iter > 0) {
			for (auto &net : nets_to_route) {
				adjust_bounding_box_2(*net);
			}
		}
		vector<pair<sink_t *, net_t *>> scheduled_sinks;
		schedule_nets_bounding_box(nets_to_route, scheduled_sinks);
		zlog_info(delta_log, "Num scheduled sink: %lu\n", scheduled_sinks.size());

		schedule_timer.stop();
		nanosecond_type schedule_time = schedule_timer.elapsed().wall;

		int time = 0;

		perf_t perf;
		perf.num_heap_pushes = 0;

		nanosecond_type route_time = 0;

		int total_num_sinks = 0;

		while (scheduled_sinks.size() > 0) {
			zlog_info(delta_log, "Time: %d\n", time);
			cpu_timer route_timer;
			route_timer.start();

			/* this for loop can be parallelized with parallel_for */
			for (const auto &sink_pair : scheduled_sinks) {
				net_t *net;
				sink_t *sink;
				std::tie(sink, net) = sink_pair;
				net_t virtual_net = *net;
				virtual_net.sinks.clear();
				virtual_net.sinks.push_back(*sink);
				virtual_net.current_source = sink->source;
				/*virtual_net.current_bounding_box = sink->current_bounding_box;*/
				/*if (time == 0) {*/
				/*trace_rip_up_net(traces[net->local_id], g, params.pres_fac);*/
				/*}*/
				trace_rip_up_segment(traces[net->local_id], g, sink->rr_node, params.pres_fac);
				/*printf("Net: %d\n", virtual_net.vpr_id);*/
				route_net(g, virtual_net, params, state, route_trees[net->local_id], traces[net->local_id], net_timing[net->vpr_id], &perf);

				net->sink_routed[sink->id] = true;
				/* save this so that we can use it in next iteration */
				sink->previous_bounding_box = sink->current_bounding_box;
				/*route_tree_merge(route_trees[net->id][0], route_trees[net->id][virtual_net.sinks[0].id]);*/
				++total_num_sinks;
			}

			/* TODO: fix the loop range and update_sink_bounding_boxes */
			for (auto &sink_pair : scheduled_sinks) {
				net_t *net = sink_pair.second;
				update_sink_bounding_boxes(*net, route_trees[net->local_id], g, params.astar_fac, perf);
			}

			route_timer.stop();
			cpu_times elapsed = route_timer.elapsed();
			route_time += elapsed.wall;
			/*zlog_info(delta_log, "Scheduling took %s\n", format(elapsed).c_str());*/
			zlog_info(delta_log, "Routing took %g\n", elapsed.wall / 1e9);

			schedule_timer.start();
		
			nets_to_route.clear();
			int current_local_id = 0;
			for (auto &net : nets) {
				if (net.has_sink) {
					nets_to_route.push_back(&net);
					net.current_local_id = current_local_id++;
				}
			}
			scheduled_sinks.clear();
			if (iter > 0) {
				for (auto &net : nets_to_route) {
					adjust_bounding_box_2(*net);
				}
			}
			schedule_nets_bounding_box(nets_to_route, scheduled_sinks);

			schedule_timer.stop();
			schedule_time += schedule_timer.elapsed().wall;

			zlog_info(delta_log, "Num scheduled sink: %lu\n", scheduled_sinks.size());

			++time;
		}

		zlog_info(delta_log, "Average concurrency: %g Num heap pushes: %d\n", (float)total_num_sinks/time, perf.num_heap_pushes);

		bool fail = false;
		for (const auto &net : nets) {
			for (int i = 0; i < net.sinks.size(); ++i) {
				if (!net.sink_routed[net.sinks[i].id]) {
					char buffer[256];
					sprintf_rr_node(net.sinks[i].rr_node, buffer);
					zlog_error(delta_log, "Error: %s of net %d is not routed\n", buffer, net.vpr_id);
					fail = true;
				}
			}
		}
		assert(!fail);

		for_all_vertices(g, [] (RRNode &v) -> void { v.properties.recalc_occ = 0; });

		for (const auto &net : nets) {
			check_route_tree(route_trees[net.local_id], net, g);
		}

		bool valid = true;
		for_all_vertices(g, [&valid] (RRNode &v) -> void {
				char buffer[256];
				sprintf_rr_node(id(v), buffer);
				if (v.properties.recalc_occ != v.properties.occ) {
					zlog_error(delta_log, "Node %s occ mismatch, recalc: %d original: %d\n", buffer, v.properties.recalc_occ, v.properties.occ);
					valid = false;
				}
				});
		assert(valid);

		for (const auto &net : nets) {
			vector<int> overused_rr_node;
			get_overused_nodes(route_trees[net.local_id], get_vertex(route_trees[net.local_id].graph, route_trees[net.local_id].root_rt_node_id), g, overused_rr_node);
			if (!overused_rr_node.empty()) {
				zlog_level(delta_log, ROUTER_V3, "Net %d bb_rank %d overused nodes:\n", net.vpr_id, net.bb_area_rank);
				for (const auto &item : overused_rr_node) {
					zlog_level(delta_log, ROUTER_V3, "%d ", item);
				}
				zlog_level(delta_log, ROUTER_V3, "\n");
			}
		}

		if (feasible_routing(g)) {
			zlog_info(delta_log, "Routed in %d iterations\n", iter+1);
			routed = true;
		} else {
			int num_overused_nodes = 0;
			for_all_vertices(g, [&num_overused_nodes] (const RRNode &v) -> void {
					if (v.properties.occ > v.properties.capacity) {
					++num_overused_nodes;
					}
					});
			zlog_info(delta_log, "Num overused nodes: %d/%d (%.2f)\n", num_overused_nodes, num_vertices(g), num_overused_nodes*100.0/num_vertices(g));

			if (iter == 0) {
				params.pres_fac = opts->initial_pres_fac;
				update_costs(g, params.pres_fac, 0);
			} else {
				params.pres_fac *= opts->pres_fac_mult;

				/* Avoid overflow for high iteration counts, even if acc_cost is big */
				params.pres_fac = std::min(params.pres_fac, static_cast<float>(HUGE_POSITIVE_FLOAT / 1e5));

				update_costs(g, params.pres_fac, opts->acc_fac);
			}

			for (auto &net : nets) { 
				reset_current_source_sink_2(net);
			}
		}

		analyze_timing(net_timing);

		iter_timer.stop();
		cpu_times elapsed = iter_timer.elapsed();

		zlog_info(delta_log, "Iteration time: %g s. Route time: %g s. Schedule time: %g s.\n", elapsed.wall / 1e9, route_time / 1e9, schedule_time / 1e9);
	}
	return routed;
}

void dump_congestion_map(const RRGraph &g, const char *filename)
{
	FILE *file = fopen(filename, "w");

	for_all_vertices(g, [&file] (const RRNode &v) -> void {
			auto &prop = v.properties;
			if (prop.occ > prop.capacity) {
			for (int x = prop.xlow; x <= prop.xhigh; ++x) {
			for (int y = prop.ylow; y <= prop.yhigh; ++y) {
				for (int i = 0; i < prop.occ-prop.capacity; ++i) {
					fprintf(file, "%d %d\n", x, y);
				}
			}
			}
			}
			});

	fclose(file);
}

void print_parameters(t_router_opts *opts)
{

}	

void dump_route(const vector<trace_t> &traces, const char *filename)
{
	FILE *routes = fopen(filename, "w");
	assert(routes);

	int inet = 0;
	char buffer[256];
	for (const auto &trace : traces) {
		fprintf(routes, "Net %d\n", inet);
		for (const auto &item : trace.segments) {
			int sink_rr_node = item.first;
			const Segment &segment = item.second;

			sprintf_rr_node(sink_rr_node, buffer);
			fprintf(routes, "Sink: %s\n", buffer);

			for (const auto &node : segment) {
				sprintf_rr_node(node, buffer);
				fprintf(routes, "\t%s\n", buffer);
			}
		}
		++inet;
	}

	fclose(routes);
}

void test_effectiveness(vector<net_t> &nets)
{
	for (int i = 0; i < nets.size(); ++i) {
		nets[i].current_local_id = i;
	}

	vector<vector<bool>> sink_scheduled(nets.size());
	for (int i = 0; i < nets.size(); ++i) {
		sink_scheduled[nets[i].current_local_id].resize(nets[i].sinks.size(), false);
	}
	vector<sink_t *> scheduled_sinks;
	choose_sinks(nets, sink_scheduled, scheduled_sinks);

	printf("scheduled sinks size: %lu\n", scheduled_sinks.size());

	vector<net_t *> nets_to_schedule;
	for (auto &net : nets) {
		nets_to_schedule.push_back(&net);
	}

	printf("chosen sink index:\n");
	for (const auto &sink : scheduled_sinks) {
		sink->net->current_sink_index = sink->id;
		printf("net %d: %d\n", sink->net->vpr_id, sink->id);
	}

	vector<vector<const net_t *>> testtest;
	schedule_nets_fast(nets_to_schedule,  testtest);
	int total = 0;
	printf("concurrency:\n");
	for (const auto &t : testtest) {
		assert(t.size() > 0);
		total += t.size();
		printf("%lu\n", t.size());
	}
	assert(total == nets.size());
	printf("num sub iters %lu average concurrency: %g\n", testtest.size(), (float)total/testtest.size());

	std::mt19937 mt(time(NULL));

	for (auto &net : nets) {
		net.current_sink_index = 0;
	}

	printf("average concurrency for random\n");
	for (int i = 0; i < 50; ++i) {
	 	for (int j = 0; j < nets.size(); ++j) {
			std::shuffle(nets[j].sinks.begin(), nets[j].sinks.end(), mt);
		}	
		vector<vector<const net_t *>> testtest;
		schedule_nets_fast(nets_to_schedule,  testtest);
		int total = 0;
		printf("concurrency:\n");
		for (const auto &t : testtest) {
			total += t.size();
			printf("%lu\n", t.size());
		}
		assert(total == nets.size());
		printf("%lu %g\n", testtest.size(), (float)total/testtest.size());
	}
}

typedef struct RouteWorkerArgs {
	net_t virtual_net;
	/*QuadTree *virtual_net_quadrant;*/
	net_t *net;
} RouteWorkerArgs;

using pending_rtree_value = pair<box, virtual_net_t *>;

struct pending_rtree_equal {
	bool operator()(const pending_rtree_value &a, const pending_rtree_value &b) const
	{
		return bg::equals(a.first, b.first) && a.second == b.second;
	}
};

using pending_rtree_t = bgi::rtree<pending_rtree_value, bgi::rstar<64>, bgi::indexable<pending_rtree_value>, pending_rtree_equal>;

void update_pending_rtree_bounding_boxes(const net_t *net, pending_rtree_t &pending_rtree)
{
	for (const auto &virtual_net : net->virtual_nets) {
		if (!virtual_net->routed) {
			zlog_level(scheduler_log, ROUTER_V3, "Removing net %d virtual net %d bb %s\n", virtual_net->sinks[0]->net->vpr_id, virtual_net->id, static_cast<const std::ostringstream &>(std::ostringstream() << bg::dsv(virtual_net->saved_scheduler_bounding_box)).str().c_str());

			if (!pending_rtree.remove(make_pair(virtual_net->saved_scheduler_bounding_box, virtual_net))) {
				zlog_error(scheduler_log, "Failed to remove net %d virtual net %d bb %s\n", virtual_net->sinks[0]->net->vpr_id, virtual_net->id, static_cast<const std::ostringstream &>(std::ostringstream() << bg::dsv(virtual_net->saved_scheduler_bounding_box)).str().c_str());
				assert(false);
			}

			zlog_level(scheduler_log, ROUTER_V3, "Reinserting net %d virtual net %d bb %s\n", virtual_net->sinks[0]->net->vpr_id, virtual_net->id, static_cast<const std::ostringstream &>(std::ostringstream() << bg::dsv(virtual_net->scheduler_bounding_box)).str().c_str());
			virtual_net->saved_scheduler_bounding_box = virtual_net->scheduler_bounding_box;

			pending_rtree.insert(make_pair(virtual_net->scheduler_bounding_box, virtual_net));
		}
	}
}

unsigned long num_leaf_node_pred_calls;
unsigned long num_internal_node_pred_calls;

template<typename Func>
void dispatch_virtual_nets(pending_rtree_t &pending_rtree, vector<virtual_net_t *> &in_flight_virtual_nets, tbb::concurrent_bounded_queue<virtual_net_t *> *update_requests, const Func &func, tbb::atomic<int> &num_dispatched_virtual_nets, sched_perf_t *perf)
{
	auto dispatch_start = std::chrono::high_resolution_clock::now();

	bool dispatched_nets = false;

	/*while (!pending_rtree.empty() && !dispatched_nets) {*/
		if (update_requests) {
			while (update_requests->size() > 0) {
				virtual_net_t *routed_virtual_net;
				update_requests->pop(routed_virtual_net); 
				assert(routed_virtual_net->routed);

				/*zlog_level(scheduler_log, ROUTER_V3, "Removing dispatched net %d virtual net %d bb %s\n", routed_virtual_net->sinks[0]->net->vpr_id, routed_virtual_net->id, static_cast<const std::ostringstream &>(std::ostringstream() << bg::dsv(routed_virtual_net->saved_scheduler_bounding_box)).str().c_str());*/
				/*printf("Removing dispatched net %d virtual net %d bb %s\n", routed_virtual_net->sinks[0]->net->vpr_id, routed_virtual_net->id, static_cast<const std::ostringstream &>(std::ostringstream() << bg::dsv(routed_virtual_net->saved_scheduler_bounding_box)).str().c_str());*/

				assert(pending_rtree.remove(make_pair(routed_virtual_net->saved_scheduler_bounding_box, routed_virtual_net)));

				auto iter = find(begin(in_flight_virtual_nets), end(in_flight_virtual_nets), routed_virtual_net);
				assert(iter != end(in_flight_virtual_nets));
				in_flight_virtual_nets.erase(iter);

				update_pending_rtree_bounding_boxes(routed_virtual_net->sinks[0]->net, pending_rtree);
			}
		}

		num_leaf_node_pred_calls = 0;
		num_internal_node_pred_calls = 0;

		for (auto iter = pending_rtree.qbegin(bgi::disjoint(&in_flight_virtual_nets)); iter != pending_rtree.qend(); ++iter) {
			/*auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::high_resolution_clock::now() - dispatch_start);*/
			/*zlog_info(scheduler_log, "Dispatch took %lld ns\n", elapsed.count());*/

#ifdef __linux__
			/*__itt_frame_begin_v3(dispatch_domain, NULL);*/
#endif
			bool not_in_flight = std::find_if(begin(in_flight_virtual_nets), end(in_flight_virtual_nets), [&iter] (virtual_net_t *in_flight) -> bool { return in_flight->sinks[0]->net == iter->second->sinks[0]->net; }) == end(in_flight_virtual_nets);

			if (not_in_flight) {
				func(iter->second);
				/*feeder.add(iter->second);*/
				/*assert(!sem_post(sargs->consume_sem));*/

				in_flight_virtual_nets.push_back(iter->second);

				/*iter = pending_rtree.qbegin(bgi::disjoint(&in_flight_virtual_nets));*/

				/*zlog_level(scheduler_log, ROUTER_V3, "Dispatched net %d virtual net %d\n", iter->second->sinks[0]->net->vpr_id, iter->second->id);*/
				/*printf("Dispatched net %d virtual net %d\n", iter->second->sinks[0]->net->vpr_id, iter->second->id);*/

				++num_dispatched_virtual_nets;

				dispatched_nets = true;
			}


#ifdef __linux__
			/*__itt_frame_end_v3(dispatch_domain, NULL);*/
#endif
		}
	/*}*/

	if (perf) {
		perf->total_dispatch_time += std::chrono::high_resolution_clock::now()-dispatch_start;
		perf->num_leaf_node_pred_calls += num_leaf_node_pred_calls;
		perf->num_internal_node_pred_calls += num_internal_node_pred_calls;
	}

	/* verification, to remove later */
	for (int i = 0; i < in_flight_virtual_nets.size(); ++i) {
		for (int j = i+1; j < in_flight_virtual_nets.size(); ++j) {
			assert(!bg::intersects(in_flight_virtual_nets[i]->scheduler_bounding_box, in_flight_virtual_nets[j]->scheduler_bounding_box));
		}
	}
}

using sch_state_lock_t = tbb::spin_mutex;

class RouteWorker2 {
	private:
		t_router_opts *opts;

		route_parameters_t &params;

		/* route state */
		vector<std::mutex> &net_locks;
		RRGraph &g;
		route_state_t *state;
		vector<route_tree_t> &route_trees;
		t_net_timing *net_timing;

		/* scheduler state */
		sch_state_lock_t &sch_state_lock;
		pending_rtree_t &pending_rtree;
		vector<virtual_net_t *> &in_flight_virtual_nets;
		tbb::concurrent_bounded_queue<virtual_net_t *> &update_requests;

		/* per thread profiling */
		perf_t *perf;
		sched_perf_t *sched_perf;

		tbb::atomic<int> &num_dispatched_virtual_nets;

	public:
		RouteWorker2(
				t_router_opts *opts,

				route_parameters_t &params,

				vector<std::mutex> &net_locks,
				RRGraph &g,
				route_state_t *state,
				vector<route_tree_t> &route_trees,
				t_net_timing *net_timing,

				sch_state_lock_t &sch_state_lock,
				pending_rtree_t &pending_rtree,
				vector<virtual_net_t *> &in_flight_virtual_nets,
				tbb::concurrent_bounded_queue<virtual_net_t *> &update_requests,

				perf_t *perf,
				sched_perf_t *sched_perf,

				tbb::atomic<int> &num_dispatched_virtual_nets
				) :
			opts(opts),

			params(params),

			net_locks(net_locks),
			g(g),
			state(state),
			route_trees(route_trees),
			net_timing(net_timing),

			sch_state_lock(sch_state_lock),
			pending_rtree(pending_rtree),
			in_flight_virtual_nets(in_flight_virtual_nets),
			update_requests(update_requests),

			perf(perf),
			sched_perf(sched_perf),

			num_dispatched_virtual_nets(num_dispatched_virtual_nets)
			{
			}

		void operator()(virtual_net_t *current_virtual_net, tbb::parallel_do_feeder<virtual_net_t *> &feeder) const
		{
			net_t *net = current_virtual_net->sinks[0]->net;

			char buffer[256];
			sprintf_rr_node(current_virtual_net->sinks[0]->rr_node, buffer);

			zlog_level(delta_log, ROUTER_V2, "Routing net %d virtual net %d\n", net->vpr_id, current_virtual_net->id);
			/*printf("sr %d %d\n", net->vpr_id, current_virtual_net->id);*/

			bool locked;
			if (net_locks[net->local_id].try_lock()) {
				locked = true;
			} else {
				zlog_error(delta_log, "Error: Multiple threads are trying to route virtual nets that belong to the same net\n");
				assert(false);
			}

			using clock = std::chrono::high_resolution_clock;

			auto rip_up_start = clock::now();

			route_tree_rip_up_marked(route_trees[net->local_id], g, params.pres_fac, false);

			perf->total_rip_up_time += clock::now()-rip_up_start;

			auto route_start = clock::now();

			route_net_2(g, net->vpr_id, current_virtual_net->source, current_virtual_net->current_sinks, params, state, route_trees[net->local_id], net_timing[net->vpr_id], perf, false);

			current_virtual_net->routed = true;

			perf->total_route_time += clock::now()-route_start;

			auto update_start = clock::now();

			vector<virtual_net_t *> unrouted;
			for (const auto &virtual_net : net->virtual_nets) {
				if (!virtual_net->routed) {
					unrouted.push_back(virtual_net);
				}
			}

			/* can be parallelized using parallel_for */
			tbb::parallel_for(tbb::blocked_range<size_t>(0, unrouted.size(), opts->grain_size),
							[&unrouted, this, &net] (const tbb::blocked_range<size_t> &r) -> void {
			/*for (const auto &virtual_net : net->virtual_nets) {*/
				/*if (!virtual_net->routed) {*/
					/*++net->num_bounding_box_updates;*/

							for (int i = r.begin(); i != r.end(); ++i) {
							virtual_net_t *virtual_net = unrouted[i];
					zlog_level(delta_log, ROUTER_V3, "Updating bounding boxes of net %d virtual net %d\n", virtual_net->sinks[0]->net->vpr_id, virtual_net->id);
					/* remove unrouted virtual net from rtree and reinsert after updating it's scheduler bounding box */
					update_virtual_net_bounding_box(*virtual_net, route_trees[net->local_id], g, params.astar_fac, perf);
					update_virtual_net_scheduler_bounding_box(*virtual_net, bg::make_inverse<box>());
				/*} else {*/
					/*zlog_level(delta_log, ROUTER_V3, "NOT updating bounding boxes of net %d virtual net %d because it is routed\n", virtual_net->sinks[0]->net->vpr_id, virtual_net->id);*/
				/*}*/
					}
			});

			perf->total_update_time += clock::now()-update_start;

			if (locked) {
				net_locks[net->local_id].unlock();
			}

			update_requests.push(current_virtual_net);

			if (!pending_rtree.empty() && !sch_state_lock.try_lock()) {
				sch_state_lock.lock();
			}
			/*printf("ss %d %d\n", net->vpr_id, current_virtual_net->id);*/
			/*if (sch_state_lock.try_lock()) {*/
			dispatch_virtual_nets(pending_rtree, in_flight_virtual_nets, &update_requests,
					[&feeder] (virtual_net_t *virtual_net) -> void {
					feeder.add(virtual_net);
					}, num_dispatched_virtual_nets, sched_perf);

			sch_state_lock.unlock();
			/*printf("ds %d %d\n", net->vpr_id, current_virtual_net->id);*/
			/*}*/
			/*printf("dr %d %d\n", net->vpr_id, current_virtual_net->id);*/
		}
};

class RouteWorker {
	private:
		std::mutex &lock;
		RRGraph &g;
		route_parameters_t &params;
		route_state_t *state;
		vector<route_tree_t> &route_trees;
		vector<trace_t> **current_traces_ptr;
		vector<trace_t> **prev_traces_ptr;
		t_net_timing *net_timing;
		perf_t &perf;
		QuadTree<virtual_net_t *> &in_flight_qt;
		/*QuadTree<sink_t *> &*/
		vector<virtual_net_t *> &unrouted_virtual_nets;
		int &num_spawned_workers;
		int &max_num_spawned_workers;

	public:
		RouteWorker(
				std::mutex &lock,
				RRGraph &g,
				route_parameters_t &params,
				route_state_t *state,
				vector<route_tree_t> &route_trees,
				vector<trace_t> **current_traces_ptr,
				vector<trace_t> **prev_traces_ptr,
				t_net_timing *net_timing,
				perf_t &perf,
				QuadTree<virtual_net_t *> &in_flight_qt,
				vector<virtual_net_t *> &unrouted_virtual_nets,
				int &num_spawned_workers,
				int &max_num_spawned_workers
				) :
			lock(lock),
			g(g),
			params(params),
			state(state),
			route_trees(route_trees),
		   	current_traces_ptr(current_traces_ptr),
		   	prev_traces_ptr(prev_traces_ptr),
		   	net_timing(net_timing),
			perf(perf),
			in_flight_qt(in_flight_qt),
			unrouted_virtual_nets(unrouted_virtual_nets),
			num_spawned_workers(num_spawned_workers),
			max_num_spawned_workers(max_num_spawned_workers)
	{
	}
		void operator()(virtual_net_t *virtual_net, tbb::parallel_do_feeder<virtual_net_t *> &feeder) const
		{
			net_t *net = virtual_net->sinks[0]->net;

			char buffer[256];
			sprintf_rr_node(virtual_net->sinks[0]->rr_node, buffer);
			zlog_info(delta_log, "Routing net %d sink %s\n", net->vpr_id, buffer);

			trace_rip_up_segment((**prev_traces_ptr)[net->local_id], g, virtual_net->sinks[0]->rr_node, params.pres_fac);
			route_net(g, net->vpr_id, virtual_net->source, virtual_net->sinks, params, state, route_trees[net->local_id], (**current_traces_ptr)[net->local_id], (**prev_traces_ptr)[net->local_id], net_timing[net->vpr_id], &perf);

			net->sink_routed[virtual_net->sinks[0]->id] = true;

			update_sink_bounding_boxes_2(*net, route_trees[net->local_id], g, params.astar_fac);
			lock.lock();
			in_flight_qt.remove({ virtual_net->sinks[0]->current_bounding_box, virtual_net });
			--num_spawned_workers;
			lock.unlock();

			/*QuadTree *virtual_net_quadrant;*/

			int i = 0;
			int num_spawned = 0;
			while (i < unrouted_virtual_nets.size() && num_spawned < 4) {
				lock.lock();

				auto virtual_net = make_pair(unrouted_virtual_nets[i]->sinks[0]->current_bounding_box, unrouted_virtual_nets[i]);
				if (!in_flight_qt.has_overlapping_boxes(virtual_net)) {
					in_flight_qt.insert(virtual_net);

					feeder.add(unrouted_virtual_nets[i]);

					sprintf_rr_node(unrouted_virtual_nets[i]->sinks[0]->rr_node, buffer);
					zlog_info(delta_log, "\tSpawning virtual net %d sink %s\n", unrouted_virtual_nets[i]->sinks[0]->net->vpr_id, buffer);
					unrouted_virtual_nets.erase(begin(unrouted_virtual_nets)+i);
					++num_spawned;
				} else {
					++i;
				}

				lock.unlock();
				++i;
			}
		}
};

FILE *sched_output_log;

void *test_scheduler_thread_tmp(void *args)
{
	SchedulerArgs *sargs = (SchedulerArgs *)args;
	int num_virtual_nets_dispatched = 0;
	while (num_virtual_nets_dispatched < sargs->virtual_nets->size()) {
		/*zlog_level(scheduler_log, ROUTER_V2, "Begin scheduler\n");*/
		virtual_net_t *vnet = &(*sargs->virtual_nets)[num_virtual_nets_dispatched];
		net_t *net = vnet->sinks[0]->net;

		net->sink_routed[vnet->sinks[0]->id] = true;

		sargs->pending_virtual_nets->push(vnet);

		vnet->routed = true;

		++num_virtual_nets_dispatched;

		assert(!sem_post(sargs->consume_sem));

		while (sem_wait(sargs->produce_sem) == -1 && errno == EINTR);

		update_sink_bounding_boxes_2(*net, (*sargs->route_trees)[net->local_id], *sargs->g, sargs->params->astar_fac);
	}

	return nullptr;
}

static void *test_scheduler_thread(void *args)
{
	SchedulerArgs *sargs = (SchedulerArgs *)args;
	int num_virtual_nets_dispatched = 0;
	vector<bool> ripped_up(sargs->num_nets, false);

	char filename[256];
	sprintf(filename, "/Volumes/DATA/scheduler_%d.log", *sargs->iteration);
	sched_output_log = fopen(filename, "w");

	zlog_level(scheduler_log, ROUTER_V1, "Scheduler thread iter %d\n", *sargs->iteration);

	for (auto &virtual_nets : *sargs->virtual_nets_by_net) {
		for (auto &vnet : virtual_nets) {
			/*zlog_level(scheduler_log, ROUTER_V2, "Begin scheduler\n");*/
			net_t *net = vnet.sinks[0]->net;

			zlog_level(scheduler_log, ROUTER_V3, "Current: (%lu sinks) net %d virtual net %d\n", net->sinks.size(), net->vpr_id, vnet.id);

			if (!ripped_up[net->local_id]) {
				route_tree_mark_nodes_to_be_ripped((*sargs->route_trees)[net->local_id], *sargs->g, 3);
				ripped_up[net->local_id] = true;
			}

			update_sink_bounding_boxes_5(virtual_nets, (*sargs->route_trees)[net->local_id], *sargs->g, sargs->params->astar_fac);

			/*net->sink_routed[vnet.sinks[0]->id] = true;*/
			/*vnet.valid = false;*/

			bool route_vnet = any_of(begin(vnet.sinks), end(vnet.sinks), [&sargs, &net] (const sink_t *sink) -> bool {
					RouteTreeNode *sink_rt_node = route_tree_get_rt_node((*sargs->route_trees)[net->local_id], sink->rr_node);
					return !sink_rt_node || sink_rt_node->properties.pending_rip_up;
					});

			if (route_vnet) {
				zlog_level(scheduler_log, ROUTER_V3, "Dispatching net %d virtual net %d\n", net->vpr_id, vnet.id);

				sargs->pending_virtual_nets->push(&vnet);

				assert(!sem_post(sargs->consume_sem));

				while (sem_wait(sargs->produce_sem) == -1 && errno == EINTR);
			} else {
				for (const auto &sink : vnet.sinks) {
					zlog_level(scheduler_log, ROUTER_V3, "Sink %d (%d) already routed\n", sink->id, sink->rr_node);
					net->sink_routed[sink->id] = true;
				}
				vnet.routed = true;
			}
		}
	}

	fclose(sched_output_log);

	return nullptr;
}

namespace boost {
	namespace geometry {
		template <>
			inline bool disjoint(const boost::geometry::model::box<boost::geometry::model::point<int, 2, boost::geometry::cs::cartesian>> &new_box,
					std::vector<virtual_net_t *> * const &in_flight_virtual_nets)
			{
				++num_leaf_node_pred_calls;
				return std::all_of(std::begin(*in_flight_virtual_nets), std::end(*in_flight_virtual_nets), [&new_box] (virtual_net_t *in_flight_virtual_net) -> bool {
						return bg::disjoint(new_box, in_flight_virtual_net->scheduler_bounding_box);
						});
			}

		template<>
			inline bool covered_by(const boost::geometry::model::box<boost::geometry::model::point<int, 2, boost::geometry::cs::cartesian>> &node_bounding_box,
					std::vector<virtual_net_t *> * const &in_flight_virtual_nets)
			{
				++num_internal_node_pred_calls;
				return std::any_of(std::begin(*in_flight_virtual_nets), std::end(*in_flight_virtual_nets), [&node_bounding_box] (virtual_net_t *in_flight_virtual_net) -> bool {
						return bg::covered_by(node_bounding_box, in_flight_virtual_net->scheduler_bounding_box);
						});
			}
	}
}

std::chrono::high_resolution_clock::time_point start_time;

static void *scheduler_thread_2(void *args)
{
	SchedulerArgs *sargs = (SchedulerArgs *)args;

	char filename[256];
	sprintf(filename, "/Volumes/DATA/scheduler_%d.log", *sargs->iteration);
	sched_output_log = fopen(filename, "w");

	zlog_level(scheduler_log, ROUTER_V1, "Scheduler thread iter %d\n", *sargs->iteration);

	sargs->perf.num_updates = 0;
	sargs->perf.num_leaf_node_pred_calls = 0;
	sargs->perf.num_internal_node_pred_calls = 0;
	sargs->perf.total_build_time = std::chrono::high_resolution_clock::duration::zero();
	sargs->perf.total_dispatch_time = std::chrono::high_resolution_clock::duration::zero();
	sargs->perf.total_rtree_update_time = std::chrono::high_resolution_clock::duration::zero();
	sargs->perf.total_wait_time = std::chrono::high_resolution_clock::duration::zero();

	int num_routed_virtual_nets = 0;
	tbb::atomic<int> num_virtual_nets = 0;

	using value = pair<box, virtual_net_t *>;
	struct equal {
		bool operator()(const value &a, const value &b) const
		{
			return bg::equals(a.first, b.first) && a.second == b.second;
		}
	};

	vector<virtual_net_t *> in_flight_virtual_nets;

#ifdef __linux__
	/*__itt_frame_begin_v3(pD, NULL);*/
#endif

	auto build_start = std::chrono::high_resolution_clock::now();

	vector<value> bulk;
	/*std::mutex bulk_lock;*/

	/*tbb::parallel_for(tbb::blocked_range<size_t>(0, sargs->nets->size(), sargs->opts->grain_size),*/
			/*[&sargs, &num_virtual_nets, &bulk, &bulk_lock] (const tbb::blocked_range<size_t> &r) -> void {*/
	/*for (auto inet = r.begin(); inet != r.end(); ++inet) {*/
		/*const net_t &net = (*sargs->nets)[inet];*/

	for (const auto &net : *sargs->nets) {
		route_tree_mark_nodes_to_be_ripped((*sargs->route_trees)[net.local_id], *sargs->g, 10000);

		for (const auto &virtual_net : net.virtual_nets) {
			update_virtual_net_current_sinks(*virtual_net, (*sargs->route_trees)[net.local_id]);
			virtual_net->routed = virtual_net->current_sinks.empty();
			if (!virtual_net->routed) {
				/* the sequence of this 2 calls are important */
				update_virtual_net_bounding_box(*virtual_net, (*sargs->route_trees)[net.local_id], *sargs->g, sargs->params->astar_fac, nullptr);
				update_virtual_net_scheduler_bounding_box(*virtual_net, (*sargs->route_trees)[net.local_id].scheduler_bounding_box);

				zlog_level(delta_log, ROUTER_V3, "\n");

				zlog_level(scheduler_log, ROUTER_V3, "Adding net %d virtual net %d bb %s\n", virtual_net->sinks[0]->net->vpr_id, virtual_net->id, static_cast<const std::ostringstream &>(std::ostringstream() << bg::dsv(virtual_net->scheduler_bounding_box)).str().c_str());

				virtual_net->saved_scheduler_bounding_box = virtual_net->scheduler_bounding_box;
				/*pending_rtree.insert(make_pair(virtual_net->scheduler_bounding_box, virtual_net));*/
				/*bulk_lock.lock();*/
				bulk.push_back(make_pair(virtual_net->scheduler_bounding_box, virtual_net));
				/*bulk_lock.unlock();*/

				++num_virtual_nets;
			} else {
				zlog_level(scheduler_log, ROUTER_V3, "NOT adding net %d virtual net %d because it is routed\n", virtual_net->sinks[0]->net->vpr_id, virtual_net->id);
			}

			zlog_level(scheduler_log, ROUTER_V3, "\n");
		}
	}

	bg::index::rtree<value, bgi::rstar<64>, bgi::indexable<value>, equal> pending_rtree(bulk);
			/*});*/

	/*FILE *sbb = fopen("scheduler_bounding_box.txt", "w");*/
	/*for (const auto &net : *sargs->nets) {*/
		/*for (const auto &virtual_net : net.virtual_nets) {*/
			/*fprintf(sbb, "%s\n", static_cast<const std::ostringstream &>(std::ostringstream() << bg::dsv(virtual_net->scheduler_bounding_box)).str().c_str());*/
		/*}*/
	/*}*/
	/*fclose(sbb);*/

#ifdef __linux__ 
	/*__itt_frame_end_v3(pD, NULL);*/
#endif

	/*bg::index::rtree<value, bgi::rstar<16>, bgi::indexable<value>, equal> pending_rtree(bulk);*/

	/*auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::high_resolution_clock::now()-build_start);*/
	sargs->perf.total_build_time += std::chrono::high_resolution_clock::now()-build_start;
	/*zlog_info(scheduler_log, "Number of virtual nets: %d\n", num_virtual_nets);*/
	/*zlog_info(scheduler_log, "Building pending rtree took %g s\n", elapsed.count() / 1e9);*/
	/*printf("Building pending rtree took %g s\n", elapsed.count() / 1e9);*/

	/*vector<bool> ripped_up(sargs->nets->size(), false);*/

	int num_dispatch_rounds = 0;

	while (num_routed_virtual_nets < num_virtual_nets) {
		/*auto pred = bg::index::satisfies([](const pair<box, virtual_net_t *> &v)->bool{return true;});*/
		/*bgi::rtree<pair<box, virtual_net_t *>, bgi::rstar<16>>::const_query_iterator iter;*/
		/*switch (in_flight_virtual_nets.size()) {*/
			/*case 0:*/
				/*iter = pending_rtree.qbegin(*/
						/*bgi::satisfies([] (const pair<box, virtual_net_t *> &v) -> bool{ return true; })*/
						/*);*/
				/*break;*/
			/*case 1:*/
				/*iter = pending_rtree.qbegin(*/
						/*bgi::disjoint(in_flight_virtual_nets[0]->bounding_box)*/
						/*);*/
				/*break;*/
			/*default:*/
				/*auto pred = bgi::disjoint(in_flight_virtual_nets[0]->bounding_box) &&*/
						/*bgi::disjoint(in_flight_virtual_nets[1]->bounding_box);*/
		 
				/*for (const auto &in_flight : in_flight_virtual_nets) {*/
					/*pred = pred && bgi::disjoint(in_flight->bounding_box);*/
				/*}*/

				/*iter = pending_rtree.qbegin(pred);*/
				/*break;*/
		/*}*/
		/*auto pred = bgi::disjoint(point(0, 0));*/
		/*for (const auto &in_flight : in_flight_virtual_nets) {*/
			/*pred = pred && bgi::disjoint(in_flight->bounding_box);*/
		/*}*/
		zlog_level(scheduler_log, ROUTER_V3, "Going to dispatching virtual net. Pending rtree size: %lu\n", pending_rtree.size());

		/*const auto &pred = bgi::satisfies([&in_flight_virtual_nets] (const pair<box, virtual_net_t *> &v) -> bool {*/
				/*bool dis = std::find_if(begin(in_flight_virtual_nets), end(in_flight_virtual_nets), [&v] (virtual_net_t *in_flight) -> bool { return in_flight->sinks[0]->net == v.second->sinks[0]->net; }) == end(in_flight_virtual_nets);*/
				/*[>bool dis = true;<]*/
				/*++num_pred_calls;*/
				/*for (int i = 0; i < in_flight_virtual_nets.size() && dis; ++i) {*/
					/*dis &= bg::disjoint(v.first, in_flight_virtual_nets[i]->scheduler_bounding_box);*/
				/*}*/
				/*return dis;*/
				/*}*/
				/*);*/

		/*vector<virtual_net_t *> dispatched_virtual_nets;*/
		auto dispatch_start = std::chrono::high_resolution_clock::now();
		num_leaf_node_pred_calls = 0;
		num_internal_node_pred_calls = 0;
		/*vector<box> boxes;*/
		int num_dispatched_nets = 0;
		for (auto iter = pending_rtree.qbegin(bgi::disjoint(&in_flight_virtual_nets)); iter != pending_rtree.qend(); ++iter) {
			/* remove the iter from pending_rtree */
			/* add the iter to in_flight_virtual_nets */
			/* loop again */
			/*auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::high_resolution_clock::now() - dispatch_start);*/
			/*zlog_info(scheduler_log, "Dispatch took %lld ns\n", elapsed.count());*/

#ifdef __linux__
			/*__itt_frame_begin_v3(dispatch_domain, NULL);*/
#endif
			bool not_in_flight = std::find_if(begin(in_flight_virtual_nets), end(in_flight_virtual_nets), [&iter] (virtual_net_t *in_flight) -> bool { return in_flight->sinks[0]->net == iter->second->sinks[0]->net; }) == end(in_flight_virtual_nets);

			if (not_in_flight) {
				sargs->pending_virtual_nets->push(iter->second);

				/*assert(!sem_post(sargs->consume_sem));*/

				in_flight_virtual_nets.push_back(iter->second);

				/*printf("Dispatching net %d virtual net %d\n", iter->second->sinks[0]->net->vpr_id, iter->second->id);*/
				/*printf("Num internal %d Num leaf %d\n", num_internal_node_pred_calls, num_leaf_node_pred_calls);*/
				++num_dispatched_nets;

				/*dispatched_virtual_nets.push_back(iter->second);*/

				/*iter = pending_rtree.qbegin(bgi::disjoint(&in_flight_virtual_nets));*/
			}

			/*zlog_level(scheduler_log, ROUTER_V3, "Dispatching net %d virtual net %d\n", iter->second->sinks[0]->net->vpr_id, iter->second->id);*/

#ifdef __linux__
			/*__itt_frame_end_v3(dispatch_domain, NULL);*/
#endif

			/*auto before_post = std::chrono::high_resolution_clock::now();*/
			/*auto post_time = std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::high_resolution_clock::now()-before_post);*/

			/*zlog_info(scheduler_log, "Posted at %lld ns in %lld\n", std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::high_resolution_clock::now()-start_time).count(), post_time.count());*/

			/*dispatch_start = std::chrono::high_resolution_clock::now();*/
		}
		/*printf("Num virtual nets: %d Num pred calls: %d\n", num_virtual_nets, num_pred_calls);*/
		sargs->perf.total_dispatch_time += std::chrono::high_resolution_clock::now()-dispatch_start;
		sargs->perf.num_leaf_node_pred_calls += num_leaf_node_pred_calls;
		sargs->perf.num_internal_node_pred_calls += num_internal_node_pred_calls;

		++num_dispatch_rounds;

		/*for (const auto &virtual_net : dispatched_virtual_nets) {*/
			/*zlog_level(scheduler_log, ROUTER_V3, "Removing dispatched net %d virtual net %d bb %s\n", virtual_net->sinks[0]->net->vpr_id, virtual_net->id, static_cast<const std::ostringstream &>(std::ostringstream() << bg::dsv(virtual_net->saved_scheduler_bounding_box)).str().c_str());*/
			/*assert(pending_rtree.remove(make_pair(virtual_net->saved_scheduler_bounding_box, virtual_net)));*/
		/*}*/

		/* verification, to remove later */
		for (int i = 0; i < in_flight_virtual_nets.size(); ++i) {
			for (int j = i+1; j < in_flight_virtual_nets.size(); ++j) {
				assert(!bg::intersects(in_flight_virtual_nets[i]->scheduler_bounding_box, in_flight_virtual_nets[j]->scheduler_bounding_box));
			}
		}

#ifdef __linux__
		/*__itt_frame_begin_v3(update_domain, NULL);*/
#endif

		/*for (int i = 0; i < num_dispatched_nets; ++i) {*/
		do {
			zlog_level(scheduler_log, ROUTER_V3, "Going to wait for net to finish routing\n");
			/* we need to wait for existing in flight nets to complete */
			auto wait_start = std::chrono::high_resolution_clock::now();

			virtual_net_t *routed_vnet;
			sargs->routed_virtual_nets->pop(routed_vnet);
			/*while (sem_wait(sargs->produce_sem) == -1 && errno == EINTR);*/

			sargs->perf.total_wait_time += std::chrono::high_resolution_clock::now()-wait_start;
			/* remove the completed net from in_flight_virtual_net */

			auto update_start = std::chrono::high_resolution_clock::now();

			/*assert(sargs->routed_virtual_nets->try_pop(routed_vnet));*/

			net_t *net = routed_vnet->sinks[0]->net;

			bool locked;
			if ((*sargs->net_locks)[net->local_id].try_lock()) {
				locked = true;
			} else {
				assert(false);
			}

			assert(routed_vnet->routed);

			zlog_level(scheduler_log, ROUTER_V3, "Removing routed net %d virtual net %d bb %s\n", routed_vnet->sinks[0]->net->vpr_id, routed_vnet->id, static_cast<const std::ostringstream &>(std::ostringstream() << bg::dsv(routed_vnet->saved_scheduler_bounding_box)).str().c_str());

			assert(pending_rtree.remove(make_pair(routed_vnet->saved_scheduler_bounding_box, routed_vnet)));

			auto iter = find(begin(in_flight_virtual_nets), end(in_flight_virtual_nets), routed_vnet);
			assert(iter != end(in_flight_virtual_nets));
			in_flight_virtual_nets.erase(iter);
			/*std::remove(begin(in_flight_virtual_nets), end(in_flight_virtual_nets), routed_vnet);*/

			for (const auto &virtual_net : net->virtual_nets) {
				if (!virtual_net->routed) {
					zlog_level(scheduler_log, ROUTER_V3, "Removing net %d virtual net %d bb %s\n", virtual_net->sinks[0]->net->vpr_id, virtual_net->id, static_cast<const std::ostringstream &>(std::ostringstream() << bg::dsv(virtual_net->saved_scheduler_bounding_box)).str().c_str());

					if (!pending_rtree.remove(make_pair(virtual_net->saved_scheduler_bounding_box, virtual_net))) {
						zlog_error(scheduler_log, "Failed to remove net %d virtual net %d bb %s\n", virtual_net->sinks[0]->net->vpr_id, virtual_net->id, static_cast<const std::ostringstream &>(std::ostringstream() << bg::dsv(virtual_net->saved_scheduler_bounding_box)).str().c_str());
						assert(false);
					}

					zlog_level(scheduler_log, ROUTER_V3, "Reinserting net %d virtual net %d bb %s\n", virtual_net->sinks[0]->net->vpr_id, virtual_net->id, static_cast<const std::ostringstream &>(std::ostringstream() << bg::dsv(virtual_net->scheduler_bounding_box)).str().c_str());
					virtual_net->saved_scheduler_bounding_box = virtual_net->scheduler_bounding_box;

					pending_rtree.insert(make_pair(virtual_net->scheduler_bounding_box, virtual_net));
				}
			}

			if (locked) {
				(*sargs->net_locks)[net->local_id].unlock();
			}
			/* do we need to remark all the unrouted virtual net's route tree again? */
			/* need to remove the unrouted virtual nets from the pending rtree and reinsert because the bounding boxes have been update */
			/*for (const auto &net : *sargs->nets) {*/
				/*if (!net.ripped_up) {*/
					/*route_tree_mark_nodes_to_be_ripped((*sargs->route_trees)[net.local_id], *sargs->g, 3);*/
					/*update_sinks(net->virtual_nets, (*sargs->route_trees)[net->local_id], *sargs->g, sargs->params->astar_fac);*/
				/*}*/
			/*}*/

			++num_routed_virtual_nets;

			sargs->perf.total_rtree_update_time += std::chrono::high_resolution_clock::now()-update_start;

			++sargs->perf.num_updates;
		} while (sargs->routed_virtual_nets->size() > 0);
		/*}*/

#ifdef __linux__ 
		/*__itt_frame_end_v3(update_domain, NULL);*/
#endif
	}

	assert(sargs->perf.num_updates == num_virtual_nets);

	printf("Num virtual nets: %d Num dispatch rounds: %d Average concurrency: %g\n", num_virtual_nets, num_dispatch_rounds, (float)num_virtual_nets/num_dispatch_rounds);

	return nullptr;
}

static void *scheduler_thread(void *args)
{
	SchedulerArgs *sargs = (SchedulerArgs *)args;
	int num_virtual_nets_dispatched = 0;
	while (num_virtual_nets_dispatched < sargs->virtual_nets->size()) {
		/*zlog_level(scheduler_log, ROUTER_V2, "Begin scheduler\n");*/

#ifdef __linux__
		__itt_frame_begin_v3(pD, NULL);
#endif
		for (auto iter = begin(*sargs->virtual_nets); iter != end(*sargs->virtual_nets); ++iter) {
			if (iter->routed) {
				continue;
			}

			if (*sargs->iteration > 0) {
				adjust_bounding_box(*iter->sinks[0]);
			}

			auto virtual_net = make_pair(iter->sinks[0]->current_bounding_box, &(*iter));

			if (!sargs->in_flight_qt->has_overlapping_boxes(virtual_net)) {
				zlog_level(scheduler_log, ROUTER_V3, "Spawning net %d sink %d %d-%d %d-%d\n", iter->sinks[0]->net->vpr_id, iter->sinks[0]->id, virtual_net.first.xmin, virtual_net.first.xmax, virtual_net.first.ymin, virtual_net.first.ymax);

				iter->sinks[0]->net->sink_routed[iter->sinks[0]->id] = true;

				sargs->in_flight_qt->insert(virtual_net);

				sargs->pending_virtual_nets->push(&(*iter));

				iter->routed = true;

				++num_virtual_nets_dispatched;

				assert(!sem_post(sargs->consume_sem));

				/*sprintf_rr_node(unrouted_virtual_nets[i]->sink->rr_node, buffer);*/
			}
		}
#ifdef __linux__ 
		__itt_frame_end_v3(pD, NULL);
#endif

		int num_routed_virtual_nets = sargs->routed_virtual_nets->size();
		/* TODO: when num_routed_virtual_nets is zero, we should wait here instead of 
		 * going into the scheduler loop again because we can be sure that no nets will be scheduled */

		/*zlog_info(delta_log, "num_routed_virtual_nets: %d\n", num_routed_virtual_nets);*/

		for (int i = 0; i < num_routed_virtual_nets; ++i) {
			while (sem_wait(sargs->produce_sem) == -1 && errno == EINTR);
			virtual_net_t *routed_vnet;
			assert(sargs->routed_virtual_nets->try_pop(routed_vnet));
			net_t *net = routed_vnet->sinks[0]->net;
			zlog_level(scheduler_log, ROUTER_V3, "Removing net %d sink %d %d-%d %d-%d\n", net->vpr_id, routed_vnet->sinks[0]->id,
					routed_vnet->sinks[0]->current_bounding_box.xmin,
					routed_vnet->sinks[0]->current_bounding_box.xmax,
					routed_vnet->sinks[0]->current_bounding_box.ymin,
					routed_vnet->sinks[0]->current_bounding_box.ymax
					);
			sargs->in_flight_qt->remove({ routed_vnet->sinks[0]->current_bounding_box, routed_vnet });

			update_sink_bounding_boxes_2(*net, (*sargs->route_trees)[net->local_id], *sargs->g, sargs->params->astar_fac);
		}
	}

	/* no more valid nets but we might still have nets being routed */
	while (!sargs->in_flight_qt->empty()) {
		/* wait for inflight nets to be routed */
		while (sem_wait(sargs->produce_sem) == -1 && errno == EINTR);
		sargs->in_flight_qt->print_items([] (const pair<bounding_box_t, virtual_net_t *> &item) -> void {
				/*zlog_level(delta_log, ROUTER_V2, "Net %d sink %d in in_flight_qt\n", item.second->sink->net->vpr_id, item.second->sink->id);*/
				});
		
		virtual_net_t *routed_vnet;
		assert(sargs->routed_virtual_nets->try_pop(routed_vnet));
		net_t *net = routed_vnet->sinks[0]->net;
		zlog_level(scheduler_log, ROUTER_V3, "Last removing net %d sink %d\n", net->vpr_id, routed_vnet->sinks[0]->id);
		sargs->in_flight_qt->remove({ routed_vnet->sinks[0]->current_bounding_box, routed_vnet });

		assert(all_of(begin(net->sink_routed), end(net->sink_routed), [] (bool item) -> bool { return item; }));
	}

	return nullptr;
}

static void *worker_thread_3(void *args)
{
	WorkerArgs *wargs = (WorkerArgs *)args;

	char tid[256];
	sprintf(tid, "%d", wargs->tid);
	assert(!zlog_put_mdc("tid", tid));

	using clock = std::chrono::high_resolution_clock;

	while (true) {
		/*while (sem_wait(wargs->consume_sem) == -1 && errno == EINTR);*/
		auto wait_start = clock::now();

		virtual_net_t *current_virtual_net = nullptr;	
		wargs->pending_virtual_nets->pop(current_virtual_net);
		/*assert(!sem_wait(wargs->consume_sem));*/

		wargs->perf.total_wait_time += clock::now()-wait_start;


		/*zlog_info(scheduler_log, "Started at %lld ns\n", std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::high_resolution_clock::now()-start_time).count());*/


		net_t *net = current_virtual_net->sinks[0]->net;

		char buffer[256];
		sprintf_rr_node(current_virtual_net->sinks[0]->rr_node, buffer);
		/*zlog_info(delta_log, "Routing net %d sink %s\n", net->vpr_id, buffer);*/

		bool locked;
		if (wargs->net_locks[net->local_id].try_lock()) {
			locked = true;
		} else {
			zlog_error(delta_log, "Error: Multiple threads are trying to route virtual nets that belong to the same net\n");
			assert(false);
		}

		zlog_level(delta_log, ROUTER_V2, "Ripping up segments for net %d virtual net %d\n", net->vpr_id, current_virtual_net->id);
		/*trace_rip_up_segment((**wargs->prev_traces_ptr)[net->local_id], wargs->g, current_virtual_net->sinks[0]->rr_node, wargs->params.pres_fac);*/
		/*route_tree_rip_up_segment_2(wargs->route_trees[net->local_id], current_virtual_net->sinks[0]->rr_node, wargs->g, wargs->params.pres_fac);*/

		auto rip_up_start = clock::now();

		route_tree_rip_up_marked(wargs->route_trees[net->local_id], wargs->g, wargs->params.pres_fac, false);

		wargs->perf.total_rip_up_time += clock::now()-rip_up_start;

		auto route_start = clock::now();

		route_net_2(wargs->g, net->vpr_id, current_virtual_net->source, current_virtual_net->current_sinks, wargs->params, wargs->state, wargs->route_trees[net->local_id], wargs->net_timing[net->vpr_id], &wargs->perf, false);

		/*for (const auto &sink : current_virtual_net->sinks) {*/
			/*net->sink_routed[sink->id] = true;*/
		/*}*/
		current_virtual_net->routed = true;

		wargs->perf.total_route_time += clock::now()-route_start;

		auto update_start = clock::now();

		for (const auto &virtual_net : net->virtual_nets) {
			if (!virtual_net->routed) {
				++net->num_bounding_box_updates;

				zlog_level(delta_log, ROUTER_V3, "Updating bounding boxes of net %d virtual net %d\n", virtual_net->sinks[0]->net->vpr_id, virtual_net->id);
				/* remove unrouted virtual net from rtree and reinsert after updating it's scheduler bounding box */
				update_virtual_net_bounding_box(*virtual_net, wargs->route_trees[net->local_id], wargs->g, wargs->params.astar_fac, &wargs->perf);

				auto scheduler_box_start = clock::now();

				update_virtual_net_scheduler_bounding_box(*virtual_net, bg::make_inverse<box>());

				wargs->perf.total_scheduler_box_time += clock::now()-scheduler_box_start;
			} else {
				zlog_level(delta_log, ROUTER_V3, "NOT updating bounding boxes of net %d virtual net %d because it is routed\n", virtual_net->sinks[0]->net->vpr_id, virtual_net->id);
			}
		}

		/*current_virtual_net->sinks[0]->previous_bounding_box = current_virtual_net->sinks[0]->current_bounding_box;*/

		wargs->perf.total_update_time += clock::now()-update_start;

		if (locked) {
			wargs->net_locks[net->local_id].unlock();
		}

		auto push_start = clock::now();

		wargs->routed_virtual_nets->push(current_virtual_net);
		zlog_level(delta_log, ROUTER_V3, "Done routing net %d virtual_net %d in %lld\n", net->vpr_id, current_virtual_net->id,
				std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::high_resolution_clock::now()-route_start).count());

		wargs->perf.total_push_time += clock::now()-push_start;

		/*assert(!sem_post(wargs->produce_sem));*/
		/*zlog_info(scheduler_log, "Done routing in %lld ns\n", std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::high_resolution_clock::now()-route_start_time).count());*/
	}
}

static void *worker_thread_2(void *args)
{
	WorkerArgs *wargs = (WorkerArgs *)args;

	char tid[256];
	sprintf(tid, "%d", wargs->tid);
	assert(!zlog_put_mdc("tid", tid));

	while (true) {
		while (sem_wait(wargs->consume_sem) == -1 && errno == EINTR);

		virtual_net_t *current_virtual_net = nullptr;	
		assert(wargs->pending_virtual_nets->try_pop(current_virtual_net));

		net_t *net = current_virtual_net->sinks[0]->net;

		char buffer[256];
		sprintf_rr_node(current_virtual_net->sinks[0]->rr_node, buffer);
		/*zlog_info(delta_log, "Routing net %d sink %s\n", net->vpr_id, buffer);*/

		bool locked;
		if (wargs->net_locks[net->local_id].try_lock()) {
			locked = true;
		} else {
			assert(false);
		}

		trace_rip_up_segment((**wargs->prev_traces_ptr)[net->local_id], wargs->g, current_virtual_net->sinks[0]->rr_node, wargs->params.pres_fac);
		vector<sink_t *> sinks = { current_virtual_net->sinks[0] };
		route_net(wargs->g, net->vpr_id, current_virtual_net->source, sinks, wargs->params, wargs->state, wargs->route_trees[net->local_id], (**wargs->current_traces_ptr)[net->local_id], (**wargs->prev_traces_ptr)[net->local_id], wargs->net_timing[net->vpr_id], nullptr);

		/*net->sink_routed[current_virtual_net->sinks[0]->id] = true;*/

		current_virtual_net->sinks[0]->previous_bounding_box = current_virtual_net->sinks[0]->current_bounding_box;

		if (locked) {
			wargs->net_locks[net->local_id].unlock();
		}

		wargs->routed_virtual_nets->push(current_virtual_net);
		assert(!sem_post(wargs->produce_sem));
		zlog_level(delta_log, ROUTER_V2, "Done routing net %d sink %d\n", net->vpr_id, current_virtual_net->sinks[0]->id);
	}
}

static void worker_thread_internal(WorkerArgs *wargs)
{
	char tid[256];
	sprintf(tid, "%d", wargs->tid);
	assert(!zlog_put_mdc("tid", tid));

	while (true) {
		virtual_net_t *current_virtual_net = nullptr;	
		wargs->lock.lock();
		/*for (auto iter = begin(wargs->virtual_nets); iter != end(wargs->virtual_nets) && !current_virtual_net; ) {*/
			/*if (!iter->valid) {*/
				/*++iter;*/
				/*continue;*/
			/*}*/

			/*if (*wargs->iteration > 0) {*/
				/*adjust_bounding_box(*iter->sinks[0]);*/
			/*}*/

			/*auto virtual_net = make_pair(iter->sinks[0]->current_bounding_box, &(*iter));*/

			/*if (!wargs->in_flight_qt.has_overlapping_boxes(virtual_net)) {*/
				/*wargs->in_flight_qt.insert(virtual_net);*/

				/*current_virtual_net = &(*iter);*/

				/*iter->valid = false;*/
				/*++wargs->num_routed_virtual_nets;*/

				/*[>sprintf_rr_node(unrouted_virtual_nets[i]->sink->rr_node, buffer);<]*/
				/*[>zlog_info(delta_log, "\tSpawning virtual net %d sink %s\n", unrouted_virtual_nets[i]->sink->net->vpr_id, buffer);<]*/
			/*}*/

			/*++iter;*/
		/*}*/
		/*if (!current_virtual_net) {	*/
			/*if (wargs->num_routed_virtual_nets == wargs->virtual_nets.size()) {*/
				/*[> wait for other threads to complete <]*/
				/*wargs->lock.unlock();*/
				/*my_pthread_barrier_wait(wargs->barrier, 0);*/
				/*break;*/
			/*} else {*/
				/*[> we didn't manage to find a non overlapping net to route, try to find again <]*/
				/*wargs->lock.unlock();*/
				/*continue;*/
			/*}*/
		/*}*/
		wargs->lock.unlock();

		net_t *net = current_virtual_net->sinks[0]->net;

		char buffer[256];
		sprintf_rr_node(current_virtual_net->sinks[0]->rr_node, buffer);
		/*zlog_info(delta_log, "Routing net %d sink %s\n", net->vpr_id, buffer);*/

		if (!wargs->net_locks[net->local_id].try_lock()) {
			assert(false);
		} else {
			wargs->net_locks[net->local_id].unlock();
		}

		trace_rip_up_segment((**wargs->prev_traces_ptr)[net->local_id], wargs->g, current_virtual_net->sinks[0]->rr_node, wargs->params.pres_fac);
		vector<sink_t *> sinks = { current_virtual_net->sinks[0] };
		route_net(wargs->g, net->vpr_id, current_virtual_net->source, sinks, wargs->params, wargs->state, wargs->route_trees[net->local_id], (**wargs->current_traces_ptr)[net->local_id], (**wargs->prev_traces_ptr)[net->local_id], wargs->net_timing[net->vpr_id], nullptr);

		net->sink_routed[current_virtual_net->sinks[0]->id] = true;

		current_virtual_net->sinks[0]->previous_bounding_box = current_virtual_net->sinks[0]->current_bounding_box;

		wargs->lock.lock();
		update_sink_bounding_boxes_2(*net, wargs->route_trees[net->local_id], wargs->g, wargs->params.astar_fac);
		wargs->in_flight_qt.remove({ current_virtual_net->sinks[0]->current_bounding_box, current_virtual_net });
		wargs->lock.unlock();
	}
}

static void *worker_thread(void *args)
{
	WorkerArgs *wargs = (WorkerArgs *)args;

	while (true) {
		zlog_info(delta_log, "Thread started\n");
		pthread_mutex_lock(wargs->start_cond_lock);
		pthread_cond_wait(wargs->start_cond, wargs->start_cond_lock);
		pthread_mutex_unlock(wargs->start_cond_lock);

		worker_thread_internal(wargs);
	}
}

FILE *current_output_log;

int zlog_custom_output(zlog_msg_t *msg)
{
	fprintf(current_output_log, "%s", msg->buf);
	fflush(current_output_log);
	return 0;
}

int zlog_sched_custom_output(zlog_msg_t *msg)
{
	fprintf(sched_output_log, "%s", msg->buf);
	fflush(sched_output_log);
	return 0;
}

void test_perf()
{
	std::mutex lock1;
	pthread_mutex_t lock2;
	tbb::spin_mutex lock3;
	tbb::queuing_mutex lock4;
	tbb::queuing_mutex::scoped_lock lock5;
	pthread_mutex_init(&lock2, NULL);

	auto time = std::chrono::high_resolution_clock::now();

	for (int i = 0; i < 1000000; ++i) {
		lock1.lock();
		lock1.unlock();
	}

	auto elapsed = std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::high_resolution_clock::now()-time);
	printf("std mutex took %lld\n", elapsed.count());

	time = std::chrono::high_resolution_clock::now();

	for (int i = 0; i < 1000000; ++i) {
		pthread_mutex_lock(&lock2);
		pthread_mutex_unlock(&lock2);
	}

	elapsed = std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::high_resolution_clock::now()-time);
	printf("pthread mutex took %lld\n", elapsed.count());

	time = std::chrono::high_resolution_clock::now();

	for (int i = 0; i < 1000000; ++i) {
		lock3.lock();
		lock3.unlock();
	}

	elapsed = std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::high_resolution_clock::now()-time);
	printf("spin_mutex took %lld\n", elapsed.count());

	time = std::chrono::high_resolution_clock::now();

	for (int i = 0; i < 1000000; ++i) {
		lock5.acquire(lock4);
		lock5.release();
	}

	elapsed = std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::high_resolution_clock::now()-time);
	printf("queuing_mutex took %lld\n", elapsed.count());

	time = std::chrono::high_resolution_clock::now();
	for (int i = 0; i < 1000000; ++i) {
		zlog_info(static_log, "zlog test string\n");
	}
	elapsed = std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::high_resolution_clock::now()-time);
	printf("zlog_info static took %lld\n", elapsed.count());

	zlog_put_mdc("iter", "0");
	zlog_set_record("custom_output", zlog_custom_output);
	time = std::chrono::high_resolution_clock::now();
	for (int i = 0; i < 1000000; ++i) {
		zlog_info(dynamic_log, "zlog test string\n");
	}
	elapsed = std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::high_resolution_clock::now()-time);
	printf("zlog_info dynamic took %lld\n", elapsed.count());

	time = std::chrono::high_resolution_clock::now();
	for (int i = 0; i < 1000000; ++i) {
		printf("zlog test string\n");
	}
	elapsed = std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::high_resolution_clock::now()-time);
	printf("printf took %lld\n", elapsed.count());
}

sem_t *create_sem(const char *name)
{
	sem_t *sem = sem_open(name, O_CREAT | O_EXCL, S_IWUSR | S_IRUSR, 0);
	if (sem == SEM_FAILED) {
		if (errno == EEXIST) {
			printf("Existing start route sem sem, recreating\n");
			assert(!sem_unlink(name));
			sem = sem_open(name, O_CREAT | O_EXCL, S_IWUSR | S_IRUSR, 0);
			assert(sem != SEM_FAILED);
		} else {
			char str[256];
			strerror_r(errno, str, sizeof(str));
			printf("failed to create start route sem: %s\n", str);
			exit(-1);
		}
	}
	return sem;
}

void test_route_tree()
{
	RRGraph g;
	add_vertex(g, 10);
	for (int i = 0; i < 10; ++i) {
		get_vertex(g, i).properties.occ = 0;
		get_vertex(g, i).properties.capacity = 1;
	}

	route_tree_t rt;
	route_tree_init(rt);
	for (int i = 0; i < 10; ++i) {
		route_tree_add_rr_node(rt, get_vertex(g, i));
	}
	route_tree_set_root(rt, 0);
	route_tree_add_edge_between_rr_node(rt, 0, 1);
	route_tree_add_edge_between_rr_node(rt, 1, 2);
	route_tree_add_edge_between_rr_node(rt, 1, 3);
	route_tree_add_edge_between_rr_node(rt, 2, 4);
	route_tree_add_edge_between_rr_node(rt, 4, 5);
	route_tree_add_edge_between_rr_node(rt, 4, 6);
	route_tree_add_edge_between_rr_node(rt, 6, 7);

	/* case */
	get_vertex(g, 0).properties.occ = 2;
	route_tree_mark_nodes_to_be_ripped(rt, g, 3);
	assert(route_tree_get_rt_node(rt, 0)->properties.pending_rip_up);
	assert(route_tree_get_rt_node(rt, 1)->properties.pending_rip_up);
	assert(route_tree_get_rt_node(rt, 2)->properties.pending_rip_up);
	assert(route_tree_get_rt_node(rt, 3)->properties.pending_rip_up);
	assert(route_tree_get_rt_node(rt, 4)->properties.pending_rip_up);
	assert(route_tree_get_rt_node(rt, 5)->properties.pending_rip_up);
	assert(route_tree_get_rt_node(rt, 6)->properties.pending_rip_up);
	assert(route_tree_get_rt_node(rt, 7)->properties.pending_rip_up);

	/* case */
	get_vertex(g, 0).properties.occ = 0;
	get_vertex(g, 2).properties.occ = 2;
	route_tree_mark_nodes_to_be_ripped(rt, g, 3);
	assert(!route_tree_get_rt_node(rt, 0)->properties.pending_rip_up);
	assert(!route_tree_get_rt_node(rt, 1)->properties.pending_rip_up);
	assert(!route_tree_get_rt_node(rt, 3)->properties.pending_rip_up);
	assert(route_tree_get_rt_node(rt, 2)->properties.pending_rip_up);
	assert(route_tree_get_rt_node(rt, 4)->properties.pending_rip_up);
	assert(route_tree_get_rt_node(rt, 5)->properties.pending_rip_up);
	assert(route_tree_get_rt_node(rt, 6)->properties.pending_rip_up);
	assert(route_tree_get_rt_node(rt, 7)->properties.pending_rip_up);
}

void test_clustering(int, int);

void create_clustered_virtual_nets(vector<net_t> &nets, int num_clusters, int sink_bb_area_threshold, vector<vector<virtual_net_t>> &virtual_nets);

namespace boost {
	namespace geometry {
		namespace detail {
			namespace dispatch {

				/*template <typename Point, typename Box>*/
					/*struct covered_by<Point, Box, point_tag, box_tag>*/
					/*{*/
						/*template <typename Strategy>*/
							/*static inline bool apply(Point const& point, Box const& box, Strategy const& strategy)*/
							/*{*/

								/*return true;*/
							/*}*/
					/*}*/
				/*bool covered_by(const box &b, const vector<box> &boxes)*/
				/*{*/
				/*return true;*/
				/*}*/
			}
		}
	}
}

void test_geometry()
{
	using polygon = bg::model::polygon<point>;
	using multi_polygon = bg::model::multi_polygon<polygon>;

	polygon poly;
	bg::append(poly.outer(), point(0, 0));
	bg::append(poly.outer(), point(0, 5));
	bg::append(poly.outer(), point(5, 5));
	bg::append(poly.outer(), point(5, 0));
	bg::append(poly.outer(), point(0, 0));

	polygon poly2;
	bg::append(poly2.outer(), point(0, 0));
	bg::append(poly2.outer(), point(0, 2));
	bg::append(poly2.outer(), point(2, 2));
	bg::append(poly2.outer(), point(2, 0));
	bg::append(poly2.outer(), point(0, 0));

	multi_polygon mpoly;
	mpoly.push_back(poly);

	/*covered_by*/

	point p(3, 3);
	/*box */
		/*vector<*/
	/*bg::covered_by(poly2, m);*/
}

void print_virtual_nets(const vector<net_t> &nets);

pending_rtree_t build_pending_rtree(const vector<net_t> &nets, const RRGraph &g, const route_parameters_t &params, vector<route_tree_t> &route_trees, sched_perf_t *perf)
{
	vector<pending_rtree_value> bulk;

	for (const auto &net : nets) {
		route_tree_mark_nodes_to_be_ripped(route_trees[net.local_id], g, 10000);

		for (const auto &virtual_net : net.virtual_nets) {
			update_virtual_net_current_sinks(*virtual_net, route_trees[net.local_id]);
			virtual_net->routed = virtual_net->current_sinks.empty();
			if (!virtual_net->routed) {
				/* the sequence of this 2 calls are important */
				update_virtual_net_bounding_box(*virtual_net, route_trees[net.local_id], g, params.astar_fac, nullptr);
				update_virtual_net_scheduler_bounding_box(*virtual_net, route_trees[net.local_id].scheduler_bounding_box);

				zlog_level(delta_log, ROUTER_V3, "\n");

				zlog_level(scheduler_log, ROUTER_V3, "Adding net %d virtual net %d bb %s\n", virtual_net->sinks[0]->net->vpr_id, virtual_net->id, static_cast<const std::ostringstream &>(std::ostringstream() << bg::dsv(virtual_net->scheduler_bounding_box)).str().c_str());

				virtual_net->saved_scheduler_bounding_box = virtual_net->scheduler_bounding_box;
				/*pending_rtree.insert(make_pair(virtual_net->scheduler_bounding_box, virtual_net));*/
				/*bulk_lock.lock();*/
				bulk.push_back(make_pair(virtual_net->scheduler_bounding_box, virtual_net));
				/*bulk_lock.unlock();*/
			} else {
				zlog_level(scheduler_log, ROUTER_V3, "NOT adding net %d virtual net %d because it is routed\n", virtual_net->sinks[0]->net->vpr_id, virtual_net->id);
			}

			zlog_level(scheduler_log, ROUTER_V3, "\n");
		}
	}

	return pending_rtree_t(bulk);
}

bool tbb_greedy_route(t_router_opts *opts)
{
	tbb::task_scheduler_init init(opts->num_threads);

	start_time = std::chrono::high_resolution_clock::now();

/*#ifdef __linux__*/
	/*printf("domain flags: %d\n", pD->flags);*/
	/*pD->flags = 1;*/
/*#endif*/
	
	init_logging();
	zlog_set_record("custom_output", zlog_custom_output);
	zlog_set_record("sched_custom_output", zlog_sched_custom_output);
	current_output_log = fopen("/Volumes/DATA/main.log", "w");

	zlog_info(delta_log, "Using boost version %d\n", BOOST_VERSION);

#if (TBB_USE_DEBUG==1)
	zlog_warn(delta_log, "Using debug version of TBB\n");
#endif

	/*test_rtree();*/

	/*test_route_tree();*/

	/*test_geometry();*/

	/*test_clustering(opts->num_threads, opts->bb_factor);*/
	/*return;*/

	/*test_perf();*/
	/*exit(-1);*/

	/*test_quadrant();*/

	/*auto time = std::chrono::high_resolution_clock::now();*/

	/*printf("test\n");*/

	/*auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::high_resolution_clock::now()-time);*/
	/*printf("printf took %g\n", elapsed.count());*/

	/*test_heap();*/

	/*exit(0);*/

	/*test_scheduler();*/

	/*test_rtree();*/

	/*test_graph();*/

	RRGraph g;
	init_graph(g);

	vector<net_t> nets;
	vector<net_t> global_nets;
	init_nets(nets, global_nets, opts->bb_factor);

	/*test_effectiveness(nets);*/

	/*exit(0);*/

	/*sort_sinks(nets);*/

	/*dump_net_bounding_boxes(nets);*/

	/* calculating net bounding box area rank */
	vector<net_t *> nets_ptr(nets.size());
	for (int i = 0; i < nets.size(); ++i) {
		nets_ptr[i] = &nets[i];
	}
	extern s_bb *route_bb;
	sort(nets_ptr.begin(), nets_ptr.end(), [] (const net_t *a, const net_t *b) -> bool {
			return get_bounding_box_area(route_bb[a->vpr_id]) > get_bounding_box_area(route_bb[b->vpr_id]);
			});
	int rank = 0;
	for (auto &net : nets_ptr) {
		net->bb_area_rank = rank++;
	}
	/* end */ 

	extern int nx;
	extern int ny;
	zlog_info(delta_log, "Num nets: %lu Num global nets: %lu nx: %d ny: %d\n", nets.size(), global_nets.size(), nx, ny);
	
	t_net_timing *net_timing = new t_net_timing[nets.size()+global_nets.size()];
	init_net_timing(nets, global_nets, net_timing);

	vector<route_tree_t> route_trees(nets.size());
	for (int i = 0; i < nets.size(); ++i) {
		route_tree_init(route_trees[i]);
	}

	vector<trace_t> prev_traces(nets.size());
	vector<trace_t> traces(nets.size());
	for (int i = 0; i < nets.size(); ++i) {
		trace_init(traces[i]);
		trace_init(prev_traces[i]);
	}

	vector<trace_t> *prev_traces_ptr = &prev_traces;
	vector<trace_t> *current_traces_ptr = &traces;

	route_state_t *state = new route_state_t[num_vertices(g)];

	for (int i = 0; i < num_vertices(g); ++i) {
		state[i].rr_node = -1;
		state[i].known_cost = std::numeric_limits<float>::max();
		state[i].cost = std::numeric_limits<float>::max();
		state[i].prev_edge = nullptr;
		state[i].upstream_R = -1;
		state[i].delay = std::numeric_limits<float>::max();
	}


	/*route_tree_t **route_trees = new route_tree_t *[nets.size()+global_nets.size()];*/
	/*for (int i = 0; i < nets.size(); ++i) {*/
		/*route_trees[nets[i].id] = new route_tree_t[nets[i].sinks.size()+1];*/
		/*for (int j = 0; j <= nets[i].sinks.size(); ++j) {*/
			/*route_tree_init(route_trees[nets[i].id][j]);*/
		/*}*/
	/*}*/
	/*for (int i = 0; i < global_nets.size(); ++i) {*/
		/*route_trees[nets[i].id] = new route_tree_t[nets[i].sinks.size()+1];*/
		/*for (int j = 0; j <= nets[i].sinks.size(); ++j) {*/
			/*route_tree_init(route_trees[nets[i].id][j]);*/
		/*}*/
	/*}*/

	route_parameters_t params;
	params.criticality_exp = opts->criticality_exp;
	params.astar_fac = opts->astar_fac;
	params.max_criticality = opts->max_criticality;
	params.bend_cost = opts->bend_cost;
	params.pres_fac = opts->first_iter_pres_fac; /* Typically 0 -> ignore cong. */

	/*int num_virtual_nets = 0;*/
	/*for (auto &net : nets) {*/
		/*for (auto &sink : net.sinks) {*/
			/*++num_virtual_nets;*/
		/*}*/
	/*}*/
	/*vector<virtual_net_t> virtual_nets(num_virtual_nets);*/
	/*num_virtual_nets = 0;*/
	/*for (auto &net : nets) {*/
		/*for (auto &sink : net.sinks) {*/
			/*virtual_nets[num_virtual_nets].valid = true;*/
			/*virtual_nets[num_virtual_nets].source = &net.source;*/
			/*virtual_nets[num_virtual_nets].sinks.push_back(&sink);*/
			/*++num_virtual_nets;*/
		/*}*/
	/*}*/
	/*sort(begin(virtual_nets), end(virtual_nets), [] (const virtual_net_t &a, const virtual_net_t &b) -> bool {*/
			/*return get_bounding_box_area(a.sinks[0]->current_bounding_box) > get_bounding_box_area(b.sinks[0]->current_bounding_box);*/
			/*});*/
	/*QuadTree<virtual_net_t *> qt(fpga_bb, 4);*/
	/*for (auto &vnet : virtual_nets) {*/
		/*qt.insert(make_pair(vnet.sinks[0]->current_bounding_box, &vnet));*/
	/*}*/
	/*zlog_info(delta_log, "Num levels in quadtree for all virtual nets: %d\n", qt.num_levels());*/
	/*qt.print_num_items(0);*/
	unsigned int max_num_sinks = 0;
	for (const auto &net : nets) {
		if (net.sinks.size() > max_num_sinks) {
			max_num_sinks = net.sinks.size();
		}
	}
	printf("max_num_sinks: %u\n", max_num_sinks);
	vector<vector<virtual_net_t>> virtual_nets_by_net;
	create_clustered_virtual_nets(nets, 5, opts->max_sink_bb_area, virtual_nets_by_net);

	/*print_virtual_nets(nets);*/

	vector<std::mutex> net_locks(nets.size());

	for_all_vertices(g, [] (RRNode &v) -> void {
			v.properties.acc_cost = 1;
			v.properties.pres_cost = 1;
			v.properties.occ = 0;
			});

	char buffer[256];

	if (current_output_log) {
		fclose(current_output_log);
	}

	std::mt19937 mt(time(NULL));

	bool routed = false;
	nanosecond_type total_route_time = 0;
	for (int iter = 0; iter < opts->max_router_iterations && !routed; ++iter) {
		sprintf(buffer, "/Volumes/DATA/iter_%d.txt", iter);
		current_output_log = fopen(buffer, "w");

		zlog_info(delta_log, "Routing iteration: %d\n", iter);

		/*perf.num_heap_pushes = 0;*/

		cpu_timer iter_timer;
		iter_timer.start();

		for (auto &net : nets) {
			update_sink_criticalities(net, net_timing[net.vpr_id], params);
		}

		for (auto &net : nets) {
			net.num_bounding_box_updates = 0;
			net.num_nearest_iters = 0;
			net.total_point_tree_size = 0;
		}

		/*for (auto &net : nets) {*/
			/*[>reset_current_source_sink_3(net);<]*/
			/*fill(begin(net.sink_routed), end(net.sink_routed), false);*/
		/*}*/
		/*for (auto &virtual_nets : virtual_nets_by_net) {*/
			/*for (auto &vnet : virtual_nets) {*/
				/*vnet.routed = false;*/
			/*}*/
		/*}*/
		/*num_routed_virtual_nets = 0;*/
		nanosecond_type route_time = 0;

		cpu_timer route_timer;

		route_timer.start();
#ifdef __linux__
		__itt_frame_begin_v3(pD, NULL);
							/*__itt_task_begin(pD, __itt_null, __itt_null, shMainTask);*/
#endif

		pending_rtree_t pending_rtree = build_pending_rtree(nets, g, params, route_trees, nullptr);

		vector<virtual_net_t *> in_flight_virtual_nets;

		sch_state_lock_t sched_state_lock;
		perf_t perf;
		sched_perf_t sched_perf;
		tbb::concurrent_bounded_queue<virtual_net_t *> update_requests;
		int num_virtual_nets = pending_rtree.size();
		tbb::atomic<int> num_dispatched_virtual_nets = 0;

		RouteWorker2 worker(
				opts,
				params,
				net_locks,
				g,
				state,
				route_trees,
				net_timing,
				sched_state_lock,
				pending_rtree,
				in_flight_virtual_nets,
				update_requests,
				&perf,
				&sched_perf,
				num_dispatched_virtual_nets);

		vector<virtual_net_t *> bootstrap;
		dispatch_virtual_nets(pending_rtree, in_flight_virtual_nets, nullptr,
				[&bootstrap] (virtual_net_t *virtual_net) -> void {
				bootstrap.push_back(virtual_net);
				}, num_dispatched_virtual_nets, &sched_perf);

		tbb::parallel_do(bootstrap, worker);

#ifdef __linux__
		__itt_frame_end_v3(pD, NULL);
							/*__itt_task_begin(pD, __itt_null, __itt_null, shMainTask);*/
#endif

		route_timer.stop();
		cpu_times elapsed = route_timer.elapsed();
		route_time = elapsed.wall;
		total_route_time += elapsed.wall;
		/*zlog_info(delta_log, "Scheduling took %s\n", format(elapsed).c_str());*/
		zlog_info(delta_log, "Routing took %g\n", elapsed.wall / 1e9);
		printf("Routing took %g\n", elapsed.wall / 1e9);
		/*zlog_info(delta_log, "Num heap pushes: %d\n", perf.num_heap_pushes);*/
		/*for (int i = 0; i < opts->num_threads; ++i) {*/
			/*printf("Thread %d num heap pushes: %lu\n", i, args[i]->perf.num_heap_pushes);*/
			/*printf("Thread %d num heap pops: %lu\n", i, args[i]->perf.num_heap_pops);*/
			/*printf("Thread %d num neighbor visits: %lu\n", i, args[i]->perf.num_neighbor_visits);*/
			/*printf("Thread %d total wait time: %g s\n", i, std::chrono::duration_cast<std::chrono::nanoseconds>(args[i]->perf.total_wait_time).count() / 1e9);*/
			/*printf("Thread %d total rip up time: %g s\n", i, std::chrono::duration_cast<std::chrono::nanoseconds>(args[i]->perf.total_rip_up_time).count() / 1e9);*/
			/*printf("Thread %d total route time: %g s\n", i, std::chrono::duration_cast<std::chrono::nanoseconds>(args[i]->perf.total_route_time).count() / 1e9);*/
			/*printf("Thread %d total update time: %g s\n", i, std::chrono::duration_cast<std::chrono::nanoseconds>(args[i]->perf.total_update_time).count() / 1e9);*/
			/*printf("\tThread %d total get nearest time: %g s (%g %%)\n", i, std::chrono::duration_cast<std::chrono::nanoseconds>(args[i]->perf.total_get_nearest_time).count() / 1e9, std::chrono::duration_cast<std::chrono::nanoseconds>(args[i]->perf.total_get_nearest_time).count()*100.0 / std::chrono::duration_cast<std::chrono::nanoseconds>(args[i]->perf.total_update_time).count());*/
			/*printf("Thread %d total push time: %g s\n", i, std::chrono::duration_cast<std::chrono::nanoseconds>(args[i]->perf.total_push_time).count() / 1e9);*/
		/*}*/
		/*printf("Scheduler total rtree build time: %g s\n", std::chrono::duration_cast<std::chrono::nanoseconds>(sargs.perf.total_build_time).count() / 1e9);*/
		/*printf("Scheduler total dispatch time: %g s\n", std::chrono::duration_cast<std::chrono::nanoseconds>(sargs.perf.total_dispatch_time).count() / 1e9);*/
		/*printf("Scheduler total rtree update time: %g s\n", std::chrono::duration_cast<std::chrono::nanoseconds>(sargs.perf.total_rtree_update_time).count() / 1e9);*/
		/*printf("Scheduler total wait time: %g s\n", std::chrono::duration_cast<std::chrono::nanoseconds>(sargs.perf.total_wait_time).count() / 1e9);*/
		/*printf("Scheduler num updates: %lu\n", sargs.perf.num_updates);*/
		/*printf("Scheduler num leaf pred: %lu\n", sargs.perf.num_leaf_node_pred_calls);*/
		/*printf("Scheduler num internal pred: %lu\n", sargs.perf.num_internal_node_pred_calls);*/

		printf("-- Sorted by num_nearest_iters --\n");
		vector<pair<unsigned long, const net_t *>> net_bb_updates;
		for (const auto &net : nets) {
			net_bb_updates.push_back(make_pair(net.num_nearest_iters, &net));
		}
		std::sort(begin(net_bb_updates), end(net_bb_updates), std::greater<pair<int, const net_t *>>());
		for (int i = 0; i < std::min(nets.size(), (unsigned long)10); ++i) {
			printf("Net %d Sinks: %lu Virtual nets: %lu BB: %d BB updates: %lu Nearest iters: %lu Point tree size: %lu\n",
					net_bb_updates[i].second->vpr_id,
					net_bb_updates[i].second->sinks.size(),
					net_bb_updates[i].second->virtual_nets.size(),
					net_bb_updates[i].second->bb_area_rank,
					net_bb_updates[i].second->num_bounding_box_updates,
					net_bb_updates[i].second->num_nearest_iters,
					net_bb_updates[i].second->total_point_tree_size);
		}

		printf("-- Sorted by num_bounding_box_updates --\n");
		net_bb_updates.clear();
		for (const auto &net : nets) {
			net_bb_updates.push_back(make_pair(net.num_bounding_box_updates, &net));
		}
		std::sort(begin(net_bb_updates), end(net_bb_updates), std::greater<pair<int, const net_t *>>());
		for (int i = 0; i < std::min(nets.size(), (unsigned long)10); ++i) {
			printf("Net %d Sinks: %lu Virtual nets: %lu BB: %d BB updates: %lu Nearest iters: %lu Point tree size: %lu\n",
					net_bb_updates[i].second->vpr_id,
					net_bb_updates[i].second->sinks.size(),
					net_bb_updates[i].second->virtual_nets.size(),
					net_bb_updates[i].second->bb_area_rank,
					net_bb_updates[i].second->num_bounding_box_updates,
					net_bb_updates[i].second->num_nearest_iters,
					net_bb_updates[i].second->total_point_tree_size);
		}

		printf("-- Sorted by total_point_tree_size --\n");
		net_bb_updates.clear();
		for (const auto &net : nets) {
			net_bb_updates.push_back(make_pair(net.total_point_tree_size, &net));
		}
		std::sort(begin(net_bb_updates), end(net_bb_updates), std::greater<pair<int, const net_t *>>());
		for (int i = 0; i < std::min(nets.size(), (unsigned long)10); ++i) {
			printf("Net %d Sinks: %lu Virtual nets: %lu BB: %d BB updates: %lu Nearest iters: %lu Point tree size: %lu\n",
					net_bb_updates[i].second->vpr_id,
					net_bb_updates[i].second->sinks.size(),
					net_bb_updates[i].second->virtual_nets.size(),
					net_bb_updates[i].second->bb_area_rank,
					net_bb_updates[i].second->num_bounding_box_updates,
					net_bb_updates[i].second->num_nearest_iters,
					net_bb_updates[i].second->total_point_tree_size);
		}

		/* checking */
		if (!pending_rtree.empty()) {
			zlog_error(delta_log, "There are still %lu in pending rtree\n", pending_rtree.size());
			assert(false);
		}
		if (!in_flight_virtual_nets.empty()) {
			zlog_error(delta_log, "There are still %lu in flight virtual nets\n", in_flight_virtual_nets.size());
			assert(false);
		}

		for (auto &virtual_nets : virtual_nets_by_net) {
			for (auto &vnet : virtual_nets) {
				assert(vnet.routed);
			}
		}

		for_all_vertices(g, [] (RRNode &v) -> void { v.properties.recalc_occ = 0; });

		for (const auto &net : nets) {
			check_route_tree(route_trees[net.local_id], net, g);
		}

		bool valid = true;
		for_all_vertices(g, [&valid, &buffer] (RRNode &v) -> void {
				sprintf_rr_node(id(v), buffer);
				if (v.properties.recalc_occ != v.properties.occ) {
					zlog_error(delta_log, "Node %s occ mismatch, recalc: %d original: %d\n", buffer, v.properties.recalc_occ, v.properties.occ);
					valid = false;
				}
				});
		assert(valid);

		for (const auto &net : nets) {
			vector<int> overused_rr_node;
			get_overused_nodes(route_trees[net.local_id], get_vertex(route_trees[net.local_id].graph, route_trees[net.local_id].root_rt_node_id), g, overused_rr_node);
			if (!overused_rr_node.empty()) {
				zlog_level(delta_log, ROUTER_V1, "Net %d bb_rank %d overused nodes:\n", net.vpr_id, net.bb_area_rank);
				for (const auto &item : overused_rr_node) {
					zlog_level(delta_log, ROUTER_V1, "%d ", item);
				}
				zlog_level(delta_log, ROUTER_V1, "\n");
			}
		}

		/*char filename[256];*/
		/*sprintf(filename, "congestion_%d.txt", iter);*/
		/*dump_congestion_map(g, filename);*/

		if (feasible_routing(g)) {
			zlog_info(delta_log, "Routed in %d iterations. Total route time: %g\n", iter+1, total_route_time / 1e9);
			printf("Routed in %d iterations. Total route time: %g\n", iter+1, total_route_time / 1e9);
			dump_route(*current_traces_ptr, "route.txt");
			routed = true;
		} else {
			int num_overused_nodes = 0;
			for_all_vertices(g, [&num_overused_nodes] (const RRNode &v) -> void {
					if (v.properties.occ > v.properties.capacity) {
					++num_overused_nodes;
					}
					});
			zlog_info(delta_log, "Num overused nodes: %d/%d (%.2f)\n", num_overused_nodes, num_vertices(g), num_overused_nodes*100.0/num_vertices(g));

			if (iter == 0) {
				params.pres_fac = opts->initial_pres_fac;
				update_costs(g, params.pres_fac, 0);
			} else {
				params.pres_fac *= opts->pres_fac_mult;

				/* Avoid overflow for high iteration counts, even if acc_cost is big */
				params.pres_fac = std::min(params.pres_fac, static_cast<float>(HUGE_POSITIVE_FLOAT / 1e5));

				update_costs(g, params.pres_fac, opts->acc_fac);
			}

			for (auto &net : nets) {
				/*zlog_level(delta_log, ROUTER_V1, "Net %d trace stats:\n", net.vpr_id);*/
				/*adjust_bb_factor((*current_traces_ptr)[net.local_id], net, g, opts->bb_expand_threshold);*/
			}

			for (auto &net : nets) {
				overused_stats((*current_traces_ptr)[net.local_id], route_trees[net.local_id], net, g);
			}

			/* swap the pointers */
			std::swap(current_traces_ptr, prev_traces_ptr);
			for (const auto &trace : *current_traces_ptr) {
				assert(trace_empty(trace));
			}
		}

		analyze_timing(net_timing);

		iter_timer.stop();
		elapsed = iter_timer.elapsed();

		zlog_info(delta_log, "Iteration time: %g s. Route time: %g s. \n", elapsed.wall / 1e9, route_time / 1e9);

		if (current_output_log && fclose(current_output_log)) {
			char str[256];
			strerror_r(errno, str, sizeof(str));
			printf("failed to close file: %s\n", str);
			assert(false);
		}
	}

	if (!routed) {
		printf("Failed to route in %d iterations. Total route time: %g\n", opts->max_router_iterations, total_route_time / 1e9);
	}

	return routed;
}

bool greedy_route(t_router_opts *opts)
{
	/*tbb::task_scheduler_init init(opts->num_threads);*/

	start_time = std::chrono::high_resolution_clock::now();

/*#ifdef __linux__*/
	/*printf("domain flags: %d\n", pD->flags);*/
	/*pD->flags = 1;*/
/*#endif*/
	
	init_logging();
	zlog_set_record("custom_output", zlog_custom_output);
	zlog_set_record("sched_custom_output", zlog_sched_custom_output);
	current_output_log = fopen("/Volumes/DATA/main.log", "w");

	zlog_info(delta_log, "Using boost version %d\n", BOOST_VERSION);

#if (TBB_USE_DEBUG==1)
	zlog_warn(delta_log, "Using debug version of TBB\n");
#endif

	/*test_rtree();*/

	/*test_route_tree();*/

	/*test_geometry();*/

	/*test_clustering(opts->num_threads, opts->bb_factor);*/
	/*return;*/

	/*test_perf();*/
	/*exit(-1);*/

	/*test_quadrant();*/

	/*auto time = std::chrono::high_resolution_clock::now();*/

	/*printf("test\n");*/

	/*auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::high_resolution_clock::now()-time);*/
	/*printf("printf took %g\n", elapsed.count());*/

	/*test_heap();*/

	/*exit(0);*/

	/*test_scheduler();*/

	/*test_rtree();*/

	/*test_graph();*/

	RRGraph g;
	init_graph(g);

	vector<net_t> nets;
	vector<net_t> global_nets;
	init_nets(nets, global_nets, opts->bb_factor);

	/*test_effectiveness(nets);*/

	/*exit(0);*/

	/*sort_sinks(nets);*/

	/*dump_net_bounding_boxes(nets);*/

	/* calculating net bounding box area rank */
	vector<net_t *> nets_ptr(nets.size());
	for (int i = 0; i < nets.size(); ++i) {
		nets_ptr[i] = &nets[i];
	}
	extern s_bb *route_bb;
	sort(nets_ptr.begin(), nets_ptr.end(), [] (const net_t *a, const net_t *b) -> bool {
			return get_bounding_box_area(route_bb[a->vpr_id]) > get_bounding_box_area(route_bb[b->vpr_id]);
			});
	int rank = 0;
	for (auto &net : nets_ptr) {
		net->bb_area_rank = rank++;
	}
	/* end */ 

	extern int nx;
	extern int ny;
	zlog_info(delta_log, "Num nets: %lu Num global nets: %lu nx: %d ny: %d\n", nets.size(), global_nets.size(), nx, ny);
	
	t_net_timing *net_timing = new t_net_timing[nets.size()+global_nets.size()];
	init_net_timing(nets, global_nets, net_timing);

	vector<route_tree_t> route_trees(nets.size());
	for (int i = 0; i < nets.size(); ++i) {
		route_tree_init(route_trees[i]);
	}

	vector<trace_t> prev_traces(nets.size());
	vector<trace_t> traces(nets.size());
	for (int i = 0; i < nets.size(); ++i) {
		trace_init(traces[i]);
		trace_init(prev_traces[i]);
	}

	vector<trace_t> *prev_traces_ptr = &prev_traces;
	vector<trace_t> *current_traces_ptr = &traces;

	route_state_t *state = new route_state_t[num_vertices(g)];

	for (int i = 0; i < num_vertices(g); ++i) {
		state[i].rr_node = -1;
		state[i].known_cost = std::numeric_limits<float>::max();
		state[i].cost = std::numeric_limits<float>::max();
		state[i].prev_edge = nullptr;
		state[i].upstream_R = -1;
		state[i].delay = std::numeric_limits<float>::max();
	}


	/*route_tree_t **route_trees = new route_tree_t *[nets.size()+global_nets.size()];*/
	/*for (int i = 0; i < nets.size(); ++i) {*/
		/*route_trees[nets[i].id] = new route_tree_t[nets[i].sinks.size()+1];*/
		/*for (int j = 0; j <= nets[i].sinks.size(); ++j) {*/
			/*route_tree_init(route_trees[nets[i].id][j]);*/
		/*}*/
	/*}*/
	/*for (int i = 0; i < global_nets.size(); ++i) {*/
		/*route_trees[nets[i].id] = new route_tree_t[nets[i].sinks.size()+1];*/
		/*for (int j = 0; j <= nets[i].sinks.size(); ++j) {*/
			/*route_tree_init(route_trees[nets[i].id][j]);*/
		/*}*/
	/*}*/

	route_parameters_t params;
	params.criticality_exp = opts->criticality_exp;
	params.astar_fac = opts->astar_fac;
	params.max_criticality = opts->max_criticality;
	params.bend_cost = opts->bend_cost;
	params.pres_fac = opts->first_iter_pres_fac; /* Typically 0 -> ignore cong. */

	std::mutex lock;
	extern int nx;
	extern int ny;
	bounding_box_t fpga_bb;
	fpga_bb.xmin = 0;
	fpga_bb.ymin = 0;
	fpga_bb.xmax = nx+2;
	fpga_bb.ymax = ny+2;
	QuadTree<virtual_net_t *> in_flight_qt(fpga_bb);

	/*int num_virtual_nets = 0;*/
	/*for (auto &net : nets) {*/
		/*for (auto &sink : net.sinks) {*/
			/*++num_virtual_nets;*/
		/*}*/
	/*}*/
	/*vector<virtual_net_t> virtual_nets(num_virtual_nets);*/
	/*num_virtual_nets = 0;*/
	/*for (auto &net : nets) {*/
		/*for (auto &sink : net.sinks) {*/
			/*virtual_nets[num_virtual_nets].valid = true;*/
			/*virtual_nets[num_virtual_nets].source = &net.source;*/
			/*virtual_nets[num_virtual_nets].sinks.push_back(&sink);*/
			/*++num_virtual_nets;*/
		/*}*/
	/*}*/
	/*sort(begin(virtual_nets), end(virtual_nets), [] (const virtual_net_t &a, const virtual_net_t &b) -> bool {*/
			/*return get_bounding_box_area(a.sinks[0]->current_bounding_box) > get_bounding_box_area(b.sinks[0]->current_bounding_box);*/
			/*});*/
	/*QuadTree<virtual_net_t *> qt(fpga_bb, 4);*/
	/*for (auto &vnet : virtual_nets) {*/
		/*qt.insert(make_pair(vnet.sinks[0]->current_bounding_box, &vnet));*/
	/*}*/
	/*zlog_info(delta_log, "Num levels in quadtree for all virtual nets: %d\n", qt.num_levels());*/
	/*qt.print_num_items(0);*/
	unsigned int max_num_sinks = 0;
	for (const auto &net : nets) {
		if (net.sinks.size() > max_num_sinks) {
			max_num_sinks = net.sinks.size();
		}
	}
	printf("max_num_sinks: %u\n", max_num_sinks);
	vector<vector<virtual_net_t>> virtual_nets_by_net;
	create_clustered_virtual_nets(nets, 5, opts->max_sink_bb_area, virtual_nets_by_net);

	/*print_virtual_nets(nets);*/

	int num_routed_virtual_nets;

	pthread_cond_t start_cond;
	pthread_mutex_t start_cond_lock;
	my_pthread_barrier_t barrier;

	/*assert(!pthread_mutex_init(&start_cond_lock, 0));*/
	/*assert(!pthread_cond_init(&start_cond, 0));*/
	/*my_pthread_barrier_init(&barrier, NULL, opts->num_threads);*/

	/*perf_t perf;*/

	vector<std::mutex> net_locks(nets.size());

	WorkerArgs **args = new WorkerArgs*[opts->num_threads];

#ifdef __linux__
	/*sem_t produce_sem_instance;*/
	/*sem_t consume_sem_instance;*/
	/*assert(!sem_init(&produce_sem_instance, 0, 0));*/
	/*assert(!sem_init(&consume_sem_instance, 0, 0));*/
	/*sem_t *produce_sem = &produce_sem_instance;*/
	/*sem_t *consume_sem = &consume_sem_instance;*/
	sem_t *produce_sem = nullptr;
	sem_t *consume_sem = nullptr;
	printf("Using sem_init\n");
#else
	/*sem_t *produce_sem = create_sem("produce_sem");*/
	/*sem_t *consume_sem = create_sem("consume_sem");*/
	sem_t *produce_sem = nullptr;
	sem_t *consume_sem = nullptr;
	printf("Using sem_open\n");
#endif

	tbb::concurrent_bounded_queue<virtual_net_t *> pending_virtual_nets;
	tbb::concurrent_bounded_queue<virtual_net_t *> routed_virtual_nets;

	/* start the threads first so that they have ample time to finish reach the blocking cond wait part */
	vector<pthread_t> threads(opts->num_threads);
	for (int i = 0; i < opts->num_threads; ++i) {
		args[i] = new WorkerArgs(
				i,
				opts->num_threads,

				&start_cond_lock,
				&start_cond,
				&barrier,

				lock,
				in_flight_qt,
				/*nullptr,*/
				num_routed_virtual_nets,

				net_locks,

				g,
				params,
				state,
				route_trees,
				&prev_traces_ptr,
				&current_traces_ptr,
				net_timing,

				consume_sem,
				produce_sem,
				
				&pending_virtual_nets,
				&routed_virtual_nets
				);

		args[i]->virtual_nets_by_net = &virtual_nets_by_net;
		pthread_create(&threads[i], NULL, worker_thread_3, args[i]);
	}

	SchedulerArgs sargs;
	sargs.virtual_nets = nullptr;
	sargs.virtual_nets_by_net = &virtual_nets_by_net;
	sargs.consume_sem = consume_sem;
	sargs.produce_sem = produce_sem;
	sargs.pending_virtual_nets = &pending_virtual_nets;
	sargs.routed_virtual_nets = &routed_virtual_nets;
	sargs.in_flight_qt = &in_flight_qt;
	sargs.route_trees = &route_trees;
	sargs.g = &g;
	sargs.params = &params;
	sargs.num_nets = nets.size();
	sargs.nets = &nets;
	sargs.net_locks = &net_locks;
	sargs.opts = opts;

	for_all_vertices(g, [] (RRNode &v) -> void {
			v.properties.acc_cost = 1;
			v.properties.pres_cost = 1;
			v.properties.occ = 0;
			});

	char buffer[256];

	if (current_output_log) {
		fclose(current_output_log);
	}

	std::mt19937 mt(time(NULL));

	bool routed = false;
	nanosecond_type total_route_time = 0;
	for (int iter = 0; iter < opts->max_router_iterations && !routed; ++iter) {
		sprintf(buffer, "/Volumes/DATA/iter_%d.txt", iter);
		current_output_log = fopen(buffer, "w");

		zlog_info(delta_log, "Routing iteration: %d\n", iter);

		for (int i = 0; i < opts->num_threads; ++i) {
			args[i]->iteration = &iter;

			args[i]->perf.num_heap_pushes = 0;
			args[i]->perf.num_heap_pops = 0;
			args[i]->perf.num_neighbor_visits = 0;
			args[i]->perf.total_wait_time = std::chrono::high_resolution_clock::duration::zero();
			args[i]->perf.total_rip_up_time = std::chrono::high_resolution_clock::duration::zero();
			args[i]->perf.total_route_time = std::chrono::high_resolution_clock::duration::zero();
			args[i]->perf.total_update_time = std::chrono::high_resolution_clock::duration::zero();
			args[i]->perf.total_push_time = std::chrono::high_resolution_clock::duration::zero();
			args[i]->perf.total_centroid_time = std::chrono::high_resolution_clock::duration::zero();
			args[i]->perf.total_get_nearest_time = std::chrono::high_resolution_clock::duration::zero();
			args[i]->perf.total_verificaion_time = std::chrono::high_resolution_clock::duration::zero();
			args[i]->perf.total_expansion_time = std::chrono::high_resolution_clock::duration::zero();
			args[i]->perf.total_scheduler_box_time = std::chrono::high_resolution_clock::duration::zero();
		}
		sargs.iteration = &iter;

		/*perf.num_heap_pushes = 0;*/

		cpu_timer iter_timer;
		iter_timer.start();

		for (auto &net : nets) {
			update_sink_criticalities(net, net_timing[net.vpr_id], params);
		}

		for (auto &net : nets) {
			net.num_bounding_box_updates = 0;
			net.num_nearest_iters = 0;
			net.total_point_tree_size = 0;
		}

		/*for (auto &net : nets) {*/
			/*[>reset_current_source_sink_3(net);<]*/
			/*fill(begin(net.sink_routed), end(net.sink_routed), false);*/
		/*}*/
		/*for (auto &virtual_nets : virtual_nets_by_net) {*/
			/*for (auto &vnet : virtual_nets) {*/
				/*vnet.routed = false;*/
			/*}*/
		/*}*/
		/*num_routed_virtual_nets = 0;*/
		if (!in_flight_qt.empty()) {
			in_flight_qt.print_items([] (const pair<bounding_box_t, virtual_net_t *> &item) -> void {
					printf("Net %d sink %d in in_flight_qt\n", item.second->sinks[0]->net->vpr_id, item.second->sinks[0]->id);
					});
			assert(false);
		}

		nanosecond_type route_time = 0;

		cpu_timer route_timer;

		route_timer.start();

		/*std::shuffle(begin(virtual_nets), end(virtual_nets), mt);*/
		scheduler_thread_2(&sargs);

		route_timer.stop();
		cpu_times elapsed = route_timer.elapsed();
		route_time = elapsed.wall;
		total_route_time += elapsed.wall;
		/*zlog_info(delta_log, "Scheduling took %s\n", format(elapsed).c_str());*/
		zlog_info(delta_log, "Routing took %g\n", elapsed.wall / 1e9);
		printf("Routing took %g\n", elapsed.wall / 1e9);
		/*zlog_info(delta_log, "Num heap pushes: %d\n", perf.num_heap_pushes);*/
		for (int i = 0; i < opts->num_threads; ++i) {
			printf("Thread %d num heap pushes: %lu\n", i, args[i]->perf.num_heap_pushes);
			printf("Thread %d num heap pops: %lu\n", i, args[i]->perf.num_heap_pops);
			printf("Thread %d num neighbor visits: %lu\n", i, args[i]->perf.num_neighbor_visits);
			printf("Thread %d total wait time: %g s\n", i, std::chrono::duration_cast<std::chrono::nanoseconds>(args[i]->perf.total_wait_time).count() / 1e9);
			printf("Thread %d total rip up time: %g s\n", i, std::chrono::duration_cast<std::chrono::nanoseconds>(args[i]->perf.total_rip_up_time).count() / 1e9);
			printf("Thread %d total route time: %g s\n", i, std::chrono::duration_cast<std::chrono::nanoseconds>(args[i]->perf.total_route_time).count() / 1e9);
			printf("Thread %d total update time: %g s\n", i, std::chrono::duration_cast<std::chrono::nanoseconds>(args[i]->perf.total_update_time).count() / 1e9);
			printf("\tThread %d total centroid time: %g s (%g %%)\n", i, std::chrono::duration_cast<std::chrono::nanoseconds>(args[i]->perf.total_centroid_time).count() / 1e9, std::chrono::duration_cast<std::chrono::nanoseconds>(args[i]->perf.total_centroid_time).count()*100.0 / std::chrono::duration_cast<std::chrono::nanoseconds>(args[i]->perf.total_update_time).count());
			printf("\tThread %d total get nearest time: %g s (%g %%)\n", i, std::chrono::duration_cast<std::chrono::nanoseconds>(args[i]->perf.total_get_nearest_time).count() / 1e9, std::chrono::duration_cast<std::chrono::nanoseconds>(args[i]->perf.total_get_nearest_time).count()*100.0 / std::chrono::duration_cast<std::chrono::nanoseconds>(args[i]->perf.total_update_time).count());
			printf("\tThread %d total verification time: %g s (%g %%)\n", i, std::chrono::duration_cast<std::chrono::nanoseconds>(args[i]->perf.total_verificaion_time).count() / 1e9, std::chrono::duration_cast<std::chrono::nanoseconds>(args[i]->perf.total_verificaion_time).count()*100.0 / std::chrono::duration_cast<std::chrono::nanoseconds>(args[i]->perf.total_update_time).count());
			printf("\tThread %d total expansion time: %g s (%g %%)\n", i, std::chrono::duration_cast<std::chrono::nanoseconds>(args[i]->perf.total_expansion_time).count() / 1e9, std::chrono::duration_cast<std::chrono::nanoseconds>(args[i]->perf.total_expansion_time).count()*100.0 / std::chrono::duration_cast<std::chrono::nanoseconds>(args[i]->perf.total_update_time).count());
			printf("\tThread %d total scheduler box time: %g s (%g %%)\n", i, std::chrono::duration_cast<std::chrono::nanoseconds>(args[i]->perf.total_scheduler_box_time).count() / 1e9, std::chrono::duration_cast<std::chrono::nanoseconds>(args[i]->perf.total_scheduler_box_time).count()*100.0 / std::chrono::duration_cast<std::chrono::nanoseconds>(args[i]->perf.total_update_time).count());
			printf("Thread %d total push time: %g s\n", i, std::chrono::duration_cast<std::chrono::nanoseconds>(args[i]->perf.total_push_time).count() / 1e9);
		}
		printf("Scheduler total rtree build time: %g s\n", std::chrono::duration_cast<std::chrono::nanoseconds>(sargs.perf.total_build_time).count() / 1e9);
		printf("Scheduler total dispatch time: %g s\n", std::chrono::duration_cast<std::chrono::nanoseconds>(sargs.perf.total_dispatch_time).count() / 1e9);
		printf("Scheduler total rtree update time: %g s\n", std::chrono::duration_cast<std::chrono::nanoseconds>(sargs.perf.total_rtree_update_time).count() / 1e9);
		printf("Scheduler total wait time: %g s\n", std::chrono::duration_cast<std::chrono::nanoseconds>(sargs.perf.total_wait_time).count() / 1e9);
		printf("Scheduler num updates: %lu\n", sargs.perf.num_updates);
		printf("Scheduler num leaf pred: %lu\n", sargs.perf.num_leaf_node_pred_calls);
		printf("Scheduler num internal pred: %lu\n", sargs.perf.num_internal_node_pred_calls);

		printf("-- Sorted by num_nearest_iters --\n");
		vector<pair<unsigned long, const net_t *>> net_bb_updates;
		for (const auto &net : nets) {
			net_bb_updates.push_back(make_pair(net.num_nearest_iters, &net));
		}
		std::sort(begin(net_bb_updates), end(net_bb_updates), std::greater<pair<int, const net_t *>>());
		for (int i = 0; i < std::min(nets.size(), (unsigned long)10); ++i) {
			printf("Net %d Sinks: %lu Virtual nets: %lu BB: %d BB updates: %lu Nearest iters: %lu Point tree size: %lu\n",
					net_bb_updates[i].second->vpr_id,
					net_bb_updates[i].second->sinks.size(),
					net_bb_updates[i].second->virtual_nets.size(),
					net_bb_updates[i].second->bb_area_rank,
					net_bb_updates[i].second->num_bounding_box_updates,
					net_bb_updates[i].second->num_nearest_iters,
					net_bb_updates[i].second->total_point_tree_size);
		}

		printf("-- Sorted by num_bounding_box_updates --\n");
		net_bb_updates.clear();
		for (const auto &net : nets) {
			net_bb_updates.push_back(make_pair(net.num_bounding_box_updates, &net));
		}
		std::sort(begin(net_bb_updates), end(net_bb_updates), std::greater<pair<int, const net_t *>>());
		for (int i = 0; i < std::min(nets.size(), (unsigned long)10); ++i) {
			printf("Net %d Sinks: %lu Virtual nets: %lu BB: %d BB updates: %lu Nearest iters: %lu Point tree size: %lu\n",
					net_bb_updates[i].second->vpr_id,
					net_bb_updates[i].second->sinks.size(),
					net_bb_updates[i].second->virtual_nets.size(),
					net_bb_updates[i].second->bb_area_rank,
					net_bb_updates[i].second->num_bounding_box_updates,
					net_bb_updates[i].second->num_nearest_iters,
					net_bb_updates[i].second->total_point_tree_size);
		}

		printf("-- Sorted by total_point_tree_size --\n");
		net_bb_updates.clear();
		for (const auto &net : nets) {
			net_bb_updates.push_back(make_pair(net.total_point_tree_size, &net));
		}
		std::sort(begin(net_bb_updates), end(net_bb_updates), std::greater<pair<int, const net_t *>>());
		for (int i = 0; i < std::min(nets.size(), (unsigned long)10); ++i) {
			printf("Net %d Sinks: %lu Virtual nets: %lu BB: %d BB updates: %lu Nearest iters: %lu Point tree size: %lu\n",
					net_bb_updates[i].second->vpr_id,
					net_bb_updates[i].second->sinks.size(),
					net_bb_updates[i].second->virtual_nets.size(),
					net_bb_updates[i].second->bb_area_rank,
					net_bb_updates[i].second->num_bounding_box_updates,
					net_bb_updates[i].second->num_nearest_iters,
					net_bb_updates[i].second->total_point_tree_size);
		}

		/* checking */
		assert(pending_virtual_nets.empty());
		assert(routed_virtual_nets.empty());

		for (auto &virtual_nets : virtual_nets_by_net) {
			for (auto &vnet : virtual_nets) {
				assert(vnet.routed);
			}
		}

		for_all_vertices(g, [] (RRNode &v) -> void { v.properties.recalc_occ = 0; });

		for (const auto &net : nets) {
			check_route_tree(route_trees[net.local_id], net, g);
		}

		bool valid = true;
		for_all_vertices(g, [&valid, &buffer] (RRNode &v) -> void {
				sprintf_rr_node(id(v), buffer);
				if (v.properties.recalc_occ != v.properties.occ) {
					zlog_error(delta_log, "Node %s occ mismatch, recalc: %d original: %d\n", buffer, v.properties.recalc_occ, v.properties.occ);
					valid = false;
				}
				});
		assert(valid);

		for (const auto &net : nets) {
			vector<int> overused_rr_node;
			get_overused_nodes(route_trees[net.local_id], get_vertex(route_trees[net.local_id].graph, route_trees[net.local_id].root_rt_node_id), g, overused_rr_node);
			if (!overused_rr_node.empty()) {
				zlog_level(delta_log, ROUTER_V1, "Net %d bb_rank %d overused nodes:\n", net.vpr_id, net.bb_area_rank);
				for (const auto &item : overused_rr_node) {
					zlog_level(delta_log, ROUTER_V1, "%d ", item);
				}
				zlog_level(delta_log, ROUTER_V1, "\n");
			}
		}

		/*char filename[256];*/
		/*sprintf(filename, "congestion_%d.txt", iter);*/
		/*dump_congestion_map(g, filename);*/

		if (feasible_routing(g)) {
			zlog_info(delta_log, "Routed in %d iterations. Total route time: %g\n", iter+1, total_route_time / 1e9);
			printf("Routed in %d iterations. Total route time: %g\n", iter+1, total_route_time / 1e9);
			dump_route(*current_traces_ptr, "route.txt");
			routed = true;
		} else {
			int num_overused_nodes = 0;
			for_all_vertices(g, [&num_overused_nodes] (const RRNode &v) -> void {
					if (v.properties.occ > v.properties.capacity) {
					++num_overused_nodes;
					}
					});
			zlog_info(delta_log, "Num overused nodes: %d/%d (%.2f)\n", num_overused_nodes, num_vertices(g), num_overused_nodes*100.0/num_vertices(g));

			if (iter == 0) {
				params.pres_fac = opts->initial_pres_fac;
				update_costs(g, params.pres_fac, 0);
			} else {
				params.pres_fac *= opts->pres_fac_mult;

				/* Avoid overflow for high iteration counts, even if acc_cost is big */
				params.pres_fac = std::min(params.pres_fac, static_cast<float>(HUGE_POSITIVE_FLOAT / 1e5));

				update_costs(g, params.pres_fac, opts->acc_fac);
			}

			for (auto &net : nets) {
				/*zlog_level(delta_log, ROUTER_V1, "Net %d trace stats:\n", net.vpr_id);*/
				/*adjust_bb_factor((*current_traces_ptr)[net.local_id], net, g, opts->bb_expand_threshold);*/
			}

			for (auto &net : nets) {
				overused_stats((*current_traces_ptr)[net.local_id], route_trees[net.local_id], net, g);
			}

			/* swap the pointers */
			std::swap(current_traces_ptr, prev_traces_ptr);
			for (const auto &trace : *current_traces_ptr) {
				assert(trace_empty(trace));
			}
		}

		analyze_timing(net_timing);

		iter_timer.stop();
		elapsed = iter_timer.elapsed();

		zlog_info(delta_log, "Iteration time: %g s. Route time: %g s. \n", elapsed.wall / 1e9, route_time / 1e9);

		if (current_output_log && fclose(current_output_log)) {
			char str[256];
			strerror_r(errno, str, sizeof(str));
			printf("failed to close file: %s\n", str);
			assert(false);
		}
	}

	if (!routed) {
		printf("Failed to route in %d iterations. Total route time: %g\n", opts->max_router_iterations, total_route_time / 1e9);
	}

	return routed;
}

bool _old_greedy_route_4(t_router_opts *opts)
{
	tbb::task_scheduler_init init(opts->num_threads);

/*#ifdef __linux__*/
	/*printf("domain flags: %d\n", pD->flags);*/
	/*pD->flags = 1;*/
/*#endif*/
	
	init_logging();
	zlog_set_record("custom_output", zlog_custom_output);
	current_output_log = fopen("/Volumes/DATA/main.log", "w");

	zlog_info(delta_log, "Using boost version %d\n", BOOST_VERSION);

#if (TBB_USE_DEBUG==1)
	zlog_warn(delta_log, "Using debug version of TBB\n");
#endif

	/*test_perf();*/
	/*exit(-1);*/

	/*test_quadrant();*/

	/*auto time = std::chrono::high_resolution_clock::now();*/

	/*printf("test\n");*/

	/*auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::high_resolution_clock::now()-time);*/
	/*printf("printf took %g\n", elapsed.count());*/

	/*test_heap();*/

	/*exit(0);*/

	/*test_scheduler();*/

	/*test_rtree();*/

	/*test_graph();*/

	RRGraph g;
	init_graph(g);

	vector<net_t> nets;
	vector<net_t> global_nets;
	init_nets(nets, global_nets, opts->bb_factor);

	/*test_effectiveness(nets);*/

	/*exit(0);*/

	sort_sinks(nets);

	/*dump_net_bounding_boxes(nets);*/

	/* calculating net bounding box area rank */
	vector<net_t *> nets_ptr(nets.size());
	for (int i = 0; i < nets.size(); ++i) {
		nets_ptr[i] = &nets[i];
	}
	extern s_bb *route_bb;
	sort(nets_ptr.begin(), nets_ptr.end(), [] (const net_t *a, const net_t *b) -> bool {
			return get_bounding_box_area(route_bb[a->vpr_id]) > get_bounding_box_area(route_bb[b->vpr_id]);
			});
	int rank = 0;
	for (auto &net : nets_ptr) {
		net->bb_area_rank = rank++;
	}
	/* end */ 

	extern int nx;
	extern int ny;
	zlog_info(delta_log, "Num nets: %lu Num global nets: %lu nx: %d ny: %d\n", nets.size(), global_nets.size(), nx, ny);
	
	t_net_timing *net_timing = new t_net_timing[nets.size()+global_nets.size()];
	init_net_timing(nets, global_nets, net_timing);

	vector<route_tree_t> route_trees(nets.size());
	for (int i = 0; i < nets.size(); ++i) {
		route_tree_init(route_trees[i]);
	}

	vector<trace_t> prev_traces(nets.size());
	vector<trace_t> traces(nets.size());
	for (int i = 0; i < nets.size(); ++i) {
		trace_init(traces[i]);
		trace_init(prev_traces[i]);
	}

	vector<trace_t> *prev_traces_ptr = &prev_traces;
	vector<trace_t> *current_traces_ptr = &traces;

	route_state_t *state = new route_state_t[num_vertices(g)];

	for (int i = 0; i < num_vertices(g); ++i) {
		state[i].rr_node = -1;
		state[i].known_cost = std::numeric_limits<float>::max();
		state[i].cost = std::numeric_limits<float>::max();
		state[i].prev_edge = nullptr;
		state[i].upstream_R = -1;
		state[i].delay = std::numeric_limits<float>::max();
	}


	/*route_tree_t **route_trees = new route_tree_t *[nets.size()+global_nets.size()];*/
	/*for (int i = 0; i < nets.size(); ++i) {*/
		/*route_trees[nets[i].id] = new route_tree_t[nets[i].sinks.size()+1];*/
		/*for (int j = 0; j <= nets[i].sinks.size(); ++j) {*/
			/*route_tree_init(route_trees[nets[i].id][j]);*/
		/*}*/
	/*}*/
	/*for (int i = 0; i < global_nets.size(); ++i) {*/
		/*route_trees[nets[i].id] = new route_tree_t[nets[i].sinks.size()+1];*/
		/*for (int j = 0; j <= nets[i].sinks.size(); ++j) {*/
			/*route_tree_init(route_trees[nets[i].id][j]);*/
		/*}*/
	/*}*/

	route_parameters_t params;
	params.criticality_exp = opts->criticality_exp;
	params.astar_fac = opts->astar_fac;
	params.max_criticality = opts->max_criticality;
	params.bend_cost = opts->bend_cost;
	params.pres_fac = opts->first_iter_pres_fac; /* Typically 0 -> ignore cong. */

	std::mutex lock;
	extern int nx;
	extern int ny;
	bounding_box_t fpga_bb;
	fpga_bb.xmin = 0;
	fpga_bb.ymin = 0;
	fpga_bb.xmax = nx+2;
	fpga_bb.ymax = ny+2;
	QuadTree<virtual_net_t *> in_flight_qt(fpga_bb);

	int num_virtual_nets = 0;
	for (auto &net : nets) {
		for (auto &sink : net.sinks) {
			++num_virtual_nets;
		}
	}
	vector<virtual_net_t> virtual_nets(num_virtual_nets);
	num_virtual_nets = 0;
	for (auto &net : nets) {
		for (auto &sink : net.sinks) {
			virtual_nets[num_virtual_nets].routed = false;
			virtual_nets[num_virtual_nets].source = &net.source;
			virtual_nets[num_virtual_nets].sinks.push_back(&sink);
			++num_virtual_nets;
		}
	}
	sort(begin(virtual_nets), end(virtual_nets), [] (const virtual_net_t &a, const virtual_net_t &b) -> bool {
			return get_bounding_box_area(a.sinks[0]->current_bounding_box) > get_bounding_box_area(b.sinks[0]->current_bounding_box);
			});
	QuadTree<virtual_net_t *> qt(fpga_bb, 4);
	for (auto &vnet : virtual_nets) {
		qt.insert(make_pair(vnet.sinks[0]->current_bounding_box, &vnet));
	}
	zlog_info(delta_log, "Num levels in quadtree for all virtual nets: %d\n", qt.num_levels());
	qt.print_num_items(0);

	int num_routed_virtual_nets;

	pthread_cond_t start_cond;
	pthread_mutex_t start_cond_lock;
	my_pthread_barrier_t barrier;

	assert(!pthread_mutex_init(&start_cond_lock, 0));
	assert(!pthread_cond_init(&start_cond, 0));
	my_pthread_barrier_init(&barrier, NULL, opts->num_threads);

	perf_t perf;

	vector<std::mutex> net_locks(nets.size());

	WorkerArgs **args = new WorkerArgs*[opts->num_threads];

	sem_t *produce_sem = create_sem("produce_sem");
	sem_t *consume_sem = create_sem("consume_sem");

	tbb::concurrent_bounded_queue<virtual_net_t *> pending_virtual_nets;
	tbb::concurrent_bounded_queue<virtual_net_t *> routed_virtual_nets;

	/* start the threads first so that they have ample time to finish reach the blocking cond wait part */
	vector<pthread_t> threads(opts->num_threads);
	for (int i = 0; i < opts->num_threads; ++i) {
		args[i] = new WorkerArgs(
				i,
				opts->num_threads,

				&start_cond_lock,
				&start_cond,
				&barrier,

				lock,
				in_flight_qt,
				/*virtual_nets,*/
				num_routed_virtual_nets,

				net_locks,

				g,
				params,
				state,
				route_trees,
				&prev_traces_ptr,
				&current_traces_ptr,
				net_timing,

				consume_sem,
				produce_sem,
				
				&pending_virtual_nets,
				&routed_virtual_nets
				);

		pthread_create(&threads[i], NULL, worker_thread_2, args[i]);
	}

	SchedulerArgs sargs;
	sargs.virtual_nets = &virtual_nets;
	sargs.consume_sem = consume_sem;
	sargs.produce_sem = produce_sem;
	sargs.pending_virtual_nets = &pending_virtual_nets;
	sargs.routed_virtual_nets = &routed_virtual_nets;
	sargs.in_flight_qt = &in_flight_qt;
	sargs.route_trees = &route_trees;
	sargs.g = &g;
	sargs.params = &params;

	for_all_vertices(g, [] (RRNode &v) -> void {
			v.properties.acc_cost = 1;
			v.properties.pres_cost = 1;
			v.properties.occ = 0;
			});

	char buffer[256];

	fclose(current_output_log);

	std::mt19937 mt(time(NULL));

	bool routed = false;
	for (int iter = 0; iter < opts->max_router_iterations && !routed; ++iter) {
		sprintf(buffer, "/Volumes/DATA/iter_%d.txt", iter);
		current_output_log = fopen(buffer, "w");
		assert(current_output_log);

		zlog_info(delta_log, "Routing iteration: %d\n", iter);

		for (int i = 0; i < opts->num_threads; ++i) {
			args[i]->iteration = &iter;
		}
		sargs.iteration = &iter;

		perf.num_heap_pushes = 0;

		cpu_timer iter_timer;
		iter_timer.start();

		for (int i = 0; i < nets.size(); ++i) {
			route_tree_clear(route_trees[i]);
		}

		for (auto &net : nets) {
			update_sink_criticalities(net, net_timing[net.vpr_id], params);
		}

		for (auto &net : nets) {
			reset_current_source_sink_3(net);
			fill(begin(net.sink_routed), end(net.sink_routed), false);
		}
		for (auto &vnet : virtual_nets) {
			vnet.routed = false;
		}
		num_routed_virtual_nets = 0;
		if (!in_flight_qt.empty()) {
			in_flight_qt.print_items([] (const pair<bounding_box_t, virtual_net_t *> &item) -> void {
					printf("Net %d sink %d in in_flight_qt\n", item.second->sinks[0]->net->vpr_id, item.second->sinks[0]->id);
					});
			assert(false);
		}

		nanosecond_type route_time = 0;

		cpu_timer route_timer;

		route_timer.start();

		/*std::shuffle(begin(virtual_nets), end(virtual_nets), mt);*/
		test_scheduler_thread_tmp(&sargs);

		route_timer.stop();
		cpu_times elapsed = route_timer.elapsed();
		route_time += elapsed.wall;
		/*zlog_info(delta_log, "Scheduling took %s\n", format(elapsed).c_str());*/
		zlog_info(delta_log, "Routing took %g\n", elapsed.wall / 1e9);
		zlog_info(delta_log, "Num heap pushes: %d\n", perf.num_heap_pushes);

		/* checking */
		for_all_vertices(g, [] (RRNode &v) -> void { v.properties.recalc_occ = 0; });

		for (const auto &net : nets) {
			check_route_tree(route_trees[net.local_id], net, g);
		}

		bool valid = true;
		for_all_vertices(g, [&valid, &buffer] (RRNode &v) -> void {
				sprintf_rr_node(id(v), buffer);
				if (v.properties.recalc_occ != v.properties.occ) {
					zlog_error(delta_log, "Node %s occ mismatch, recalc: %d original: %d\n", buffer, v.properties.recalc_occ, v.properties.occ);
					valid = false;
				}
				});
		assert(valid);

		for (const auto &net : nets) {
			vector<int> overused_rr_node;
			get_overused_nodes(route_trees[net.local_id], get_vertex(route_trees[net.local_id].graph, route_trees[net.local_id].root_rt_node_id), g, overused_rr_node);
			if (!overused_rr_node.empty()) {
				zlog_level(delta_log, ROUTER_V1, "Net %d bb_rank %d overused nodes:\n", net.vpr_id, net.bb_area_rank);
				for (const auto &item : overused_rr_node) {
					zlog_level(delta_log, ROUTER_V1, "%d ", item);
				}
				zlog_level(delta_log, ROUTER_V1, "\n");
			}
		}

		for (const auto &vnet : virtual_nets) {
			assert(vnet.routed);
		}

		char filename[256];
		sprintf(filename, "congestion_%d.txt", iter);
		dump_congestion_map(g, filename);

		if (feasible_routing(g)) {
			zlog_info(delta_log, "Routed in %d iterations\n", iter+1);
			dump_route(*current_traces_ptr, "route.txt");
			routed = true;
		} else {
			int num_overused_nodes = 0;
			for_all_vertices(g, [&num_overused_nodes] (const RRNode &v) -> void {
					if (v.properties.occ > v.properties.capacity) {
					++num_overused_nodes;
					}
					});
			zlog_info(delta_log, "Num overused nodes: %d/%d (%.2f)\n", num_overused_nodes, num_vertices(g), num_overused_nodes*100.0/num_vertices(g));

			if (iter == 0) {
				params.pres_fac = opts->initial_pres_fac;
				update_costs(g, params.pres_fac, 0);
			} else {
				params.pres_fac *= opts->pres_fac_mult;

				/* Avoid overflow for high iteration counts, even if acc_cost is big */
				params.pres_fac = std::min(params.pres_fac, static_cast<float>(HUGE_POSITIVE_FLOAT / 1e5));

				update_costs(g, params.pres_fac, opts->acc_fac);
			}

			for (auto &net : nets) {
				/*zlog_level(delta_log, ROUTER_V1, "Net %d trace stats:\n", net.vpr_id);*/
				adjust_bb_factor((*current_traces_ptr)[net.local_id], net, g, opts->bb_expand_threshold);
			}

			for (auto &net : nets) {
				overused_stats((*current_traces_ptr)[net.local_id], route_trees[net.local_id], net, g);
			}

			/* swap the pointers */
			std::swap(current_traces_ptr, prev_traces_ptr);
			for (const auto &trace : *current_traces_ptr) {
				assert(trace_empty(trace));
			}
		}

		analyze_timing(net_timing);

		iter_timer.stop();
		elapsed = iter_timer.elapsed();

		zlog_info(delta_log, "Iteration time: %g s. Route time: %g s. \n", elapsed.wall / 1e9, route_time / 1e9);

		if (fclose(current_output_log)) {
			char str[256];
			strerror_r(errno, str, sizeof(str));
			printf("failed to close file: %s\n", str);
			assert(false);
		}
	}
	return routed;
}

bool _old_greedy_route_3(t_router_opts *opts)
{
	tbb::task_scheduler_init init(opts->num_threads);

/*#ifdef __linux__*/
	/*printf("domain flags: %d\n", pD->flags);*/
	/*pD->flags = 1;*/
/*#endif*/
	
	init_logging();

	zlog_info(delta_log, "Using boost version %d\n", BOOST_VERSION);

#if (TBB_USE_DEBUG==1)
	zlog_warn(delta_log, "Using debug version of TBB\n");
#endif

	test_perf();

	exit(0);

	/*test_quadrant();*/

	/*auto time = std::chrono::high_resolution_clock::now();*/

	/*printf("test\n");*/

	/*auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::high_resolution_clock::now()-time);*/
	/*printf("printf took %g\n", elapsed.count());*/

	/*test_heap();*/

	/*exit(0);*/

	/*test_scheduler();*/

	/*test_rtree();*/

	/*test_graph();*/

	RRGraph g;
	init_graph(g);

	vector<net_t> nets;
	vector<net_t> global_nets;
	init_nets(nets, global_nets, opts->bb_factor);

	/*test_effectiveness(nets);*/

	/*exit(0);*/

	sort_sinks(nets);

	/*dump_net_bounding_boxes(nets);*/

	/* calculating net bounding box area rank */
	vector<net_t *> nets_ptr(nets.size());
	for (int i = 0; i < nets.size(); ++i) {
		nets_ptr[i] = &nets[i];
	}
	extern s_bb *route_bb;
	sort(nets_ptr.begin(), nets_ptr.end(), [] (const net_t *a, const net_t *b) -> bool {
			return get_bounding_box_area(route_bb[a->vpr_id]) > get_bounding_box_area(route_bb[b->vpr_id]);
			});
	int rank = 0;
	for (auto &net : nets_ptr) {
		net->bb_area_rank = rank++;
	}
	/* end */ 

	extern int nx;
	extern int ny;
	zlog_info(delta_log, "Num nets: %lu Num global nets: %lu nx: %d ny: %d\n", nets.size(), global_nets.size(), nx, ny);
	
	t_net_timing *net_timing = new t_net_timing[nets.size()+global_nets.size()];
	init_net_timing(nets, global_nets, net_timing);

	vector<route_tree_t> route_trees(nets.size());
	for (int i = 0; i < nets.size(); ++i) {
		route_tree_init(route_trees[i]);
	}

	vector<trace_t> prev_traces(nets.size());
	vector<trace_t> traces(nets.size());
	for (int i = 0; i < nets.size(); ++i) {
		trace_init(traces[i]);
		trace_init(prev_traces[i]);
	}

	vector<trace_t> *prev_traces_ptr = &prev_traces;
	vector<trace_t> *current_traces_ptr = &traces;

	route_state_t *state = new route_state_t[num_vertices(g)];

	for (int i = 0; i < num_vertices(g); ++i) {
		state[i].rr_node = -1;
		state[i].known_cost = std::numeric_limits<float>::max();
		state[i].cost = std::numeric_limits<float>::max();
		state[i].prev_edge = nullptr;
		state[i].upstream_R = -1;
		state[i].delay = std::numeric_limits<float>::max();
	}


	/*route_tree_t **route_trees = new route_tree_t *[nets.size()+global_nets.size()];*/
	/*for (int i = 0; i < nets.size(); ++i) {*/
		/*route_trees[nets[i].id] = new route_tree_t[nets[i].sinks.size()+1];*/
		/*for (int j = 0; j <= nets[i].sinks.size(); ++j) {*/
			/*route_tree_init(route_trees[nets[i].id][j]);*/
		/*}*/
	/*}*/
	/*for (int i = 0; i < global_nets.size(); ++i) {*/
		/*route_trees[nets[i].id] = new route_tree_t[nets[i].sinks.size()+1];*/
		/*for (int j = 0; j <= nets[i].sinks.size(); ++j) {*/
			/*route_tree_init(route_trees[nets[i].id][j]);*/
		/*}*/
	/*}*/

	route_parameters_t params;
	params.criticality_exp = opts->criticality_exp;
	params.astar_fac = opts->astar_fac;
	params.max_criticality = opts->max_criticality;
	params.bend_cost = opts->bend_cost;
	params.pres_fac = opts->first_iter_pres_fac; /* Typically 0 -> ignore cong. */

	std::mutex lock;
	extern int nx;
	extern int ny;
	bounding_box_t fpga_bb;
	fpga_bb.xmin = 0;
	fpga_bb.ymin = 0;
	fpga_bb.xmax = nx+2;
	fpga_bb.ymax = ny+2;
	QuadTree<virtual_net_t *> in_flight_qt(fpga_bb);

	int num_virtual_nets = 0;
	for (auto &net : nets) {
		for (auto &sink : net.sinks) {
			++num_virtual_nets;
		}
	}
	vector<virtual_net_t> virtual_nets(num_virtual_nets);
	num_virtual_nets = 0;
	for (auto &net : nets) {
		for (auto &sink : net.sinks) {
			virtual_nets[num_virtual_nets].routed = false;
			virtual_nets[num_virtual_nets].source = &net.source;
			virtual_nets[num_virtual_nets].sinks.push_back(&sink);
			++num_virtual_nets;
		}
	}
	int num_routed_virtual_nets;

	pthread_cond_t start_cond;
	pthread_mutex_t start_cond_lock;
	my_pthread_barrier_t barrier;

	assert(!pthread_mutex_init(&start_cond_lock, 0));
	assert(!pthread_cond_init(&start_cond, 0));
	my_pthread_barrier_init(&barrier, NULL, opts->num_threads);

	perf_t perf;

	vector<std::mutex> net_locks(nets.size());

	WorkerArgs **args = new WorkerArgs*[opts->num_threads];

	sem_t produce_sem, consume_sem;

	/* start the threads first so that they have ample time to finish reach the blocking cond wait part */
	vector<pthread_t> threads(opts->num_threads);
	for (int i = 0; i < opts->num_threads; ++i) {
		/*args[i] = new WorkerArgs(*/
				/*i,*/
				/*opts->num_threads,*/

				/*&start_cond_lock,*/
				/*&start_cond,*/
				/*&barrier,*/

				/*lock,*/
				/*in_flight_qt,*/
				/*virtual_nets,*/
				/*num_routed_virtual_nets,*/

				/*net_locks,*/

				/*g,*/
				/*params,*/
				/*state,*/
				/*route_trees,*/
				/*&prev_traces_ptr,*/
				/*&current_traces_ptr,*/
				/*net_timing,*/
				/*perf,*/

				/*&consume_sem,*/
				/*&produce_sem*/
				/*);*/

		if (i > 0) {
			pthread_create(&threads[i], NULL, worker_thread, args[i]);
		}
	}

	for_all_vertices(g, [] (RRNode &v) -> void {
			v.properties.acc_cost = 1;
			v.properties.pres_cost = 1;
			v.properties.occ = 0;
			});

	char buffer[256];

	bool routed = false;
	for (int iter = 0; iter < opts->max_router_iterations && !routed; ++iter) {
		zlog_info(delta_log, "Routing iteration: %d\n", iter);

		for (int i = 0; i < opts->num_threads; ++i) {
			args[i]->iteration = &iter;
		}

		perf.num_heap_pushes = 0;

		cpu_timer iter_timer;
		iter_timer.start();

		for (int i = 0; i < nets.size(); ++i) {
			route_tree_clear(route_trees[i]);
		}

		for (auto &net : nets) {
			update_sink_criticalities(net, net_timing[net.vpr_id], params);
		}

		for (auto &net : nets) {
			reset_current_source_sink_3(net);
			fill(begin(net.sink_routed), end(net.sink_routed), false);
		}
		for (auto &vnet : virtual_nets) {
			vnet.routed = false;
		}
		num_routed_virtual_nets = 0;
		assert(in_flight_qt.empty());

		nanosecond_type route_time = 0;

		cpu_timer route_timer;

		route_timer.start();

		zlog_info(delta_log, "Going to wake all threads\n");
		/* start all worker threads */
        pthread_cond_broadcast(&start_cond);

		worker_thread_internal(args[0]);

		route_timer.stop();
		cpu_times elapsed = route_timer.elapsed();
		route_time += elapsed.wall;
		/*zlog_info(delta_log, "Scheduling took %s\n", format(elapsed).c_str());*/
		zlog_info(delta_log, "Routing took %g\n", elapsed.wall / 1e9);
		zlog_info(delta_log, "Num heap pushes: %d\n", perf.num_heap_pushes);

		/* checking */
		for_all_vertices(g, [] (RRNode &v) -> void { v.properties.recalc_occ = 0; });

		for (const auto &net : nets) {
			check_route_tree(route_trees[net.local_id], net, g);
		}

		bool valid = true;
		for_all_vertices(g, [&valid, &buffer] (RRNode &v) -> void {
				sprintf_rr_node(id(v), buffer);
				if (v.properties.recalc_occ != v.properties.occ) {
					zlog_error(delta_log, "Node %s occ mismatch, recalc: %d original: %d\n", buffer, v.properties.recalc_occ, v.properties.occ);
					valid = false;
				}
				});
		assert(valid);

		for (const auto &net : nets) {
			vector<int> overused_rr_node;
			get_overused_nodes(route_trees[net.local_id], get_vertex(route_trees[net.local_id].graph, route_trees[net.local_id].root_rt_node_id), g, overused_rr_node);
			if (!overused_rr_node.empty()) {
				zlog_level(delta_log, ROUTER_V1, "Net %d bb_rank %d overused nodes:\n", net.vpr_id, net.bb_area_rank);
				for (const auto &item : overused_rr_node) {
					zlog_level(delta_log, ROUTER_V1, "%d ", item);
				}
				zlog_level(delta_log, ROUTER_V1, "\n");
			}
		}

		for (const auto &vnet : virtual_nets) {
			assert(vnet.routed);
		}

		char filename[256];
		sprintf(filename, "congestion_%d.txt", iter);
		dump_congestion_map(g, filename);

		for (auto &net : nets) {
			/*zlog_level(delta_log, ROUTER_V1, "Net %d trace stats:\n", net.vpr_id);*/
			adjust_bb_factor((*current_traces_ptr)[net.local_id], net, g, opts->bb_expand_threshold);
		}

		if (feasible_routing(g)) {
			zlog_info(delta_log, "Routed in %d iterations\n", iter+1);
			dump_route(*current_traces_ptr, "route.txt");
			routed = true;
		} else {
			int num_overused_nodes = 0;
			for_all_vertices(g, [&num_overused_nodes] (const RRNode &v) -> void {
					if (v.properties.occ > v.properties.capacity) {
					++num_overused_nodes;
					}
					});
			zlog_info(delta_log, "Num overused nodes: %d/%d (%.2f)\n", num_overused_nodes, num_vertices(g), num_overused_nodes*100.0/num_vertices(g));

			if (iter == 0) {
				params.pres_fac = opts->initial_pres_fac;
				update_costs(g, params.pres_fac, 0);
			} else {
				params.pres_fac *= opts->pres_fac_mult;

				/* Avoid overflow for high iteration counts, even if acc_cost is big */
				params.pres_fac = std::min(params.pres_fac, static_cast<float>(HUGE_POSITIVE_FLOAT / 1e5));

				update_costs(g, params.pres_fac, opts->acc_fac);
			}

			/* swap the pointers */
			std::swap(current_traces_ptr, prev_traces_ptr);
			for (const auto &trace : *current_traces_ptr) {
				assert(trace_empty(trace));
			}
		}

		analyze_timing(net_timing);

		iter_timer.stop();
		elapsed = iter_timer.elapsed();

		zlog_info(delta_log, "Iteration time: %g s. Route time: %g s. \n", elapsed.wall / 1e9, route_time / 1e9);
	}
	return routed;
}

bool _old_greedy_route_2(t_router_opts *opts)
{
	tbb::task_scheduler_init init(opts->num_threads);

/*#ifdef __linux__*/
	/*printf("domain flags: %d\n", pD->flags);*/
	/*pD->flags = 1;*/
/*#endif*/
	
	init_logging();

	zlog_info(delta_log, "Using boost version %d\n", BOOST_VERSION);

#if (TBB_USE_DEBUG==1)
	zlog_warn(delta_log, "Using debug version of TBB\n");
#endif

	test_quadrant();

	auto time = std::chrono::high_resolution_clock::now();

	printf("test\n");

	auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::high_resolution_clock::now()-time);
	printf("printf took %lld\n", elapsed.count());

	/*test_heap();*/

	/*exit(0);*/

	/*test_scheduler();*/

	test_rtree();

	/*test_graph();*/

	RRGraph g;
	init_graph(g);

	vector<net_t> nets;
	vector<net_t> global_nets;
	init_nets(nets, global_nets, opts->bb_factor);

	/*test_effectiveness(nets);*/

	/*exit(0);*/

	sort_sinks(nets);

	/*dump_net_bounding_boxes(nets);*/

	/* calculating net bounding box area rank */
	vector<net_t *> nets_ptr(nets.size());
	for (int i = 0; i < nets.size(); ++i) {
		nets_ptr[i] = &nets[i];
	}
	extern s_bb *route_bb;
	sort(nets_ptr.begin(), nets_ptr.end(), [] (const net_t *a, const net_t *b) -> bool {
			return get_bounding_box_area(route_bb[a->vpr_id]) > get_bounding_box_area(route_bb[b->vpr_id]);
			});
	int rank = 0;
	for (auto &net : nets_ptr) {
		net->bb_area_rank = rank++;
	}
	/* end */ 

	extern int nx;
	extern int ny;
	zlog_info(delta_log, "Num nets: %lu Num global nets: %lu nx: %d ny: %d\n", nets.size(), global_nets.size(), nx, ny);
	
	t_net_timing *net_timing = new t_net_timing[nets.size()+global_nets.size()];
	init_net_timing(nets, global_nets, net_timing);

	vector<route_tree_t> route_trees(nets.size());
	for (int i = 0; i < nets.size(); ++i) {
		route_tree_init(route_trees[i]);
	}

	vector<trace_t> prev_traces(nets.size());
	vector<trace_t> traces(nets.size());
	for (int i = 0; i < nets.size(); ++i) {
		trace_init(traces[i]);
		trace_init(prev_traces[i]);
	}

	vector<trace_t> *prev_traces_ptr = &prev_traces;
	vector<trace_t> *current_traces_ptr = &traces;

	route_state_t *state = new route_state_t[num_vertices(g)];

	for (int i = 0; i < num_vertices(g); ++i) {
		state[i].rr_node = -1;
		state[i].known_cost = std::numeric_limits<float>::max();
		state[i].cost = std::numeric_limits<float>::max();
		state[i].prev_edge = nullptr;
		state[i].upstream_R = -1;
		state[i].delay = std::numeric_limits<float>::max();
	}


	/*route_tree_t **route_trees = new route_tree_t *[nets.size()+global_nets.size()];*/
	/*for (int i = 0; i < nets.size(); ++i) {*/
		/*route_trees[nets[i].id] = new route_tree_t[nets[i].sinks.size()+1];*/
		/*for (int j = 0; j <= nets[i].sinks.size(); ++j) {*/
			/*route_tree_init(route_trees[nets[i].id][j]);*/
		/*}*/
	/*}*/
	/*for (int i = 0; i < global_nets.size(); ++i) {*/
		/*route_trees[nets[i].id] = new route_tree_t[nets[i].sinks.size()+1];*/
		/*for (int j = 0; j <= nets[i].sinks.size(); ++j) {*/
			/*route_tree_init(route_trees[nets[i].id][j]);*/
		/*}*/
	/*}*/

	route_parameters_t params;
	params.criticality_exp = opts->criticality_exp;
	params.astar_fac = opts->astar_fac;
	params.max_criticality = opts->max_criticality;
	params.bend_cost = opts->bend_cost;
	params.pres_fac = opts->first_iter_pres_fac; /* Typically 0 -> ignore cong. */

	for_all_vertices(g, [] (RRNode &v) -> void {
			v.properties.acc_cost = 1;
			v.properties.pres_cost = 1;
			v.properties.occ = 0;
			});

	char buffer[256];

	bool routed = false;
	for (int iter = 0; iter < opts->max_router_iterations && !routed; ++iter) {
		zlog_info(delta_log, "Routing iteration: %d\n", iter);

		sprintf(buffer, "concurrency_dump_%d.txt", iter);
		FILE *concurrency_dump = fopen(buffer, "w");
		assert(concurrency_dump);

		cpu_timer iter_timer;
		iter_timer.start();

		for (int i = 0; i < nets.size(); ++i) {
			route_tree_clear(route_trees[i]);
		}

		for (auto &net : nets) {
			update_sink_criticalities(net, net_timing[net.vpr_id], params);
		}

		cpu_timer schedule_timer;
		schedule_timer.start();

		vector<net_t *> nets_to_route;
		int current_local_id = 0;
		for (auto &net : nets) {
			if (net.has_sink) {
				nets_to_route.push_back(&net);
				net.current_local_id = current_local_id++;
			}
		}
		if (iter > 0) {
			for (auto &net : nets_to_route) {
				adjust_bounding_box(*net);
			}
		}
		vector<net_t *> scheduled_nets;
		schedule_nets_faster(nets_to_route, scheduled_nets);

		schedule_timer.stop();
		nanosecond_type schedule_time = schedule_timer.elapsed().wall;
		zlog_info(delta_log, "Scheduling took %g\n", schedule_time / 1e9);

		zlog_info(delta_log, "Initial concurrency: %lu\n", scheduled_nets.size());

		/*zlog_level(delta_log, ROUTER_V3, "-- Start concurrency stats --\n");*/
		/*int temp = 0;*/
		/*for (const auto &nets : scheduled_nets_at_time) {*/
			/*zlog_level(delta_log, ROUTER_V3, "Concurrency at time %d: %d\n", temp, nets.size());*/
			/*++temp;*/
		/*}*/
		/*zlog_level(delta_log, ROUTER_V3, "Average concurrency: %g\n", (float)nets_to_route.size()/scheduled_nets_at_time.size());*/
		/*zlog_level(delta_log, ROUTER_V3, "-- End concurrency stats --\n");*/

		int isink = 0;

		perf_t perf;
		perf.num_heap_pushes = 0;

		nanosecond_type route_time = 0;

		cpu_timer route_timer;

		route_timer.start();

		std::mutex lock;
		extern int nx;
		extern int ny;
		bounding_box_t fpga_bb;
		fpga_bb.xmin = 0;
		fpga_bb.ymin = 0;
		fpga_bb.xmax = nx+2;
		fpga_bb.ymax = ny+2;
		QuadTree<virtual_net_t *> in_flight_qt(fpga_bb);
		for (auto &net : nets) { 
			fill(begin(net.sink_routed), end(net.sink_routed), false);
		}

		vector<virtual_net_t *> virtual_nets;

		for (const auto &net : scheduled_nets) {
			virtual_net_t *virtual_net = new virtual_net_t;

			virtual_net->source = &net->source;
			virtual_net->sinks[0] = &net->sinks[net->current_sink_index];
			virtual_nets.push_back(virtual_net);

			in_flight_qt.insert({ virtual_net->sinks[0]->current_bounding_box, virtual_net });
		}

		vector<virtual_net_t *> unrouted_virtual_nets;
		for (auto &net : nets) {
			for (auto &sink : net.sinks) {
				auto iter = find_if(begin(virtual_nets), end(virtual_nets), [&sink] (const virtual_net_t *v) -> bool { return v->sinks[0] == &sink; });
				if (iter == end(virtual_nets)) {
					virtual_net_t *unrouted_virtual_net = new virtual_net_t;
					unrouted_virtual_net->source = &net.source;
					unrouted_virtual_net->sinks.push_back(&sink);
					unrouted_virtual_nets.push_back(unrouted_virtual_net);
				} 
			}
		}

		int num_spawned_workers = virtual_nets.size();
		int max_num_spawned_workers = opts->num_threads;

		tbb::parallel_do(begin(virtual_nets), end(virtual_nets), RouteWorker(
					lock,
					g,
					params,
					state,
					route_trees,
					&current_traces_ptr,
					&prev_traces_ptr,
					net_timing,
					perf,
					in_flight_qt,
					unrouted_virtual_nets,
					num_spawned_workers,
					max_num_spawned_workers
					));

		route_timer.stop();
		cpu_times elapsed = route_timer.elapsed();
		route_time += elapsed.wall;
		/*zlog_info(delta_log, "Scheduling took %s\n", format(elapsed).c_str());*/
		zlog_info(delta_log, "Routing took %g\n", elapsed.wall / 1e9);

		schedule_timer.start();

		zlog_info(delta_log, "Num heap pushes: %d\n", perf.num_heap_pushes);

		for_all_vertices(g, [] (RRNode &v) -> void { v.properties.recalc_occ = 0; });

		for (const auto &net : nets) {
			check_route_tree(route_trees[net.local_id], net, g);
		}

		bool valid = true;
		for_all_vertices(g, [&valid, &buffer] (RRNode &v) -> void {
				sprintf_rr_node(id(v), buffer);
				if (v.properties.recalc_occ != v.properties.occ) {
					zlog_error(delta_log, "Node %s occ mismatch, recalc: %d original: %d\n", buffer, v.properties.recalc_occ, v.properties.occ);
					valid = false;
				}
				});
		assert(valid);

		for (const auto &net : nets) {
			vector<int> overused_rr_node;
			get_overused_nodes(route_trees[net.local_id], get_vertex(route_trees[net.local_id].graph, route_trees[net.local_id].root_rt_node_id), g, overused_rr_node);
			if (!overused_rr_node.empty()) {
				zlog_level(delta_log, ROUTER_V1, "Net %d bb_rank %d overused nodes:\n", net.vpr_id, net.bb_area_rank);
				for (const auto &item : overused_rr_node) {
					zlog_level(delta_log, ROUTER_V1, "%d ", item);
				}
				zlog_level(delta_log, ROUTER_V1, "\n");
			}
		}

		char filename[256];
		sprintf(filename, "congestion_%d.txt", iter);
		dump_congestion_map(g, filename);

		for (auto &net : nets) {
			/*zlog_level(delta_log, ROUTER_V1, "Net %d trace stats:\n", net.vpr_id);*/
			adjust_bb_factor((*current_traces_ptr)[net.local_id], net, g, opts->bb_expand_threshold);
		}

		if (feasible_routing(g)) {
			zlog_info(delta_log, "Routed in %d iterations\n", iter+1);
			dump_route(*current_traces_ptr, "route.txt");
			routed = true;
		} else {
			int num_overused_nodes = 0;
			for_all_vertices(g, [&num_overused_nodes] (const RRNode &v) -> void {
					if (v.properties.occ > v.properties.capacity) {
					++num_overused_nodes;
					}
					});
			zlog_info(delta_log, "Num overused nodes: %d/%d (%.2f)\n", num_overused_nodes, num_vertices(g), num_overused_nodes*100.0/num_vertices(g));

			if (iter == 0) {
				params.pres_fac = opts->initial_pres_fac;
				update_costs(g, params.pres_fac, 0);
			} else {
				params.pres_fac *= opts->pres_fac_mult;

				/* Avoid overflow for high iteration counts, even if acc_cost is big */
				params.pres_fac = std::min(params.pres_fac, static_cast<float>(HUGE_POSITIVE_FLOAT / 1e5));

				update_costs(g, params.pres_fac, opts->acc_fac);
			}

			/* swap the pointers */
			std::swap(current_traces_ptr, prev_traces_ptr);
			for (const auto &trace : *current_traces_ptr) {
				assert(trace_empty(trace));
			}
		}

		analyze_timing(net_timing);

		iter_timer.stop();
		elapsed = iter_timer.elapsed();

		zlog_info(delta_log, "Iteration time: %g s. Route time: %g s. Schedule time: %g s.\n", elapsed.wall / 1e9, route_time / 1e9, schedule_time / 1e9);
	}
	return routed;
}

bool partitioning_route(t_router_opts *opts)
{
	tbb::task_scheduler_init init(opts->num_threads);

/*#ifdef __linux__*/
	/*printf("domain flags: %d\n", pD->flags);*/
	/*pD->flags = 1;*/
/*#endif*/
	
	init_logging();

	zlog_info(delta_log, "Using boost version %d\n", BOOST_VERSION);

#if (TBB_USE_DEBUG==1)
	zlog_warn(delta_log, "Using debug version of TBB\n");
#endif

	test_quadrant();

	auto time = std::chrono::high_resolution_clock::now();

	printf("test\n");

	auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::high_resolution_clock::now()-time);
	printf("printf took %lld\n", elapsed.count());

	/*test_heap();*/

	/*exit(0);*/

	/*test_scheduler();*/

	test_rtree();

	/*test_graph();*/

	RRGraph g;
	init_graph(g);

	vector<net_t> nets;
	vector<net_t> global_nets;
	init_nets(nets, global_nets, opts->bb_factor);

	/*test_effectiveness(nets);*/

	/*exit(0);*/

	sort_sinks(nets);

	test_quadtree(nets);

	/*dump_net_bounding_boxes(nets);*/

	/* calculating net bounding box area rank */
	vector<net_t *> nets_ptr(nets.size());
	for (int i = 0; i < nets.size(); ++i) {
		nets_ptr[i] = &nets[i];
	}
	extern s_bb *route_bb;
	sort(nets_ptr.begin(), nets_ptr.end(), [] (const net_t *a, const net_t *b) -> bool {
			return get_bounding_box_area(route_bb[a->vpr_id]) > get_bounding_box_area(route_bb[b->vpr_id]);
			});
	int rank = 0;
	for (auto &net : nets_ptr) {
		net->bb_area_rank = rank++;
	}
	/* end */ 

	extern int nx;
	extern int ny;
	zlog_info(delta_log, "Num nets: %lu Num global nets: %lu nx: %d ny: %d\n", nets.size(), global_nets.size(), nx, ny);
	
	t_net_timing *net_timing = new t_net_timing[nets.size()+global_nets.size()];
	init_net_timing(nets, global_nets, net_timing);

	vector<route_tree_t> route_trees(nets.size());
	for (int i = 0; i < nets.size(); ++i) {
		route_tree_init(route_trees[i]);
	}

	vector<trace_t> prev_traces(nets.size());
	vector<trace_t> traces(nets.size());
	for (int i = 0; i < nets.size(); ++i) {
		trace_init(traces[i]);
		trace_init(prev_traces[i]);
	}

	vector<trace_t> *prev_traces_ptr = &prev_traces;
	vector<trace_t> *current_traces_ptr = &traces;

	route_state_t *state = new route_state_t[num_vertices(g)];

	for (int i = 0; i < num_vertices(g); ++i) {
		state[i].rr_node = -1;
		state[i].known_cost = std::numeric_limits<float>::max();
		state[i].cost = std::numeric_limits<float>::max();
		state[i].prev_edge = nullptr;
		state[i].upstream_R = -1;
		state[i].delay = std::numeric_limits<float>::max();
	}


	/*route_tree_t **route_trees = new route_tree_t *[nets.size()+global_nets.size()];*/
	/*for (int i = 0; i < nets.size(); ++i) {*/
		/*route_trees[nets[i].id] = new route_tree_t[nets[i].sinks.size()+1];*/
		/*for (int j = 0; j <= nets[i].sinks.size(); ++j) {*/
			/*route_tree_init(route_trees[nets[i].id][j]);*/
		/*}*/
	/*}*/
	/*for (int i = 0; i < global_nets.size(); ++i) {*/
		/*route_trees[nets[i].id] = new route_tree_t[nets[i].sinks.size()+1];*/
		/*for (int j = 0; j <= nets[i].sinks.size(); ++j) {*/
			/*route_tree_init(route_trees[nets[i].id][j]);*/
		/*}*/
	/*}*/

	route_parameters_t params;
	params.criticality_exp = opts->criticality_exp;
	params.astar_fac = opts->astar_fac;
	params.max_criticality = opts->max_criticality;
	params.bend_cost = opts->bend_cost;
	params.pres_fac = opts->first_iter_pres_fac; /* Typically 0 -> ignore cong. */

	for_all_vertices(g, [] (RRNode &v) -> void {
			v.properties.acc_cost = 1;
			v.properties.pres_cost = 1;
			v.properties.occ = 0;
			});

	char buffer[256];

	bool routed = false;
	for (int iter = 0; iter < opts->max_router_iterations && !routed; ++iter) {
		zlog_info(delta_log, "Routing iteration: %d\n", iter);

		sprintf(buffer, "concurrency_dump_%d.txt", iter);
		FILE *concurrency_dump = fopen(buffer, "w");
		assert(concurrency_dump);

		cpu_timer iter_timer;
		iter_timer.start();

		for (int i = 0; i < nets.size(); ++i) {
			route_tree_clear(route_trees[i]);
		}

		for (auto &net : nets) {
			update_sink_criticalities(net, net_timing[net.vpr_id], params);
		}

		cpu_timer schedule_timer;
		schedule_timer.start();

		vector<net_t *> nets_to_route;
		int current_local_id = 0;
		for (auto &net : nets) {
			if (net.has_sink) {
				nets_to_route.push_back(&net);
				net.current_local_id = current_local_id++;
			}
		}
		vector<vector<const net_t *>> scheduled_nets_at_time;
		if (iter > 0) {
			for (auto &net : nets_to_route) {
				adjust_bounding_box(*net);
			}
		}
		switch (opts->scheduler) {
			case SchedulerType::IND:
				schedule_nets_ind(nets_to_route, scheduled_nets_at_time);
				break;
			case SchedulerType::FAST:
				schedule_nets_fast(nets_to_route, scheduled_nets_at_time);
				break;
			default:
				zlog_error(delta_log, "Error: Unknown scheduler\n");
				assert(false);
				break;
		}

		schedule_timer.stop();
		nanosecond_type schedule_time = schedule_timer.elapsed().wall;
		zlog_info(delta_log, "Scheduling took %g\n", schedule_time / 1e9);

		zlog_info(delta_log, "Average concurrency: %g\n", (float)nets_to_route.size()/scheduled_nets_at_time.size());

		/*zlog_level(delta_log, ROUTER_V3, "-- Start concurrency stats --\n");*/
		/*int temp = 0;*/
		/*for (const auto &nets : scheduled_nets_at_time) {*/
			/*zlog_level(delta_log, ROUTER_V3, "Concurrency at time %d: %d\n", temp, nets.size());*/
			/*++temp;*/
		/*}*/
		/*zlog_level(delta_log, ROUTER_V3, "Average concurrency: %g\n", (float)nets_to_route.size()/scheduled_nets_at_time.size());*/
		/*zlog_level(delta_log, ROUTER_V3, "-- End concurrency stats --\n");*/

		int isink = 0;

		perf_t perf;
		perf.num_heap_pushes = 0;

		nanosecond_type route_time = 0;

		while (scheduled_nets_at_time.size() > 0) {
			zlog_info(delta_log, "isink: %d\n", isink);
			tbb::atomic<int> num_nets_routed = 0;
			cpu_timer route_timer;
			route_timer.start();

			int time = 0;
			for (const auto &nets_at_time : scheduled_nets_at_time) {
				/* this for loop can be parallelized with parallel_for */
				if (true) {
					zlog_info(delta_log, "num nets at time %d: %lu\n", time, nets_at_time.size());
#ifdef __linux__
					__itt_frame_begin_v3(pD, NULL);
							/*__itt_task_begin(pD, __itt_null, __itt_null, shMainTask);*/
#endif
					extern std::chrono::time_point<std::chrono::high_resolution_clock> program_start;
					auto now = std::chrono::high_resolution_clock::now();
					auto frame_start = std::chrono::duration_cast<std::chrono::microseconds>(now-program_start);
					zlog_info(delta_log, "Frame started at %lld\n", frame_start.count());
					tbb::parallel_for(tbb::blocked_range<size_t>(0,nets_at_time.size(),opts->grain_size),
#ifdef __linux__
							[&nets_at_time, &g, &params, &state, &route_trees, &prev_traces_ptr, &current_traces_ptr, &net_timing, &num_nets_routed, &time, &perf] (const tbb::blocked_range<size_t> &r) -> void {
#else
							[&nets_at_time, &g, &params, &state, &route_trees, &prev_traces_ptr, &current_traces_ptr, &net_timing, &num_nets_routed, &time, &perf] (const tbb::blocked_range<size_t> &r) -> void {
#endif
#ifdef __linux__
					/*__itt_frame_begin_v3(pD, NULL);*/
							__itt_task_begin(pD, __itt_null, __itt_null, shMyTask);
#endif
							extern std::chrono::time_point<std::chrono::high_resolution_clock> program_start;
							auto now = std::chrono::high_resolution_clock::now();
							auto elapsed = std::chrono::duration_cast<std::chrono::microseconds>(now-program_start);
							uint64_t tid;
#ifdef __linux__
							tid = syscall(SYS_gettid);
#else
							pthread_threadid_np(pthread_self(), &tid);
#endif
							zlog_info(delta_log, "[%lld] parallel_for %d at %lld us %lu-%lu\n", tid, time, elapsed.count(), r.begin(), r.end()-1);
						for (auto inet = r.begin(); inet != r.end(); ++inet) {
							const net_t *net = nets_at_time[inet];
							net_t virtual_net = *net;
							virtual_net.sinks.clear();
							virtual_net.sinks.push_back(net->sinks[net->current_sink_index]);
							/*if (time == 0) {*/
							/*trace_rip_up_net(traces[net->local_id], g, params.pres_fac);*/
							/*}*/
							trace_rip_up_segment((*prev_traces_ptr)[net->local_id], g, net->sinks[net->current_sink_index].rr_node, params.pres_fac);
							/*printf("Net: %d\n", virtual_net.vpr_id);*/
							route_net(g, virtual_net, params, state, route_trees[net->local_id], (*current_traces_ptr)[net->local_id], net_timing[net->vpr_id], &perf);
							/*route_tree_merge(route_trees[net->id][0], route_trees[net->id][virtual_net.sinks[0].id]);*/
							++num_nets_routed;
						}
#ifdef __linux__ 
					/*__itt_frame_end_v3(pD, NULL);*/
					__itt_task_end(pD);
#endif
							});
#ifdef __linux__ 
					__itt_frame_end_v3(pD, NULL);
					/*__itt_task_end(pD);*/
#endif
				} else {
					for (const auto &net : nets_at_time) {
						net_t virtual_net = *net;
						virtual_net.sinks.clear();
						virtual_net.sinks.push_back(net->sinks[net->current_sink_index]);
						/*if (time == 0) {*/
						/*trace_rip_up_net(traces[net->local_id], g, params.pres_fac);*/
						/*}*/
						trace_rip_up_segment((*prev_traces_ptr)[net->local_id], g, net->sinks[net->current_sink_index].rr_node, params.pres_fac);
						/*printf("Net: %d\n", virtual_net.vpr_id);*/
						route_net(g, virtual_net, params, state, route_trees[net->local_id], (*current_traces_ptr)[net->local_id], net_timing[net->vpr_id], &perf);
						/*route_tree_merge(route_trees[net->id][0], route_trees[net->id][virtual_net.sinks[0].id]);*/
						++num_nets_routed;
					}
				}

				fprintf(concurrency_dump, "%lu\n", nets_at_time.size());
			}

			/* HACK: use this as a marker */
			fprintf(concurrency_dump, "32\n");

			if (isink == 0) {
				assert(num_nets_routed == nets.size());
			}
			assert(num_nets_routed == nets_to_route.size());

			for (auto &net : nets_to_route) {
				set_to_next_sink(*net, route_trees[net->local_id], g, params.astar_fac, perf);
			}

			route_timer.stop();
			cpu_times elapsed = route_timer.elapsed();
			route_time += elapsed.wall;
			/*zlog_info(delta_log, "Scheduling took %s\n", format(elapsed).c_str());*/
			zlog_info(delta_log, "Routing took %g\n", elapsed.wall / 1e9);

			schedule_timer.start();
		
			nets_to_route.clear();
			int current_local_id = 0;
			for (auto &net : nets) {
				if (net.has_sink) {
					nets_to_route.push_back(&net);
					net.current_local_id = current_local_id++;
				}
			}
			scheduled_nets_at_time.clear();
			if (iter > 0) {
				for (auto &net : nets_to_route) {
					adjust_bounding_box(*net);
				}
			}
			switch (opts->scheduler) {
				case SchedulerType::IND:
					schedule_nets_ind(nets_to_route, scheduled_nets_at_time);
					break;
				case SchedulerType::FAST:
					schedule_nets_fast(nets_to_route, scheduled_nets_at_time);
					break;
				default:
					zlog_error(delta_log, "Error: Unknown scheduler\n");
					assert(false);
					break;
			}

			schedule_timer.stop();
			schedule_time += schedule_timer.elapsed().wall;
			zlog_info(delta_log, "Scheduling took %g\n", schedule_timer.elapsed().wall / 1e9);
			zlog_info(delta_log, "Average concurrency: %g\n", (float)nets_to_route.size()/scheduled_nets_at_time.size());

			++isink;
		}

		fclose(concurrency_dump);

		zlog_info(delta_log, "Num heap pushes: %d\n", perf.num_heap_pushes);

		for_all_vertices(g, [] (RRNode &v) -> void { v.properties.recalc_occ = 0; });

		for (const auto &net : nets) {
			check_route_tree(route_trees[net.local_id], net, g);
		}

		bool valid = true;
		for_all_vertices(g, [&valid, &buffer] (RRNode &v) -> void {
				sprintf_rr_node(id(v), buffer);
				if (v.properties.recalc_occ != v.properties.occ) {
					zlog_error(delta_log, "Node %s occ mismatch, recalc: %d original: %d\n", buffer, v.properties.recalc_occ, v.properties.occ);
					valid = false;
				}
				});
		assert(valid);

		for (const auto &net : nets) {
			vector<int> overused_rr_node;
			get_overused_nodes(route_trees[net.local_id], get_vertex(route_trees[net.local_id].graph, route_trees[net.local_id].root_rt_node_id), g, overused_rr_node);
			if (!overused_rr_node.empty()) {
				zlog_level(delta_log, ROUTER_V1, "Net %d bb_rank %d overused nodes:\n", net.vpr_id, net.bb_area_rank);
				for (const auto &item : overused_rr_node) {
					zlog_level(delta_log, ROUTER_V1, "%d ", item);
				}
				zlog_level(delta_log, ROUTER_V1, "\n");
			}
		}

		char filename[256];
		sprintf(filename, "congestion_%d.txt", iter);
		dump_congestion_map(g, filename);

		for (auto &net : nets) {
			/*zlog_level(delta_log, ROUTER_V1, "Net %d trace stats:\n", net.vpr_id);*/
			adjust_bb_factor((*current_traces_ptr)[net.local_id], net, g, opts->bb_expand_threshold);
		}

		if (feasible_routing(g)) {
			zlog_info(delta_log, "Routed in %d iterations\n", iter+1);
			dump_route(*current_traces_ptr, "route.txt");
			routed = true;
		} else {
			int num_overused_nodes = 0;
			for_all_vertices(g, [&num_overused_nodes] (const RRNode &v) -> void {
					if (v.properties.occ > v.properties.capacity) {
					++num_overused_nodes;
					}
					});
			zlog_info(delta_log, "Num overused nodes: %d/%d (%.2f)\n", num_overused_nodes, num_vertices(g), num_overused_nodes*100.0/num_vertices(g));

			if (iter == 0) {
				params.pres_fac = opts->initial_pres_fac;
				update_costs(g, params.pres_fac, 0);
			} else {
				params.pres_fac *= opts->pres_fac_mult;

				/* Avoid overflow for high iteration counts, even if acc_cost is big */
				params.pres_fac = std::min(params.pres_fac, static_cast<float>(HUGE_POSITIVE_FLOAT / 1e5));

				update_costs(g, params.pres_fac, opts->acc_fac);
			}

			for (auto &net : nets) { 
				reset_current_source_sink(net);
			}

			/* swap the pointers */
			std::swap(current_traces_ptr, prev_traces_ptr);
			for (const auto &trace : *current_traces_ptr) {
				assert(trace_empty(trace));
			}
		}

		analyze_timing(net_timing);

		iter_timer.stop();
		cpu_times elapsed = iter_timer.elapsed();

		zlog_info(delta_log, "Iteration time: %g s. Route time: %g s. Schedule time: %g s.\n", elapsed.wall / 1e9, route_time / 1e9, schedule_time / 1e9);
	}
	return routed;
}

bool hybrid_route(t_router_opts *opts)
{
	const int max_iter = 50;

	test_interval();

	init_logging();

	RRGraph g;
	vector<net_t> nets;
	vector<net_t> global_nets;
	init_graph(g);
	init_nets(nets, global_nets, opts->bb_factor);
	
	const int num_parts = 4;

	graph_t<int, int> dep_g;
	/*load_overlapping_nets(nets);*/

	/*const char *metis_filename = "metis.graph";*/
	/*write_metis_file(dep_g, metis_filename);*/

	/*run_metis(metis_filename, num_parts);*/

	/*char part_filename[256];*/
	/*sprintf(part_filename, "%s.part.%d", metis_filename, num_parts);*/
	/*vector<partition_t> partitions(num_parts);*/
	/*init_nets_partition(nets, part_filename, num_parts, partitions);*/

	/*schedule_nets_greedy(nets);*/

	std::sort(nets.begin(), nets.end(), [] (const net_t &a, const net_t &b) -> bool {
		return a.sinks.size() > b.sinks.size(); }
		);

	t_net_timing *net_timing = new t_net_timing[nets.size()];
	init_net_timing(nets, global_nets, net_timing);

	exit(-1);

	route_tree_t *route_trees = new route_tree_t[nets.size()];

	route_parameters_t params;
	params.criticality_exp = opts->criticality_exp;
	params.astar_fac = opts->astar_fac;
	params.max_criticality = opts->max_criticality;
	params.bend_cost = opts->bend_cost;

	for_all_vertices(g, [] (RRNode &v) -> void {
			v.properties.acc_cost = 1;
			v.properties.pres_cost = 1;
			v.properties.occ = 0;
			});

	params.pres_fac = opts->first_iter_pres_fac; /* Typically 0 -> ignore cong. */

	bool routed = false;
	for (int iter = 0; iter < max_iter && !routed; ++iter) {
		zlog_info(delta_log, "Routing iteration: %d\n", iter);
		fine_grained_route_nets(g, nets, params, route_trees, net_timing);
		/*if (use_net_level_parallelism()) {*/
			/*[>delta_stepping_route_nets();<]*/
		/*} else {*/
			/*fine_grained_route_nets();*/
		/*}*/

		if (feasible_routing(g)) {
			zlog_info(delta_log, "Routed in %d iterations\n", iter);
			routed = true;
		} else {
			int num_overused_nodes = 0;
			for_all_vertices(g, [&num_overused_nodes] (const RRNode &v) -> void {
					if (v.properties.occ > v.properties.capacity) {
					++num_overused_nodes;
					}
					});
			zlog_info(delta_log, "Num overused nodes: %d/%d (%.2f)\n", num_overused_nodes, num_vertices(g), num_overused_nodes*100.0/num_vertices(g));

			if (iter == 0) {
				params.pres_fac = opts->initial_pres_fac;
				update_costs(g, params.pres_fac, 0);
			} else {
				params.pres_fac *= opts->pres_fac_mult;

				/* Avoid overflow for high iteration counts, even if acc_cost is big */
				params.pres_fac = std::min(params.pres_fac, static_cast<float>(HUGE_POSITIVE_FLOAT / 1e5));

				update_costs(g, params.pres_fac, opts->acc_fac);
			}
		}

		analyze_timing(net_timing);
	}
	return routed;
}
