#include "pch.h"

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
#include "cluster.h"
#include "args.h"
#include "init.h"
#include "func.h"
#include "partition.h"

#ifdef __linux__
extern __itt_domain* pD;
extern __itt_domain* dispatch_domain;
extern __itt_domain* update_domain;
extern __itt_string_handle *shMyTask;
extern __itt_string_handle *shMainTask;
#endif

bool locking_route(t_router_opts *opts)
{
	tbb::task_scheduler_init init(opts->num_threads);

	using clock = std::chrono::high_resolution_clock;

	extern clock::time_point start_time;

	extern FILE *current_output_log;

	start_time = clock::now();

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

	RRGraph g;
	init_graph(g);

	vector<net_t> nets;
	vector<net_t> global_nets;
	init_nets(nets, global_nets, opts->bb_factor);

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

	vector<route_state_t *> states(opts->num_threads);

	for (int i = 0; i < opts->num_threads; ++i) {
		states[i] = new route_state_t[num_vertices(g)];
		for (int j = 0; j < num_vertices(g); ++j) {
			states[i][j].rr_node = -1;
			states[i][j].known_cost = std::numeric_limits<float>::max();
			states[i][j].cost = std::numeric_limits<float>::max();
			states[i][j].prev_edge = nullptr;
			states[i][j].upstream_R = -1;
			states[i][j].delay = std::numeric_limits<float>::max();
		}
	}

	route_parameters_t params;
	params.criticality_exp = opts->criticality_exp;
	params.astar_fac = opts->astar_fac;
	params.max_criticality = opts->max_criticality;
	params.bend_cost = opts->bend_cost;
	params.pres_fac = opts->first_iter_pres_fac; /* Typically 0 -> ignore cong. */

	//vector<vector<virtual_net_t>> virtual_nets_by_net;
	//create_clustered_virtual_nets(nets, 5, opts->max_sink_bb_area, virtual_nets_by_net);

	for_all_vertices(g, [] (RRNode &v) -> void {
			v.properties.acc_cost = 1;
			v.properties.pres_cost = 1;
			v.properties.occ = 0;
			});

	char buffer[256];

	if (current_output_log) {
		fclose(current_output_log);
	}

	vector<pair<box, net_t *>> nets_to_route;
	vector<pair<box, net_t *>> nets_to_partition;
	vector<vector<int>> overlaps;
	vector<vector<int>> partitions;
	vector<bool> has_interpartition_overlap;
	
	for (auto &net : nets) {
		box b = bg::make_inverse<box>();

		bg::expand(b, point(net.source.x, net.source.y));
		for (const auto &sink : net.sinks) {
			bg::expand(b, point(sink.x, sink.y));
		}
		bg::subtract_value(b.min_corner(), 1+opts->bb_factor);
		bg::add_value(b.max_corner(), opts->bb_factor);

		nets_to_route.push_back(make_pair(b, &net));
	}
	std::sort(begin(nets_to_route), end(nets_to_route), [] (const pair<box, net_t *> &a, const pair<box, net_t *> &b) -> bool {
			return a.second->sinks.size() > b.second->sinks.size();
			});

	bool routed = false;
	unsigned long prev_num_overused_nodes = std::numeric_limits<unsigned long>::max();
	bool use_partitioned = false;

	vector<perf_t> perfs(opts->num_threads);
	vector<lock_perf_t> lock_perfs(opts->num_threads);
	clock::duration total_route_time = clock::duration::zero();
	const int rip_up_all_period = 10;
	int use_partitioned_iter = -1;
	
	for (int iter = 0; iter < opts->max_router_iterations && !routed; ++iter) {
		sprintf(buffer, "/Volumes/DATA/iter_%d.txt", iter);
		current_output_log = fopen(buffer, "w");

		zlog_info(delta_log, "Routing iteration: %d\n", iter);
		bool rip_up_all = (use_partitioned_iter % rip_up_all_period) == 0;
		printf("Routing iteration: %d Use partitioned: %d Rip up all: %d\n", iter, use_partitioned ? 1 : 0, rip_up_all ? 1 : 0);

		auto iter_start = clock::now();

		for (int i = 0; i < opts->num_threads; ++i) {
			lock_perfs[i].num_lock_tries = 0;
			lock_perfs[i].num_lock_waits = 0;
		}

		for (auto &net : nets) {
			update_sink_criticalities(net, net_timing[net.vpr_id], params);
		}

		for (auto &net : nets) {
			net.num_bounding_box_updates = 0;
			net.num_nearest_iters = 0;
			net.total_point_tree_size = 0;
		}

		auto route_start = clock::now();
#ifdef __linux__ 
		__itt_frame_begin_v3(pD, NULL);
#endif
		//tbb::enumerable_thread_specific<state_t *> state_tls;
		//
		tbb::atomic<int> net_index = 0;
		vector<tbb::spin_mutex> debug_lock(opts->num_threads);

		/*std::shuffle(begin(virtual_nets), end(virtual_nets), mt);*/
		if (!use_partitioned) {
			tbb::parallel_for(tbb::blocked_range<int>(0, opts->num_threads, 1),
					[&] (const tbb::blocked_range<int> &range) {

					assert(range.end()-range.begin() == 1);

					int tid = range.begin();

					//assert(debug_lock[tid].try_lock());
					
					perf_t local_perf;
					local_perf.num_heap_pushes = 0;
					local_perf.num_heap_pops = 0;
					local_perf.num_neighbor_visits = 0;

					lock_perf_t local_lock_perf;
					local_lock_perf.num_lock_tries = 0;
					local_lock_perf.num_lock_waits = 0;

					int i;
					while ((i = net_index++) < nets_to_route.size()) {
						net_t *net = nets_to_route[i].second;

						route_tree_mark_congested_nodes_to_be_ripped(route_trees[net->local_id], g);
						route_tree_rip_up_marked(route_trees[net->local_id], g, params.pres_fac, true, &local_lock_perf);

						vector<sink_t *> sinks;	
						sinks.reserve(net->sinks.size());
						for (auto &sink : net->sinks) {
							const RouteTreeNode *sink_rt_node = route_tree_get_rt_node(route_trees[net->local_id], sink.rr_node);
							if (!sink_rt_node)  {
								sinks.push_back(&sink);
							} else {
								assert(!sink_rt_node->properties.pending_rip_up);
							}
						}
						if (!sinks.empty()) {
							route_net_2(g, net->vpr_id, &net->source, sinks, params, states[tid], route_trees[net->local_id], net_timing[net->vpr_id], true, &local_perf, &local_lock_perf);
						}
					}

					perfs[tid] = local_perf;
					lock_perfs[tid] = local_lock_perf; 

					//debug_lock[tid].unlock();
					});
		} else {
			tbb::parallel_for(tbb::blocked_range<int>(0, opts->num_threads, 1),
					[&] (const tbb::blocked_range<int> &range) {

					assert(range.end()-range.begin() == 1);

					int tid = range.begin();

					perf_t local_perf;
					local_perf.num_heap_pushes = 0;
					local_perf.num_heap_pops = 0;
					local_perf.num_neighbor_visits = 0;

					lock_perf_t local_lock_perf;
					local_lock_perf.num_lock_tries = 0;
					local_lock_perf.num_lock_waits = 0;

					//assert(debug_lock[tid].try_lock());

					for (int i = 0; i < partitions[tid].size(); ++i) {
						net_t *net = nets_to_partition[partitions[tid][i]].second;

						if (rip_up_all) {
							route_tree_mark_all_nodes_to_be_ripped(route_trees[net->local_id], g);
						} else {
							route_tree_mark_congested_nodes_to_be_ripped(route_trees[net->local_id], g);
						}
						route_tree_rip_up_marked(route_trees[net->local_id], g, params.pres_fac, true, &local_lock_perf);

						vector<sink_t *> sinks;	
						sinks.reserve(net->sinks.size());
						for (auto &sink : net->sinks) {
							const RouteTreeNode *sink_rt_node = route_tree_get_rt_node(route_trees[net->local_id], sink.rr_node);
							if (!sink_rt_node)  {
								sinks.push_back(&sink);
							} else {
								assert(!sink_rt_node->properties.pending_rip_up);
							}
						}
						if (!sinks.empty()) {
							route_net_2(g, net->vpr_id, &net->source, sinks, params, states[tid], route_trees[net->local_id], net_timing[net->vpr_id], true, &local_perf, &local_lock_perf);
						}
					}

					perfs[tid] = local_perf;
					lock_perfs[tid] = local_lock_perf; 

					//debug_lock[tid].unlock();
				});
			++use_partitioned_iter;
		}

#ifdef __linux__ 
		__itt_frame_end_v3(pD, NULL);
#endif
		
		auto route_time = clock::now()-route_start;
		total_route_time += route_time;

		zlog_info(delta_log, "Routing took %g\n", std::chrono::duration_cast<std::chrono::nanoseconds>(route_time).count() / 1e9);
		printf("Routing took %g\n", std::chrono::duration_cast<std::chrono::nanoseconds>(route_time).count() / 1e9);
		for (int i = 0; i < opts->num_threads; ++i) {
			printf("Thread %d num lock waits/tries: %lu/%lu (%g)\n", i, lock_perfs[i].num_lock_waits, lock_perfs[i].num_lock_tries, (float)lock_perfs[i].num_lock_waits*100/lock_perfs[i].num_lock_tries);
		}

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

		/*char filename[256];*/
		/*sprintf(filename, "congestion_%d.txt", iter);*/
		/*dump_congestion_map(g, filename);*/

		if (feasible_routing(g)) {
			zlog_info(delta_log, "Routed in %d iterations. Total route time: %g\n", iter+1, std::chrono::duration_cast<std::chrono::nanoseconds>(total_route_time).count() / 1e9);
			printf("Routed in %d iterations. Total route time: %g\n", iter+1, std::chrono::duration_cast<std::chrono::nanoseconds>(total_route_time).count() / 1e9);
			//dump_route(*current_traces_ptr, "route.txt");
			routed = true;
		} else {
			unsigned long num_overused_nodes = 0;
			for_all_vertices(g, [&num_overused_nodes] (const RRNode &v) -> void {
					if (v.properties.occ > v.properties.capacity) {
					++num_overused_nodes;
					}
					});
			zlog_info(delta_log, "Num overused nodes: %lu/%d (%.2f)\n", num_overused_nodes, num_vertices(g), num_overused_nodes*100.0/num_vertices(g));
			printf("Num overused nodes: %lu/%d (%.2f)\n", num_overused_nodes, num_vertices(g), num_overused_nodes*100.0/num_vertices(g));

			if (!use_partitioned && num_overused_nodes >= prev_num_overused_nodes && opts->num_threads > 1) {
				nets_to_partition.clear();
				tbb::spin_mutex lock;
				tbb::parallel_for(tbb::blocked_range<int>(0, nets_to_route.size(), 1024),
						[&] (const tbb::blocked_range<int> &range) {

						for (int i = range.begin(); i != range.end(); ++i) {
							net_t *net = nets_to_route[i].second;

							bool marked = route_tree_mark_congested_nodes_to_be_ripped(route_trees[net->local_id], g);
							if (marked) {
								lock.lock();
								nets_to_partition.push_back(nets_to_route[i]);
								lock.unlock();
							}
						}
					});

				overlaps.clear();
				partitions.clear();
				has_interpartition_overlap.clear();
				assert(!nets_to_partition.empty());
				int num_partitions = std::min((unsigned long)opts->num_threads, nets_to_partition.size());
				partition_nets(nets_to_partition, num_partitions, 1000, overlaps, partitions, has_interpartition_overlap);
				use_partitioned = true;
				for (int i = 0; i < opts->num_threads; ++i) {
					printf("Partition %d size %lu\n", i, partitions[i].size());
				}

				use_partitioned_iter = 0;
			} else {
				/* check whether we need to do this */
				//use_partitioned = false;
			}

			prev_num_overused_nodes = num_overused_nodes;

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
				//overused_stats((*current_traces_ptr)[net.local_id], route_trees[net.local_id], net, g);
			}
		}

		analyze_timing(net_timing);

		auto iter_time = clock::now()-iter_start;

		zlog_info(delta_log, "Iteration time: %g s. Route time: %g s. \n", std::chrono::duration_cast<std::chrono::nanoseconds>(iter_time).count() / 1e9, route_time / 1e9);

		if (current_output_log && fclose(current_output_log)) {
			char str[256];
			strerror_r(errno, str, sizeof(str));
			printf("failed to close file: %s\n", str);
			assert(false);
		}
	}

	if (!routed) {
		printf("Failed to route in %d iterations. Total route time: %g\n", opts->max_router_iterations, std::chrono::duration_cast<std::chrono::nanoseconds>(total_route_time).count() / 1e9);
	}

	return routed;
}

bool locking_route__1(t_router_opts *opts)
{
	tbb::task_scheduler_init init(opts->num_threads);

	using clock = std::chrono::high_resolution_clock;

	extern clock::time_point start_time;

	extern FILE *current_output_log;

	start_time = clock::now();

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

	RRGraph g;
	init_graph(g);

	vector<net_t> nets;
	vector<net_t> global_nets;
	init_nets(nets, global_nets, opts->bb_factor);

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

	vector<route_state_t *> states(opts->num_threads);

	for (int i = 0; i < opts->num_threads; ++i) {
		states[i] = new route_state_t[num_vertices(g)];
		for (int j = 0; j < num_vertices(g); ++j) {
			states[i][j].rr_node = -1;
			states[i][j].known_cost = std::numeric_limits<float>::max();
			states[i][j].cost = std::numeric_limits<float>::max();
			states[i][j].prev_edge = nullptr;
			states[i][j].upstream_R = -1;
			states[i][j].delay = std::numeric_limits<float>::max();
		}
	}

	route_parameters_t params;
	params.criticality_exp = opts->criticality_exp;
	params.astar_fac = opts->astar_fac;
	params.max_criticality = opts->max_criticality;
	params.bend_cost = opts->bend_cost;
	params.pres_fac = opts->first_iter_pres_fac; /* Typically 0 -> ignore cong. */

	//vector<vector<virtual_net_t>> virtual_nets_by_net;
	//create_clustered_virtual_nets(nets, 5, opts->max_sink_bb_area, virtual_nets_by_net);

	for_all_vertices(g, [] (RRNode &v) -> void {
			v.properties.acc_cost = 1;
			v.properties.pres_cost = 1;
			v.properties.occ = 0;
			});

	char buffer[256];

	if (current_output_log) {
		fclose(current_output_log);
	}

	vector<pair<box, net_t *>> nets_to_route;
	vector<pair<box, net_t *>> nets_to_partition;
	vector<vector<int>> overlaps;
	vector<vector<int>> partitions;
	vector<bool> has_interpartition_overlap;
	
	for (auto &net : nets) {
		box b = bg::make_inverse<box>();

		bg::expand(b, point(net.source.x, net.source.y));
		for (const auto &sink : net.sinks) {
			bg::expand(b, point(sink.x, sink.y));
		}
		bg::subtract_value(b.min_corner(), 1+opts->bb_factor);
		bg::add_value(b.max_corner(), opts->bb_factor);

		nets_to_route.push_back(make_pair(b, &net));
	}
	std::sort(begin(nets_to_route), end(nets_to_route), [] (const pair<box, net_t *> &a, const pair<box, net_t *> &b) -> bool {
			return a.second->sinks.size() > b.second->sinks.size();
			});

	bool routed = false;
	unsigned long prev_num_overused_nodes = std::numeric_limits<unsigned long>::max();
	bool use_partitioned = false;

	vector<perf_t> perfs(opts->num_threads);
	vector<lock_perf_t> lock_perfs(opts->num_threads);
	clock::duration total_route_time = clock::duration::zero();
	const int rip_up_all_period = 10;
	
	for (int iter = 0; iter < opts->max_router_iterations && !routed; ++iter) {
		sprintf(buffer, "/Volumes/DATA/iter_%d.txt", iter);
		current_output_log = fopen(buffer, "w");

		zlog_info(delta_log, "Routing iteration: %d\n", iter);
		bool rip_up_all = (iter % rip_up_all_period) == 0 && !use_partitioned;
		printf("Routing iteration: %d Use partitioned: %d Rip up all: %d\n", iter, use_partitioned ? 1 : 0, rip_up_all ? 1 : 0);

		auto iter_start = clock::now();

		for (int i = 0; i < opts->num_threads; ++i) {
			lock_perfs[i].num_lock_tries = 0;
			lock_perfs[i].num_lock_waits = 0;
		}

		for (auto &net : nets) {
			update_sink_criticalities(net, net_timing[net.vpr_id], params);
		}

		for (auto &net : nets) {
			net.num_bounding_box_updates = 0;
			net.num_nearest_iters = 0;
			net.total_point_tree_size = 0;
		}

		auto route_start = clock::now();
#ifdef __linux__ 
		__itt_frame_begin_v3(pD, NULL);
#endif
		//tbb::enumerable_thread_specific<state_t *> state_tls;
		//
		tbb::atomic<int> net_index = 0;
		vector<tbb::spin_mutex> debug_lock(opts->num_threads);

		/*std::shuffle(begin(virtual_nets), end(virtual_nets), mt);*/
		if (!use_partitioned) {
			tbb::parallel_for(tbb::blocked_range<int>(0, opts->num_threads, 1),
					[&] (const tbb::blocked_range<int> &range) {

					assert(range.end()-range.begin() == 1);

					int tid = range.begin();

					assert(debug_lock[tid].try_lock());

					int i;
					while ((i = net_index++) < nets_to_route.size()) {
						net_t *net = nets_to_route[i].second;

						if (rip_up_all) {
							route_tree_mark_all_nodes_to_be_ripped(route_trees[net->local_id], g);
						} else {
							route_tree_mark_congested_nodes_to_be_ripped(route_trees[net->local_id], g);
						}
						route_tree_rip_up_marked(route_trees[net->local_id], g, params.pres_fac, true, &lock_perfs[tid]);

						vector<sink_t *> sinks;	
						sinks.reserve(net->sinks.size());
						for (auto &sink : net->sinks) {
							const RouteTreeNode *sink_rt_node = route_tree_get_rt_node(route_trees[net->local_id], sink.rr_node);
							if (!sink_rt_node)  {
								sinks.push_back(&sink);
							} else {
								assert(!sink_rt_node->properties.pending_rip_up);
							}
						}
						if (!sinks.empty()) {
							route_net_2(g, net->vpr_id, &net->source, sinks, params, states[tid], route_trees[net->local_id], net_timing[net->vpr_id], true, &perfs[tid], &lock_perfs[tid]);
						}
					}

					debug_lock[tid].unlock();
					});
		} else {
			tbb::parallel_for(tbb::blocked_range<int>(0, opts->num_threads, 1),
					[&] (const tbb::blocked_range<int> &range) {

					assert(range.end()-range.begin() == 1);

					int tid = range.begin();

					assert(debug_lock[tid].try_lock());

					for (int i = 0; i < partitions[tid].size(); ++i) {
						net_t *net = nets_to_partition[partitions[tid][i]].second;

						route_tree_mark_all_nodes_to_be_ripped(route_trees[net->local_id], g);
						route_tree_rip_up_marked(route_trees[net->local_id], g, params.pres_fac, true, &lock_perfs[tid]);

						vector<sink_t *> sinks;	
						sinks.reserve(net->sinks.size());
						for (auto &sink : net->sinks) {
							const RouteTreeNode *sink_rt_node = route_tree_get_rt_node(route_trees[net->local_id], sink.rr_node);
							if (!sink_rt_node)  {
								sinks.push_back(&sink);
							} else {
								assert(!sink_rt_node->properties.pending_rip_up);
							}
						}
						if (!sinks.empty()) {
							route_net_2(g, net->vpr_id, &net->source, sinks, params, states[tid], route_trees[net->local_id], net_timing[net->vpr_id], true, &perfs[tid], &lock_perfs[tid]);
						}
					}
					debug_lock[tid].unlock();
				});
		}

#ifdef __linux__ 
		__itt_frame_end_v3(pD, NULL);
#endif
		
		auto route_time = clock::now()-route_start;
		total_route_time += route_time;

		zlog_info(delta_log, "Routing took %g\n", std::chrono::duration_cast<std::chrono::nanoseconds>(route_time).count() / 1e9);
		printf("Routing took %g\n", std::chrono::duration_cast<std::chrono::nanoseconds>(route_time).count() / 1e9);
		for (int i = 0; i < opts->num_threads; ++i) {
			printf("Thread %d num lock waits/tries: %lu/%lu (%g)\n", i, lock_perfs[i].num_lock_waits, lock_perfs[i].num_lock_tries, (float)lock_perfs[i].num_lock_waits*100/lock_perfs[i].num_lock_tries);
		}

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

		/*char filename[256];*/
		/*sprintf(filename, "congestion_%d.txt", iter);*/
		/*dump_congestion_map(g, filename);*/

		if (feasible_routing(g)) {
			zlog_info(delta_log, "Routed in %d iterations. Total route time: %g\n", iter+1, std::chrono::duration_cast<std::chrono::nanoseconds>(total_route_time).count() / 1e9);
			printf("Routed in %d iterations. Total route time: %g\n", iter+1, std::chrono::duration_cast<std::chrono::nanoseconds>(total_route_time).count() / 1e9);
			//dump_route(*current_traces_ptr, "route.txt");
			routed = true;
		} else {
			unsigned long num_overused_nodes = 0;
			for_all_vertices(g, [&num_overused_nodes] (const RRNode &v) -> void {
					if (v.properties.occ > v.properties.capacity) {
					++num_overused_nodes;
					}
					});
			zlog_info(delta_log, "Num overused nodes: %lu/%d (%.2f)\n", num_overused_nodes, num_vertices(g), num_overused_nodes*100.0/num_vertices(g));
			printf("Num overused nodes: %lu/%d (%.2f)\n", num_overused_nodes, num_vertices(g), num_overused_nodes*100.0/num_vertices(g));

			if ((!rip_up_all && num_overused_nodes > prev_num_overused_nodes && opts->num_threads > 1)) {
				nets_to_partition.clear();
				tbb::spin_mutex lock;
				tbb::parallel_for(tbb::blocked_range<int>(0, nets_to_route.size(), 1024),
						[&] (const tbb::blocked_range<int> &range) {

						for (int i = range.begin(); i != range.end(); ++i) {
							net_t *net = nets_to_route[i].second;

							bool marked = route_tree_mark_congested_nodes_to_be_ripped(route_trees[net->local_id], g);
							if (marked) {
								lock.lock();
								nets_to_partition.push_back(nets_to_route[i]);
								lock.unlock();
							}
						}
					});

				overlaps.clear();
				partitions.clear();
				has_interpartition_overlap.clear();
				partition_nets(nets_to_partition, opts->num_threads, 1000, overlaps, partitions, has_interpartition_overlap);
				use_partitioned = true;
				for (int i = 0; i < opts->num_threads; ++i) {
					printf("Partition %d size %d\n", i, partitions[i].size());
				}
			} else {
				/* check whether we need to do this */
				//use_partitioned = false;
			}

			prev_num_overused_nodes = num_overused_nodes;

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
				//overused_stats((*current_traces_ptr)[net.local_id], route_trees[net.local_id], net, g);
			}
		}

		analyze_timing(net_timing);

		auto iter_time = clock::now()-iter_start;

		zlog_info(delta_log, "Iteration time: %g s. Route time: %g s. \n", std::chrono::duration_cast<std::chrono::nanoseconds>(iter_time).count() / 1e9, route_time / 1e9);

		if (current_output_log && fclose(current_output_log)) {
			char str[256];
			strerror_r(errno, str, sizeof(str));
			printf("failed to close file: %s\n", str);
			assert(false);
		}
	}

	if (!routed) {
		printf("Failed to route in %d iterations. Total route time: %g\n", opts->max_router_iterations, std::chrono::duration_cast<std::chrono::nanoseconds>(total_route_time).count() / 1e9);
	}

	return routed;
}

bool locking_route_0(t_router_opts *opts)
{
	tbb::task_scheduler_init init(opts->num_threads);

	using clock = std::chrono::high_resolution_clock;

	extern clock::time_point start_time;

	extern FILE *current_output_log;

	start_time = clock::now();

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

	RRGraph g;
	init_graph(g);

	vector<net_t> nets;
	vector<net_t> global_nets;
	init_nets(nets, global_nets, opts->bb_factor);

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

	vector<route_state_t *> states(opts->num_threads);

	for (int i = 0; i < opts->num_threads; ++i) {
		states[i] = new route_state_t[num_vertices(g)];
		for (int j = 0; j < num_vertices(g); ++j) {
			states[i][j].rr_node = -1;
			states[i][j].known_cost = std::numeric_limits<float>::max();
			states[i][j].cost = std::numeric_limits<float>::max();
			states[i][j].prev_edge = nullptr;
			states[i][j].upstream_R = -1;
			states[i][j].delay = std::numeric_limits<float>::max();
		}
	}

	route_parameters_t params;
	params.criticality_exp = opts->criticality_exp;
	params.astar_fac = opts->astar_fac;
	params.max_criticality = opts->max_criticality;
	params.bend_cost = opts->bend_cost;
	params.pres_fac = opts->first_iter_pres_fac; /* Typically 0 -> ignore cong. */

	//vector<vector<virtual_net_t>> virtual_nets_by_net;
	//create_clustered_virtual_nets(nets, 5, opts->max_sink_bb_area, virtual_nets_by_net);

	for_all_vertices(g, [] (RRNode &v) -> void {
			v.properties.acc_cost = 1;
			v.properties.pres_cost = 1;
			v.properties.occ = 0;
			});

	char buffer[256];

	if (current_output_log) {
		fclose(current_output_log);
	}

	vector<pair<box, net_t *>> nets_to_partition;
	vector<vector<int>> overlaps;
	vector<vector<int>> partitions;
	vector<bool> has_interpartition_overlap;
	for (auto &net : nets) {
		box b = bg::make_inverse<box>();

		bg::expand(b, point(net.source.x, net.source.y));
		for (const auto &sink : net.sinks) {
			bg::expand(b, point(sink.x, sink.y));
		}
		bg::subtract_value(b.min_corner(), 1+opts->bb_factor);
		bg::add_value(b.max_corner(), opts->bb_factor);

		nets_to_partition.push_back(make_pair(b, &net));
	}
	std::sort(begin(nets_to_partition), end(nets_to_partition), [] (const pair<box, net_t *> &a, const pair<box, net_t *> &b) -> bool {
			return a.second->sinks.size() > b.second->sinks.size();
			});

	//if (opts->num_threads == 1) {
		//partitions.resize(1);
		//for (int i = 0; i < nets.size(); ++i) {
			//partitions[0].push_back(i);
		//}
	//} else {
		//partition_nets(nets_to_partition, opts->num_threads, overlaps, partitions, has_interpartition_overlap);
	//}

	std::mt19937 mt(time(NULL));

	bool routed = false;
	clock::duration total_route_time = clock::duration::zero();

	tbb::concurrent_queue<net_t *> nets_to_route;
	vector<perf_t> perfs(opts->num_threads);
	vector<lock_perf_t> lock_perfs(opts->num_threads);

	for (int iter = 0; iter < opts->max_router_iterations && !routed; ++iter) {
		sprintf(buffer, "/Volumes/DATA/iter_%d.txt", iter);
		current_output_log = fopen(buffer, "w");

		zlog_info(delta_log, "Routing iteration: %d\n", iter);

		auto iter_start = clock::now();

		for (int i = 0; i < opts->num_threads; ++i) {
			lock_perfs[i].num_lock_tries = 0;
			lock_perfs[i].num_lock_waits = 0;
		}

		for (auto &net : nets) {
			update_sink_criticalities(net, net_timing[net.vpr_id], params);
		}

		for (auto &net : nets) {
			net.num_bounding_box_updates = 0;
			net.num_nearest_iters = 0;
			net.total_point_tree_size = 0;
		}

		auto route_start = clock::now();
#ifdef __linux__ 
		__itt_frame_begin_v3(pD, NULL);
#endif
		//tbb::enumerable_thread_specific<state_t *> state_tls;
		//
		tbb::atomic<int> net_index = 0;
		vector<tbb::spin_mutex> debug_lock(opts->num_threads);

		/*std::shuffle(begin(virtual_nets), end(virtual_nets), mt);*/
		tbb::parallel_for(tbb::blocked_range<int>(0, opts->num_threads, 1),
				[&] (const tbb::blocked_range<int> &range) {

				assert(range.end()-range.begin() == 1);

				int tid = range.begin();

				assert(debug_lock[tid].try_lock());

				int i;
				while ((i = net_index++) < nets_to_partition.size()) {
					net_t *net = nets_to_partition[i].second;

					route_tree_mark_nodes_to_be_ripped(route_trees[net->local_id], g, 10000);
					route_tree_rip_up_marked(route_trees[net->local_id], g, params.pres_fac, true, &lock_perfs[tid]);

					vector<sink_t *> sinks;	
					sinks.reserve(net->sinks.size());
					for (auto &sink : net->sinks) {
						const RouteTreeNode *sink_rt_node = route_tree_get_rt_node(route_trees[net->local_id], sink.rr_node);
						if (!sink_rt_node)  {
							sinks.push_back(&sink);
						} else {
							assert(!sink_rt_node->properties.pending_rip_up);
						}
					}
					if (!sinks.empty()) {
						route_net_2(g, net->vpr_id, &net->source, sinks, params, states[tid], route_trees[net->local_id], net_timing[net->vpr_id], true, &perfs[tid], &lock_perfs[tid]);
					}
				}

				debug_lock[tid].unlock();
				});

#ifdef __linux__ 
		__itt_frame_end_v3(pD, NULL);
#endif
		
		auto route_time = clock::now()-route_start;
		total_route_time += route_time;

		zlog_info(delta_log, "Routing took %g\n", std::chrono::duration_cast<std::chrono::nanoseconds>(route_time).count() / 1e9);
		printf("Routing took %g\n", std::chrono::duration_cast<std::chrono::nanoseconds>(route_time).count() / 1e9);
		for (int i = 0; i < opts->num_threads; ++i) {
			printf("Thread %d num lock waits/tries: %lu/%lu (%g)\n", i, lock_perfs[i].num_lock_waits, lock_perfs[i].num_lock_tries, (float)lock_perfs[i].num_lock_waits*100/lock_perfs[i].num_lock_tries);
		}

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

		/*char filename[256];*/
		/*sprintf(filename, "congestion_%d.txt", iter);*/
		/*dump_congestion_map(g, filename);*/

		if (feasible_routing(g)) {
			zlog_info(delta_log, "Routed in %d iterations. Total route time: %g\n", iter+1, std::chrono::duration_cast<std::chrono::nanoseconds>(total_route_time).count() / 1e9);
			printf("Routed in %d iterations. Total route time: %g\n", iter+1, std::chrono::duration_cast<std::chrono::nanoseconds>(total_route_time).count() / 1e9);
			//dump_route(*current_traces_ptr, "route.txt");
			routed = true;
		} else {
			int num_overused_nodes = 0;
			for_all_vertices(g, [&num_overused_nodes] (const RRNode &v) -> void {
					if (v.properties.occ > v.properties.capacity) {
					++num_overused_nodes;
					}
					});
			zlog_info(delta_log, "Num overused nodes: %d/%d (%.2f)\n", num_overused_nodes, num_vertices(g), num_overused_nodes*100.0/num_vertices(g));
			printf("Num overused nodes: %d/%d (%.2f)\n", num_overused_nodes, num_vertices(g), num_overused_nodes*100.0/num_vertices(g));

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
				//overused_stats((*current_traces_ptr)[net.local_id], route_trees[net.local_id], net, g);
			}
		}

		analyze_timing(net_timing);

		auto iter_time = clock::now()-iter_start;

		zlog_info(delta_log, "Iteration time: %g s. Route time: %g s. \n", std::chrono::duration_cast<std::chrono::nanoseconds>(iter_time).count() / 1e9, route_time / 1e9);

		if (current_output_log && fclose(current_output_log)) {
			char str[256];
			strerror_r(errno, str, sizeof(str));
			printf("failed to close file: %s\n", str);
			assert(false);
		}
	}

	if (!routed) {
		printf("Failed to route in %d iterations. Total route time: %g\n", opts->max_router_iterations, std::chrono::duration_cast<std::chrono::nanoseconds>(total_route_time).count() / 1e9);
	}

	return routed;
}

bool locking_route_2(t_router_opts *opts)
{
	tbb::task_scheduler_init init(opts->num_threads);

	using clock = std::chrono::high_resolution_clock;

	extern clock::time_point start_time;

	extern FILE *current_output_log;

	start_time = clock::now();

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

	RRGraph g;
	init_graph(g);

	vector<net_t> nets;
	vector<net_t> global_nets;
	init_nets(nets, global_nets, opts->bb_factor);

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

	vector<route_state_t *> states(opts->num_threads);

	for (int i = 0; i < opts->num_threads; ++i) {
		states[i] = new route_state_t[num_vertices(g)];
		for (int j = 0; j < num_vertices(g); ++j) {
			states[i][j].rr_node = -1;
			states[i][j].known_cost = std::numeric_limits<float>::max();
			states[i][j].cost = std::numeric_limits<float>::max();
			states[i][j].prev_edge = nullptr;
			states[i][j].upstream_R = -1;
			states[i][j].delay = std::numeric_limits<float>::max();
		}
	}

	route_parameters_t params;
	params.criticality_exp = opts->criticality_exp;
	params.astar_fac = opts->astar_fac;
	params.max_criticality = opts->max_criticality;
	params.bend_cost = opts->bend_cost;
	params.pres_fac = opts->first_iter_pres_fac; /* Typically 0 -> ignore cong. */

	//vector<vector<virtual_net_t>> virtual_nets_by_net;
	//create_clustered_virtual_nets(nets, 5, opts->max_sink_bb_area, virtual_nets_by_net);

	for_all_vertices(g, [] (RRNode &v) -> void {
			v.properties.acc_cost = 1;
			v.properties.pres_cost = 1;
			v.properties.occ = 0;
			});

	char buffer[256];

	if (current_output_log) {
		fclose(current_output_log);
	}

	vector<pair<box, net_t *>> nets_to_partition;
	vector<vector<int>> overlaps;
	vector<vector<int>> partitions;
	vector<bool> has_interpartition_overlap;
	for (auto &net : nets) {
		box b = bg::make_inverse<box>();

		bg::expand(b, point(net.source.x, net.source.y));
		for (const auto &sink : net.sinks) {
			bg::expand(b, point(sink.x, sink.y));
		}
		bg::subtract_value(b.min_corner(), 1+opts->bb_factor);
		bg::add_value(b.max_corner(), opts->bb_factor);

		nets_to_partition.push_back(make_pair(b, &net));
	}
	std::sort(begin(nets_to_partition), end(nets_to_partition), [] (const pair<box, net_t *> &a, const pair<box, net_t *> &b) -> bool {
			return a.second->sinks.size() > b.second->sinks.size();
			});

	//if (opts->num_threads == 1) {
		//partitions.resize(1);
		//for (int i = 0; i < nets.size(); ++i) {
			//partitions[0].push_back(i);
		//}
	//} else {
		//partition_nets(nets_to_partition, opts->num_threads, overlaps, partitions, has_interpartition_overlap);
	//}

	std::mt19937 mt(time(NULL));

	bool routed = false;
	clock::duration total_route_time = clock::duration::zero();

	tbb::concurrent_queue<net_t *> nets_to_route;

	for (int iter = 0; iter < opts->max_router_iterations && !routed; ++iter) {
		sprintf(buffer, "/Volumes/DATA/iter_%d.txt", iter);
		current_output_log = fopen(buffer, "w");

		zlog_info(delta_log, "Routing iteration: %d\n", iter);

		/*perf.num_heap_pushes = 0;*/

		auto iter_start = clock::now();

		for (auto &net : nets) {
			update_sink_criticalities(net, net_timing[net.vpr_id], params);
		}

		for (auto &net : nets) {
			net.num_bounding_box_updates = 0;
			net.num_nearest_iters = 0;
			net.total_point_tree_size = 0;
		}

		auto route_start = clock::now();
#ifdef __linux__ 
		__itt_frame_begin_v3(pD, NULL);
#endif
		//tbb::enumerable_thread_specific<state_t *> state_tls;
		//
		tbb::atomic<int> net_index = 0;

		/*std::shuffle(begin(virtual_nets), end(virtual_nets), mt);*/
		tbb::parallel_for(tbb::blocked_range<int>(0, opts->num_threads, 1),
				[&] (const tbb::blocked_range<int> &range) {

				assert(range.end()-range.begin() == 1);

				int tid = range.begin();

				int i;
				while ((i = net_index++) < nets_to_partition.size()) {
					net_t *net = nets_to_partition[i].second;

					route_tree_mark_nodes_to_be_ripped(route_trees[net->local_id], g, 10000);
					route_tree_rip_up_marked(route_trees[net->local_id], g, params.pres_fac, true, nullptr);

					vector<sink_t *> sinks;	
					sinks.reserve(net->sinks.size());
					for (auto &sink : net->sinks) {
						const RouteTreeNode *sink_rt_node = route_tree_get_rt_node(route_trees[net->local_id], sink.rr_node);
						if (!sink_rt_node)  {
							sinks.push_back(&sink);
						} else {
							assert(!sink_rt_node->properties.pending_rip_up);
						}
					}
					if (!sinks.empty()) {
						route_net_2(g, net->vpr_id, &net->source, sinks, params, states[tid], route_trees[net->local_id], net_timing[net->vpr_id], true, nullptr, nullptr);
					}
				}
				});

#ifdef __linux__ 
		__itt_frame_end_v3(pD, NULL);
#endif
		
		auto route_time = clock::now()-route_start;
		total_route_time += route_time;
		/*zlog_info(delta_log, "Scheduling took %s\n", format(elapsed).c_str());*/
		zlog_info(delta_log, "Routing took %g\n", std::chrono::duration_cast<std::chrono::nanoseconds>(route_time).count() / 1e9);
		printf("Routing took %g\n", std::chrono::duration_cast<std::chrono::nanoseconds>(route_time).count() / 1e9);
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

		/*char filename[256];*/
		/*sprintf(filename, "congestion_%d.txt", iter);*/
		/*dump_congestion_map(g, filename);*/

		if (feasible_routing(g)) {
			zlog_info(delta_log, "Routed in %d iterations. Total route time: %g\n", iter+1, std::chrono::duration_cast<std::chrono::nanoseconds>(total_route_time).count() / 1e9);
			printf("Routed in %d iterations. Total route time: %g\n", iter+1, std::chrono::duration_cast<std::chrono::nanoseconds>(total_route_time).count() / 1e9);
			//dump_route(*current_traces_ptr, "route.txt");
			routed = true;
		} else {
			int num_overused_nodes = 0;
			for_all_vertices(g, [&num_overused_nodes] (const RRNode &v) -> void {
					if (v.properties.occ > v.properties.capacity) {
					++num_overused_nodes;
					}
					});
			zlog_info(delta_log, "Num overused nodes: %d/%d (%.2f)\n", num_overused_nodes, num_vertices(g), num_overused_nodes*100.0/num_vertices(g));
			printf("Num overused nodes: %d/%d (%.2f)\n", num_overused_nodes, num_vertices(g), num_overused_nodes*100.0/num_vertices(g));

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
				//overused_stats((*current_traces_ptr)[net.local_id], route_trees[net.local_id], net, g);
			}
		}

		analyze_timing(net_timing);

		auto iter_time = clock::now()-iter_start;

		zlog_info(delta_log, "Iteration time: %g s. Route time: %g s. \n", std::chrono::duration_cast<std::chrono::nanoseconds>(iter_time).count() / 1e9, route_time / 1e9);

		if (current_output_log && fclose(current_output_log)) {
			char str[256];
			strerror_r(errno, str, sizeof(str));
			printf("failed to close file: %s\n", str);
			assert(false);
		}
	}

	if (!routed) {
		printf("Failed to route in %d iterations. Total route time: %g\n", opts->max_router_iterations, std::chrono::duration_cast<std::chrono::nanoseconds>(total_route_time).count() / 1e9);
	}

	return routed;
}

bool locking_route_3(t_router_opts *opts)
{
	tbb::task_scheduler_init init(opts->num_threads);

	using clock = std::chrono::high_resolution_clock;

	extern clock::time_point start_time;

	extern FILE *current_output_log;

	start_time = clock::now();

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

	RRGraph g;
	init_graph(g);

	vector<net_t> nets;
	vector<net_t> global_nets;
	init_nets(nets, global_nets, opts->bb_factor);

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

	vector<route_state_t *> states(opts->num_threads);

	for (int i = 0; i < opts->num_threads; ++i) {
		states[i] = new route_state_t[num_vertices(g)];
		for (int j = 0; j < num_vertices(g); ++j) {
			states[i][j].rr_node = -1;
			states[i][j].known_cost = std::numeric_limits<float>::max();
			states[i][j].cost = std::numeric_limits<float>::max();
			states[i][j].prev_edge = nullptr;
			states[i][j].upstream_R = -1;
			states[i][j].delay = std::numeric_limits<float>::max();
		}
	}

	route_parameters_t params;
	params.criticality_exp = opts->criticality_exp;
	params.astar_fac = opts->astar_fac;
	params.max_criticality = opts->max_criticality;
	params.bend_cost = opts->bend_cost;
	params.pres_fac = opts->first_iter_pres_fac; /* Typically 0 -> ignore cong. */

	//vector<vector<virtual_net_t>> virtual_nets_by_net;
	//create_clustered_virtual_nets(nets, 5, opts->max_sink_bb_area, virtual_nets_by_net);

	for_all_vertices(g, [] (RRNode &v) -> void {
			v.properties.acc_cost = 1;
			v.properties.pres_cost = 1;
			v.properties.occ = 0;
			});

	char buffer[256];

	if (current_output_log) {
		fclose(current_output_log);
	}

	vector<pair<box, net_t *>> nets_to_partition;
	vector<vector<int>> overlaps;
	vector<vector<int>> partitions;
	vector<bool> has_interpartition_overlap;
	for (auto &net : nets) {
		box b = bg::make_inverse<box>();

		bg::expand(b, point(net.source.x, net.source.y));
		for (const auto &sink : net.sinks) {
			bg::expand(b, point(sink.x, sink.y));
		}
		bg::subtract_value(b.min_corner(), 1+opts->bb_factor);
		bg::add_value(b.max_corner(), opts->bb_factor);

		nets_to_partition.push_back(make_pair(b, &net));
	}

	if (opts->num_threads == 1) {
		partitions.resize(1);
		for (int i = 0; i < nets.size(); ++i) {
			partitions[0].push_back(i);
		}
	} else {
		partition_nets(nets_to_partition, opts->num_threads, 1, overlaps, partitions, has_interpartition_overlap);
	}

	std::mt19937 mt(time(NULL));

	bool routed = false;
	clock::duration total_route_time = clock::duration::zero();

	for (int iter = 0; iter < opts->max_router_iterations && !routed; ++iter) {
		sprintf(buffer, "/Volumes/DATA/iter_%d.txt", iter);
		current_output_log = fopen(buffer, "w");

		zlog_info(delta_log, "Routing iteration: %d\n", iter);

		/*perf.num_heap_pushes = 0;*/

		auto iter_start = clock::now();

		for (auto &net : nets) {
			update_sink_criticalities(net, net_timing[net.vpr_id], params);
		}

		for (auto &net : nets) {
			net.num_bounding_box_updates = 0;
			net.num_nearest_iters = 0;
			net.total_point_tree_size = 0;
		}

		auto route_start = clock::now();
#ifdef __linux__ 
		__itt_frame_begin_v3(pD, NULL);
#endif

		/*std::shuffle(begin(virtual_nets), end(virtual_nets), mt);*/
		tbb::parallel_for(tbb::blocked_range<int>(0, opts->num_threads, 1),
				[&] (const tbb::blocked_range<int> &range) {
				assert(range.end()-range.begin() == 1);

				int tid = range.begin();
				for (int i = 0; i < partitions[tid].size(); ++i) {
					net_t *net = nets_to_partition[partitions[tid][i]].second;

					route_tree_mark_nodes_to_be_ripped(route_trees[net->local_id], g, 5);
					route_tree_rip_up_marked(route_trees[net->local_id], g, params.pres_fac, true, nullptr);

					vector<sink_t *> sinks;	
					sinks.reserve(net->sinks.size());
					for (auto &sink : net->sinks) {
						const RouteTreeNode *sink_rt_node = route_tree_get_rt_node(route_trees[net->local_id], sink.rr_node);
						if (!sink_rt_node)  {
							sinks.push_back(&sink);
						} else {
							assert(!sink_rt_node->properties.pending_rip_up);
						}
					}
					if (!sinks.empty()) {
						route_net_2(g, net->vpr_id, &net->source, sinks, params, states[tid], route_trees[net->local_id], net_timing[net->vpr_id], true, nullptr, nullptr);
					}
				}
				});

#ifdef __linux__ 
		__itt_frame_end_v3(pD, NULL);
#endif
		
		auto route_time = clock::now()-route_start;
		total_route_time += route_time;
		/*zlog_info(delta_log, "Scheduling took %s\n", format(elapsed).c_str());*/
		zlog_info(delta_log, "Routing took %g\n", std::chrono::duration_cast<std::chrono::nanoseconds>(route_time).count() / 1e9);
		printf("Routing took %g\n", std::chrono::duration_cast<std::chrono::nanoseconds>(route_time).count() / 1e9);
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

		/*char filename[256];*/
		/*sprintf(filename, "congestion_%d.txt", iter);*/
		/*dump_congestion_map(g, filename);*/

		if (feasible_routing(g)) {
			zlog_info(delta_log, "Routed in %d iterations. Total route time: %g\n", iter+1, std::chrono::duration_cast<std::chrono::nanoseconds>(total_route_time).count() / 1e9);
			printf("Routed in %d iterations. Total route time: %g\n", iter+1, std::chrono::duration_cast<std::chrono::nanoseconds>(total_route_time).count() / 1e9);
			//dump_route(*current_traces_ptr, "route.txt");
			routed = true;
		} else {
			int num_overused_nodes = 0;
			for_all_vertices(g, [&num_overused_nodes] (const RRNode &v) -> void {
					if (v.properties.occ > v.properties.capacity) {
					++num_overused_nodes;
					}
					});
			zlog_info(delta_log, "Num overused nodes: %d/%d (%.2f)\n", num_overused_nodes, num_vertices(g), num_overused_nodes*100.0/num_vertices(g));
			printf("Num overused nodes: %d/%d (%.2f)\n", num_overused_nodes, num_vertices(g), num_overused_nodes*100.0/num_vertices(g));

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
				//overused_stats((*current_traces_ptr)[net.local_id], route_trees[net.local_id], net, g);
			}
		}

		analyze_timing(net_timing);

		auto iter_time = clock::now()-iter_start;

		zlog_info(delta_log, "Iteration time: %g s. Route time: %g s. \n", std::chrono::duration_cast<std::chrono::nanoseconds>(iter_time).count() / 1e9, route_time / 1e9);

		if (current_output_log && fclose(current_output_log)) {
			char str[256];
			strerror_r(errno, str, sizeof(str));
			printf("failed to close file: %s\n", str);
			assert(false);
		}
	}

	if (!routed) {
		printf("Failed to route in %d iterations. Total route time: %g\n", opts->max_router_iterations, std::chrono::duration_cast<std::chrono::nanoseconds>(total_route_time).count() / 1e9);
	}

	return routed;
}
