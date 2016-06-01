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
#include "congestion.h"
#include "router.h"
#include "partition.h"

#ifdef LMAO
extern __itt_domain* pD;
extern __itt_domain* dispatch_domain;
extern __itt_domain* update_domain;
extern __itt_string_handle *shMyTask;
extern __itt_string_handle *shMainTask;
#endif

using std::chrono::duration_cast;
using std::chrono::nanoseconds;

extern vector<vector<FILE *>> delta_log_files;

static int delta_log_output(zlog_msg_t *msg)
{
	int tid;
	int iter;

	if (sscanf(msg->path, "iter_%d_tid_%d.log", &iter, &tid) != 2) {
		return 0;
	}

	concurrent_log_impl(msg, delta_log_files, iter, tid);

	return 0;
}

void dump_all_net_bounding_boxes(const char *circuit_name, const vector<net_t> &nets)
{
	char filename[256];
	sprintf(filename, "%s_all_net_heatmap.txt", circuit_name);
	FILE *file = fopen(filename, "w");

	fprintf(file, "%lu 0\n", nets.size());
	extern int nx;
	extern int ny;
	fprintf(file, "%d %d\n", nx+2, ny+2);
	for (const auto &net : nets) {
		for (int x = net.bounding_box.xmin; x <= net.bounding_box.xmax; ++x) {
			for (int y = net.bounding_box.ymin; y <= net.bounding_box.ymax; ++y) {
				fprintf(file, "%d %d\n", x, y);
			}
		}
	}
	fclose(file);
}

void dump_all_net_bounding_boxes_area(const char *circuit_name, const vector<net_t> &nets)
{
	char filename[256];
	sprintf(filename, "%s_all_net_bb_area.txt", circuit_name);
	FILE *file = fopen(filename, "w");

	for (const auto &net : nets) {
		fprintf(file, "%d\n", abs(net.bounding_box.xmin-net.bounding_box.xmax)*abs(net.bounding_box.ymin-net.bounding_box.ymax));
	}
	fclose(file);
}

void dump_rr_graph_occ(const congestion_locked_t *congestion, int num_vertices, const char *filename)
{
	FILE *file = fopen(filename, "w");
	for (int i = 0; i < num_vertices; ++i) {
		fprintf(file, "%d\n", congestion[i].cong.occ);
	}
	fclose(file);
}

void bounding_box_overlap_stats(const vector<net_t> &nets)
{
	vector<pair<box, const net_t *>> nets_with_box;
	vector<float> ubvecs = { 1.01, 10, 100, 1000, 10000 };
	
	vector<vector<int>> overlaps;
	vector<vector<int>> partitions;
	vector<bool> has_interpartition_overlap;
	for (auto &net : nets) {
		box b(point(net.bounding_box.xmin, net.bounding_box.ymin),
				point(net.bounding_box.xmax, net.bounding_box.ymax));
		nets_with_box.push_back(make_pair(b, &net));
	}
	for (const auto &ubvec : ubvecs) {
		overlaps.clear();
		partitions.clear();
		has_interpartition_overlap.clear();
		partition_nets(nets_with_box, 8, ubvec, overlaps, partitions, has_interpartition_overlap);
	}
}

bool locking_route(t_router_opts *opts, int run);

bool locking_route_driver(t_router_opts *opts)
{
	int num_routed = 0;
	for (int i = 0; i < opts->num_runs; ++i) {
		bool routed = locking_route(opts, i);
		if (routed) {
			++num_routed;
		}
	}
	printf("%d/%d routed\n", num_routed, opts->num_runs);
	return false;
}

bool locking_route(t_router_opts *opts, int run)
{
	tbb::task_scheduler_init init(opts->num_threads);

	using clock = std::chrono::high_resolution_clock;

	//extern clock::time_point start_time;

	//start_time = clock::now();

/*#ifdef __linux__*/
	/*printf("domain flags: %d\n", pD->flags);*/
	/*pD->flags = 1;*/
/*#endif*/
	
	init_logging();
	zlog_set_record("custom_output", delta_log_output);
    delta_log_files.resize(opts->max_router_iterations);
    for (int i = 0; i < opts->max_router_iterations; ++i) {
        delta_log_files[i].resize(opts->num_threads, nullptr);
    }

	zlog_put_mdc("iter", "0");
	zlog_put_mdc("tid", "0");

	zlog_info(delta_log, "Using boost version %d\n", BOOST_VERSION);

#if (TBB_USE_DEBUG==1)
	zlog_warn(delta_log, "Using debug version of TBB\n");
#endif

	//box a(point(0, 0), point(5,5));
	//box b(point(5, 0), point(10,5));
	//box inter;
	//bg::intersection(a, b, inter);
	//bg::add_value(inter.max_corner(), 1);
	//int lol_area = bg::area(inter);
	//printf("we're here\n");

	RRGraph g;
	init_graph(g);

	vector<net_t> nets;
	vector<net_t> global_nets;
	init_nets(nets, global_nets, opts->bb_factor, opts->large_bb);

	printf("Num nets: %lu\n", nets.size());
	extern char *s_circuit_name;
	printf("Circuit name: %s\n", s_circuit_name);
	//dump_all_net_bounding_boxes(s_circuit_name, nets);
	//dump_all_net_bounding_boxes_area(s_circuit_name, nets);
	//bounding_box_overlap_stats(nets);
	//exit(0);

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
			states[i][j].rr_node = RRGraph::null_vertex();
			states[i][j].known_cost = std::numeric_limits<float>::max();
			states[i][j].cost = std::numeric_limits<float>::max();
			states[i][j].prev_edge = RRGraph::null_edge();
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
	congestion_locked_t *congestion = new congestion_locked_t[num_vertices(g)];
	for (int i = 0; i < num_vertices(g); ++i) {
		congestion[i].cong.acc_cost = 1;
		congestion[i].cong.pres_cost = 1;
		congestion[i].cong.occ = 0;
	}

	char buffer[256];

	vector<pair<box, net_t *>> nets_to_route;
	vector<pair<box, net_t *>> nets_to_partition;
	//vector<vector<int>> overlaps;
	vector<vector<pair<int, int>>> overlaps;
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
	bool initialized_initial_num_overused_nodes = false;
	unsigned long initial_num_overused_nodes = 0;
	bool use_partitioned = false;

	vector<perf_t> perfs(opts->num_threads);
	vector<lock_perf_t> lock_perfs(opts->num_threads);
	vector<clock::time_point> greedy_end_time(opts->num_threads);
	vector<clock::time_point> partitioned_end_time(opts->num_threads);
	clock::duration total_greedy_route_time = clock::duration::zero();
	clock::duration total_partitioned_route_time = clock::duration::zero();
	clock::duration total_update_cost_time = clock::duration::zero();
	clock::duration total_partitioning_time = clock::duration::zero();
	clock::duration total_analyze_timing_time = clock::duration::zero();
	clock::duration total_iter_time = clock::duration::zero();

	//int greedy_rip_up_all_period = 3;
	int partitioned_iter = -1;
	int iter;
	int next_greedy_rip_up_iter = 0;
	float crit_path_delay;

	for (iter = 0; iter < opts->max_router_iterations && !routed; ++iter) {
		clock::duration greedy_route_time = clock::duration::zero();
		clock::duration partitioned_route_time = clock::duration::zero();
		clock::duration update_cost_time = clock::duration::zero();
		clock::duration partitioning_time = clock::duration::zero();
		clock::duration analyze_timing_time = clock::duration::zero();
		clock::duration iter_time = clock::duration::zero();

		zlog_info(delta_log, "Routing iteration: %d\n", iter);
		bool partitioned_rip_up_all = opts->rip_up_always || ((partitioned_iter % opts->rip_up_period) == 0);
		bool greedy_rip_up_all = opts->rip_up_always;//(next_greedy_rip_up_iter == iter);
		printf("Routing iteration: %d Use partitioned: %d Greedy rip up all: %d Partitioned rip up all: %d\n", iter, use_partitioned ? 1 : 0, greedy_rip_up_all ? 1 : 0, partitioned_rip_up_all ? 1 : 0);

		for (auto &net : nets) {
			net.num_bounding_box_updates = 0;
			net.num_nearest_iters = 0;
			net.total_point_tree_size = 0;
		}

		auto iter_start = clock::now();

		//auto route_start = clock::now();

#ifdef LMAO 
		__itt_frame_begin_v3(pD, NULL);
#endif
		//tbb::enumerable_thread_specific<state_t *> state_tls;
		//
		tbb::atomic<int> net_index = 0;
		vector<tbb::spin_mutex> debug_lock(opts->num_threads);
		tbb::atomic<int> num_nets_routed = 0;

		/*std::shuffle(begin(virtual_nets), end(virtual_nets), mt);*/
		if (!use_partitioned) {
			auto greedy_route_start = clock::now();

			tbb::parallel_for(tbb::blocked_range<int>(0, opts->num_threads, 1),
					[&] (const tbb::blocked_range<int> &range) {

					assert(range.end()-range.begin() == 1);

					int tid = range.begin();

					char local_buffer[256];
					sprintf(local_buffer, "%d", iter);
					zlog_put_mdc("iter", local_buffer);

					sprintf(local_buffer, "%d", tid);
					zlog_put_mdc("tid", local_buffer);

					//assert(debug_lock[tid].try_lock());
					
					perf_t local_perf;
					local_perf.num_heap_pushes = 0;
					local_perf.num_heap_pops = 0;
					local_perf.num_neighbor_visits = 0;

					lock_perf_t local_lock_perf;
					local_lock_perf.num_lock_tries = 0;
					local_lock_perf.num_lock_waits = 0;
					local_lock_perf.total_wait_time = clock::duration::zero();

					if (opts->work_conserving) {
						int i;
						while ((i = net_index++) < nets_to_route.size()) {
							net_t *net = nets_to_route[i].second;

							update_sink_criticalities(*net, net_timing[net->vpr_id], params);

							//auto rip_up_start = clock::now();
							zlog_level(delta_log, ROUTER_V1, "Ripping up net %d\n", net->vpr_id);
							if (greedy_rip_up_all) {
								route_tree_mark_all_nodes_to_be_ripped(route_trees[net->local_id], g);
							} else {
								route_tree_mark_congested_nodes_to_be_ripped(route_trees[net->local_id], g, congestion);
							}
							route_tree_rip_up_marked(route_trees[net->local_id], g, congestion, params.pres_fac, true, &local_lock_perf);

							//local_perf.total_rip_up_time += clock::now()-rip_up_start;

							//auto route_start = clock::now();

							vector<sink_t *> sinks;	
							sinks.reserve(net->sinks.size());
							for (auto &sink : net->sinks) {
								RouteTreeNode sink_rt_node = route_tree_get_rt_node(route_trees[net->local_id], sink.rr_node);
								if (sink_rt_node == RouteTree::null_vertex())  {
									sinks.push_back(&sink);
								} else {
									const auto &sink_rt_node_p = get_vertex_props(route_trees[net->local_id].graph, sink_rt_node);
									assert(!sink_rt_node_p.pending_rip_up);
								}
							}

							if (!sinks.empty()) {
								route_net_with_fine_grain_lock(g, net->vpr_id, &net->source, sinks, params, states[tid], congestion, route_trees[net->local_id], net_timing[net->vpr_id], true, &local_perf, &local_lock_perf);

								++num_nets_routed;
							}

							//local_perf.total_route_time += clock::now()-rip_up_start;
						}
					} else {
						int i = tid;
						while (i < nets_to_route.size()) {
							//while ((i = net_index++) < nets_to_route.size()) {
							net_t *net = nets_to_route[i].second;

							update_sink_criticalities(*net, net_timing[net->vpr_id], params);

							//auto rip_up_start = clock::now();
							zlog_level(delta_log, ROUTER_V1, "Ripping up net %d\n", net->vpr_id);
							if (greedy_rip_up_all) {
								route_tree_mark_all_nodes_to_be_ripped(route_trees[net->local_id], g);
							} else {
								route_tree_mark_congested_nodes_to_be_ripped(route_trees[net->local_id], g, congestion);
							}
							route_tree_rip_up_marked(route_trees[net->local_id], g, congestion, params.pres_fac, true, &local_lock_perf);

							//local_perf.total_rip_up_time += clock::now()-rip_up_start;

							//auto route_start = clock::now();

							vector<sink_t *> sinks;	
							sinks.reserve(net->sinks.size());
							for (auto &sink : net->sinks) {
								RouteTreeNode sink_rt_node = route_tree_get_rt_node(route_trees[net->local_id], sink.rr_node);
								if (sink_rt_node == RouteTree::null_vertex())  {
									sinks.push_back(&sink);
								} else {
									const auto &sink_rt_node_p = get_vertex_props(route_trees[net->local_id].graph, sink_rt_node);
									assert(!sink_rt_node_p.pending_rip_up);
								}
							}

							if (!sinks.empty()) {
								route_net_with_fine_grain_lock(g, net->vpr_id, &net->source, sinks, params, states[tid], congestion, route_trees[net->local_id], net_timing[net->vpr_id], true, &local_perf, &local_lock_perf);

								++num_nets_routed;
							}

							//local_perf.total_route_time += clock::now()-rip_up_start;
							i += opts->num_threads;
						}
					}

					greedy_end_time[tid] = clock::now();

					perfs[tid] = local_perf;
					lock_perfs[tid] = local_lock_perf; 

					//debug_lock[tid].unlock();
					});

			greedy_route_time = clock::now()-greedy_route_start;
		} else {
			auto partitioned_route_start = clock::now();

			if (false) {
				for (const auto &n : nets_to_partition) {
					net_t *net = n.second;
					if (partitioned_rip_up_all) {
						route_tree_mark_all_nodes_to_be_ripped(route_trees[net->local_id], g);
					} else {
						route_tree_mark_congested_nodes_to_be_ripped(route_trees[net->local_id], g, congestion);
					}
					route_tree_rip_up_marked(route_trees[net->local_id], g, congestion, params.pres_fac, true, nullptr);

					vector<sink_t *> sinks;	
					sinks.reserve(net->sinks.size());
					for (auto &sink : net->sinks) {
						RouteTreeNode sink_rt_node = route_tree_get_rt_node(route_trees[net->local_id], sink.rr_node);
						if (sink_rt_node == RouteTree::null_vertex())  {
							sinks.push_back(&sink);
						} else {
							const auto &sink_rt_node_p = get_vertex_props(route_trees[net->local_id].graph, sink_rt_node);
							assert(!sink_rt_node_p.pending_rip_up);
						}
					}
					if (!sinks.empty()) {
						route_net_with_fine_grain_lock(g, net->vpr_id, &net->source, sinks, params, states[0], congestion, route_trees[net->local_id], net_timing[net->vpr_id], true, nullptr, nullptr);
					}
				}
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
						local_lock_perf.total_wait_time = clock::duration::zero();

						//assert(debug_lock[tid].try_lock());

						for (int i = 0; i < partitions[tid].size(); ++i) {
							net_t *net = nets_to_partition[partitions[tid][i]].second;

							update_sink_criticalities(*net, net_timing[net->vpr_id], params);

							if (partitioned_rip_up_all) {
								route_tree_mark_all_nodes_to_be_ripped(route_trees[net->local_id], g);
							} else {
								route_tree_mark_congested_nodes_to_be_ripped(route_trees[net->local_id], g, congestion);
							}
							route_tree_rip_up_marked(route_trees[net->local_id], g, congestion, params.pres_fac, true, &local_lock_perf);

							vector<sink_t *> sinks;	
							sinks.reserve(net->sinks.size());
							for (auto &sink : net->sinks) {
								RouteTreeNode sink_rt_node = route_tree_get_rt_node(route_trees[net->local_id], sink.rr_node);
								if (sink_rt_node == RouteTree::null_vertex()) {
									sinks.push_back(&sink);
								} else {
									const auto &sink_rt_node_p = get_vertex_props(route_trees[net->local_id].graph, sink_rt_node);
									assert(!sink_rt_node_p.pending_rip_up);
								}
							}
							if (!sinks.empty()) {
								route_net_with_fine_grain_lock(g, net->vpr_id, &net->source, sinks, params, states[tid], congestion, route_trees[net->local_id], net_timing[net->vpr_id], true, &local_perf, &local_lock_perf);

								++num_nets_routed;
							}
						}

						partitioned_end_time[tid] = clock::now();

						perfs[tid] = local_perf;
						lock_perfs[tid] = local_lock_perf; 

						//debug_lock[tid].unlock();
					});
			}

			++partitioned_iter;

			partitioned_route_time = clock::now()-partitioned_route_start;
		}

		//if (greedy_rip_up_all) {
			//next_greedy_rip_up_iter += greedy_rip_up_all_period;
			//++greedy_rip_up_all_period;
			//prev_num_overused_nodes = std::numeric_limits<unsigned long>::max();
		//}

#ifdef LMAO 
		__itt_frame_end_v3(pD, NULL);
#endif
		
		//route_time = clock::now()-route_start;

		iter_time = clock::now()-iter_start;

		/* checking */
		for (int i = 0; i < num_vertices(g); ++i) {
			congestion[i].cong.recalc_occ = 0; 
		}

		for (const auto &net : nets) {
			check_route_tree(route_trees[net.local_id], net, g);
			recalculate_occ(route_trees[net.local_id], g, congestion);
		}

		bool valid = true;
		for (int i = 0; i < num_vertices(g); ++i) {
			sprintf_rr_node(i, buffer);
			if (congestion[i].cong.recalc_occ != congestion[i].cong.occ) {
				zlog_error(delta_log, "Node %s occ mismatch, recalc: %d original: %d\n", buffer, congestion[i].cong.recalc_occ, congestion[i].cong.occ);
				valid = false;
			}
		}
		assert(valid);

		for (const auto &net : nets) {
			vector<int> overused_rr_node;
			assert(route_trees[net.local_id].root_rt_nodes.size() == 1);
			get_overused_nodes(route_trees[net.local_id], route_trees[net.local_id].root_rt_nodes[0], g, congestion, overused_rr_node);
			if (!overused_rr_node.empty()) {
				zlog_level(delta_log, ROUTER_V1, "Net %d bb_rank %d overused nodes:\n", net.vpr_id, net.bb_area_rank);
				for (const auto &item : overused_rr_node) {
					zlog_level(delta_log, ROUTER_V1, "%d ", item);
				}
				zlog_level(delta_log, ROUTER_V1, "\n");
			}
		}

		iter_start = clock::now();

		if (feasible_routing(g, congestion)) {
			//dump_route(*current_traces_ptr, "route.txt");
			routed = true;
		} else {
			unsigned long num_overused_nodes = 0;
			for (int i = 0; i < num_vertices(g); ++i) {
				if (congestion[i].cong.occ > get_vertex_props(g, i).capacity) {
					++num_overused_nodes;
				}
			}
			zlog_info(delta_log, "Num overused nodes: %lu/%d (%.2f)\n", num_overused_nodes, num_vertices(g), num_overused_nodes*100.0/num_vertices(g));
			printf("Num overused nodes: %lu/%d (%.2f)\n", num_overused_nodes, num_vertices(g), num_overused_nodes*100.0/num_vertices(g));

			if (iter == 0) {
				initial_num_overused_nodes = num_overused_nodes;
			}

			if (!use_partitioned && num_overused_nodes >= prev_num_overused_nodes && iter > 5 && opts->num_threads > 1) {
			//if (!use_partitioned && iter > 0 && (float)num_overused_nodes/initial_num_overused_nodes < opts->transition_threshold && opts->num_threads > 1) {
				auto partitioning_start = clock::now();

				nets_to_partition.clear();
				tbb::spin_mutex lock;
				tbb::parallel_for(tbb::blocked_range<int>(0, nets_to_route.size(), 1024),
						[&] (const tbb::blocked_range<int> &range) {

						for (int i = range.begin(); i != range.end(); ++i) {
							net_t *net = nets_to_route[i].second;

							bool marked = route_tree_mark_congested_nodes_to_be_ripped(route_trees[net->local_id], g, congestion);
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
				int num_partitions = opts->num_threads;
				//int num_partitions = std::min((unsigned long)opts->num_threads, nets_to_partition.size());
				partition_nets_overlap_area_metric(nets_to_partition, num_partitions, 1000000, overlaps, partitions, has_interpartition_overlap);
				use_partitioned = true;
				for (int i = 0; i < opts->num_threads; ++i) {
					printf("Partition %d size %lu\n", i, partitions[i].size());
				}

				partitioned_iter = 0;

				partitioning_time = clock::now()-partitioning_start;
			} else {
				/* check whether we need to do this */
				//use_partitioned = false;
			}

			prev_num_overused_nodes = num_overused_nodes;

			auto update_cost_start = clock::now();

			if (iter == 0) {
				params.pres_fac = opts->initial_pres_fac;
				update_costs(g, congestion, params.pres_fac, 0);
			} else {
				params.pres_fac *= opts->pres_fac_mult;

				/* Avoid overflow for high iteration counts, even if acc_cost is big */
				params.pres_fac = std::min(params.pres_fac, static_cast<float>(HUGE_POSITIVE_FLOAT / 1e5));

				update_costs(g, congestion, params.pres_fac, opts->acc_fac);
			}

			update_cost_time = clock::now()-update_cost_start;
		}

		auto analyze_timing_start = clock::now();

		crit_path_delay = analyze_timing(net_timing);

		analyze_timing_time = clock::now()-analyze_timing_start;

		iter_time += clock::now()-iter_start;

		//total_route_time += route_time;
		total_greedy_route_time += greedy_route_time;
		total_partitioned_route_time += partitioned_route_time;
		total_update_cost_time += update_cost_time;
		total_partitioning_time += partitioning_time;
		total_analyze_timing_time += analyze_timing_time;
		total_iter_time += iter_time;

		printf("Iteration time: %g s.\n", duration_cast<nanoseconds>(iter_time).count() / 1e9);
			printf("\tRoute time: %g s.\n", duration_cast<nanoseconds>(greedy_route_time+partitioned_route_time).count() / 1e9);
				printf("\t\tGreedy route time: %g s.\n", duration_cast<nanoseconds>(greedy_route_time).count() / 1e9);
				printf("\t\tPartitioned route time: %g s.\n", duration_cast<nanoseconds>(partitioned_route_time).count() / 1e9);
			printf("\tUpdate cost time: %g s.\n", duration_cast<nanoseconds>(update_cost_time).count() / 1e9);
			printf("\tPartitioning time: %g s.\n", duration_cast<nanoseconds>(partitioning_time).count() / 1e9);
			printf("\tAnalyze timing time: %g s.\n", duration_cast<nanoseconds>(analyze_timing_time).count() / 1e9);
		printf("Critical path: %g ns\n", crit_path_delay);

		printf("Num nets routed: %d\n", num_nets_routed);

		unsigned long total_num_heap_pushes = 0;
		unsigned long total_num_heap_pops = 0;
		unsigned long total_num_neighbor_visits = 0;

		for (int i = 0; i < opts->num_threads; ++i) {
			total_num_heap_pushes += perfs[i].num_heap_pushes;
			total_num_heap_pops += perfs[i].num_heap_pops;
			total_num_neighbor_visits += perfs[i].num_neighbor_visits;
		}

		for (int i = 0; i < opts->num_threads; ++i) {
			printf("Thread %d num lock waits/tries: %lu/%lu (%g)\n", i, lock_perfs[i].num_lock_waits, lock_perfs[i].num_lock_tries, (float)lock_perfs[i].num_lock_waits*100/lock_perfs[i].num_lock_tries);
			printf("Thread %d total wait time: %g (%g)\n", i, duration_cast<nanoseconds>(lock_perfs[i].total_wait_time).count() / 1e9, 100.0*lock_perfs[i].total_wait_time/(greedy_route_time+partitioned_route_time));
			printf("Thread %d num_heap_pushes: %lu\n", i, perfs[i].num_heap_pushes);
			printf("Thread %d num_heap_pops: %lu\n", i, perfs[i].num_heap_pops);
			printf("Thread %d num_neighbor_visits: %lu\n", i, perfs[i].num_neighbor_visits);
		}

		printf("total_num_heap_pushes: %lu\n", total_num_heap_pushes);
		printf("total_num_heap_pops: %lu\n", total_num_heap_pops); 
		printf("total_num_neighbor_visits: %lu\n", total_num_neighbor_visits);

		clock::time_point greedy_earliest_end_time = *std::min_element(begin(greedy_end_time), end(greedy_end_time));
		clock::time_point partitioned_earliest_end_time = *std::min_element(begin(partitioned_end_time), end(partitioned_end_time));

		for (int i = 0; i < opts->num_threads; ++i) {
			printf("Thread %d greedy wait time %g (%g)\n", i, duration_cast<nanoseconds>(greedy_end_time[i]-greedy_earliest_end_time).count() / 1e9, 100.0*(greedy_end_time[i]-greedy_earliest_end_time)/greedy_route_time);
		}
		for (int i = 0; i < opts->num_threads; ++i) {
			if (partitioned_route_time > clock::duration::zero()) {
				printf("Thread %d partitioned wait time %g (%g)\n", i, duration_cast<nanoseconds>(partitioned_end_time[i]-partitioned_earliest_end_time).count() / 1e9, 100.0*(partitioned_end_time[i]-partitioned_earliest_end_time)/partitioned_route_time);
			}
		}
	}

	sprintf(buffer, "%s_run_%d_rr_graph_occ.txt", s_circuit_name, run);
	dump_rr_graph_occ(congestion, num_vertices(g), buffer);

	if (routed) {
		printf("Routed in %d iterations. Total iteration time: %g\n", iter, duration_cast<nanoseconds>(total_iter_time).count() / 1e9);
			printf("\tTotal route time: %g s.\n", duration_cast<nanoseconds>(total_greedy_route_time+total_partitioned_route_time).count() / 1e9);
				printf("\t\tTotal greedy route time: %g s.\n", duration_cast<nanoseconds>(total_greedy_route_time).count() / 1e9);
				printf("\t\tTotal partitioned route time: %g s.\n", duration_cast<nanoseconds>(total_partitioned_route_time).count() / 1e9);
			printf("\tTotal update cost time: %g s.\n", duration_cast<nanoseconds>(total_update_cost_time).count() / 1e9);
			printf("\tTotal partitioning time: %g s.\n", duration_cast<nanoseconds>(total_partitioning_time).count() / 1e9);
			printf("\tTotal analyze timing time: %g s.\n", duration_cast<nanoseconds>(total_analyze_timing_time).count() / 1e9);

		printf("Final critical path: %g ns\n", crit_path_delay);
	} else {
		printf("Failed to route in %d iterations. Total iteration time: %g\n", opts->max_router_iterations, duration_cast<nanoseconds>(total_iter_time).count() / 1e9);
			printf("\tTotal route time: %g s.\n", duration_cast<nanoseconds>(total_greedy_route_time+total_partitioned_route_time).count() / 1e9);
				printf("\t\tTotal greedy route time: %g s.\n", duration_cast<nanoseconds>(total_greedy_route_time).count() / 1e9);
				printf("\t\tTotal partitioned route time: %g s.\n", duration_cast<nanoseconds>(total_partitioned_route_time).count() / 1e9);
			printf("\tTotal update cost time: %g s.\n", duration_cast<nanoseconds>(total_update_cost_time).count() / 1e9);
			printf("\tTotal partitioning time: %g s.\n", duration_cast<nanoseconds>(total_partitioning_time).count() / 1e9);
			printf("\tTotal analyze timing time: %g s.\n", duration_cast<nanoseconds>(total_analyze_timing_time).count() / 1e9);
	}

	delete_graph(g);
	delete_net_timing(nets, global_nets, net_timing);	
	delete [] net_timing;
	for (int i = 0; i < opts->num_threads; ++i) {
		delete [] states[i];
	}

	return routed;
}

bool locking_route__1(t_router_opts *opts)
{
	tbb::task_scheduler_init init(opts->num_threads);

	using clock = std::chrono::high_resolution_clock;

	//extern clock::time_point start_time;

	//extern FILE *current_output_log;

	//start_time = clock::now();

/*#ifdef __linux__*/
	/*printf("domain flags: %d\n", pD->flags);*/
	/*pD->flags = 1;*/
/*#endif*/
	
	init_logging();
	//zlog_set_record("custom_output", zlog_custom_output);
	//zlog_set_record("sched_custom_output", zlog_sched_custom_output);
	//current_output_log = fopen("/Volumes/DATA/main.log", "w");

	zlog_info(delta_log, "Using boost version %d\n", BOOST_VERSION);

#if (TBB_USE_DEBUG==1)
	zlog_warn(delta_log, "Using debug version of TBB\n");
#endif

	RRGraph g;
	init_graph(g);

	vector<net_t> nets;
	vector<net_t> global_nets;
	init_nets(nets, global_nets, opts->bb_factor, opts->large_bb);

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
			states[i][j].rr_node = RRGraph::null_vertex();
			states[i][j].known_cost = std::numeric_limits<float>::max();
			states[i][j].cost = std::numeric_limits<float>::max();
			states[i][j].prev_edge = RRGraph::null_edge();
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

	//for_all_vertices(g, [] (RRNode &v) -> void {
			//v.properties.acc_cost = 1;
			//v.properties.pres_cost = 1;
			//v.properties.occ = 0;
			//});

	char buffer[256];

	//if (current_output_log) {
		//fclose(current_output_log);
	//}

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
		//current_output_log = fopen(buffer, "w");

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
#ifdef LMAO 
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
							//route_tree_mark_congested_nodes_to_be_ripped(route_trees[net->local_id], g, nullptr);
						}
						route_tree_rip_up_marked(route_trees[net->local_id], g, nullptr, params.pres_fac, true, &lock_perfs[tid]);

						vector<sink_t *> sinks;	
						sinks.reserve(net->sinks.size());
						for (auto &sink : net->sinks) {
							RouteTreeNode sink_rt_node = route_tree_get_rt_node(route_trees[net->local_id], sink.rr_node);
							if (!sink_rt_node)  {
								sinks.push_back(&sink);
							} else {
								const auto &sink_rt_node_p = get_vertex_props(route_trees[net->local_id].graph, sink_rt_node);
								assert(!sink_rt_node_p.pending_rip_up);
							}
						}
						if (!sinks.empty()) {
							//route_net_2(g, net->vpr_id, &net->source, sinks, params, states[tid], route_trees[net->local_id], net_timing[net->vpr_id], true, &perfs[tid], &lock_perfs[tid]);
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
						route_tree_rip_up_marked(route_trees[net->local_id], g, nullptr, params.pres_fac, true, &lock_perfs[tid]);

						vector<sink_t *> sinks;	
						sinks.reserve(net->sinks.size());
						for (auto &sink : net->sinks) {
							RouteTreeNode sink_rt_node = route_tree_get_rt_node(route_trees[net->local_id], sink.rr_node);
							if (!sink_rt_node)  {
								sinks.push_back(&sink);
							} else {
								const auto &sink_rt_node_p = get_vertex_props(route_trees[net->local_id].graph, sink_rt_node);
								assert(!sink_rt_node_p.pending_rip_up);
							}
						}
						if (!sinks.empty()) {
							//route_net_2(g, net->vpr_id, &net->source, sinks, params, states[tid], route_trees[net->local_id], net_timing[net->vpr_id], true, &perfs[tid], &lock_perfs[tid]);
						}
					}
					debug_lock[tid].unlock();
				});
		}

#ifdef LMAO 
		__itt_frame_end_v3(pD, NULL);
#endif
		
		auto route_time = clock::now()-route_start;
		total_route_time += route_time;

		zlog_info(delta_log, "Routing took %g\n", duration_cast<nanoseconds>(route_time).count() / 1e9);
		printf("Routing took %g\n", duration_cast<nanoseconds>(route_time).count() / 1e9);
		for (int i = 0; i < opts->num_threads; ++i) {
			printf("Thread %d num lock waits/tries: %lu/%lu (%g)\n", i, lock_perfs[i].num_lock_waits, lock_perfs[i].num_lock_tries, (float)lock_perfs[i].num_lock_waits*100/lock_perfs[i].num_lock_tries);
		}

		/* checking */
		//for_all_vertices(g, [] (RRNode &v) -> void { v.properties.recalc_occ = 0; });

		for (const auto &net : nets) {
			check_route_tree(route_trees[net.local_id], net, g);
		}

		bool valid = true;
		//for_all_vertices(g, [&valid, &buffer] (RRNode &v) -> void {
				//sprintf_rr_node(id(v), buffer);
				//if (v.properties.recalc_occ != v.properties.occ) {
					//zlog_error(delta_log, "Node %s occ mismatch, recalc: %d original: %d\n", buffer, v.properties.recalc_occ, v.properties.occ);
					//valid = false;
				//}
				//});
		assert(valid);

		for (const auto &net : nets) {
			vector<int> overused_rr_node;
			assert(route_trees[net.local_id].root_rt_nodes.size() == 1);
			//get_overused_nodes(route_trees[net.local_id], route_trees[net.local_id].root_rt_nodes[0], g, congestion, overused_rr_node);
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

		if (true/*feasible_routing(g, nullptr)*/) {
			zlog_info(delta_log, "Routed in %d iterations. Total route time: %g\n", iter+1, duration_cast<nanoseconds>(total_route_time).count() / 1e9);
			printf("Routed in %d iterations. Total route time: %g\n", iter+1, duration_cast<nanoseconds>(total_route_time).count() / 1e9);
			//dump_route(*current_traces_ptr, "route.txt");
			routed = true;
		} else {
			unsigned long num_overused_nodes = 0;
			//for_all_vertices(g, [&num_overused_nodes] (const RRNode &v) -> void {
					//if (v.properties.occ > v.properties.capacity) {
					//++num_overused_nodes;
					//}
					//});
			zlog_info(delta_log, "Num overused nodes: %lu/%d (%.2f)\n", num_overused_nodes, num_vertices(g), num_overused_nodes*100.0/num_vertices(g));
			printf("Num overused nodes: %lu/%d (%.2f)\n", num_overused_nodes, num_vertices(g), num_overused_nodes*100.0/num_vertices(g));

			if ((!rip_up_all && num_overused_nodes > prev_num_overused_nodes && opts->num_threads > 1)) {
				nets_to_partition.clear();
				tbb::spin_mutex lock;
				tbb::parallel_for(tbb::blocked_range<int>(0, nets_to_route.size(), 1024),
						[&] (const tbb::blocked_range<int> &range) {

						for (int i = range.begin(); i != range.end(); ++i) {
							net_t *net = nets_to_route[i].second;

							//bool marked = route_tree_mark_congested_nodes_to_be_ripped(route_trees[net->local_id], g, nullptr);
							//if (marked) {
								//lock.lock();
								//nets_to_partition.push_back(nets_to_route[i]);
								//lock.unlock();
							//}
						}
					});

				overlaps.clear();
				partitions.clear();
				has_interpartition_overlap.clear();
				partition_nets(nets_to_partition, opts->num_threads, 1000, overlaps, partitions, has_interpartition_overlap);
				use_partitioned = true;
				for (int i = 0; i < opts->num_threads; ++i) {
					printf("Partition %d size %lu\n", i, partitions[i].size());
				}
			} else {
				/* check whether we need to do this */
				//use_partitioned = false;
			}

			prev_num_overused_nodes = num_overused_nodes;

			if (iter == 0) {
				params.pres_fac = opts->initial_pres_fac;
				//update_costs(g, params.pres_fac, 0);
			} else {
				params.pres_fac *= opts->pres_fac_mult;

				/* Avoid overflow for high iteration counts, even if acc_cost is big */
				params.pres_fac = std::min(params.pres_fac, static_cast<float>(HUGE_POSITIVE_FLOAT / 1e5));

				//update_costs(g, params.pres_fac, opts->acc_fac);
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

		zlog_info(delta_log, "Iteration time: %g s. Route time: %g s. \n", duration_cast<nanoseconds>(iter_time).count() / 1e9, route_time / 1e9);

		//if (current_output_log && fclose(current_output_log)) {
			//char str[256];
			//strerror_r(errno, str, sizeof(str));
			//printf("failed to close file: %s\n", str);
			//assert(false);
		//}
	}

	if (!routed) {
		printf("Failed to route in %d iterations. Total route time: %g\n", opts->max_router_iterations, duration_cast<nanoseconds>(total_route_time).count() / 1e9);
	}

	return routed;
}

bool locking_route_0(t_router_opts *opts)
{
	tbb::task_scheduler_init init(opts->num_threads);

	using clock = std::chrono::high_resolution_clock;

	//extern clock::time_point start_time;

	//extern FILE *current_output_log;

	//start_time = clock::now();

/*#ifdef __linux__*/
	/*printf("domain flags: %d\n", pD->flags);*/
	/*pD->flags = 1;*/
/*#endif*/
	
	init_logging();
	//zlog_set_record("custom_output", zlog_custom_output);
	//zlog_set_record("sched_custom_output", zlog_sched_custom_output);
	//current_output_log = fopen("/Volumes/DATA/main.log", "w");

	zlog_info(delta_log, "Using boost version %d\n", BOOST_VERSION);

#if (TBB_USE_DEBUG==1)
	zlog_warn(delta_log, "Using debug version of TBB\n");
#endif

	RRGraph g;
	init_graph(g);

	vector<net_t> nets;
	vector<net_t> global_nets;
	init_nets(nets, global_nets, opts->bb_factor, opts->large_bb);

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
			states[i][j].rr_node = RRGraph::null_vertex();
			states[i][j].known_cost = std::numeric_limits<float>::max();
			states[i][j].cost = std::numeric_limits<float>::max();
			states[i][j].prev_edge = RRGraph::null_edge();
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

	//for_all_vertices(g, [] (RRNode &v) -> void {
			//v.properties.acc_cost = 1;
			//v.properties.pres_cost = 1;
			//v.properties.occ = 0;
			//});

	char buffer[256];

	//if (current_output_log) {
		//fclose(current_output_log);
	//}

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
		//current_output_log = fopen(buffer, "w");

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
#ifdef LMAO 
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

					//route_tree_mark_nodes_to_be_ripped(route_trees[net->local_id], g, nullptr, 10000);
					route_tree_rip_up_marked(route_trees[net->local_id], g, nullptr, params.pres_fac, true, &lock_perfs[tid]);

					vector<sink_t *> sinks;	
					sinks.reserve(net->sinks.size());
					for (auto &sink : net->sinks) {
						RouteTreeNode sink_rt_node = route_tree_get_rt_node(route_trees[net->local_id], sink.rr_node);
						if (!sink_rt_node)  {
							sinks.push_back(&sink);
						} else {
							const auto &sink_rt_node_p = get_vertex_props(route_trees[net->local_id].graph, sink_rt_node);
							assert(!sink_rt_node_p.pending_rip_up);
						}
					}
					if (!sinks.empty()) {
						//route_net_2(g, net->vpr_id, &net->source, sinks, params, states[tid], route_trees[net->local_id], net_timing[net->vpr_id], true, &perfs[tid], &lock_perfs[tid]);
					}
				}

				debug_lock[tid].unlock();
				});

#ifdef LMAO 
		__itt_frame_end_v3(pD, NULL);
#endif
		
		auto route_time = clock::now()-route_start;
		total_route_time += route_time;

		zlog_info(delta_log, "Routing took %g\n", duration_cast<nanoseconds>(route_time).count() / 1e9);
		printf("Routing took %g\n", duration_cast<nanoseconds>(route_time).count() / 1e9);
		for (int i = 0; i < opts->num_threads; ++i) {
			printf("Thread %d num lock waits/tries: %lu/%lu (%g)\n", i, lock_perfs[i].num_lock_waits, lock_perfs[i].num_lock_tries, (float)lock_perfs[i].num_lock_waits*100/lock_perfs[i].num_lock_tries);
		}

		/* checking */
		//for_all_vertices(g, [] (RRNode &v) -> void { v.properties.recalc_occ = 0; });

		for (const auto &net : nets) {
			check_route_tree(route_trees[net.local_id], net, g);
		}

		bool valid = true;
		//for_all_vertices(g, [&valid, &buffer] (RRNode &v) -> void {
				//sprintf_rr_node(id(v), buffer);
				//if (v.properties.recalc_occ != v.properties.occ) {
					//zlog_error(delta_log, "Node %s occ mismatch, recalc: %d original: %d\n", buffer, v.properties.recalc_occ, v.properties.occ);
					//valid = false;
				//}
				//});
		assert(valid);

		for (const auto &net : nets) {
			vector<int> overused_rr_node;
			//get_overused_nodes(route_trees[net.local_id], route_trees[net.local_id].root_rt_nodes[0], g, congestion, overused_rr_node);
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

		if (true /*feasible_routing(g, nullptr)*/) {
			zlog_info(delta_log, "Routed in %d iterations. Total route time: %g\n", iter+1, duration_cast<nanoseconds>(total_route_time).count() / 1e9);
			printf("Routed in %d iterations. Total route time: %g\n", iter+1, duration_cast<nanoseconds>(total_route_time).count() / 1e9);
			//dump_route(*current_traces_ptr, "route.txt");
			routed = true;
		} else {
			int num_overused_nodes = 0;
			//for_all_vertices(g, [&num_overused_nodes] (const RRNode &v) -> void {
					//if (v.properties.occ > v.properties.capacity) {
					//++num_overused_nodes;
					//}
					//});
			zlog_info(delta_log, "Num overused nodes: %d/%d (%.2f)\n", num_overused_nodes, num_vertices(g), num_overused_nodes*100.0/num_vertices(g));
			printf("Num overused nodes: %d/%d (%.2f)\n", num_overused_nodes, num_vertices(g), num_overused_nodes*100.0/num_vertices(g));

			if (iter == 0) {
				params.pres_fac = opts->initial_pres_fac;
				//update_costs(g, params.pres_fac, 0);
			} else {
				params.pres_fac *= opts->pres_fac_mult;

				/* Avoid overflow for high iteration counts, even if acc_cost is big */
				params.pres_fac = std::min(params.pres_fac, static_cast<float>(HUGE_POSITIVE_FLOAT / 1e5));

				//update_costs(g, params.pres_fac, opts->acc_fac);
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

		zlog_info(delta_log, "Iteration time: %g s. Route time: %g s. \n", duration_cast<nanoseconds>(iter_time).count() / 1e9, route_time / 1e9);

		//if (current_output_log && fclose(current_output_log)) {
			//char str[256];
			//strerror_r(errno, str, sizeof(str));
			//printf("failed to close file: %s\n", str);
			//assert(false);
		//}
	}

	if (!routed) {
		printf("Failed to route in %d iterations. Total route time: %g\n", opts->max_router_iterations, duration_cast<nanoseconds>(total_route_time).count() / 1e9);
	}

	return routed;
}

bool locking_route_2(t_router_opts *opts)
{
	tbb::task_scheduler_init init(opts->num_threads);

	using clock = std::chrono::high_resolution_clock;

	//extern clock::time_point start_time;

	//extern FILE *current_output_log;

	//start_time = clock::now();

/*#ifdef __linux__*/
	/*printf("domain flags: %d\n", pD->flags);*/
	/*pD->flags = 1;*/
/*#endif*/
	
	init_logging();
	//zlog_set_record("custom_output", zlog_custom_output);
	//zlog_set_record("sched_custom_output", zlog_sched_custom_output);
	//current_output_log = fopen("/Volumes/DATA/main.log", "w");

	zlog_info(delta_log, "Using boost version %d\n", BOOST_VERSION);

#if (TBB_USE_DEBUG==1)
	zlog_warn(delta_log, "Using debug version of TBB\n");
#endif

	RRGraph g;
	init_graph(g);

	vector<net_t> nets;
	vector<net_t> global_nets;
	init_nets(nets, global_nets, opts->bb_factor, opts->large_bb);

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
			states[i][j].rr_node = RRGraph::null_vertex();
			states[i][j].known_cost = std::numeric_limits<float>::max();
			states[i][j].cost = std::numeric_limits<float>::max();
			states[i][j].prev_edge = RRGraph::null_edge();
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

	//for_all_vertices(g, [] (RRNode &v) -> void {
			//v.properties.acc_cost = 1;
			//v.properties.pres_cost = 1;
			//v.properties.occ = 0;
			//});

	char buffer[256];

	//if (current_output_log) {
		//fclose(current_output_log);
	//}

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
		//current_output_log = fopen(buffer, "w");

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
#ifdef LMAO 
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

					//route_tree_mark_nodes_to_be_ripped(route_trees[net->local_id], g, nullptr, 10000);
					route_tree_rip_up_marked(route_trees[net->local_id], g, nullptr, params.pres_fac, true, nullptr);

					vector<sink_t *> sinks;	
					sinks.reserve(net->sinks.size());
					for (auto &sink : net->sinks) {
						RouteTreeNode sink_rt_node = route_tree_get_rt_node(route_trees[net->local_id], sink.rr_node);
						if (!sink_rt_node)  {
							sinks.push_back(&sink);
						} else {
							const auto &sink_rt_node_p = get_vertex_props(route_trees[net->local_id].graph, sink_rt_node);
							assert(!sink_rt_node_p.pending_rip_up);
						}
					}
					if (!sinks.empty()) {
						//route_net_2(g, net->vpr_id, &net->source, sinks, params, states[tid], route_trees[net->local_id], net_timing[net->vpr_id], true, nullptr, nullptr);
					}
				}
				});

#ifdef LMAO 
		__itt_frame_end_v3(pD, NULL);
#endif
		
		auto route_time = clock::now()-route_start;
		total_route_time += route_time;
		/*zlog_info(delta_log, "Scheduling took %s\n", format(elapsed).c_str());*/
		zlog_info(delta_log, "Routing took %g\n", duration_cast<nanoseconds>(route_time).count() / 1e9);
		printf("Routing took %g\n", duration_cast<nanoseconds>(route_time).count() / 1e9);
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
		//for_all_vertices(g, [] (RRNode &v) -> void { v.properties.recalc_occ = 0; });

		for (const auto &net : nets) {
			check_route_tree(route_trees[net.local_id], net, g);
		}

		bool valid = true;
		//for_all_vertices(g, [&valid, &buffer] (RRNode &v) -> void {
				//sprintf_rr_node(id(v), buffer);
				//if (v.properties.recalc_occ != v.properties.occ) {
					//zlog_error(delta_log, "Node %s occ mismatch, recalc: %d original: %d\n", buffer, v.properties.recalc_occ, v.properties.occ);
					//valid = false;
				//}
				//});
		assert(valid);

		for (const auto &net : nets) {
			vector<int> overused_rr_node;
			//get_overused_nodes(route_trees[net.local_id], route_trees[net.local_id].root_rt_nodes[0], g, congestion, overused_rr_node);
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

		if (true/*feasible_routing(g, nullptr)*/) {
			zlog_info(delta_log, "Routed in %d iterations. Total route time: %g\n", iter+1, duration_cast<nanoseconds>(total_route_time).count() / 1e9);
			printf("Routed in %d iterations. Total route time: %g\n", iter+1, duration_cast<nanoseconds>(total_route_time).count() / 1e9);
			//dump_route(*current_traces_ptr, "route.txt");
			routed = true;
		} else {
			int num_overused_nodes = 0;
			//for_all_vertices(g, [&num_overused_nodes] (const RRNode &v) -> void {
					//if (v.properties.occ > v.properties.capacity) {
					//++num_overused_nodes;
					//}
					//});
			zlog_info(delta_log, "Num overused nodes: %d/%d (%.2f)\n", num_overused_nodes, num_vertices(g), num_overused_nodes*100.0/num_vertices(g));
			printf("Num overused nodes: %d/%d (%.2f)\n", num_overused_nodes, num_vertices(g), num_overused_nodes*100.0/num_vertices(g));

			if (iter == 0) {
				params.pres_fac = opts->initial_pres_fac;
				//update_costs(g, params.pres_fac, 0);
			} else {
				params.pres_fac *= opts->pres_fac_mult;

				/* Avoid overflow for high iteration counts, even if acc_cost is big */
				params.pres_fac = std::min(params.pres_fac, static_cast<float>(HUGE_POSITIVE_FLOAT / 1e5));

				//update_costs(g, params.pres_fac, opts->acc_fac);
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

		zlog_info(delta_log, "Iteration time: %g s. Route time: %g s. \n", duration_cast<nanoseconds>(iter_time).count() / 1e9, route_time / 1e9);

		//if (current_output_log && fclose(current_output_log)) {
			//char str[256];
			//strerror_r(errno, str, sizeof(str));
			//printf("failed to close file: %s\n", str);
			//assert(false);
		//}
	}

	if (!routed) {
		printf("Failed to route in %d iterations. Total route time: %g\n", opts->max_router_iterations, duration_cast<nanoseconds>(total_route_time).count() / 1e9);
	}

	return routed;
}

bool locking_route_3(t_router_opts *opts)
{
	tbb::task_scheduler_init init(opts->num_threads);

	using clock = std::chrono::high_resolution_clock;

	//extern clock::time_point start_time;

	//extern FILE *current_output_log;

	//start_time = clock::now();

/*#ifdef __linux__*/
	/*printf("domain flags: %d\n", pD->flags);*/
	/*pD->flags = 1;*/
/*#endif*/
	
	init_logging();
	//zlog_set_record("custom_output", zlog_custom_output);
	//zlog_set_record("sched_custom_output", zlog_sched_custom_output);
	//current_output_log = fopen("/Volumes/DATA/main.log", "w");

	zlog_info(delta_log, "Using boost version %d\n", BOOST_VERSION);

#if (TBB_USE_DEBUG==1)
	zlog_warn(delta_log, "Using debug version of TBB\n");
#endif

	RRGraph g;
	init_graph(g);

	vector<net_t> nets;
	vector<net_t> global_nets;
	init_nets(nets, global_nets, opts->bb_factor, opts->large_bb);

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
			states[i][j].rr_node = RRGraph::null_vertex();
			states[i][j].known_cost = std::numeric_limits<float>::max();
			states[i][j].cost = std::numeric_limits<float>::max();
			states[i][j].prev_edge = RRGraph::null_edge();
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

	//for_all_vertices(g, [] (RRNode &v) -> void {
			//v.properties.acc_cost = 1;
			//v.properties.pres_cost = 1;
			//v.properties.occ = 0;
			//});

	char buffer[256];

	//if (current_output_log) {
		//fclose(current_output_log);
	//}

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
		//current_output_log = fopen(buffer, "w");

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
#ifdef LMAO 
		__itt_frame_begin_v3(pD, NULL);
#endif

		/*std::shuffle(begin(virtual_nets), end(virtual_nets), mt);*/
		tbb::parallel_for(tbb::blocked_range<int>(0, opts->num_threads, 1),
				[&] (const tbb::blocked_range<int> &range) {
				assert(range.end()-range.begin() == 1);

				int tid = range.begin();
				for (int i = 0; i < partitions[tid].size(); ++i) {
					net_t *net = nets_to_partition[partitions[tid][i]].second;

					//route_tree_mark_nodes_to_be_ripped(route_trees[net->local_id], g, nullptr, 5);
					route_tree_rip_up_marked(route_trees[net->local_id], g, nullptr, params.pres_fac, true, nullptr);

					vector<sink_t *> sinks;	
					sinks.reserve(net->sinks.size());
					for (auto &sink : net->sinks) {
						RouteTreeNode sink_rt_node = route_tree_get_rt_node(route_trees[net->local_id], sink.rr_node);
						if (!sink_rt_node)  {
							sinks.push_back(&sink);
						} else {
							const auto &sink_rt_node_p = get_vertex_props(route_trees[net->local_id].graph, sink_rt_node);
							assert(!sink_rt_node_p.pending_rip_up);
						}
					}
					if (!sinks.empty()) {
						//route_net_2(g, net->vpr_id, &net->source, sinks, params, states[tid], route_trees[net->local_id], net_timing[net->vpr_id], true, nullptr, nullptr);
					}
				}
				});

#ifdef LMAO 
		__itt_frame_end_v3(pD, NULL);
#endif
		
		auto route_time = clock::now()-route_start;
		total_route_time += route_time;
		/*zlog_info(delta_log, "Scheduling took %s\n", format(elapsed).c_str());*/
		zlog_info(delta_log, "Routing took %g\n", duration_cast<nanoseconds>(route_time).count() / 1e9);
		printf("Routing took %g\n", duration_cast<nanoseconds>(route_time).count() / 1e9);
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
		//for_all_vertices(g, [] (RRNode &v) -> void { v.properties.recalc_occ = 0; });

		for (const auto &net : nets) {
			check_route_tree(route_trees[net.local_id], net, g);
		}

		bool valid = true;
		//for_all_vertices(g, [&valid, &buffer] (RRNode &v) -> void {
				//sprintf_rr_node(id(v), buffer);
				//if (v.properties.recalc_occ != v.properties.occ) {
					//zlog_error(delta_log, "Node %s occ mismatch, recalc: %d original: %d\n", buffer, v.properties.recalc_occ, v.properties.occ);
					//valid = false;
				//}
				//});
		assert(valid);

		for (const auto &net : nets) {
			vector<int> overused_rr_node;
			//get_overused_nodes(route_trees[net.local_id], route_trees[net.local_id].root_rt_nodes[0], g, congestion, overused_rr_node);
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

		if (true/*feasible_routing(g, nullptr)*/) {
			zlog_info(delta_log, "Routed in %d iterations. Total route time: %g\n", iter+1, duration_cast<nanoseconds>(total_route_time).count() / 1e9);
			printf("Routed in %d iterations. Total route time: %g\n", iter+1, duration_cast<nanoseconds>(total_route_time).count() / 1e9);
			//dump_route(*current_traces_ptr, "route.txt");
			routed = true;
		} else {
			int num_overused_nodes = 0;
			//for_all_vertices(g, [&num_overused_nodes] (const RRNode &v) -> void {
					//if (v.properties.occ > v.properties.capacity) {
					//++num_overused_nodes;
					//}
					//});
			zlog_info(delta_log, "Num overused nodes: %d/%d (%.2f)\n", num_overused_nodes, num_vertices(g), num_overused_nodes*100.0/num_vertices(g));
			printf("Num overused nodes: %d/%d (%.2f)\n", num_overused_nodes, num_vertices(g), num_overused_nodes*100.0/num_vertices(g));

			if (iter == 0) {
				params.pres_fac = opts->initial_pres_fac;
				//update_costs(g, params.pres_fac, 0);
			} else {
				params.pres_fac *= opts->pres_fac_mult;

				/* Avoid overflow for high iteration counts, even if acc_cost is big */
				params.pres_fac = std::min(params.pres_fac, static_cast<float>(HUGE_POSITIVE_FLOAT / 1e5));

				//update_costs(g, params.pres_fac, opts->acc_fac);
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

		zlog_info(delta_log, "Iteration time: %g s. Route time: %g s. \n", duration_cast<nanoseconds>(iter_time).count() / 1e9, route_time / 1e9);

		//if (current_output_log && fclose(current_output_log)) {
			//char str[256];
			//strerror_r(errno, str, sizeof(str));
			//printf("failed to close file: %s\n", str);
			//assert(false);
		//}
	}

	if (!routed) {
		printf("Failed to route in %d iterations. Total route time: %g\n", opts->max_router_iterations, duration_cast<nanoseconds>(total_route_time).count() / 1e9);
	}

	return routed;
}