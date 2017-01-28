#include "pch.h"

#include "vpr_types.h"
#include "path_delay.h"
#include "rr_graph.h"

#include "log.h"
#include "barrier.h"
#include "graph.h"
#include "fast_graph.h"
#include "filtered_graph.h"
#include "route.h"
#include "route_tree.h"
#include "trace.h"
#include "old_misr.h"
#include "geometry.h"
#include "quadtree.h"
#include "utility.h"
#include "net_cluster.h"
#include "args.h"
#include "init.h"
#include "router.h"
#include "congestion.h"
#include "metis_partitioner.h"
#include "fm.h"

bool mpi_spatial_route(t_router_opts *opts, struct s_det_routing_arch det_routing_arch, t_direct_inf *directs, int num_directs, t_segment_inf *segment_inf, t_timing_inf timing_inf)
{
    using std::chrono::duration_cast;
    using std::chrono::nanoseconds;
    using clock = std::chrono::high_resolution_clock;

    init_logging();
    zlog_set_record("custom_output", delta_log_output);
    zlog_set_record("missing_edge", missing_edge_log_output);
    zlog_set_record("ss", ss_log_output);
    delta_log_files.resize(opts->max_router_iterations);
    for (int i = 0; i < opts->max_router_iterations; ++i) {
        delta_log_files[i].resize(opts->num_threads, nullptr);
    }
    missing_edge_log_files.resize(opts->max_router_iterations);
    for (int i = 0; i < opts->max_router_iterations; ++i) {
        missing_edge_log_files[i].resize(opts->num_threads, nullptr);
    }

    //fclose(fopen(LOG_PATH_PREFIX"iter__tid_.log", "w"));

    test_fast_graph();
    test_topo();
    test_fm();
    test_filter_graph();
    test_partition_graph();
    test_connected_components();

    vector<RRGraph *> graphs;
    rr_graph_partitioner partitioner;
    partitioner.init(opts, &det_routing_arch, directs, num_directs, segment_inf, &timing_inf);
    partitioner.partition(opts->num_threads, graphs);
    //for (const auto &g : graphs) {
        //routability(*g);
    //}
    //
    RRGraph combined_g;
    add_vertex(combined_g, num_vertices(partitioner.orig_g));

    for (const auto &g : graphs) {
        for (const auto &e : get_edges(*g)) {
            int from = get_source(*g, e);
            int to = get_target(*g, e);
            const auto &from_ver = get_vertex_props(*g, from);
            const auto &to_ver = get_vertex_props(*g, to);

            if (is_channel(from_ver) && is_channel(to_ver)) {
                assert(!has_edge(combined_g, from, to));
                add_edge(combined_g, from, to);
            }
        }
    }

    for (const auto &e : get_edges(partitioner.orig_g)) {
        int from = get_source(partitioner.orig_g, e);
        int to = get_target(partitioner.orig_g, e);
        const auto &from_ver = get_vertex_props(partitioner.orig_g, from);
        const auto &to_ver = get_vertex_props(partitioner.orig_g, to);
        if (!is_channel(from_ver) || !is_channel(to_ver)) {
            assert(!has_edge(combined_g, from, to));
            add_edge(combined_g, from, to);
        }
    }

    printf("Combined/Orig graph has %d/%d (%g) edges.\n", num_edges(combined_g), num_edges(partitioner.orig_g), 100.0*num_edges(combined_g)/num_edges(partitioner.orig_g));

    //goto lol;

    //extern int num_types;
    //extern struct s_type_descriptor *type_descriptors;
    //extern int nx, ny;
    //extern struct s_grid_tile **grid;

    //free_rr_graph();

    //int warnings;

    //build_rr_graph(GRAPH_UNIDIR, num_types, type_descriptors, nx, ny, grid,
            //opts->fixed_channel_width, NULL, det_routing_arch.switch_block_type,
            //det_routing_arch.Fs, det_routing_arch.num_segment,
            //det_routing_arch.num_switch, segment_inf,
            //det_routing_arch.global_route_switch,
            //det_routing_arch.delayless_switch, timing_inf,
            //det_routing_arch.wire_to_ipin_switch, opts->base_cost_type,
            //directs, num_directs, FALSE,
            //&warnings);

    //RRGraph channel_with_interior_g;
    //init_channel_only_graph(channel_with_interior_g);

    //dump_rr_graph(channel_with_interior_g, "/Volumes/DATA/rr_graph_channel_with_interior.txt");

    //RRGraph orig_g;
    //init_graph(orig_g);

    //dump_rr_graph(orig_g, "/Volumes/DATA/rr_graph.txt");

    //free_rr_graph();
    //for (int i = 0; i < det_routing_arch.num_segment; ++i) {
        //for (int j = 0; j < segment_inf[i].sb_len; ++j) {
            //if (j != 0 && j != segment_inf[i].sb_len-1) {
                //segment_inf[i].sb[j] = FALSE;
            //}
        //}
    //}
    //build_rr_graph(GRAPH_UNIDIR, num_types, type_descriptors, nx, ny, grid,
            //opts->fixed_channel_width, NULL, det_routing_arch.switch_block_type,
            //det_routing_arch.Fs, det_routing_arch.num_segment,
            //det_routing_arch.num_switch, segment_inf,
            //det_routing_arch.global_route_switch,
            //det_routing_arch.delayless_switch, timing_inf,
            //det_routing_arch.wire_to_ipin_switch, opts->base_cost_type,
            //directs, num_directs, FALSE,
            //&warnings);

    //RRGraph channel_without_interior_g;
    //init_channel_only_graph(channel_without_interior_g);

    //dump_rr_graph(channel_without_interior_g, "/Volumes/DATA/rr_graph_channel_without_interior.txt");

    //vector<vector<vector<int>>> all_partition_components(opts->num_threads);
    ////for (int x = 0; x < nx+1; ++x) {
        ////for (int y = 0; y < ny+1; ++y) {
            ////init_partitioned_graph_5(opts->num_threads, channel_with_interior_g, channel_without_interior_g, orig_g, graphs, x, y);
        ////}
    ////}
    ////init_partitioned_graph_4(opts->num_threads, channel_with_interior_g, channel_without_interior_g, orig_g, graphs, 14, 14);
    //init_partitioned_graph_4(opts->num_threads, channel_with_interior_g, channel_without_interior_g, orig_g, graphs, all_partition_components);

    //for (int i = 0; i < graphs.size(); ++i) {
        //char filename[256];

        //sprintf(filename, "/Volumes/DATA/rr_graph_%d.txt", i);
        //dump_rr_graph(*graphs[i], filename);

        //sprintf(filename, "/Volumes/DATA/rr_graph_%d_edges.txt", i);
        //dump_edges(*graphs[i], filename);
    //}

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

	vector<route_state_t *> states;
	congestion_t *congestion;
	vector<route_tree_t> route_trees;
	t_net_timing *net_timing;
	init_route_structs(partitioner.orig_g, nets, global_nets, opts->num_threads, states, &congestion, route_trees, &net_timing);

    route_parameters_t params;
    params.criticality_exp = opts->criticality_exp;
    params.astar_fac = opts->astar_fac;
    params.max_criticality = opts->max_criticality;
    params.bend_cost = opts->bend_cost;
	params.pres_fac = opts->first_iter_pres_fac; /* Typically 0 -> ignore cong. */

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

    int iter;
    float crit_path_delay;
	int current_level = 0;

    for (iter = 0; iter < opts->max_router_iterations && !routed; ++iter) {
        clock::duration greedy_route_time = clock::duration::zero();
        clock::duration partitioned_route_time = clock::duration::zero();
        clock::duration update_cost_time = clock::duration::zero();
        clock::duration partitioning_time = clock::duration::zero();
        clock::duration analyze_timing_time = clock::duration::zero();
        clock::duration iter_time = clock::duration::zero();

        //for (int i = 0; i < opts->num_threads; ++i) {
            //sprintf(buffer, LOG_PATH_PREFIX"iter_%d_tid_%d.log", iter, i);
            //fclose(fopen(buffer, "w"));
        //}

        zlog_info(delta_log, "Routing iteration: %d\n", iter);
        printf("Routing iteration: %d\n", iter);
        
        auto iter_start = clock::now();

        //auto route_start = clock::now();

        for (auto &net : nets) {
            update_sink_criticalities(net, net_timing[net.vpr_id], params);
        }

        //tbb::enumerable_thread_specific<state_t *> state_tls;
        //
        tbb::atomic<int> net_index = 0;
        vector<tbb::spin_mutex> debug_lock(opts->num_threads);
        vector<int> thread_num_nets_routed(opts->num_threads);
        vector<int> thread_num_nets_to_route(opts->num_threads);
        vector<int> thread_num_sinks_routed(opts->num_threads);
        vector<int> thread_num_sinks_to_route(opts->num_threads);
        vector<int> thread_bfs_num_sinks_routed(opts->num_threads);
		vector<int> thread_num_interpartition_sinks(opts->num_threads);
        vector<int> graph_used_by_net(nets.size(), -1);
        vector<vector<pseudo_net_t *>> partition_pseudo_nets_0(opts->num_threads);
        vector<vector<pseudo_net_t *>> partition_pseudo_nets_1(opts->num_threads);
        vector<vector<pseudo_net_t *>> *current_partition_pseudo_nets = &partition_pseudo_nets_0;
        vector<vector<pseudo_net_t *>> *new_partition_pseudo_nets = &partition_pseudo_nets_1;
        vector<tbb::spin_mutex> pseudo_nets_locks(opts->num_threads);
        vector<vector<pseudo_net_t *>> net_pseudo_nets(nets.size());
        tbb::atomic<bool> has_interpartition_sinks = false;
        vector<vector<interpartition_sink_t>> interpartition_sinks(nets.size()); 
        vector<int> net_next_pid(nets.size(), -1);
        vector<int> net_initial_pid(nets.size(), -1);

        auto greedy_route_start = clock::now();
		int num_threads = pow(2, partitioner.num_levels-current_level-1);

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

                int local_num_nets_to_route = 0;
                int local_num_nets_routed = 0;
                int local_num_sinks_routed = 0;
                int local_num_sinks_to_route = 0;
                int local_bfs_num_sinks_routed = 0;
				int local_num_interpartition_sinks = 0;

				if (tid < num_threads) {
					int i = tid;
					//while ((i = net_index++) < nets_to_route.size()) {
					while (i < nets_to_route.size()) {
						net_t *net = nets_to_route[i].second;

						//auto rip_up_start = clock::now();
						//if (greedy_rip_up_all) {
						route_tree_mark_all_nodes_to_be_ripped(route_trees[net->local_id], *graphs[tid]);
						//} else {
						//route_tree_mark_congested_nodes_to_be_ripped(route_trees[net->local_id], *graphs[tid], congestion);
						//}
						route_tree_rip_up_marked(route_trees[net->local_id], *graphs[tid], congestion, params.pres_fac, true, &local_lock_perf);

						//local_perf.total_rip_up_time += clock::now()-rip_up_start;

						//auto route_start = clock::now();

						vector<sink_t *> sinks;	
						get_sinks_to_route(net, route_trees[net->local_id], sinks);

						if (!sinks.empty()) {
							route_net_4(partitioner.orig_g, partitioner.result_pid_by_level[current_level], tid, net->vpr_id, &net->source, sinks, params, states[tid], congestion, route_trees[net->local_id], net_timing[net->vpr_id], interpartition_sinks[net->local_id], true, &local_perf, &local_lock_perf);

							if (!has_interpartition_sinks) {
								has_interpartition_sinks = !interpartition_sinks[net->local_id].empty();
							}

							local_num_interpartition_sinks += interpartition_sinks[net->local_id].size();
							++local_num_nets_routed;
							++local_num_nets_to_route;
							local_num_sinks_to_route += sinks.size();

							local_num_sinks_routed += sinks.size();
						}

						//local_perf.total_route_time += clock::now()-rip_up_start;
						i += num_threads;
					}
				}

                greedy_end_time[tid] = clock::now();

                perfs[tid] = local_perf;
                lock_perfs[tid] = local_lock_perf; 
                thread_num_nets_routed[tid] = local_num_nets_routed;
                thread_num_nets_to_route[tid] = local_num_nets_to_route;
                thread_num_sinks_routed[tid] = local_num_sinks_routed;
                thread_num_sinks_to_route[tid] = local_num_sinks_to_route;
                thread_bfs_num_sinks_routed[tid] = local_bfs_num_sinks_routed;
                thread_num_interpartition_sinks[tid] = local_num_interpartition_sinks;

                //debug_lock[tid].unlock();
                });

        greedy_route_time = clock::now()-greedy_route_start;

        auto partitioned_route_start = clock::now();

		//if (has_interpartition_sinks) {
			//for (const auto &net_interpartition_sinks : interpartition_sinks) {
				//for (const auto &is : net_interpartition_sinks) {
					//vector<RRNode> added_nodes;
					//for (const auto &node : is.path) {
						//if (node.update_cost) {
							//added_nodes.push_back(node.rr_node_id);
						//}
					//}
					//update_one_cost(partitioner.orig_g, congestion, added_nodes.begin(), added_nodes.end(), 1, params.pres_fac, false, nullptr);
				//}
			//}
		//}

        partitioned_route_time = clock::now()-partitioned_route_start;

        //if (greedy_rip_up_all) {
            //next_greedy_rip_up_iter += greedy_rip_up_all_period;
            //++greedy_rip_up_all_period;
            //prev_num_overused_nodes = std::numeric_limits<unsigned long>::max();
        //}

        //route_time = clock::now()-route_start;

        iter_time = clock::now()-iter_start;

        int total_num_nets_routed = 0;
        int total_num_nets_to_route = 0;
        int total_num_sinks_routed = 0;
        int total_num_sinks_to_route = 0;
        int total_bfs_num_sinks_routed = 0;
		int total_num_interpartition_sinks = 0;
        for (int i = 0; i < opts->num_threads; ++i) {
            total_num_nets_to_route += thread_num_nets_to_route[i];
            total_num_nets_routed += thread_num_nets_routed[i];
            total_num_sinks_routed += thread_num_sinks_routed[i];
            total_num_sinks_to_route += thread_num_sinks_to_route[i];
            total_bfs_num_sinks_routed += thread_bfs_num_sinks_routed[i];
			total_num_interpartition_sinks += thread_num_interpartition_sinks[i];
        }

		printf("num nets routed: ");
		for (int i = 0; i < opts->num_threads; ++i) {
			printf("%d/%d (%g), ", thread_num_nets_routed[i], thread_num_nets_to_route[i], thread_num_nets_routed[i]*100.0/thread_num_nets_to_route[i]);
		}
		printf("\n");

		printf("num sinks routed: ");
		for (int i = 0; i < opts->num_threads; ++i) {
			printf("%d/%d (%g), ", thread_num_sinks_routed[i], thread_num_sinks_to_route[i], thread_num_sinks_routed[i]*100.0/thread_num_sinks_to_route[i]);
		}
		printf("\n");

        printf("Total num nets routed: %d/%d (%g)\n", total_num_nets_routed, total_num_nets_to_route, total_num_nets_routed*100.0/total_num_nets_to_route);
        printf("Total num interpartition sinks routed: %d/%d (%g)\n", total_num_interpartition_sinks, total_num_sinks_to_route, total_num_interpartition_sinks*100.0/total_num_sinks_to_route);
        printf("Total num sinks routed: %d/%d (%g)\n", total_bfs_num_sinks_routed+total_num_sinks_routed, total_num_sinks_to_route, (total_bfs_num_sinks_routed+total_num_sinks_routed)*100.0/total_num_sinks_to_route);

		/* checking */
        for (int i = 0; i < num_vertices(*graphs[0]); ++i) {
            congestion[i].recalc_occ = 0; 
        }

        for (const auto &net : nets) {
            check_route_tree(route_trees[net.local_id], net, *graphs[0]);
            recalculate_occ(route_trees[net.local_id], *graphs[0], congestion);
        }

        bool valid = true;
        for (int i = 0; i < num_vertices(*graphs[0]); ++i) {
            sprintf_rr_node(i, buffer);
            if (congestion[i].recalc_occ != congestion[i].occ) {
                zlog_error(delta_log, "Node %s occ mismatch, recalc: %d original: %d\n", buffer, congestion[i].recalc_occ, congestion[i].occ);
                valid = false;
            }
        }
        assert(valid);

        int overused_total_bb_rank = 0;
        int num_congested_nets = 0;
        for (const auto &net : nets) {
            vector<int> overused_rr_node;
            assert(route_trees[net.local_id].root_rt_nodes.size() == 1);
            get_overused_nodes(route_trees[net.local_id], route_trees[net.local_id].root_rt_nodes[0], *graphs[0], congestion, overused_rr_node);
            if (!overused_rr_node.empty()) {
                zlog_level(delta_log, ROUTER_V1, "Net %d bb_rank %d has %lu overused nodes:\n", net.vpr_id, net.bb_area_rank, overused_rr_node.size());
                for (const auto &item : overused_rr_node) {
                    zlog_level(delta_log, ROUTER_V1, "%d ", item);
                }
                zlog_level(delta_log, ROUTER_V1, "\n");
                overused_total_bb_rank += net.bb_area_rank;
                ++num_congested_nets;
            }
        }

        zlog_level(delta_log, ROUTER_V1, "Average overused net bb rank: %g\n", (float)overused_total_bb_rank/num_congested_nets);
        zlog_level(delta_log, ROUTER_V1, "Num congested nets: %d\n", num_congested_nets);

        iter_start = clock::now();

        if (feasible_routing(partitioner.orig_g, congestion)) {
            //dump_route(*current_traces_ptr, "route.txt");
            routed = true;
        } else {
            unsigned long num_overused_nodes = 0;
            for (int i = 0; i < num_vertices(*graphs[0]); ++i) {
                if (congestion[i].occ > get_vertex_props(*graphs[0], i).capacity) {
                    ++num_overused_nodes;
                }
            }
            zlog_info(delta_log, "Num overused nodes: %lu/%d (%.2f)\n", num_overused_nodes, num_vertices(*graphs[0]), num_overused_nodes*100.0/num_vertices(*graphs[0]));
            printf("Num overused nodes: %lu/%d (%.2f)\n", num_overused_nodes, num_vertices(*graphs[0]), num_overused_nodes*100.0/num_vertices(*graphs[0]));

			//if ((float)total_num_interpartition_sinks/total_num_sinks_to_route > 0.1) {
				//++current_level;
				//assert(current_level < partitioner.result_pid_by_level.size());
			//}

            auto update_cost_start = clock::now();

            if (iter == 0) {
                params.pres_fac = opts->initial_pres_fac;
                update_costs(*graphs[0], congestion, params.pres_fac, 0);
            } else {
                params.pres_fac *= opts->pres_fac_mult;

				/* Avoid overflow for high iteration counts, even if acc_cost is big */
                params.pres_fac = std::min(params.pres_fac, static_cast<float>(HUGE_POSITIVE_FLOAT / 1e5));

                update_costs(*graphs[0], congestion, params.pres_fac, opts->acc_fac);
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
                printf("\t\tGreedy route time: %g s (%g).\n", duration_cast<nanoseconds>(greedy_route_time).count() / 1e9, 100.0*duration_cast<nanoseconds>(greedy_route_time).count()/duration_cast<nanoseconds>(greedy_route_time+partitioned_route_time).count());
                printf("\t\tPartitioned route time: %g s (%g).\n", duration_cast<nanoseconds>(partitioned_route_time).count() / 1e9, 100.0*duration_cast<nanoseconds>(partitioned_route_time).count()/duration_cast<nanoseconds>(greedy_route_time+partitioned_route_time).count());
            printf("\tUpdate cost time: %g s.\n", duration_cast<nanoseconds>(update_cost_time).count() / 1e9);
            printf("\tPartitioning time: %g s.\n", duration_cast<nanoseconds>(partitioning_time).count() / 1e9);
            printf("\tAnalyze timing time: %g s.\n", duration_cast<nanoseconds>(analyze_timing_time).count() / 1e9);
        printf("Critical path: %g ns\n", crit_path_delay);

        unsigned long total_num_heap_pushes = 0;
        unsigned long total_num_heap_pops = 0;
        unsigned long total_num_neighbor_visits = 0;

        for (int i = 0; i < opts->num_threads; ++i) {
            total_num_heap_pushes += perfs[i].num_heap_pushes;
            total_num_heap_pops += perfs[i].num_heap_pops;
            total_num_neighbor_visits += perfs[i].num_neighbor_visits;
        }

        printf("num lock waits/tries: ");
        for (int i = 0; i < opts->num_threads; ++i) {
			printf("%lu/%lu (%g) ", lock_perfs[i].num_lock_waits, lock_perfs[i].num_lock_tries, (float)lock_perfs[i].num_lock_waits*100/lock_perfs[i].num_lock_tries); 
		}
		printf("\n");

		printf("total wait time: ");
		for (int i = 0; i < opts->num_threads; ++i) {
			printf("%g (%g) ", duration_cast<nanoseconds>(lock_perfs[i].total_wait_time).count() / 1e9, 100.0*lock_perfs[i].total_wait_time/(greedy_route_time+partitioned_route_time));
		}
		printf("\n");

		printf("num_heap_pushes: ");
		for (int i = 0; i < opts->num_threads; ++i) {
			printf("%lu ", perfs[i].num_heap_pushes);
		}
		printf("\n");

		printf("num_heap_pops: ");
		for (int i = 0; i < opts->num_threads; ++i) {
			printf("%lu ", perfs[i].num_heap_pops);
		}
		printf("\n");

		printf("num_neighbor_visits: ");
		for (int i = 0; i < opts->num_threads; ++i) {
			printf("%lu ", perfs[i].num_neighbor_visits);
		}
		printf("\n");

        printf("total_num_heap_pushes: %lu\n", total_num_heap_pushes);
        printf("total_num_heap_pops: %lu\n", total_num_heap_pops); 
        printf("total_num_neighbor_visits: %lu\n", total_num_neighbor_visits);

        clock::time_point greedy_earliest_end_time = *std::min_element(begin(greedy_end_time), begin(greedy_end_time)+num_threads);

        printf("greedy wait time: ");
        for (int i = 0; i < opts->num_threads; ++i) {
			printf("%g (%g), ", duration_cast<nanoseconds>(greedy_end_time[i]-greedy_earliest_end_time).count() / 1e9, 100.0*(greedy_end_time[i]-greedy_earliest_end_time)/greedy_route_time);
        }
		printf("\n");
    }

    //sprintf(buffer, "%s_run_%d_rr_graph_occ.txt", s_circuit_name, run);
    //dump_rr_graph_occ(congestion, num_vertices(g), buffer);

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

    //delete_graph(g);
    delete_net_timing(nets, global_nets, net_timing);	
    delete [] congestion;
    delete [] net_timing;
    for (int i = 0; i < opts->num_threads; ++i) {
        delete [] states[i];
    }

    return routed;
}
