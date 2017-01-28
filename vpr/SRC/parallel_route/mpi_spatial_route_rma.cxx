#include "pch.h"

#include <mpi.h>

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
#include "rr_graph_partitioner.h"

void init_datatypes();
int delta_log_output(zlog_msg_t *msg);
int ss_log_output(zlog_msg_t *msg);
int missing_edge_log_output(zlog_msg_t *msg);

void sync(congestion_t *congestion, const RRGraph &g, float pres_fac, int this_pid, int num_procs, MPI_Comm comm);
void sync_recalc_occ(congestion_t *congestion, int num_vertices, int procid, int num_procs, MPI_Comm comm);
void sync_nets(vector<net_t> &nets, vector<net_t> &global_nets, int procid, MPI_Comm comm);
void sync_net_delay(const vector<pair<box, net_t *>> &nets_to_route, int procid, int num_procs, int initial_num_procs, int *recvcounts, int *displs, int current_level, MPI_Comm comm, t_net_timing *net_timing);
void free_circuit();
void init_displ(int num_procs, int current_level, const vector<pair<box, net_t *>> &nets_to_route, int initial_num_procs, int **recvcounts, int **displs);
void get_sinks_to_route(net_t *net, const route_tree_t &rt, const vector<sink_t *> &unroutable_sinks, vector<sink_t *> &sinks_to_route);
void send_route_tree(net_t *net, const RRGraph &g, const vector<vector<sink_t *>> &routed_sinks, const vector<route_tree_t> &route_trees, int to_procid, MPI_Comm comm);
void recv_route_tree(net_t *net, const RRGraph &g, vector<vector<sink_t *>> &routed_sinks, route_state_t *states, congestion_t *congestion, float pres_fac, vector<route_tree_t> &route_trees, t_net_timing *net_timing, int from_procid, MPI_Comm comm);
void init_route_structs(const RRGraph &g, const vector<net_t> &nets, const vector<net_t> &global_nets, route_state_t **states, congestion_t **congestion, vector<route_tree_t> &route_trees, t_net_timing **net_timing);
void init_route_structs_mpi(const RRGraph &g, const vector<net_t> &nets, const vector<net_t> &global_nets, route_state_t **states, congestion_t **congestion, MPI_Win *win, vector<route_tree_t> &route_trees, t_net_timing **net_timing);
void get_sinks_to_route(net_t *net, const route_tree_t &rt, vector<sink_t *> &sinks);

extern vector<vector<FILE *>> delta_log_files;
extern vector<vector<FILE *>> missing_edge_log_files;

bool mpi_spatial_route_rma(t_router_opts *opts, struct s_det_routing_arch det_routing_arch, t_direct_inf *directs, int num_directs, t_segment_inf *segment_inf, t_timing_inf timing_inf)
{
    using std::chrono::duration_cast;
    using std::chrono::nanoseconds;
    using clock = std::chrono::high_resolution_clock;

	int num_procs, procid;
    MPI_Comm_size(MPI_COMM_WORLD, &num_procs);
    MPI_Comm_rank(MPI_COMM_WORLD, &procid);
	init_congestion_mpi_datatype();
	init_datatypes();

    init_logging();
    zlog_set_record("custom_output", delta_log_output);
    zlog_set_record("missing_edge", missing_edge_log_output);
    zlog_set_record("ss", ss_log_output);
    delta_log_files.resize(opts->max_router_iterations);
    for (int i = 0; i < opts->max_router_iterations; ++i) {
        delta_log_files[i].resize(num_procs, nullptr);
    }
    missing_edge_log_files.resize(opts->max_router_iterations);
    for (int i = 0; i < opts->max_router_iterations; ++i) {
        missing_edge_log_files[i].resize(num_procs, nullptr);
    }

    //fclose(fopen(LOG_PATH_PREFIX"iter__tid_.log", "w"));

    //test_fast_graph();
    //test_topo();
    //test_fm();
    //test_filter_graph();
    //test_partition_graph();
    //test_connected_components();

    vector<RRGraph *> graphs;
    rr_graph_partitioner partitioner;
    partitioner.init(opts, &det_routing_arch, directs, num_directs, segment_inf, &timing_inf);
    partitioner.partition(num_procs, graphs);
    //for (const auto &g : graphs) {
        //routability(*g);
    //}
    //
    //RRGraph combined_g;
    //add_vertex(combined_g, num_vertices(partitioner.orig_g));

    //for (const auto &g : graphs) {
        //for (const auto &e : get_edges(*g)) {
            //int from = get_source(*g, e);
            //int to = get_target(*g, e);
            //const auto &from_ver = get_vertex_props(*g, from);
            //const auto &to_ver = get_vertex_props(*g, to);

            //if (is_channel(from_ver) && is_channel(to_ver)) {
                //assert(!has_edge(combined_g, from, to));
                //add_edge(combined_g, from, to);
            //}
        //}
    //}

    //for (const auto &e : get_edges(partitioner.orig_g)) {
        //int from = get_source(partitioner.orig_g, e);
        //int to = get_target(partitioner.orig_g, e);
        //const auto &from_ver = get_vertex_props(partitioner.orig_g, from);
        //const auto &to_ver = get_vertex_props(partitioner.orig_g, to);
        //if (!is_channel(from_ver) || !is_channel(to_ver)) {
            //assert(!has_edge(combined_g, from, to));
            //add_edge(combined_g, from, to);
        //}
    //}

    //printf("Combined/Orig graph has %d/%d (%g) edges.\n", num_edges(combined_g), num_edges(partitioner.orig_g), 100.0*num_edges(combined_g)/num_edges(partitioner.orig_g));

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
    init_nets(nets, global_nets, opts->bb_factor, opts->large_bb);	

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

	route_state_t *states;
	congestion_t *congestion;
	vector<route_tree_t> route_trees;
	t_net_timing *net_timing;
	MPI_Win win;
	init_route_structs_mpi(partitioner.orig_g, nets, global_nets, &states, &congestion, &win, route_trees, &net_timing);

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

    //vector<clock::time_point> greedy_end_time(num_procs);
    clock::duration total_greedy_route_time = clock::duration::zero();
    clock::duration total_update_cost_time = clock::duration::zero();
    clock::duration total_analyze_timing_time = clock::duration::zero();
    clock::duration total_iter_time = clock::duration::zero();

    int iter;
    float crit_path_delay;
	int current_level = 0;

    for (iter = 0; iter < opts->max_router_iterations && !routed; ++iter) {
        clock::duration greedy_route_time = clock::duration::zero();
        clock::duration update_cost_time = clock::duration::zero();
        clock::duration partitioning_time = clock::duration::zero();
        clock::duration analyze_timing_time = clock::duration::zero();
        clock::duration iter_time = clock::duration::zero();

        //for (int i = 0; i < opts->num_threads; ++i) {
            //sprintf(buffer, LOG_PATH_PREFIX"iter_%d_tid_%d.log", iter, i);
            //fclose(fopen(buffer, "w"));
        //}

		MPI_Barrier(MPI_COMM_WORLD);
        
        auto iter_start = clock::now();

        //auto route_start = clock::now();

        //tbb::enumerable_thread_specific<state_t *> state_tls;
        //
		vector<perf_t> perfs(num_procs);
        vector<int> thread_num_nets_routed(num_procs, 0);
        vector<int> thread_num_nets_to_route(num_procs, 0);
        vector<int> thread_num_sinks_routed(num_procs, 0);
        vector<int> thread_num_sinks_to_route(num_procs, 0);
		vector<int> thread_num_interpartition_sinks(num_procs, 0);
        bool has_interpartition_sinks = false;
        vector<vector<interpartition_sink_t>> interpartition_sinks(nets.size()); 

		sprintf(buffer, "%d", iter);
		zlog_put_mdc("iter", buffer);

		sprintf(buffer, "%d", procid);
		zlog_put_mdc("tid", buffer);

		zlog_error(delta_log, "We're fucking here\n");

        zlog_info(delta_log, "Routing iteration: %d\n", iter);
        printf("Routing iteration: %d\n", iter);

        auto greedy_route_start = clock::now();

		perfs[procid].num_heap_pushes = 0;
		perfs[procid].num_heap_pops = 0;
		perfs[procid].num_neighbor_visits = 0;

		int i = procid;
		while (i < nets_to_route.size()) {
			net_t *net = nets_to_route[i].second;

            update_sink_criticalities(*net, net_timing[net->vpr_id], params);

			//printf("Routing net %d\n", net->vpr_id);

			//auto rip_up_start = clock::now();
			//if (greedy_rip_up_all) {
			route_tree_mark_all_nodes_to_be_ripped(route_trees[net->local_id], partitioner.orig_g);
			//} else {
			//route_tree_mark_congested_nodes_to_be_ripped(route_trees[net->local_id], *graphs[procid], congestion);
			//}
			route_tree_rip_up_marked_mpi_rma(route_trees[net->local_id], partitioner.orig_g, partitioner.result_pid_by_level[current_level], procid, congestion, win, params.pres_fac);

			//local_perf.total_rip_up_time += clock::now()-rip_up_start;

			//auto route_start = clock::now();

			vector<sink_t *> sinks;	
			get_sinks_to_route(net, route_trees[net->local_id], sinks);

			if (!sinks.empty()) {
				route_net_mpi_rma(partitioner.orig_g, partitioner.result_pid_by_level[current_level], procid, win, net->vpr_id, &net->source, sinks, params, states, congestion, route_trees[net->local_id], net_timing[net->vpr_id], interpartition_sinks[net->local_id], &perfs[procid]);

				if (!has_interpartition_sinks) {
					has_interpartition_sinks = !interpartition_sinks[net->local_id].empty();
				}

				thread_num_interpartition_sinks[procid] += interpartition_sinks[net->local_id].size();

				++thread_num_nets_routed[procid];
				++thread_num_nets_to_route[procid];

				thread_num_sinks_to_route[procid] += sinks.size();
				thread_num_sinks_routed[procid] += sinks.size();
			}

			//local_perf.total_route_time += clock::now()-rip_up_start;
			i += num_procs;
		}

		//greedy_end_time = clock::now();
		MPI_Barrier(MPI_COMM_WORLD);

        greedy_route_time = clock::now()-greedy_route_start;

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

        //if (greedy_rip_up_all) {
            //next_greedy_rip_up_iter += greedy_rip_up_all_period;
            //++greedy_rip_up_all_period;
            //prev_num_overused_nodes = std::numeric_limits<unsigned long>::max();
        //}

        //route_time = clock::now()-route_start;

        iter_time = clock::now()-iter_start;

		int total_num_sinks_to_route;
		if (procid == 0) {
			int total_num_nets_to_route = thread_num_nets_to_route[0];
			int total_num_nets_routed = thread_num_nets_routed[0];
			int total_num_sinks_routed = thread_num_sinks_routed[0];
			total_num_sinks_to_route = thread_num_sinks_to_route[0];
			int total_num_interpartition_sinks = thread_num_interpartition_sinks[0];

			for (int i = 1; i < num_procs; ++i) {
				MPI_Recv(&thread_num_nets_to_route[i], 1, MPI_INT, i, i, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
				MPI_Recv(&thread_num_nets_routed[i], 1, MPI_INT, i, i, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
				MPI_Recv(&thread_num_sinks_to_route[i], 1, MPI_INT, i, i, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
				MPI_Recv(&thread_num_sinks_routed[i], 1, MPI_INT, i, i, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
				MPI_Recv(&thread_num_interpartition_sinks[i], 1, MPI_INT, i, i, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

				total_num_nets_to_route += thread_num_nets_to_route[i];
				total_num_nets_routed += thread_num_nets_routed[i];
				total_num_sinks_routed += thread_num_sinks_routed[i];
				total_num_sinks_to_route += thread_num_sinks_to_route[i];
				total_num_interpartition_sinks += thread_num_interpartition_sinks[i];
			}
			assert(total_num_nets_to_route == total_num_nets_routed);
			assert(total_num_nets_to_route == nets.size());
			assert(total_num_sinks_to_route == total_num_sinks_routed);

			printf("num nets routed: ");
			for (int i = 0; i < num_procs; ++i) {
				printf("%d/%d (%g), ", thread_num_nets_routed[i], thread_num_nets_to_route[i], thread_num_nets_routed[i]*100.0/thread_num_nets_to_route[i]);
			}
			printf("\n");

			printf("num sinks routed: ");
			for (int i = 0; i < num_procs; ++i) {
				printf("%d/%d (%g), ", thread_num_sinks_routed[i], thread_num_sinks_to_route[i], thread_num_sinks_routed[i]*100.0/thread_num_sinks_to_route[i]);
			}
			printf("\n");

			printf("Total num nets routed: %d/%d (%g)\n", total_num_nets_routed, total_num_nets_to_route, total_num_nets_routed*100.0/total_num_nets_to_route);
			printf("Total num sinks routed: %d/%d (%g)\n", total_num_sinks_routed, total_num_sinks_to_route, (total_num_sinks_routed)*100.0/total_num_sinks_to_route);
			printf("Total num interpartition sinks routed: %d/%d (%g)\n", total_num_interpartition_sinks, total_num_sinks_to_route, total_num_interpartition_sinks*100.0/total_num_sinks_to_route);
		} else {
			MPI_Send(&thread_num_nets_to_route[procid], 1, MPI_INT, 0, procid, MPI_COMM_WORLD);
			MPI_Send(&thread_num_nets_routed[procid], 1, MPI_INT, 0, procid, MPI_COMM_WORLD);
			MPI_Send(&thread_num_sinks_to_route[procid], 1, MPI_INT, 0, procid, MPI_COMM_WORLD);
			MPI_Send(&thread_num_sinks_routed[procid], 1, MPI_INT, 0, procid, MPI_COMM_WORLD);
			MPI_Send(&thread_num_interpartition_sinks[procid], 1, MPI_INT, 0, procid, MPI_COMM_WORLD);
		}

		sprintf(buffer, "/Volumes/DATA/occ_before_%d.txt", procid);
		FILE *occ = fopen(buffer, "w");
		for (const auto &v : get_vertices(partitioner.orig_g)) {
			fprintf(occ, "%d\n", congestion[v].occ);
		}
		fclose(occ);

		//MPI_Win_fence(0, win);
		for (const auto &rr_node : get_vertices(partitioner.orig_g)) {
			int rr_node_pid = partitioner.result_pid_by_level[current_level][rr_node];
			int from_pid = rr_node_pid == -1 ? 0 : rr_node_pid;

			if (from_pid != procid) {
				congestion[rr_node].occ = std::numeric_limits<int>::min();
				assert(MPI_Win_lock(MPI_LOCK_SHARED, from_pid, 0, win) == MPI_SUCCESS);
				assert(MPI_Get(&congestion[rr_node].occ, 1, get_occ_dt(),
							from_pid,
							get_occ_disp(rr_node), 1, get_occ_dt(),
							win) == MPI_SUCCESS);
				assert(MPI_Win_unlock(from_pid, win) == MPI_SUCCESS);
				assert(congestion[rr_node].occ != std::numeric_limits<int>::min());
			}
		}
		//MPI_Win_fence(0, win);

		sprintf(buffer, "/Volumes/DATA/occ_after_%d.txt", procid);
		occ = fopen(buffer, "w");
		for (const auto &v : get_vertices(partitioner.orig_g)) {
			fprintf(occ, "%d\n", congestion[v].occ);
		}
		fclose(occ);

		/* checking */
        for (int i = 0; i < num_vertices(partitioner.orig_g); ++i) {
            congestion[i].recalc_occ = 0; 
        }

        for (int i = 0; i < nets_to_route.size(); ++i) {
			if (i % num_procs == procid) {
				net_t *net = nets_to_route[i].second;

				check_route_tree(route_trees[net->local_id], *net, partitioner.orig_g);
				recalculate_occ(route_trees[net->local_id], partitioner.orig_g, congestion);
			}
        }

		int *recalc_occ = new int[num_vertices(partitioner.orig_g)];
		for (int i = 0; i < num_vertices(partitioner.orig_g); ++i) {
			recalc_occ[i] = congestion[i].recalc_occ;	
		}
		if (procid == 0) {
			int *recv_recalc_occ = new int[num_vertices(partitioner.orig_g)];
			MPI_Reduce(recalc_occ, recv_recalc_occ, num_vertices(partitioner.orig_g), MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);
			for (int i = 0; i < num_vertices(partitioner.orig_g); ++i) {
				congestion[i].recalc_occ = recv_recalc_occ[i];
			}
			delete [] recv_recalc_occ;
		} else {
			MPI_Reduce(recalc_occ, nullptr, num_vertices(partitioner.orig_g), MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);
		}
		delete [] recalc_occ;
		
		//MPI_Op recalc_sum_op;
		//MPI_Op_create((MPI_User_function *)recalc_sum, 1, &recalc_sum_op);	
		//if (procid == 0) {
			//MPI_Reduce(MPI_IN_PLACE, congestion, num_vertices(partitioner.orig_g), recalc_occ_dt, recalc_sum_op, 0, MPI_COMM_WORLD);
		//} else {
			//MPI_Reduce(congestion, nullptr, num_vertices(partitioner.orig_g), recalc_occ_dt, recalc_sum_op, 0, MPI_COMM_WORLD);
		//}
		//MPI_Op_free(&recalc_sum_op);

		//if (procid == 0) {
			//int num_recvs = 0;
			//for (int i = 0; i < num_vertices(*graphs[0]); ++i) {
				////int rr_node_pid = partitioner.result_pid_by_level[current_level][i];
				////int from_pid = rr_node_pid == -1 ? 0 : rr_node_pid;

				////if (from_pid != procid) {
				//for (int from_pid = 1; from_pid < num_procs; ++from_pid) {
					//int recalc_occ;
					//MPI_Recv(&recalc_occ, 1, MPI_INT, from_pid, i, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
					//zlog_level(delta_log, ROUTER_V3, "%d Recvd %d recalc_occ %d\n", from_pid, i, recalc_occ);
					//congestion[i].recalc_occ += recalc_occ;

					//++num_recvs;
				//}
			//}
			//printf("Num recvs: %d\n", num_recvs);
		//} else {
			//int num_sends = 0;

			//for (int i = 0; i < num_vertices(*graphs[0]); ++i) {
				////int pid = partitioner.result_pid_by_level[current_level][i];
				////if (pid == procid) {
					//zlog_level(delta_log, ROUTER_V3, "Sending %d recalc_occ %d\n", i, congestion[i].recalc_occ);
					//MPI_Send(&congestion[i].recalc_occ, 1, MPI_INT, 0, i, MPI_COMM_WORLD);

					//++num_sends;
				////}
			//}

			//printf("[%d] Num sends: %d\n", procid, num_sends);
		//}

		if (procid == 0) {
			bool valid = true;
//#undef sprintf_rr_node
			for (int i = 0; i < num_vertices(partitioner.orig_g); ++i) {
				sprintf_rr_node(i, buffer);
				if (congestion[i].recalc_occ != congestion[i].occ) {
					zlog_error(delta_log, "Node %s occ mismatch, recalc: %d original: %d\n", buffer, congestion[i].recalc_occ, congestion[i].occ);
					valid = false;
				}
			}
//#define sprintf_rr_node(...)
			assert(valid);
		}

        int overused_total_bb_rank = 0;
        int num_congested_nets = 0;

        for (int i = 0; i < nets_to_route.size(); ++i) {
			if (i % num_procs == procid) {
				net_t *net = nets_to_route[i].second;

				vector<int> overused_rr_node;
				assert(route_trees[net->local_id].root_rt_nodes.size() == 1);
				get_overused_nodes(route_trees[net->local_id], route_trees[net->local_id].root_rt_nodes[0], partitioner.orig_g, congestion, overused_rr_node);
				if (!overused_rr_node.empty()) {
					zlog_level(delta_log, ROUTER_V1, "Net %d bb_rank %d has %lu overused nodes:\n", net->vpr_id, net->bb_area_rank, overused_rr_node.size());
					for (const auto &item : overused_rr_node) {
						zlog_level(delta_log, ROUTER_V1, "%d ", item);
					}
					zlog_level(delta_log, ROUTER_V1, "\n");
					overused_total_bb_rank += net->bb_area_rank;
					++num_congested_nets;
				}
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
			vector<int> overused_nodes_by_type(NUM_RR_TYPES, 0);
            for (int i = 0; i < num_vertices(partitioner.orig_g); ++i) {
                if (congestion[i].occ > get_vertex_props(partitioner.orig_g, i).capacity) {
					const auto &v_p = get_vertex_props(partitioner.orig_g, i);
					++overused_nodes_by_type[v_p.type];

                    ++num_overused_nodes;
                }
            }
            zlog_info(delta_log, "Num overused nodes: %lu/%d (%.2f)\n", num_overused_nodes, num_vertices(partitioner.orig_g), num_overused_nodes*100.0/num_vertices(partitioner.orig_g));

			if (procid == 0) {
				printf("Num overused nodes: %lu/%d (%.2f)\n", num_overused_nodes, num_vertices(partitioner.orig_g), num_overused_nodes*100.0/num_vertices(partitioner.orig_g));
				static const char *name_type[] = { "SOURCE", "SINK", "IPIN", "OPIN",
					"CHANX", "CHANY", "INTRA_CLUSTER_EDGE" };
				for (int i = 0; i < overused_nodes_by_type.size(); ++i) {
					printf("\t%s: %d (%g)\n", name_type[i], overused_nodes_by_type[i], overused_nodes_by_type[i]*100.0/num_overused_nodes);
				}
			}

			//if ((float)total_num_interpartition_sinks/total_num_sinks_to_route > 0.1) {
				//++current_level;
				//assert(current_level < partitioner.result_pid_by_level.size());
			//}

            auto update_cost_start = clock::now();

            if (iter == 0) {
                params.pres_fac = opts->initial_pres_fac;
                //update_costs_mpi(partitioner.orig_g[0], partitioner.result_pid_by_level[current_level], procid, congestion, win, params.pres_fac, 0);
				update_costs(partitioner.orig_g, congestion, params.pres_fac, 0);
            } else {
                params.pres_fac *= opts->pres_fac_mult;

				/* Avoid overflow for high iteration counts, even if acc_cost is big */
                params.pres_fac = std::min(params.pres_fac, static_cast<float>(HUGE_POSITIVE_FLOAT / 1e5));

                //update_costs_mpi(partitioner.orig_g[0], partitioner.result_pid_by_level[current_level], procid, congestion, win, params.pres_fac, opts->acc_fac);
				update_costs(partitioner.orig_g, congestion, params.pres_fac, opts->acc_fac);
            }

            update_cost_time = clock::now()-update_cost_start;
        }

        auto analyze_timing_start = clock::now();

		//MPI_Allgather(, , ,
				//, , ,
				//MPI_COMM_WORLD);
		int num_delays = 0;
		for (int i = procid; i < nets_to_route.size(); i += num_procs) {
			net_t *net = nets_to_route[i].second;
			num_delays += net->sinks.size();
		}

		int *recvcounts = new int[num_procs];
		int *displs = new int[num_procs];

		for (int i = 0; i < num_procs; ++i) {
			recvcounts[i] = 0;
		}
		for (int i = 0; i < nets_to_route.size(); ++i) {
			net_t *net = nets_to_route[i].second;
			int to_pid = i % num_procs;
			recvcounts[to_pid] += net->sinks.size();
		}
		for (int i = 0; i < num_procs; ++i) {
			displs[i] = 0;
			for (int j = 0; j < i; ++j) {
				displs[i] += recvcounts[j];
			}
		}

		float *delays = new float[num_delays];
		int idx = 0;
		for (int i = procid; i < nets_to_route.size(); i += num_procs) {
			net_t *net = nets_to_route[i].second;
			zlog_level(delta_log, ROUTER_V3, "Net %d send delays:\n", net->vpr_id);
			for (int j = 1; j <= net->sinks.size(); ++j) {
				delays[idx] = net_timing[net->vpr_id].delay[j];
				zlog_level(delta_log, ROUTER_V3, "\t%g\n", delays[idx]);
				++idx;
			}
		}

		assert(idx == num_delays);
		assert(idx == recvcounts[procid]);
		int temp_sum = 0;
		for (int i = 0; i < num_procs; ++i) {
			temp_sum += recvcounts[i];
		}

		float *all_delays;

		if (procid == 0) {
			assert(temp_sum == total_num_sinks_to_route);

			all_delays = new float[total_num_sinks_to_route];
		} else {
			all_delays = nullptr;
		}
		
		MPI_Gatherv(delays, recvcounts[procid], MPI_FLOAT, all_delays, recvcounts, displs, MPI_FLOAT, 0, MPI_COMM_WORLD);

		float *all_crits = nullptr;

		if (procid == 0) {
			for (int i = 0; i < num_procs; ++i) {
				idx = 0;
				for (int j = i; j < nets_to_route.size(); j += num_procs) {
					net_t *net = nets_to_route[j].second;
					zlog_level(delta_log, ROUTER_V3, "Net %d recv delays:\n", net->vpr_id);
					for (int k = 1; k <= net->sinks.size(); ++k) {
						net_timing[net->vpr_id].delay[k] = all_delays[displs[i] + idx];
						zlog_level(delta_log, ROUTER_V3, "\t%g\n", net_timing[net->vpr_id].delay[k]);
						++idx;
					}
				}
				assert(idx == recvcounts[i]);
			}

			crit_path_delay = analyze_timing(net_timing);

			all_crits = all_delays;

			idx = 0;
			for (int i = 0; i < num_procs; ++i) {
				idx = 0;
				for (int j = i; j < nets_to_route.size(); j += num_procs) {
					net_t *net = nets_to_route[j].second;
					zlog_level(delta_log, ROUTER_V3, "Net %d send crits:\n", net->vpr_id);
					for (int k = 1; k <= net->sinks.size(); ++k) {
						all_crits[displs[i] + idx] = net_timing[net->vpr_id].timing_criticality[k];
						zlog_level(delta_log, ROUTER_V3, "\t%g\n", all_crits[displs[i] + idx]);
						++idx;
					}
				}
				assert(idx == recvcounts[i]);
			}
		} else {
			all_crits = nullptr;
		}

		float *crits = delays;

		MPI_Scatterv(all_crits, recvcounts, displs, MPI_FLOAT, crits, recvcounts[procid], MPI_FLOAT, 0, MPI_COMM_WORLD);

		idx = 0;
		for (int i = procid; i < nets_to_route.size(); i += num_procs) {
			net_t *net = nets_to_route[i].second;
			zlog_level(delta_log, ROUTER_V3, "Net %d recv crits:\n", net->vpr_id);
			for (int j = 1; j <= net->sinks.size(); ++j) {
				net_timing[net->vpr_id].timing_criticality[j] = crits[idx];
				zlog_level(delta_log, ROUTER_V3, "\t%g\n", net_timing[net->vpr_id].timing_criticality[j]);
				++idx;
			}
		}

		delete [] delays;
		delete [] all_delays;

        analyze_timing_time = clock::now()-analyze_timing_start;

        iter_time += clock::now()-iter_start;

        //total_route_time += route_time;
        total_greedy_route_time += greedy_route_time;
        total_update_cost_time += update_cost_time;
        total_analyze_timing_time += analyze_timing_time;
        total_iter_time += iter_time;

		if (procid == 0) {
			printf("Iteration time: %g s.\n", duration_cast<nanoseconds>(iter_time).count() / 1e9);
			printf("\tRoute time: %g s.\n", duration_cast<nanoseconds>(greedy_route_time).count() / 1e9);
			printf("\tUpdate cost time: %g s.\n", duration_cast<nanoseconds>(update_cost_time).count() / 1e9);
			printf("\tAnalyze timing time: %g s.\n", duration_cast<nanoseconds>(analyze_timing_time).count() / 1e9);
			printf("Critical path: %g ns\n", crit_path_delay);

			printf("num_heap_pushes: ");
			for (int i = 0; i < num_procs; ++i) {
				printf("%lu ", perfs[i].num_heap_pushes);
			}
			printf("\n");

			printf("num_heap_pops: ");
			for (int i = 0; i < num_procs; ++i) {
				printf("%lu ", perfs[i].num_heap_pops);
			}
			printf("\n");

			printf("num_neighbor_visits: ");
			for (int i = 0; i < num_procs; ++i) {
				printf("%lu ", perfs[i].num_neighbor_visits);
			}
			printf("\n");

			unsigned long total_num_heap_pushes = 0;
			unsigned long total_num_heap_pops = 0;
			unsigned long total_num_neighbor_visits = 0;

			for (int i = 0; i < num_procs; ++i) {
				total_num_heap_pushes += perfs[i].num_heap_pushes;
				total_num_heap_pops += perfs[i].num_heap_pops;
				total_num_neighbor_visits += perfs[i].num_neighbor_visits;
			}

			printf("total_num_heap_pushes: %lu\n", total_num_heap_pushes);
			printf("total_num_heap_pops: %lu\n", total_num_heap_pops); 
			printf("total_num_neighbor_visits: %lu\n", total_num_neighbor_visits);
		}

        //clock::time_point greedy_earliest_end_time = *std::min_element(begin(greedy_end_time), begin(greedy_end_time)+num_procs);

        //printf("greedy wait time: ");
        //for (int i = 0; i < num_procs; ++i) {
			//printf("%g (%g), ", duration_cast<nanoseconds>(greedy_end_time[i]-greedy_earliest_end_time).count() / 1e9, 100.0*(greedy_end_time[i]-greedy_earliest_end_time)/greedy_route_time);
        //}
		//printf("\n");
    }

    //sprintf(buffer, "%s_run_%d_rr_graph_occ.txt", s_circuit_name, run);
    //dump_rr_graph_occ(congestion, num_vertices(g), buffer);

    if (routed) {
        printf("Routed in %d iterations. Total iteration time: %g\n", iter, duration_cast<nanoseconds>(total_iter_time).count() / 1e9);
            printf("\tTotal route time: %g s.\n", duration_cast<nanoseconds>(total_greedy_route_time).count() / 1e9);
                printf("\t\tTotal greedy route time: %g s.\n", duration_cast<nanoseconds>(total_greedy_route_time).count() / 1e9);
            printf("\tTotal update cost time: %g s.\n", duration_cast<nanoseconds>(total_update_cost_time).count() / 1e9);
            printf("\tTotal analyze timing time: %g s.\n", duration_cast<nanoseconds>(total_analyze_timing_time).count() / 1e9);

        printf("Final critical path: %g ns\n", crit_path_delay);
    } else {
        printf("Failed to route in %d iterations. Total iteration time: %g\n", opts->max_router_iterations, duration_cast<nanoseconds>(total_iter_time).count() / 1e9);
            printf("\tTotal route time: %g s.\n", duration_cast<nanoseconds>(total_greedy_route_time).count() / 1e9);
                printf("\t\tTotal greedy route time: %g s.\n", duration_cast<nanoseconds>(total_greedy_route_time).count() / 1e9);
            printf("\tTotal update cost time: %g s.\n", duration_cast<nanoseconds>(total_update_cost_time).count() / 1e9);
            printf("\tTotal analyze timing time: %g s.\n", duration_cast<nanoseconds>(total_analyze_timing_time).count() / 1e9);
    }

	MPI_Finalize();
	exit(0);

    //delete_graph(g);
    delete_net_timing(nets, global_nets, net_timing);	
    delete [] congestion;
    delete [] net_timing;
    delete [] states;

    return routed;
}
