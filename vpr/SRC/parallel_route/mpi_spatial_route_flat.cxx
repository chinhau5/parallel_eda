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
#include "scheduler.h"
#include "geometry.h"
#include "quadtree.h"
#include "utility.h"
#include "cluster.h"
#include "args.h"
#include "init.h"
#include "router.h"
#include "congestion.h"
#include "partition.h"
#include "fm.h"
#include "rr_graph_partitioner.h"

void init_datatypes();
int delta_log_output(zlog_msg_t *msg);
int ss_log_output(zlog_msg_t *msg);
int missing_edge_log_output(zlog_msg_t *msg);

void sync(congestion_t *congestion, const RRGraph &g, float pres_fac, int this_pid, int num_procs, MPI_Comm comm, mpi_perf_t *mpi_perf);
void sync_recalc_occ(congestion_t *congestion, int num_vertices, int procid, int num_procs, MPI_Comm comm);
void sync_nets(vector<net_t> &nets, vector<net_t> &global_nets, int procid, MPI_Comm comm);
void sync_net_delay(const vector<pair<box, net_t *>> &nets_to_route, int procid, int num_procs, int initial_num_procs, int *recvcounts, int *displs, int current_level, MPI_Comm comm, t_net_timing *net_timing);
void free_circuit();
void init_displ(int num_procs, int current_level, const vector<pair<box, net_t *>> &nets_to_route, int initial_num_procs, int **recvcounts, int **displs);
void get_sinks_to_route(net_t *net, const route_tree_t &rt, const vector<sink_t *> &unroutable_sinks, vector<sink_t *> &sinks_to_route);
void send_route_tree(net_t *net, const RRGraph &g, const vector<vector<sink_t *>> &routed_sinks, const vector<route_tree_t> &route_trees, int to_procid, MPI_Comm comm);
void recv_route_tree(net_t *net, const RRGraph &g, vector<vector<sink_t *>> &routed_sinks, route_state_t *states, vector<route_tree_t> &route_trees, t_net_timing *net_timing, int from_procid, MPI_Comm comm);
void init_route_structs(const RRGraph &g, const vector<net_t> &nets, const vector<net_t> &global_nets, route_state_t **states, congestion_t **congestion, vector<route_tree_t> &route_trees, t_net_timing **net_timing);

extern vector<vector<FILE *>> delta_log_files;
extern vector<vector<FILE *>> missing_edge_log_files;

bool mpi_spatial_route_flat(t_router_opts *opts, struct s_det_routing_arch det_routing_arch, t_direct_inf *directs, int num_directs, t_segment_inf *segment_inf, t_timing_inf timing_inf)
{
    using std::chrono::duration_cast;
    using std::chrono::nanoseconds;
    using clock = std::chrono::high_resolution_clock;

	int initial_num_procs, prev_num_procs, num_procs, initial_procid, prev_procid, procid;
	MPI_Comm cur_comm = MPI_COMM_WORLD;
	MPI_Comm prev_comm = cur_comm;

    MPI_Comm_size(cur_comm, &initial_num_procs);
    MPI_Comm_rank(cur_comm, &initial_procid);

	printf("[%d] Initializing router\n", initial_procid);

	num_procs = initial_num_procs;
	prev_num_procs = num_procs;

	procid = initial_procid;
	prev_procid = procid;

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
    partitioner.partition_without_ipin(1, graphs);

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
	//
	zlog_put_mdc("iter", "0");

	char buffer[256];
	sprintf(buffer, "%d", procid);
	zlog_put_mdc("tid", buffer);

    vector<net_t> nets;
    vector<net_t> global_nets;
    //init_nets(nets, global_nets, opts->bb_factor, partitioner.sink_in_nodes);	
	if (procid == 0) {
		printf("[%d] initializing nets\n", procid);
		init_nets(nets, global_nets, opts->bb_factor, opts->large_bb);	

		//printf("Num_interpatition_nets [%d nets]: ", nets.size());
		//for (int i = 2; i <= 16; i *= 2) {
			//printf("%d ", get_num_interpartition_nets(nets, i));
		//}
		//printf("\n");
		//exit(0);
	}
	printf("[%d] syncing nets\n", procid);
	sync_nets(nets, global_nets, procid, cur_comm);

	printf("[%d] done initializing nets\n", procid);

    vector<net_t *> nets_ptr(nets.size());
    for (int i = 0; i < nets.size(); ++i) {
        nets_ptr[i] = &nets[i];
    }
    extern s_bb *route_bb;
    sort(nets_ptr.begin(), nets_ptr.end(), [] (const net_t *a, const net_t *b) -> bool {
            return get_bounding_box_area(a->bounding_box) > get_bounding_box_area(b->bounding_box);
            });
    int rank = 0;
    for (auto &net : nets_ptr) {
        net->bb_area_rank = rank++;
    }

	route_state_t *states;
	congestion_t *congestion;
	vector<route_tree_t> route_trees;
	t_net_timing *net_timing;
	init_route_structs(partitioner.orig_g, nets, global_nets, &states, &congestion, route_trees, &net_timing);

    route_parameters_t params;
    params.criticality_exp = opts->criticality_exp;
    params.astar_fac = opts->astar_fac;
    params.max_criticality = opts->max_criticality;
    params.bend_cost = opts->bend_cost;
	params.pres_fac = opts->first_iter_pres_fac; /* Typically 0 -> ignore cong. */

	if (procid != 0) {
		free_circuit();
	}
	free_rr_graph();

    vector<pair<box, net_t *>> nets_to_route;
    
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
            return a.second->sinks.size()*get_bounding_box_area(a.second->bounding_box) > b.second->sinks.size()*get_bounding_box_area(b.second->bounding_box);
            });

    bool routed = false;
	bool idling = false;

    int iter;
    float crit_path_delay;
	int current_level = 0;

	vector<vector<sink_t *>> unroutable_sinks(nets.size());
	//vector<vector<sink_t *>> fixed_sinks(nets.size());
	vector<vector<sink_t *>> routed_sinks(nets.size());
	bool has_unroutable_sinks = false;
	unsigned long prev_num_overused_nodes = std::numeric_limits<unsigned long>::max();

    clock::duration total_greedy_route_time = clock::duration::zero();
    clock::duration total_update_cost_time = clock::duration::zero();
    clock::duration total_analyze_timing_time = clock::duration::zero();
    clock::duration total_iter_time = clock::duration::zero();
    clock::duration total_wait_time = clock::duration::zero();
    clock::duration total_combine_time = clock::duration::zero();

    for (iter = 0; iter < opts->max_router_iterations && !routed && !idling; ++iter) {
        clock::duration greedy_route_time = clock::duration::zero();
        clock::duration update_cost_time = clock::duration::zero();
        clock::duration analyze_timing_time = clock::duration::zero();
        clock::duration iter_time = clock::duration::zero();
        clock::duration wait_time = clock::duration::zero();
        clock::duration combine_time = clock::duration::zero();

        //for (int i = 0; i < opts->num_threads; ++i) {
            //sprintf(buffer, LOG_PATH_PREFIX"iter_%d_tid_%d.log", iter, i);
            //fclose(fopen(buffer, "w"));
        //}

        //auto route_start = clock::now();

        //tbb::enumerable_thread_specific<state_t *> state_tls;
        //
		perf_t perf;
		mpi_perf_t mpi_perf;
        int thread_num_nets_routed = 0;
        int thread_num_nets_to_route = 0;
        int thread_num_sinks_routed = 0;
        int thread_num_sinks_to_route = 0;
		//vector<vector<RRNode>> net_sinks(nets.size());

		sprintf(buffer, "%d", iter);
		zlog_put_mdc("iter", buffer);

		sprintf(buffer, "%d", procid);
		zlog_put_mdc("tid", buffer);

        zlog_info(delta_log, "Routing iteration: %d\n", iter);
        printf("Routing iteration: %d\n", iter);

		vector<ongoing_transaction_t> transactions;

		MPI_Barrier(cur_comm);
        
        auto iter_start = clock::now();

        auto greedy_route_start = clock::now();

		perf.num_heap_pushes = 0;
		perf.num_heap_pops = 0;
		perf.num_neighbor_visits = 0;

		mpi_perf.total_broadcast_time = clock::duration::zero();
		mpi_perf.total_sync_time = clock::duration::zero();
		mpi_perf.total_probe_time = clock::duration::zero();

		for (int i = procid*pow(2, current_level); i < nets_to_route.size(); i += initial_num_procs) {
			for (int j = 0; j < pow(2, current_level) && i+j < nets_to_route.size(); ++j) {
				zlog_level(delta_log, ROUTER_V1, "Net index: %d\n", i+j);

				net_t *net = nets_to_route[i+j].second;

				update_sink_criticalities(*net, net_timing[net->vpr_id], params);

				//printf("Routing net %d\n", net->vpr_id);
				//vector<sink_t *> temp_routed_sinks = routed_sinks[net->local_id];

				//vector<RRNode> sinks_to_mark;
				//for (const auto &sink : temp_routed_sinks) {
					//bool fixed = find(begin(fixed_sinks[net->local_id]), end(fixed_sinks[net->local_id]), sink) != end(fixed_sinks[net->local_id]);
					//if (!fixed) {
						//sinks_to_mark.push_back(sink->rr_node);
					//} else {
						//routed_sinks[net->local_id].push_back(sink);
					//}
				//}

				//route_tree_mark_paths_to_be_ripped(route_trees[net->local_id], partitioner.orig_g, partitioner.result_pid_by_level[current_level], procid, sinks_to_mark);
				if (opts->rip_up_always && current_level == 0) {
					route_tree_mark_all_nodes_to_be_ripped(route_trees[net->local_id], partitioner.orig_g);
				} else {
					route_tree_mark_congested_nodes_to_be_ripped(route_trees[net->local_id], partitioner.orig_g, congestion);
				}

				route_tree_rip_up_marked_mpi_send_recv(route_trees[net->local_id], partitioner.orig_g, congestion, params.pres_fac, procid, num_procs, cur_comm, transactions);

				auto sync_start	= clock::now();
				sync(congestion, partitioner.orig_g, params.pres_fac, procid, num_procs, cur_comm, &mpi_perf);
				mpi_perf.total_sync_time += clock::now()-sync_start;

				routed_sinks[net->local_id].clear();
				for (auto &sink : net->sinks) {
					if (route_tree_get_rt_node(route_trees[net->local_id], sink.rr_node) != RouteTree::null_vertex()) {
						routed_sinks[net->local_id].push_back(&sink);
					}
				}

				vector<sink_t *> sinks;	
				get_sinks_to_route(net, route_trees[net->local_id], unroutable_sinks[net->local_id], sinks);

				//local_perf.total_rip_up_time += clock::now()-rip_up_start;

				//auto route_start = clock::now();

				if (!sinks.empty()) {
					int previous_num_unroutable_sinks = unroutable_sinks[net->local_id].size();

					int previous_num_routed_sinks = routed_sinks[net->local_id].size();

					route_net_mpi_send_recv(partitioner.orig_g, procid, net->vpr_id, &net->source, sinks, params, states, congestion, num_procs, cur_comm, transactions, route_trees[net->local_id], net_timing[net->vpr_id], routed_sinks[net->local_id], unroutable_sinks[net->local_id], &perf, &mpi_perf);

					assert(routed_sinks[net->local_id].size() + unroutable_sinks[net->local_id].size() == net->sinks.size());

					if (!unroutable_sinks[net->local_id].empty() && previous_num_unroutable_sinks > 0) {
						assert(previous_num_unroutable_sinks == unroutable_sinks[net->local_id].size());
					}

					if (!has_unroutable_sinks) {
						has_unroutable_sinks = !unroutable_sinks[net->local_id].empty();
					}

					++thread_num_nets_routed;
					++thread_num_nets_to_route;

					thread_num_sinks_to_route += sinks.size();
					thread_num_sinks_routed += routed_sinks[net->local_id].size() - previous_num_routed_sinks;
				} else {
					assert(routed_sinks[net->local_id].size() + unroutable_sinks[net->local_id].size() == net->sinks.size());
					zlog_level(delta_log, ROUTER_V2, "Net %d has no sinks in the current iteration because there are %lu/%lu non-routable/all sinks\n", net->vpr_id, unroutable_sinks[net->local_id].size(), net->sinks.size());
				}

				//local_perf.total_route_time += clock::now()-rip_up_start;
			}
		}

		auto wait_start = clock::now();
		MPI_Barrier(cur_comm);
		wait_time = clock::now()-wait_start;

		auto sync_start	= clock::now();
		sync(congestion, partitioner.orig_g, params.pres_fac, procid, num_procs, cur_comm, &mpi_perf);
		mpi_perf.total_sync_time += clock::now()-sync_start;

		MPI_Barrier(cur_comm);

		//greedy_end_time = clock::now();

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
			int total_num_nets_to_route = thread_num_nets_to_route;
			int total_num_nets_routed = thread_num_nets_routed;
			int total_num_sinks_routed = thread_num_sinks_routed;
			int total_num_sinks_to_route = thread_num_sinks_to_route;

			printf("num nets routed: %d/%d (%g) ", thread_num_nets_routed, thread_num_nets_to_route, thread_num_nets_routed*100.0/thread_num_nets_to_route);
			for (int i = 1; i < num_procs; ++i) {
				int tmp_thread_num_nets_to_route;
				int tmp_thread_num_nets_routed;
				MPI_Recv(&tmp_thread_num_nets_to_route, 1, MPI_INT, i, i, cur_comm, MPI_STATUS_IGNORE);
				MPI_Recv(&tmp_thread_num_nets_routed, 1, MPI_INT, i, i, cur_comm, MPI_STATUS_IGNORE);
				total_num_nets_to_route += tmp_thread_num_nets_to_route;
				total_num_nets_routed += tmp_thread_num_nets_routed;

				printf("%d/%d (%g) ", tmp_thread_num_nets_routed, tmp_thread_num_nets_to_route, tmp_thread_num_nets_routed*100.0/tmp_thread_num_nets_to_route);
			}
			printf("\n");

			assert(total_num_nets_to_route == total_num_nets_routed);

			printf("num sinks routed: %d/%d (%g) ", thread_num_sinks_routed, thread_num_sinks_to_route, thread_num_sinks_routed*100.0/thread_num_sinks_to_route);
			for (int i = 1; i < num_procs; ++i) {
				int tmp_thread_num_sinks_to_route;
				int tmp_thread_num_sinks_routed;
				MPI_Recv(&tmp_thread_num_sinks_to_route, 1, MPI_INT, i, i, cur_comm, MPI_STATUS_IGNORE);
				MPI_Recv(&tmp_thread_num_sinks_routed, 1, MPI_INT, i, i, cur_comm, MPI_STATUS_IGNORE);

				total_num_sinks_to_route += tmp_thread_num_sinks_to_route;
				total_num_sinks_routed += tmp_thread_num_sinks_routed;

				printf("%d/%d (%g) ", tmp_thread_num_sinks_routed, tmp_thread_num_sinks_to_route, tmp_thread_num_sinks_routed*100.0/tmp_thread_num_sinks_to_route);
			}
			printf("\n");

			printf("Total num nets routed: %d/%d (%g)\n", total_num_nets_routed, total_num_nets_to_route, total_num_nets_routed*100.0/total_num_nets_to_route);
			printf("Total num sinks routed: %d/%d (%g)\n", total_num_sinks_routed, total_num_sinks_to_route, (total_num_sinks_routed)*100.0/total_num_sinks_to_route);

			unsigned long total_num_heap_pushes = perf.num_heap_pushes;
			unsigned long total_num_heap_pops = perf.num_heap_pops;
			unsigned long total_num_neighbor_visits = perf.num_neighbor_visits;

			printf("num_heap_pushes: %lu ", perf.num_heap_pushes);
			for (int i = 1; i < num_procs; ++i) {
				unsigned long tmp_num_heap_pushes;
				MPI_Recv(&tmp_num_heap_pushes, 1, MPI_UNSIGNED_LONG, i, i, cur_comm, MPI_STATUS_IGNORE);
				total_num_heap_pushes += tmp_num_heap_pushes;
				printf("%lu ", tmp_num_heap_pushes);
			}
			printf("\n");

			printf("num_heap_pops: %lu ", perf.num_heap_pops);
			for (int i = 1; i < num_procs; ++i) {
				unsigned long tmp_num_heap_pops;
				MPI_Recv(&tmp_num_heap_pops, 1, MPI_UNSIGNED_LONG, i, i, cur_comm, MPI_STATUS_IGNORE);
				total_num_heap_pops += tmp_num_heap_pops;
				printf("%lu ", tmp_num_heap_pops);
			}
			printf("\n");

			printf("num_neighbor_visits: %lu ", perf.num_neighbor_visits);
			for (int i = 1; i < num_procs; ++i) {
				unsigned long tmp_num_neighbor_visits;
				MPI_Recv(&tmp_num_neighbor_visits, 1, MPI_UNSIGNED_LONG, i, i, cur_comm, MPI_STATUS_IGNORE);
				total_num_neighbor_visits += tmp_num_neighbor_visits;
				printf("%lu ", tmp_num_neighbor_visits);
			}
			printf("\n");

			printf("total_num_heap_pushes: %lu\n", total_num_heap_pushes);
			printf("total_num_heap_pops: %lu\n", total_num_heap_pops); 
			printf("total_num_neighbor_visits: %lu\n", total_num_neighbor_visits);
		} else {
			MPI_Send(&thread_num_nets_to_route, 1, MPI_INT, 0, procid, cur_comm);
			MPI_Send(&thread_num_nets_routed, 1, MPI_INT, 0, procid, cur_comm);
			MPI_Send(&thread_num_sinks_to_route, 1, MPI_INT, 0, procid, cur_comm);
			MPI_Send(&thread_num_sinks_routed, 1, MPI_INT, 0, procid, cur_comm);

			MPI_Send(&perf.num_heap_pushes, 1, MPI_UNSIGNED_LONG, 0, procid, cur_comm);
			MPI_Send(&perf.num_heap_pops, 1, MPI_UNSIGNED_LONG, 0, procid, cur_comm);
			MPI_Send(&perf.num_neighbor_visits, 1, MPI_UNSIGNED_LONG, 0, procid, cur_comm);
		}

		/* checking */
        for (int i = 0; i < num_vertices(partitioner.orig_g); ++i) {
            congestion[i].recalc_occ = 0; 
        }

		for (int i = procid*pow(2, current_level); i < nets_to_route.size(); i += initial_num_procs) {
			for (int j = 0; j < pow(2, current_level) && i+j < nets_to_route.size(); ++j) {
				net_t *net = nets_to_route[i+j].second;

				check_route_tree(route_trees[net->local_id], *net, routed_sinks[net->local_id], partitioner.orig_g);
				recalculate_occ(route_trees[net->local_id], partitioner.orig_g, congestion);
			}
        }

		//int *recalc_occ = new int[num_vertices(partitioner.orig_g)];
		//for (int i = 0; i < num_vertices(partitioner.orig_g); ++i) {
			//recalc_occ[i] = congestion[i].recalc_occ;	
		//}
		//if (procid == 0) {
			//int *recv_recalc_occ = new int[num_vertices(partitioner.orig_g)];
			//MPI_Reduce(recalc_occ, recv_recalc_occ, num_vertices(partitioner.orig_g), MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);
			//for (int i = 0; i < num_vertices(partitioner.orig_g); ++i) {
				//congestion[i].recalc_occ = recv_recalc_occ[i];
			//}
			//delete [] recv_recalc_occ;
		//} else {
			//MPI_Reduce(recalc_occ, nullptr, num_vertices(partitioner.orig_g), MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);
		//}
		//delete [] recalc_occ;
		
		//MPI_Op recalc_sum_op;
		//MPI_Op_create((MPI_User_function *)recalc_sum, 1, &recalc_sum_op);	
		//if (procid == 0) {
			//MPI_Reduce(MPI_IN_PLACE, congestion, num_vertices(partitioner.orig_g), recalc_occ_dt, recalc_sum_op, 0, MPI_COMM_WORLD);
		//} else {
			//MPI_Reduce(congestion, nullptr, num_vertices(partitioner.orig_g), recalc_occ_dt, recalc_sum_op, 0, MPI_COMM_WORLD);
		//}
		//MPI_Op_free(&recalc_sum_op);

		sync_recalc_occ(congestion, num_vertices(partitioner.orig_g), procid, num_procs, cur_comm);

		if (procid == 0) {
			bool valid = true;
			for (int i = 0; i < num_vertices(partitioner.orig_g); ++i) {
				sprintf_rr_node_impl(i, buffer);
				if (congestion[i].recalc_occ != congestion[i].occ) {
					printf("Node %s occ mismatch, recalc: %d original: %d\n", buffer, congestion[i].recalc_occ, congestion[i].occ);
					valid = false;
				}
			}
			assert(valid);
		}

        int overused_total_bb_rank = 0;
        int num_congested_nets = 0;

		for (int i = procid*pow(2, current_level); i < nets_to_route.size(); i += initial_num_procs) {
			for (int j = 0; j < pow(2, current_level) && i+j < nets_to_route.size(); ++j) {
				net_t *net = nets_to_route[i+j].second;

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

        auto analyze_timing_start = clock::now();

		int *recvcounts = nullptr;
		int *displs = nullptr;

		init_displ(num_procs, current_level, nets_to_route, initial_num_procs, &recvcounts, &displs);

		int num_crits = 0;
		for (int i = 0; i < num_procs; ++i) {
			num_crits += recvcounts[i];
		}

		sync_net_delay(nets_to_route, procid, num_procs, initial_num_procs, recvcounts, displs, current_level, cur_comm, net_timing);

		float *all_crits;
		int idx;

		if (procid == 0) {
			crit_path_delay = analyze_timing(net_timing);

			all_crits = new float[num_crits];

			idx = 0;
			for (int pi = 0; pi < num_procs; ++pi) {
				idx = 0;
				for (int i = pi*pow(2, current_level); i < nets_to_route.size(); i += initial_num_procs) {
					for (int j = 0; j < pow(2, current_level) && i+j < nets_to_route.size(); ++j) {
						net_t *net = nets_to_route[i+j].second;
						zlog_level(delta_log, ROUTER_V3, "Net %d send crits:\n", net->vpr_id);
						for (int k = 1; k <= net->sinks.size(); ++k) {
							all_crits[displs[pi] + idx] = net_timing[net->vpr_id].timing_criticality[k];
							zlog_level(delta_log, ROUTER_V3, "\t%g\n", all_crits[displs[pi] + idx]);
							++idx;
						}
					}
				}
				assert(idx == recvcounts[pi]);
			}
		} else {
			all_crits = nullptr;
		}

		float *crits = new float[recvcounts[procid]];

		MPI_Scatterv(all_crits, recvcounts, displs, MPI_FLOAT, crits, recvcounts[procid], MPI_FLOAT, 0, cur_comm);

		idx = 0;
		for (int i = procid*pow(2, current_level); i < nets_to_route.size(); i += initial_num_procs) {
			for (int j = 0; j < pow(2, current_level) && i+j < nets_to_route.size(); ++j) {
				net_t *net = nets_to_route[i+j].second;
				zlog_level(delta_log, ROUTER_V3, "Net %d recv crits:\n", net->vpr_id);

				for (int s = 1; s <= net->sinks.size(); ++s) {
					net_timing[net->vpr_id].timing_criticality[s] = crits[idx];
					zlog_level(delta_log, ROUTER_V3, "\t%g\n", net_timing[net->vpr_id].timing_criticality[s]);
					++idx;
				}
			}
		}

		if (procid == 0) {
			delete [] all_crits;
		}
		delete [] crits;

        analyze_timing_time = clock::now()-analyze_timing_start;

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

        total_greedy_route_time += greedy_route_time;
		total_wait_time += wait_time;
        total_analyze_timing_time += analyze_timing_time;
        total_update_cost_time += update_cost_time;

		int m_routed = (feasible_routing(partitioner.orig_g, congestion) && !has_unroutable_sinks) ? 1 : 0;
		int reduced_routed;
		MPI_Allreduce(&m_routed, &reduced_routed, 1, MPI_INT, MPI_LAND, cur_comm);

        zlog_level(delta_log, ROUTER_V1, "m_routed: %d reduced_routed: %d\n", m_routed, reduced_routed);

        if (reduced_routed) {
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

			static const char *name_type[] = { "SOURCE", "SINK", "IPIN", "OPIN",
				"CHANX", "CHANY", "INTRA_CLUSTER_EDGE" };
            zlog_info(delta_log, "Num overused nodes: %lu/%d (%.2f)\n", num_overused_nodes, num_vertices(partitioner.orig_g), num_overused_nodes*100.0/num_vertices(partitioner.orig_g));
			for (int i = 0; i < overused_nodes_by_type.size(); ++i) {
				zlog_info(delta_log, "\t%s: %d (%g)\n", name_type[i], overused_nodes_by_type[i], overused_nodes_by_type[i]*100.0/num_overused_nodes);
			}

			int *all_overused_nodes_by_type = new int[num_procs*NUM_RR_TYPES];
			int *overused_nodes_by_type_send = new int[NUM_RR_TYPES];
			for (int i = 0; i < overused_nodes_by_type.size(); ++i) {
				overused_nodes_by_type_send[i] = overused_nodes_by_type[i];
			}

			MPI_Gather(overused_nodes_by_type_send, NUM_RR_TYPES, MPI_INT, all_overused_nodes_by_type, NUM_RR_TYPES, MPI_INT, 0, cur_comm);

			unsigned long *all_num_overused_nodes = new unsigned long[num_procs];
			MPI_Gather(&num_overused_nodes, 1, MPI_UNSIGNED_LONG, all_num_overused_nodes, 1, MPI_UNSIGNED_LONG, 0, cur_comm);

			if (procid == 0) {
				printf("Num overused nodes: ");
				for (int i = 0; i < num_procs; ++i) {
					printf("%lu/%d (%.2f) ", all_num_overused_nodes[i], num_vertices(partitioner.orig_g), all_num_overused_nodes[i]*100.0/num_vertices(partitioner.orig_g));
				}
				printf("\n");

				for (int i = 0; i < overused_nodes_by_type.size(); ++i) {
					printf("\t%s: ", name_type[i]);
					for (int j = 0; j < num_procs; ++j) {
						printf("%d (%g) ", all_overused_nodes_by_type[j*NUM_RR_TYPES+i], all_overused_nodes_by_type[j*NUM_RR_TYPES+i]*100.0/all_num_overused_nodes[j]);
					}
					printf("\n");
				}
			} 

			delete [] all_overused_nodes_by_type;
			delete [] overused_nodes_by_type_send;
			delete [] all_num_overused_nodes;

			int not_decreasing = (num_overused_nodes >= prev_num_overused_nodes && iter > 10) ? 1 : 0;
			//int not_decreasing = current_level+1 < partitioner.result_pid_by_level.size(); [> testing <]
			//int not_decreasing = current_level+1 <= std::log2(initial_num_procs); [> testing <]
			int reduced_not_decreasing;
			MPI_Allreduce(&not_decreasing, &reduced_not_decreasing, 1, MPI_INT, MPI_LOR, cur_comm);

			prev_num_overused_nodes = num_overused_nodes;

			zlog_level(delta_log, ROUTER_V1, "not_decreasing: %d reduced_not_decreasing: %d\n", not_decreasing, reduced_not_decreasing);

			if (reduced_not_decreasing && initial_num_procs > 1 && current_level < (int)std::log2(initial_num_procs)) {
				/* need to send route tree over */

				auto combine_start = clock::now();

				if (procid % 2 == 0) {
					/*receiver*/
					for (int i = (procid+1)*pow(2, current_level); i < nets_to_route.size(); i += initial_num_procs) {
						for (int j = 0; j < pow(2, current_level) && i+j < nets_to_route.size(); ++j) {
							net_t *net = nets_to_route[i+j].second;

							zlog_level(delta_log, ROUTER_V3, "Recving net index %d from %d\n", i+j, procid+1);

							recv_route_tree(net, partitioner.orig_g, routed_sinks, states, route_trees, net_timing, procid+1, cur_comm);
						}
					}
				} else {
					/*sender*/
					for (int i = procid*pow(2, current_level); i < nets_to_route.size(); i += initial_num_procs) {
						for (int j = 0; j < pow(2, current_level) && i+j < nets_to_route.size(); ++j) {
							net_t *net = nets_to_route[i+j].second;

							zlog_level(delta_log, ROUTER_V3, "Sending net index %d from %d\n", i+j, procid);

							send_route_tree(net, partitioner.orig_g, routed_sinks, route_trees, procid-1, cur_comm);
						}
					}
				}

				assert(num_procs % 2 == 0);
				num_procs /= 2;
				MPI_Comm new_comm;
				MPI_Comm_split(cur_comm, procid%2, procid, &new_comm);

				combine_time = clock::now()-combine_start;

				idling = procid % 2 == 1;

				if (!idling) {
					cur_comm = new_comm;
					MPI_Comm_rank(cur_comm, &procid);

					++current_level;

					//assert(current_level < partitioner.result_pid_by_level.size());

					printf("[%d] Transitioned to level %d at iteration %d\n", procid, current_level, iter);
					zlog_level(delta_log, ROUTER_V1, "Transitioned to level %d at iteration %d\n", current_level, iter);
					zlog_level(delta_log, ROUTER_V1, "New pid %d for initial pid %d\n", procid, initial_procid);

					for (int i = procid*pow(2, current_level); i < nets_to_route.size(); i += initial_num_procs) {
						for (int j = 0; j < pow(2, current_level) && i+j < nets_to_route.size(); ++j) {
							net_t *net = nets_to_route[i+j].second;

							unroutable_sinks[net->local_id].clear();

							//for (const auto &rs : routed_sinks[net->local_id]) {
							//assert(find(begin(fixed_sinks[net->local_id]), end(fixed_sinks[net->local_id]), rs) == end(fixed_sinks[net->local_id]));

							//zlog_level(delta_log, ROUTER_V3, "Fixing net %d sink %d\n", net->vpr_id, rs->id);

							//fixed_sinks[net->local_id].push_back(rs);
							//}
						}
					}

					prev_num_overused_nodes = std::numeric_limits<unsigned long>::max();
					has_unroutable_sinks = false;

					for (int i = 0; i < num_vertices(partitioner.orig_g); ++i) {
						congestion[i].recalc_occ = 0; 
					}

					for (int i = procid*pow(2, current_level); i < nets_to_route.size(); i += initial_num_procs) {
						for (int j = 0; j < pow(2, current_level) && i+j < nets_to_route.size(); ++j) {
							net_t *net = nets_to_route[i+j].second;

							if (!routed_sinks[net->local_id].empty()) {
								zlog_level(delta_log, ROUTER_V3, "Checking net index %d\n", i+j);

								check_route_tree(route_trees[net->local_id], *net, routed_sinks[net->local_id], partitioner.orig_g);
								recalculate_occ(route_trees[net->local_id], partitioner.orig_g, congestion);
							} else {
								zlog_level(delta_log, ROUTER_V3, "Not checking net index %d because of empty route tree\n", i+j);
							}
						}
					}

					sync_recalc_occ(congestion, num_vertices(partitioner.orig_g),  procid, num_procs, cur_comm);

					if (procid == 0) {
						bool valid = true;
						for (int i = 0; i < num_vertices(partitioner.orig_g); ++i) {
							sprintf_rr_node_impl(i, buffer);
							if (congestion[i].recalc_occ != congestion[i].occ) {
								printf("Node %s occ mismatch, recalc: %d original: %d\n", buffer, congestion[i].recalc_occ, congestion[i].occ);
								valid = false;
							}
						}
						assert(valid);
					}
				}
			} 

			//if ((float)total_num_interpartition_sinks/total_num_sinks_to_route > 0.1) {
				//++current_level;
				//assert(current_level < partitioner.result_pid_by_level.size());
			//}

        }

        iter_time += clock::now()-iter_start;

        total_iter_time += iter_time;
		total_combine_time += combine_time;

		if (prev_procid == 0) {
			printf("Iteration time: %g s.\n", duration_cast<nanoseconds>(iter_time).count() / 1e9);
			printf("\tRoute time: %g s.\n", duration_cast<nanoseconds>(greedy_route_time).count() / 1e9);
			float f_sync_time = duration_cast<nanoseconds>(mpi_perf.total_sync_time).count() / 1e9;
			printf("\t\tSync time: %g (%g) ", f_sync_time, f_sync_time * 100.0 / (duration_cast<nanoseconds>(greedy_route_time).count() / 1e9));
			for (int i = 1; i < prev_num_procs; ++i) {
				MPI_Recv(&f_sync_time, 1, MPI_FLOAT, i, i, prev_comm, MPI_STATUS_IGNORE);
				printf("%g (%g) ", f_sync_time, f_sync_time * 100 / (duration_cast<nanoseconds>(greedy_route_time).count() / 1e9));
			}
			printf("\n");

			float f_probe_time = duration_cast<nanoseconds>(mpi_perf.total_probe_time).count() / 1e9;
			printf("\t\t\tProbe time: %g (%g) ", f_probe_time, f_probe_time * 100.0 / (duration_cast<nanoseconds>(greedy_route_time).count() / 1e9));
			for (int i = 1; i < prev_num_procs; ++i) {
				MPI_Recv(&f_probe_time, 1, MPI_FLOAT, i, i, prev_comm, MPI_STATUS_IGNORE);
				printf("%g (%g) ", f_probe_time, f_probe_time * 100 / (duration_cast<nanoseconds>(greedy_route_time).count() / 1e9));
			}
			printf("\n");

			float f_broadcast_time = duration_cast<nanoseconds>(mpi_perf.total_broadcast_time).count() / 1e9;
			printf("\t\tBroadcast time: %g (%g) ", f_broadcast_time, f_broadcast_time * 100.0 / (duration_cast<nanoseconds>(greedy_route_time).count() / 1e9));
			for (int i = 1; i < prev_num_procs; ++i) {
				MPI_Recv(&f_broadcast_time, 1, MPI_FLOAT, i, i, prev_comm, MPI_STATUS_IGNORE);
				printf("%g (%g) ", f_broadcast_time, f_broadcast_time * 100 / (duration_cast<nanoseconds>(greedy_route_time).count() / 1e9));
			}
			printf("\n");

			printf("\tWait time: %g (%g) ", duration_cast<nanoseconds>(wait_time).count() / 1e9, duration_cast<nanoseconds>(wait_time).count() * 100.0 / duration_cast<nanoseconds>(greedy_route_time).count());
			for (int i = 1; i < prev_num_procs; ++i) {
				float f_wait_time;
				MPI_Recv(&f_wait_time, 1, MPI_FLOAT, i, i, prev_comm, MPI_STATUS_IGNORE);
				printf("%g (%g) ", f_wait_time, f_wait_time * 100 / (duration_cast<nanoseconds>(greedy_route_time).count() / 1e9));
			}
			printf("\n");

			printf("\tUpdate cost time: %g s.\n", duration_cast<nanoseconds>(update_cost_time).count() / 1e9);
			printf("\tAnalyze timing time: %g s.\n", duration_cast<nanoseconds>(analyze_timing_time).count() / 1e9);
			printf("\tCombine time: %g s.\n", duration_cast<nanoseconds>(combine_time).count() / 1e9);
			printf("Critical path: %g ns\n", crit_path_delay);
		} else {
			float f_sync_time = duration_cast<nanoseconds>(mpi_perf.total_sync_time).count() / 1e9;
			MPI_Send(&f_sync_time, 1, MPI_FLOAT, 0, prev_procid, prev_comm);

			float f_probe_time = duration_cast<nanoseconds>(mpi_perf.total_probe_time).count() / 1e9;
			MPI_Send(&f_probe_time, 1, MPI_FLOAT, 0, prev_procid, prev_comm);

			float f_broadcast_time = duration_cast<nanoseconds>(mpi_perf.total_broadcast_time).count() / 1e9;
			MPI_Send(&f_broadcast_time, 1, MPI_FLOAT, 0, prev_procid, prev_comm);

			float f_wait_time = duration_cast<nanoseconds>(wait_time).count() / 1e9;
			MPI_Send(&f_wait_time, 1, MPI_FLOAT, 0, prev_procid, prev_comm);
		}

		MPI_Barrier(prev_comm);

		prev_num_procs = num_procs;
		prev_comm = cur_comm;
		prev_procid = procid;

        //clock::time_point greedy_earliest_end_time = *std::min_element(begin(greedy_end_time), begin(greedy_end_time)+num_procs);

        //printf("greedy wait time: ");
        //for (int i = 0; i < num_procs; ++i) {
			//printf("%g (%g), ", duration_cast<nanoseconds>(greedy_end_time[i]-greedy_earliest_end_time).count() / 1e9, 100.0*(greedy_end_time[i]-greedy_earliest_end_time)/greedy_route_time);
        //}
		//printf("\n");
    }
	
done:
	MPI_Barrier(MPI_COMM_WORLD);
    //sprintf(buffer, "%s_run_%d_rr_graph_occ.txt", s_circuit_name, run);
    //dump_rr_graph_occ(congestion, num_vertices(g), buffer);
	if (initial_procid == 0) {
		if (routed) {
			printf("Routed in %d iterations. Total iteration time: %g\n", iter, duration_cast<nanoseconds>(total_iter_time).count() / 1e9);
			printf("\tTotal route time: %g s.\n", duration_cast<nanoseconds>(total_greedy_route_time).count() / 1e9);
			printf("\tTotal wait time: %g (%g) ", duration_cast<nanoseconds>(total_wait_time).count() / 1e9, duration_cast<nanoseconds>(total_wait_time).count() * 100.0 / duration_cast<nanoseconds>(total_greedy_route_time).count());
			for (int i = 1; i < initial_num_procs; ++i) {
				float f_total_wait_time;
				MPI_Recv(&f_total_wait_time, 1, MPI_FLOAT, i, i, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
				printf("%g (%g) ", f_total_wait_time, f_total_wait_time * 100.0 / (duration_cast<nanoseconds>(total_greedy_route_time).count() / 1e9));
			}
			printf("\n");
			printf("\tTotal update cost time: %g s.\n", duration_cast<nanoseconds>(total_update_cost_time).count() / 1e9);
			printf("\tTotal analyze timing time: %g s.\n", duration_cast<nanoseconds>(total_analyze_timing_time).count() / 1e9);
			printf("\tTotal combine time: %g s.\n", duration_cast<nanoseconds>(total_combine_time).count() / 1e9);

			printf("Final critical path: %g ns\n", crit_path_delay);
		} else {
			printf("Failed to route in %d iterations. Total iteration time: %g\n", opts->max_router_iterations, duration_cast<nanoseconds>(total_iter_time).count() / 1e9);
			printf("\tTotal route time: %g s.\n", duration_cast<nanoseconds>(total_greedy_route_time).count() / 1e9);
			printf("\tTotal wait time: %g (%g) ", duration_cast<nanoseconds>(total_wait_time).count() / 1e9, duration_cast<nanoseconds>(total_wait_time).count() * 100.0 / duration_cast<nanoseconds>(total_greedy_route_time).count());
			for (int i = 1; i < initial_num_procs; ++i) {
				float f_total_wait_time;
				MPI_Recv(&f_total_wait_time, 1, MPI_FLOAT, i, i, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
				printf("%g (%g) ", f_total_wait_time, f_total_wait_time * 100.0 / (duration_cast<nanoseconds>(total_greedy_route_time).count() / 1e9));
			}
			printf("\n");
			printf("\tTotal update cost time: %g s.\n", duration_cast<nanoseconds>(total_update_cost_time).count() / 1e9);
			printf("\tTotal analyze timing time: %g s.\n", duration_cast<nanoseconds>(total_analyze_timing_time).count() / 1e9);
			printf("\tTotal combine time: %g s.\n", duration_cast<nanoseconds>(total_combine_time).count() / 1e9);
		}
	} else {
		float f_total_wait_time = duration_cast<nanoseconds>(total_wait_time).count() / 1e9;
		MPI_Send(&f_total_wait_time, 1, MPI_FLOAT, 0, initial_procid, MPI_COMM_WORLD);
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
