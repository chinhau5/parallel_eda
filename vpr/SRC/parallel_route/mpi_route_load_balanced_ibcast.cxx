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
#include "route_net_mpi_ibcast.h"
#include "clock.h"
#include "queue.h"

void init_datatypes();

void sync_recalc_occ(congestion_t *congestion, int num_vertices, int procid, int num_procs, MPI_Comm comm);
void sync_nets(vector<net_t> &nets, vector<net_t> &global_nets, int procid, MPI_Comm comm);
void sync_net_delay(const vector<vector<net_t *>> &partitions, int procid, int num_procs, int *recvcounts, int *displs, int current_level, MPI_Comm comm, t_net_timing *net_timing);
void free_circuit();
void init_displ(const vector<vector<net_t*>> &partitions, int **recvcounts, int **displs);
void init_displ_nets(const vector<vector<net_t*>> &partitions, int **recvcounts, int **displs);
void get_sinks_to_route(net_t *net, const route_tree_t &rt, vector<sink_t *> &sinks_to_route);
void send_route_tree(const net_t *net, const vector<route_tree_t> &route_trees, int to_procid, MPI_Comm comm);
void recv_route_tree(const net_t *net, const RRGraph &g, vector<route_tree_t> &route_trees, t_net_timing *net_timing, int from_procid, MPI_Comm comm);
void init_route_structs(const RRGraph &g, const vector<net_t> &nets, const vector<net_t> &global_nets, route_state_t **states, congestion_t **congestion, vector<route_tree_t> &route_trees, t_net_timing **net_timing);
void broadcast_rip_up_all_ibcast(int net_id, mpi_context_t *mpi);
void progress_ibcast(const vector<net_t> &nets, congestion_t *congestion, const RRGraph &g, float pres_fac, vector<vector<RRNode>> &net_route_trees, mpi_context_t *mpi);
void progress_active_ibcast(const vector<net_t> &nets, congestion_t *congestion, const RRGraph &g, float pres_fac, vector<vector<RRNode>> &net_route_trees, mpi_context_t *mpi);
void progress_ibcast_blocking(const vector<net_t> &nets, congestion_t *congestion, const RRGraph &g, float pres_fac, vector<vector<RRNode>> &net_route_trees, mpi_context_t *mpi);
void route_net_mpi_ibcast(const RRGraph &g, int vpr_id, int net_id, const source_t *source, const vector<sink_t *> &sinks, const t_router_opts *params, float pres_fac, const vector<net_t> &nets, const vector<vector<net_t *>> &partition_nets, int current_net_index, route_state_t *state, congestion_t *congestion, route_tree_t &rt, t_net_timing &net_timing, vector<vector<RRNode>> &net_route_trees, mpi_context_t *mpi, perf_t *perf, mpi_perf_t *mpi_perf, bool delayed_progress);
int get_ibcast_buffer(mpi_context_t *mpi);
void bootstrap_ibcast(mpi_context_t *mpi);
void broadcast_as_root_large(int data_index, int count, mpi_context_t *mpi);
bool all_done_ibcast(mpi_context_t *mpi);
void dump_rr_graph(const RRGraph &g, const char *filename);

void write_header(int *packet, IbcastPacketID packet_id, unsigned int payload_size, int net_id);

struct net_graph_vertex_prop {
	int weight;
};

struct net_graph_edge_prop {
	int weight;
};

static void partition(vector<net_t *> &nets, bool load_balanced, int num_partitions, int current_level, int initial_comm_size, vector<int> &net_partition_id, vector<vector<net_t *>> &partition_nets)
{
	std::fill(begin(net_partition_id), end(net_partition_id), -1);

	partition_nets.resize(num_partitions);
	for (int i = 0; i < num_partitions; ++i) {
		partition_nets[i].clear();
	}

	if (load_balanced) {
		fast_graph_t<net_graph_vertex_prop, net_graph_edge_prop> net_g;

		add_vertex(net_g, nets.size());
		for (int i = 0; i < nets.size(); ++i) {
			get_vertex_props(net_g, nets[i]->local_id).weight = nets[i]->sinks.size();
		}

		assert(false);

		//for (int i = 0; i < nets_to_route.size(); ++i) {
			//for (int j = i + 1; j < nets_to_route.size(); ++j) {
				//if (bg::intersects(nets_to_route[i].first, nets_to_route[j].first)) {
					//add_edge(net_g, nets_to_route[i].second->local_id, nets_to_route[j].second->local_id);
				//}
			//}
		//}

		assert(nets.size() == num_vertices(net_g));

		if (num_partitions > 1) {
			partition_graph(net_g, num_partitions, 1, net_partition_id);
		} else {
			std::fill(begin(net_partition_id), end(net_partition_id), 0);
		}

		for (int i = 0; i < nets.size(); ++i) {
			int id = nets[i]->local_id;

			assert(net_partition_id[id] < partition_nets.size() && net_partition_id[id] >= 0);
			
			partition_nets[net_partition_id[id]].push_back(nets[i]);
		}
	} else {
		vector<net_t *> sorted_nets = nets;

		std::sort(begin(sorted_nets), end(sorted_nets), [] (const net_t *a, const net_t *b) -> bool {
				return a->sinks.size() > b->sinks.size();
				});

		for (int rank = 0; rank < num_partitions; ++rank) {
			for (int i = rank*pow(2, current_level); i < sorted_nets.size(); i += initial_comm_size) {
				for (int j = 0; j < pow(2, current_level) && i+j < sorted_nets.size(); ++j) {
					zlog_level(delta_log, ROUTER_V3, "Rank %d net index %d num sinks %d\n", rank, i+j, sorted_nets[i+j]->sinks.size());

					assert(sorted_nets[i+j]->local_id < net_partition_id.size());

					net_partition_id[sorted_nets[i+j]->local_id] = rank;

					partition_nets[rank].push_back(sorted_nets[i+j]);
				}
			}
		}
	}
}

static void move_route_tree(const net_t *net, int old_id, int new_id, const RRGraph &g, vector<route_tree_t> &route_trees, vector<vector<RRNode>> &net_route_trees, t_net_timing *net_timing, mpi_context_t &mpi)
{
	if (old_id == new_id) {
		return;
	}

	if (mpi.rank == old_id) { 
		zlog_level(delta_log, ROUTER_V3, "Sending net %d from %d\n", net->vpr_id, mpi.rank);

		assert(net_route_trees[net->local_id].empty());
		assert(!route_tree_empty(route_trees[net->local_id]));

		send_route_tree(net, route_trees, new_id, mpi.comm);

		for (const auto &rt_node : route_tree_get_nodes(route_trees[net->local_id])) {
			const auto &rt_node_p = get_vertex_props(route_trees[net->local_id].graph, rt_node);
			net_route_trees[net->local_id].push_back(rt_node_p.rr_node);
		}

		route_tree_clear(route_trees[net->local_id]);
	} else if (mpi.rank == new_id) {
		zlog_level(delta_log, ROUTER_V3, "Recving net %d from %d\n", net->vpr_id, old_id);

		assert(!net_route_trees[net->local_id].empty());
		assert(route_tree_empty(route_trees[net->local_id]));

		recv_route_tree(net, g, route_trees, net_timing, old_id, mpi.comm);

		net_route_trees[net->local_id].clear();
	} else {
		/* the route tree is not owned by the sending or recving process */
		assert(!net_route_trees[net->local_id].empty());
		assert(route_tree_empty(route_trees[net->local_id]));
		//if (mpi.rank < mpi.comm_size/2) {
		//printf("rank %d < comm size %d old %d new %d\n", mpi.rank, mpi.comm_size/2, old_net_partition_id[net.local_id], net_partition_id[net.local_id]);
		//assert(false);
		//}
		//assert(mpi.rank >= mpi.comm_size/2);
	}
}

static void move_route_trees(const vector<net_t *> &nets, const vector<int> &old_net_partition_id, const vector<int> &net_partition_id, const RRGraph &g, vector<route_tree_t> &route_trees, vector<vector<RRNode>> &net_route_trees, t_net_timing *net_timing, mpi_context_t &mpi)
{
	for (const auto &net : nets) {
		assert(old_net_partition_id[net->local_id] != -1);
		assert(net_partition_id[net->local_id] != -1);

		if (old_net_partition_id[net->local_id] != net_partition_id[net->local_id]) {
			move_route_tree(net, old_net_partition_id[net->local_id], net_partition_id[net->local_id], g, route_trees, net_route_trees, net_timing, mpi);
		} 
	}

	//num_recvs_called.resize(mpi.comm_size / 2);
	//for (int i = 0; i < num_recvs_called.size(); ++i) {
	//num_recvs_called[i] = 0;
	//}

	//zlog_info(delta_log, "num_recvs_required before:\n");
	//for (int i = 0; i < num_recvs_required.size(); ++i) {
	//zlog_info(delta_log, "%d ", num_recvs_required[i]);
	//}
	//zlog_info(delta_log, "\n");

	//for (int i = 0; i < mpi.comm_size; i += 2) {
	//num_recvs_required[i/2] += num_recvs_required[i+1];
	//}
	//num_recvs_required.resize(mpi.comm_size / 2);
	//orig_num_recvs_required = num_recvs_required;

	//zlog_info(delta_log, "num_recvs_required after:\n");
	//for (int i = 0; i < num_recvs_required.size(); ++i) {
	//zlog_info(delta_log, "%d ", num_recvs_required[i]);
	//}
	//zlog_info(delta_log, "\n");
	//
}

static void repartition(vector<net_t *> &nets, bool load_balanced, int num_partitions, int current_level, int initial_comm_size, vector<int> &net_partition_id, vector<vector<net_t *>> &partition_nets, const RRGraph &g, vector<route_tree_t> &route_trees, vector<vector<RRNode>> &net_route_trees, t_net_timing *net_timing, mpi_context_t &mpi)
{
	vector<int> old_net_partition_id = net_partition_id;

	partition(nets, load_balanced, num_partitions, current_level, initial_comm_size, net_partition_id, partition_nets);

	move_route_trees(nets, old_net_partition_id, net_partition_id, g, route_trees, net_route_trees, net_timing, mpi);
}

void test_queue()
{
	queue_t<int> q;
	q_init(&q, 5);
	assert(q_empty(&q));
	assert(q_size(&q) == 0);

	assert(q_push(&q, 5));
	assert(!q_empty(&q));
	assert(q_size(&q) == 1);
	assert(q_front(&q) == 5);

	assert(q_push(&q, 6));
	assert(q_front(&q) == 5);
	assert(q_size(&q) == 2);

	assert(q_pop(&q));
	assert(q_front(&q) == 6);

	assert(q_pop(&q));
	assert(q_size(&q) == 0);

	assert(!q_pop(&q));

	for (int i = 0; i < 5; ++i) {
		assert(q_size(&q) == i);
		assert(q_push(&q, i));
	}
	assert(!q_push(&q, 100));
	assert(q_size(&q) == 5);
}

int get_max_num_out_edges(const RRGraph &g)
{
	int max_num_edges = std::numeric_limits<int>::min();
	for (const auto &n : get_vertices(g)) {
		max_num_edges = std::max(max_num_edges, num_out_edges(g, n));
	} 
	return max_num_edges;
}

bool mpi_route_load_balanced_ibcast(t_router_opts *opts, struct s_det_routing_arch det_routing_arch, t_direct_inf *directs, int num_directs, t_segment_inf *segment_inf, t_timing_inf timing_inf)
{
    using std::chrono::duration_cast;
    using std::chrono::nanoseconds;
	using clock = myclock;

	int initial_comm_size, prev_comm_size, initial_rank, prev_rank;
	MPI_Comm prev_comm;

	mpi_context_t mpi;

	mpi.comm = MPI_COMM_WORLD;
    MPI_Comm_size(mpi.comm, &initial_comm_size);
    MPI_Comm_rank(mpi.comm, &initial_rank);
	mpi.comm_size = initial_comm_size;
	mpi.rank = initial_rank;

	prev_comm = mpi.comm;
	prev_comm_size = mpi.comm_size;
	prev_rank = mpi.rank;

	printf("[%d] Initializing router\n", initial_rank);

	auto init_dt_start = clock::now();

	init_congestion_mpi_datatype();
	init_datatypes();

	auto init_dt_time = clock::now()-init_dt_start;
	printf("[%d] init dt time: %g\n", initial_rank, duration_cast<nanoseconds>(init_dt_time).count() / 1e9);

    init_logging();
    zlog_set_record("custom_output", concurrent_log_impl);

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

	char buffer[256];
	//sprintf(buffer, "rr_graph_%d.txt", mpi.rank);
	//dump_rr_graph(partitioner.orig_g, buffer);

	int max_num_edges = get_max_num_out_edges(partitioner.orig_g);
	printf("max edges: %d\n", max_num_edges);
	mpi.bit_width = ceil(std::log2(max_num_edges+1));

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

	sprintf(buffer, "%d", mpi.rank);
	zlog_put_mdc("tid", buffer);

    vector<net_t> nets;
    vector<net_t> global_nets;
    //init_nets(nets, global_nets, opts->bb_factor, partitioner.sink_in_nodes);	
	if (mpi.rank == 0) {
		printf("[%d] initializing nets\n", mpi.rank);
		init_nets(nets, global_nets, opts->bb_factor, opts->large_bb);	

		//printf("Num_interpatition_nets [%d nets]: ", nets.size());
		//for (int i = 2; i <= 16; i *= 2) {
			//printf("%d ", get_num_interpartition_nets(nets, i));
		//}
		//printf("\n");
		//exit(0);
	}
	printf("[%d] syncing nets\n", mpi.rank);
	sync_nets(nets, global_nets, mpi.rank, mpi.comm);

	printf("[%d] done initializing nets\n", mpi.rank);

    vector<net_t *> nets_ptr(nets.size());
	vector<net_t *> non_congested_nets_ptr;
    for (int i = 0; i < nets.size(); ++i) {
        nets_ptr[i] = &nets[i];
    }
    //extern s_bb *route_bb;
    //sort(nets_ptr.begin(), nets_ptr.end(), [] (const net_t *a, const net_t *b) -> bool {
            //return get_bounding_box_area(a->bounding_box) > get_bounding_box_area(b->bounding_box);
            //});
    //int rank = 0;
    //for (auto &net : nets_ptr) {
        //net->bb_area_rank = rank++;
    //}

	vector<int> num_sinks(nets.size());
	for (int i = 0; i < nets.size(); ++i) {
		num_sinks[i] = nets[i].sinks.size();
	}
	sort(begin(num_sinks), end(num_sinks));
	float total_num_sinks = 0;
	for (const auto &n : num_sinks) {
		total_num_sinks += n;
	}
	printf("Num sinks, min: %d max: %d mean: %g\n", num_sinks.front(), num_sinks.back(), total_num_sinks / nets.size());

	/* per process */
	route_state_t *states;

	/* synchronized while routing */
	congestion_t *congestion;

	/* explicitly synchronized after routing */
	vector<route_tree_t> route_trees;
	vector<vector<RRNode>> net_route_trees(nets.size());
	t_net_timing *net_timing;
	vector<vector<clock::duration>> net_route_time(opts->max_router_iterations, vector<clock::duration>(nets.size(), clock::duration::zero()));
	vector<vector<clock::duration>> net_mpi_time(opts->max_router_iterations, vector<clock::duration>(nets.size(), clock::duration::zero()));

	init_route_structs(partitioner.orig_g, nets, global_nets, &states, &congestion, route_trees, &net_timing);

    route_parameters_t params;
    params.criticality_exp = opts->criticality_exp;
    params.astar_fac = opts->astar_fac;
    params.max_criticality = opts->max_criticality;
    params.bend_cost = opts->bend_cost;

	float pres_fac = opts->first_iter_pres_fac; /* Typically 0 -> ignore cong. */

	if (mpi.rank != 0) {
		free_circuit();
	}
	free_rr_graph();

    //vector<pair<box, net_t *>> nets_to_route;
    
    //for (auto &net : nets) {
        //box b = bg::make_inverse<box>();

        //bg::expand(b, point(net.source.x, net.source.y));
        //for (const auto &sink : net.sinks) {
            //bg::expand(b, point(sink.x, sink.y));
        //}
        //bg::subtract_value(b.min_corner(), 1+opts->bb_factor);
        //bg::add_value(b.max_corner(), opts->bb_factor);

        //nets_to_route.push_back(make_pair(b, &net));
    //}
    //std::sort(begin(nets_to_route), end(nets_to_route), [] (const pair<box, net_t *> &a, const pair<box, net_t *> &b) -> bool {
            ////return a.second->sinks.size()*get_bounding_box_area(a.second->bounding_box) > b.second->sinks.size()*get_bounding_box_area(b.second->bounding_box);
			//return a.second->sinks.size() > b.second->sinks.size();
            //});

	int current_level = 0;

	vector<int> net_partition_id(nets.size());
	vector<vector<net_t *>> partition_nets;

	partition(nets_ptr, opts->load_balanced, mpi.comm_size, current_level, initial_comm_size, net_partition_id, partition_nets);

	//vector<int> num_recvs_called(mpi.comm_size);
	//vector<int> num_recvs_required(mpi.comm_size, 0);
	//vector<int> orig_num_recvs_required;

	//for (int i = 0; i < num_recvs_required.size(); ++i) {
		//for (int j = i; j < nets_to_route.size(); j += mpi.comm_size) {
			//num_recvs_required[i] += nets_to_route[j].second->sinks.size();
		//}
	//}
	//orig_num_recvs_required = num_recvs_required;
	//

    bool routed = false;
	bool idling = false;

    int iter;
    float crit_path_delay;

	//vector<vector<sink_t *>> fixed_sinks(nets.size());
	unsigned long prev_num_overused_nodes = std::numeric_limits<unsigned long>::max();

    clock::duration total_actual_route_time = clock::duration::zero();
    clock::duration total_route_time = clock::duration::zero();
    clock::duration total_update_cost_time = clock::duration::zero();
    clock::duration total_analyze_timing_time = clock::duration::zero();
    clock::duration total_iter_time = clock::duration::zero();
    clock::duration total_wait_time = clock::duration::zero();
    clock::duration total_combine_time = clock::duration::zero();
    clock::duration total_last_sync_time = clock::duration::zero();
    clock::duration total_last_progress_time = clock::duration::zero();
	 //mpi 
    clock::duration total_sync_time = clock::duration::zero();
    clock::duration total_broadcast_time = clock::duration::zero();
    clock::duration total_send_testsome_time = clock::duration::zero();
	int total_calls = 0;

	int max_num_broadcasts;

	mpi.ibcast_comm.resize(mpi.comm_size);
	for (int i = 0; i < mpi.comm_size; ++i) {
		int error = MPI_Comm_dup(mpi.comm, &mpi.ibcast_comm[i]);
		assert(error == MPI_SUCCESS);
	}

	mpi.max_ibcast_count.resize(mpi.comm_size);
	std::fill(begin(mpi.max_ibcast_count), end(mpi.max_ibcast_count), 32);

	for (int i = 0; i < mpi.comm_size; ++i) {
		zlog_level(delta_log, ROUTER_V3, "Preallocating buffer size %d for rank %d\n", mpi.max_ibcast_count[i], i);
		mpi.pending_send_data_nbc.emplace_back(new vector<int>(mpi.max_ibcast_count[i]));
		mpi.pending_send_data_ref_count.emplace_back(0);
	}

	mpi.pending_send_req.resize(mpi.comm_size, MPI_REQUEST_NULL);
	mpi.pending_req_meta.resize(mpi.comm_size, { -1, -1, -1 });
	mpi.num_pending_reqs = 0;
	mpi.received_last_update.resize(mpi.comm_size);

    for (iter = 0; iter < opts->max_router_iterations && !routed && !idling; ++iter) {
        clock::duration actual_route_time = clock::duration::zero();
        clock::duration route_time = clock::duration::zero();
        clock::duration update_cost_time = clock::duration::zero();
        clock::duration analyze_timing_time = clock::duration::zero();
        clock::duration iter_time = clock::duration::zero();
        clock::duration wait_time = clock::duration::zero();
        clock::duration combine_time = clock::duration::zero();
        clock::duration last_sync_time = clock::duration::zero();
        clock::duration last_progress_time = clock::duration::zero();

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
		//
		mpi.max_send_data_size = 0;
		mpi.total_send_count = 0;
		mpi.total_send_data_size = 0;
		mpi.max_req_buffer_size = 0;

		std::fill(begin(mpi.received_last_update), end(mpi.received_last_update), false);

		sprintf(buffer, "%d", iter);
		zlog_put_mdc("iter", buffer);

		sprintf(buffer, "%d", mpi.rank);
		zlog_put_mdc("tid", buffer);

		//extern map<string, FILE *> log_files;

		//char fname_buf[256];
		//sprintf(fname_buf, "iter_%d_tid_%d.log", iter, mpi.rank);

		//auto file_iter = log_files.find(fname_buf);
		//assert(file_iter != log_files.end());
		//pair<string, FILE *> file = *file_iter;
		//int fd = fileno(file.second);
		//assert(fd != -1);

		//assert(dup2(fd, 1) != -1);

        zlog_info(delta_log, "Routing iteration: %d\n", iter);
		if (mpi.rank == 0) {
			printf("Routing iteration: %d\n", iter);
		}

		MPI_Barrier(mpi.comm);
        
        auto iter_start = clock::now();

        auto route_start = clock::now();

		perf.num_heap_pushes = 0;
		perf.num_heap_pops = 0;
		perf.num_neighbor_visits = 0;

		mpi_perf.total_sync_time = clock::duration::zero();
		mpi_perf.total_broadcast_time = clock::duration::zero();
		mpi_perf.total_send_testsome_time = clock::duration::zero();
		mpi_perf.total_calls = 0;

		//for (int i = 0; i < mpi.comm_size; ++i) {
			//num_recvs_called[i] = 0;
		//}

		//if (iter > 0) {
			//for (int i = 0; i < mpi.comm_size; ++i) {
				//num_recvs_required[i] = orig_num_recvs_required[i] * 2;
			//}
		//}
		//

		if (iter == 1 && opts->load_balanced) {
			clock::rep min_net_route_time = std::numeric_limits<clock::rep>::max();
			for (const auto &net : partition_nets[mpi.rank]) {
				min_net_route_time = std::min(net_route_time[iter-1][net->local_id].count(), min_net_route_time);
			}

			unsigned long min_net_route_time_ul = static_cast<unsigned long>(min_net_route_time);
			assert(min_net_route_time_ul == min_net_route_time);

			unsigned long reduced_min_net_route_time;

			MPI_Allreduce(&min_net_route_time_ul, &reduced_min_net_route_time, 1, MPI_UNSIGNED_LONG, MPI_MIN, mpi.comm);

			//printf("min net route time: %ld\n", min_net_route_time);
			//printf("min net route time ul: %lu\n", min_net_route_time_ul);
			//printf("min net route time reduced: %lu\n", reduced_min_net_route_time);

			vector<int> i_net_route_time;
			for (const auto &net : partition_nets[mpi.rank]) {
				i_net_route_time.push_back(net_route_time[iter-1][net->local_id].count() / reduced_min_net_route_time);
			}
			
			vector<int> all_net_route_time(nets.size());

			int *recvcounts_nets = nullptr;
			int *displs_nets = nullptr;

			init_displ_nets(partition_nets, &recvcounts_nets, &displs_nets);

			MPI_Allgatherv(i_net_route_time.data(), recvcounts_nets[mpi.rank], MPI_INT, all_net_route_time.data(), recvcounts_nets, displs_nets, MPI_INT, mpi.comm);

			assert(partition_nets.size() == mpi.comm_size);

			//for (int i = 0; i < partition_nets.size(); ++i) {
				//for (int j = 0; j < partition_nets[i].size(); ++j) {
					//assert(recvcounts_nets[i] == partition_nets[i].size());
					//get_vertex_props(net_g, partition_nets[i][j]->local_id).weight = all_net_route_time[displs_nets[i]+j];
				//}
			//}

			//if (mpi.rank == 0) {
				//FILE *file = fopen("normalized_route_time.txt", "w");

				//fprintf(file, "min %lu\n", reduced_min_net_route_time);

				//for (const auto &v : get_vertices(net_g)) {
					//fprintf(file, "%d %d\n", v, get_vertex_props(net_g, v).weight);
				//}
				//fclose(file);
			//}

			repartition(nets_ptr, opts->load_balanced, mpi.comm_size, current_level, initial_comm_size, net_partition_id, partition_nets, partitioner.orig_g, route_trees, net_route_trees, net_timing, mpi);

			//for (int i = 0; i < mpi.comm_size; ++i) {
				//int total_weight = 0;
				//for (const auto &net : partition_nets[i]) {
					//total_weight += get_vertex_props(net_g, net->local_id).weight;
				//}
				////printf("part %d weight %d\n", i, total_weight);
			//}

		}

		bootstrap_ibcast(&mpi);

		int current_net_index = 0;

		for (auto &net : partition_nets[mpi.rank]) {
			auto old_broadcast_time = mpi_perf.total_broadcast_time;
			auto old_progress_time = mpi_perf.total_send_testsome_time;

			auto net_route_start = clock::now();

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

			//route_tree_mark_paths_to_be_ripped(route_trees[net->local_id], partitioner.orig_g, partitioner.result_pid_by_level[current_level], mpi.rank, sinks_to_mark);

			if (iter > 0) {
				int route_tree_size = route_tree_num_nodes(route_trees[net->local_id]);
				int num_marked = 0;

				if (opts->rip_up_always) {
					route_tree_mark_all_nodes_to_be_ripped(route_trees[net->local_id], partitioner.orig_g);
					num_marked = route_tree_size;
				} else {
					num_marked = route_tree_mark_congested_nodes_to_be_ripped(route_trees[net->local_id], partitioner.orig_g, congestion);
				}

				int data_index = get_ibcast_buffer(&mpi);
				vector<int> &data = *mpi.pending_send_data_nbc[data_index];

				int num_nodes_ripped = 0;

				route_tree_rip_up_marked_mpi_collective(route_trees[net->local_id], partitioner.orig_g, congestion, pres_fac, [&] (const RRNode &rr_node) -> void {
						update_one_cost_internal(rr_node, partitioner.orig_g, congestion, -1, pres_fac); 

						if (num_marked < route_tree_size) {
							if (data.size() <= 2+num_nodes_ripped) {
								data.resize((2+num_nodes_ripped)*2);
							}

							data[2+num_nodes_ripped] = rr_node;
						}

						++num_nodes_ripped;
						}
						);

				assert(num_nodes_ripped == num_marked);

				zlog_level(delta_log, ROUTER_V3, "Ripped up net %d route tree\n", net->vpr_id);

				if (mpi.comm_size > 1) {
					int count;
					if (num_nodes_ripped == route_tree_size) {
						write_header(&data[0], IbcastPacketID::RIP_UP_ALL, 0, net->local_id);
						count = 2;
					} else {
						write_header(&data[0], IbcastPacketID::RIP_UP, num_nodes_ripped, net->local_id);
						count = 2+num_nodes_ripped;
					} 

					auto broadcast_start = clock::now();

					broadcast_as_root_large(data_index, count, &mpi);

					mpi_perf.total_broadcast_time += clock::now()-broadcast_start;

					auto progress_start = clock::now();

					progress_ibcast(nets, congestion, partitioner.orig_g, pres_fac, net_route_trees, &mpi);

					mpi_perf.total_send_testsome_time += clock::now()-progress_start;
				} else {
					assert(mpi.pending_send_data_ref_count[data_index] == 0);
					mpi.free_send_data_index.push(data_index);
				}
			} 

			vector<sink_t *> sinks;	
			get_sinks_to_route(net, route_trees[net->local_id], sinks);

			//local_perf.total_rip_up_time += clock::now()-rip_up_start;

			//auto route_start = clock::now();

			if (!sinks.empty()) {
				route_net_mpi_ibcast(partitioner.orig_g, net->vpr_id, net->local_id, &net->source, sinks, opts, pres_fac, nets, partition_nets, current_net_index, states, congestion, route_trees[net->local_id], net_timing[net->vpr_id], net_route_trees, &mpi, &perf, &mpi_perf, false);

				++thread_num_nets_routed;
				++thread_num_nets_to_route;

				thread_num_sinks_to_route += sinks.size();
				thread_num_sinks_routed += sinks.size();
			} else {
				zlog_level(delta_log, ROUTER_V2, "Net %d has no sinks in the current iteration\n", net->vpr_id);
			}

			auto end = clock::now();
			net_route_time[iter][net->local_id] = end-net_route_start;
			net_mpi_time[iter][net->local_id] = mpi_perf.total_broadcast_time - old_broadcast_time + mpi_perf.total_send_testsome_time - old_progress_time;

			//local_perf.total_route_time += clock::now()-rip_up_start;
			++current_net_index;
		}

		actual_route_time = clock::now() - route_start;

		int num_last_syncs = 0;

		if (mpi.comm_size > 1) {
			auto last_sync_start = clock::now();

			int data_index = get_ibcast_buffer(&mpi);

			int *data = mpi.pending_send_data_nbc[data_index]->data();
			write_header(data, IbcastPacketID::TRAILER, 0, 0);

			zlog_level(delta_log, ROUTER_V2, "Sent trailer packet\n");

			broadcast_as_root_large(data_index, 2, &mpi);

			last_sync_time += clock::now()-last_sync_start;

			auto last_progress_start = clock::now();
			while (!all_done_ibcast(&mpi) || mpi.num_pending_reqs > 0) {
				progress_ibcast(nets, congestion, partitioner.orig_g, pres_fac, net_route_trees, &mpi);
				++num_last_syncs;
			}
			//progress_ibcast_blocking(nets, congestion, partitioner.orig_g, pres_fac, net_route_trees, &mpi);
			last_progress_time += clock::now()-last_progress_start;
		}

		/* checking */
		for (int i = 0; i < mpi.pending_send_req.size(); ++i) {
			assert(mpi.pending_send_req[i] == MPI_REQUEST_NULL);
		}
		assert(mpi.num_pending_reqs == 0);
		for (int i = 0; i < mpi.pending_send_data_ref_count.size(); ++i) {
			assert(mpi.pending_send_data_ref_count[i] == 0);
		}
		zlog_level(delta_log, ROUTER_V2, "free data %d free req %d\n", mpi.free_send_data_index.size(), mpi.free_send_req_index.size());
		zlog_level(delta_log, ROUTER_V2, "data %d req %d\n", mpi.pending_send_data_nbc.size(), mpi.pending_send_req.size());
		assert(mpi.free_send_data_index.size() == mpi.pending_send_data_nbc.size()-mpi.comm_size);
		assert(mpi.free_send_req_index.size() == mpi.pending_send_req.size()-mpi.comm_size);

		//for (int pid = 0; pid < mpi.comm_size; ++pid) {
			//if (pid == mpi.rank) {
				//assert(mpi.pending_recv_req_flat[pid] == MPI_REQUEST_NULL);
			//} else {
				//assert(mpi.pending_recv_req_flat[pid] != MPI_REQUEST_NULL);
				//MPI_Cancel(&mpi.pending_recv_req_flat[pid]);
			//}
		//}
		//int error = MPI_Waitall(mpi.comm_size, mpi.pending_recv_req_flat.data(), MPI_STATUSES_IGNORE);
		//assert(error == MPI_SUCCESS);

		//vector<MPI_Request> requests;
		//for (int pid = 0; pid < mpi.comm_size; ++pid) {
			//if (pid != mpi.rank) {
				//for (auto &recv: pending_recvs[pid]) {
					//assert(recv.req != MPI_REQUEST_NULL);
					//MPI_Cancel(&recv.req);
					//assert(recv.req != MPI_REQUEST_NULL);
					//requests.push_back(recv.req);
				//}
			//}
		//}

		//MPI_Waitall(requests.size(), requests.data(), MPI_STATUSES_IGNORE);

		auto wait_start = clock::now();
		MPI_Barrier(mpi.comm);
		wait_time = clock::now()-wait_start;

		//assert(pending_recvs.empty());

		//MPI_Barrier(mpi.comm);

		//greedy_end_time = clock::now();

        route_time = clock::now()-route_start;

		//if (has_interpartition_sinks) {
			//for (const auto &net_interpartition_sinks : interpartition_sinks) {
				//for (const auto &is : net_interpartition_sinks) {
					//vector<RRNode> added_nodes;
					//for (const auto &node : is.path) {
						//if (node.update_cost) {
							//added_nodes.push_back(node.rr_node_id);
						//}
					//}
					//update_one_cost(partitioner.orig_g, congestion, added_nodes.begin(), added_nodes.end(), 1, pres_fac, false, nullptr);
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

		sprintf(buffer, "net_route_time_iter_%d_%d.txt", iter, mpi.rank);
		FILE *net_route_time_f = fopen(buffer, "w");
		for (const auto &net : partition_nets[mpi.rank]) {
			float route_time = duration_cast<nanoseconds>(net_route_time[iter][net->local_id]).count() / 1e9;
			float mpi_time = duration_cast<nanoseconds>(net_mpi_time[iter][net->local_id]).count() / 1e9;
			fprintf(net_route_time_f, "%d %lu %g %g %ld %g %ld %g\n", net->local_id, net->sinks.size(), bg::area(box(point(net->bounding_box.xmin, net->bounding_box.ymin), point(net->bounding_box.xmax, net->bounding_box.ymax))), route_time, net_route_time[iter][net->local_id].count(), mpi_time, net_mpi_time[iter][net->local_id].count(), mpi_time/route_time*100);
		}
		fclose(net_route_time_f);

		total_sync_time += mpi_perf.total_sync_time;
		total_broadcast_time += mpi_perf.total_broadcast_time;
		total_send_testsome_time += mpi_perf.total_send_testsome_time;
		total_last_progress_time += last_progress_time;
		total_last_sync_time += last_sync_time;
		total_calls += mpi_perf.total_calls;

		int total_num_sinks_to_route;
		if (mpi.rank == 0) {
			int total_num_nets_to_route = thread_num_nets_to_route;
			int total_num_nets_routed = thread_num_nets_routed;
			int total_num_sinks_routed = thread_num_sinks_routed;
			int total_num_sinks_to_route = thread_num_sinks_to_route;

			printf("num nets routed: %d/%d (%g) ", thread_num_nets_routed, thread_num_nets_to_route, thread_num_nets_routed*100.0/thread_num_nets_to_route);
			for (int i = 1; i < mpi.comm_size; ++i) {
				int tmp_thread_num_nets_to_route;
				int tmp_thread_num_nets_routed;
				MPI_Recv(&tmp_thread_num_nets_to_route, 1, MPI_INT, i, i, mpi.comm, MPI_STATUS_IGNORE);
				MPI_Recv(&tmp_thread_num_nets_routed, 1, MPI_INT, i, i, mpi.comm, MPI_STATUS_IGNORE);
				total_num_nets_to_route += tmp_thread_num_nets_to_route;
				total_num_nets_routed += tmp_thread_num_nets_routed;

				printf("%d/%d (%g) ", tmp_thread_num_nets_routed, tmp_thread_num_nets_to_route, tmp_thread_num_nets_routed*100.0/tmp_thread_num_nets_to_route);
			}
			printf("\n");

			assert(total_num_nets_to_route == total_num_nets_routed);

			printf("num sinks routed: %d/%d (%g) ", thread_num_sinks_routed, thread_num_sinks_to_route, thread_num_sinks_routed*100.0/thread_num_sinks_to_route);
			for (int i = 1; i < mpi.comm_size; ++i) {
				int tmp_thread_num_sinks_to_route;
				int tmp_thread_num_sinks_routed;
				MPI_Recv(&tmp_thread_num_sinks_to_route, 1, MPI_INT, i, i, mpi.comm, MPI_STATUS_IGNORE);
				MPI_Recv(&tmp_thread_num_sinks_routed, 1, MPI_INT, i, i, mpi.comm, MPI_STATUS_IGNORE);

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
			for (int i = 1; i < mpi.comm_size; ++i) {
				unsigned long tmp_num_heap_pushes;
				MPI_Recv(&tmp_num_heap_pushes, 1, MPI_UNSIGNED_LONG, i, i, mpi.comm, MPI_STATUS_IGNORE);
				total_num_heap_pushes += tmp_num_heap_pushes;
				printf("%lu ", tmp_num_heap_pushes);
			}
			printf("\n");

			printf("num_heap_pops: %lu ", perf.num_heap_pops);
			for (int i = 1; i < mpi.comm_size; ++i) {
				unsigned long tmp_num_heap_pops;
				MPI_Recv(&tmp_num_heap_pops, 1, MPI_UNSIGNED_LONG, i, i, mpi.comm, MPI_STATUS_IGNORE);
				total_num_heap_pops += tmp_num_heap_pops;
				printf("%lu ", tmp_num_heap_pops);
			}
			printf("\n");

			printf("num_neighbor_visits: %lu ", perf.num_neighbor_visits);
			for (int i = 1; i < mpi.comm_size; ++i) {
				unsigned long tmp_num_neighbor_visits;
				MPI_Recv(&tmp_num_neighbor_visits, 1, MPI_UNSIGNED_LONG, i, i, mpi.comm, MPI_STATUS_IGNORE);
				total_num_neighbor_visits += tmp_num_neighbor_visits;
				printf("%lu ", tmp_num_neighbor_visits);
			}
			printf("\n");

			printf("total_num_heap_pushes: %lu\n", total_num_heap_pushes);
			printf("total_num_heap_pops: %lu\n", total_num_heap_pops); 
			printf("total_num_neighbor_visits: %lu\n", total_num_neighbor_visits);

			printf("max send data size: %d ", mpi.max_send_data_size*sizeof(int));
			for (int i = 1; i < mpi.comm_size; ++i) {
				int tmp_max_send_data_size;
				MPI_Recv(&tmp_max_send_data_size, 1, MPI_INT, i, i, mpi.comm, MPI_STATUS_IGNORE);
				printf("%d ", tmp_max_send_data_size*sizeof(int));
			}
			printf("\n");
			printf("total send size/count/average: %lu/%lu/%g ", mpi.total_send_data_size*sizeof(int), mpi.total_send_count, (float)mpi.total_send_data_size*sizeof(int)/mpi.total_send_count);
			for (int i = 1; i < mpi.comm_size; ++i) {
				unsigned long tmp_total_send_data_size;
				unsigned long tmp_total_send_count;
				MPI_Recv(&tmp_total_send_data_size, 1, MPI_UNSIGNED_LONG, i, i, mpi.comm, MPI_STATUS_IGNORE);
				MPI_Recv(&tmp_total_send_count, 1, MPI_UNSIGNED_LONG, i, i, mpi.comm, MPI_STATUS_IGNORE);
				printf("%lu/%lu/%g ", tmp_total_send_data_size*sizeof(int), tmp_total_send_count, (float)tmp_total_send_data_size*sizeof(int)/tmp_total_send_count);
			}
			printf("\n");
			printf("max req buffer size: %d ", mpi.max_req_buffer_size);
			for (int i = 1; i < mpi.comm_size; ++i) {
				int tmp_max_req_buffer_size;
				MPI_Recv(&tmp_max_req_buffer_size, 1, MPI_INT, i, i, mpi.comm, MPI_STATUS_IGNORE);
				printf("%d ", tmp_max_req_buffer_size);
			}
			printf("\n");
		} else {
			MPI_Send(&thread_num_nets_to_route, 1, MPI_INT, 0, mpi.rank, mpi.comm);
			MPI_Send(&thread_num_nets_routed, 1, MPI_INT, 0, mpi.rank, mpi.comm);
			MPI_Send(&thread_num_sinks_to_route, 1, MPI_INT, 0, mpi.rank, mpi.comm);
			MPI_Send(&thread_num_sinks_routed, 1, MPI_INT, 0, mpi.rank, mpi.comm);

			MPI_Send(&perf.num_heap_pushes, 1, MPI_UNSIGNED_LONG, 0, mpi.rank, mpi.comm);
			MPI_Send(&perf.num_heap_pops, 1, MPI_UNSIGNED_LONG, 0, mpi.rank, mpi.comm);
			MPI_Send(&perf.num_neighbor_visits, 1, MPI_UNSIGNED_LONG, 0, mpi.rank, mpi.comm);

			MPI_Send(&mpi.max_send_data_size, 1, MPI_INT, 0, mpi.rank, mpi.comm);
			MPI_Send(&mpi.total_send_data_size, 1, MPI_UNSIGNED_LONG, 0, mpi.rank, mpi.comm);
			MPI_Send(&mpi.total_send_count, 1, MPI_UNSIGNED_LONG, 0, mpi.rank, mpi.comm);
			MPI_Send(&mpi.max_req_buffer_size, 1, MPI_INT, 0, mpi.rank, mpi.comm);
		}

		/* checking */
        for (int i = 0; i < num_vertices(partitioner.orig_g); ++i) {
            congestion[i].recalc_occ = 0; 
        }

		if (mpi.rank == 0) {
			for (const auto &net : non_congested_nets_ptr) {
				zlog_level(delta_log, ROUTER_V3, "Checking non-congested net vpr id %d\n", net->vpr_id);

				check_route_tree(route_trees[net->local_id], *net, partitioner.orig_g);
				recalculate_occ(route_trees[net->local_id], partitioner.orig_g, congestion);
			}
		}

		for (auto &net : partition_nets[mpi.rank]) {
			check_route_tree(route_trees[net->local_id], *net, partitioner.orig_g);
			recalculate_occ(route_trees[net->local_id], partitioner.orig_g, congestion);
		}

		sync_recalc_occ(congestion, num_vertices(partitioner.orig_g), mpi.rank, mpi.comm_size, mpi.comm);

		bool valid = true;
		for (int i = 0; i < num_vertices(partitioner.orig_g); ++i) {
			if (congestion[i].recalc_occ != congestion[i].occ) {
				sprintf_rr_node_impl(i, buffer);
				printf("[%d] Node %s occ mismatch, recalc: %d original: %d\n", mpi.rank, buffer, congestion[i].recalc_occ, congestion[i].occ);
				valid = false;
			}
		}
		assert(valid);

        int overused_total_bb_rank = 0;
        int num_congested_nets = 0;

		for (auto &net : partition_nets[mpi.rank]) {
			vector<int> overused_rr_node;
			assert(route_trees[net->local_id].root_rt_nodes.size() == 1);
			get_overused_nodes(route_trees[net->local_id], route_trees[net->local_id].root_rt_nodes[0], partitioner.orig_g, congestion, overused_rr_node);
			if (!overused_rr_node.empty()) {
				//zlog_level(delta_log, ROUTER_V1, "Net %d bb_rank %d has %lu overused nodes:\n", net->vpr_id, net->bb_area_rank, overused_rr_node.size());
				for (const auto &item : overused_rr_node) {
					zlog_level(delta_log, ROUTER_V1, "%d ", item);
				}
				zlog_level(delta_log, ROUTER_V1, "\n");
				//overused_total_bb_rank += net->bb_area_rank;
				++num_congested_nets;
			}
		}

        zlog_level(delta_log, ROUTER_V1, "Average overused net bb rank: %g\n", (float)overused_total_bb_rank/num_congested_nets);
        zlog_level(delta_log, ROUTER_V1, "Num congested nets: %d\n", num_congested_nets);

        iter_start = clock::now();

        auto analyze_timing_start = clock::now();

		int *recvcounts = nullptr;
		int *displs = nullptr;

		init_displ(partition_nets, &recvcounts, &displs);

		sync_net_delay(partition_nets, mpi.rank, mpi.comm_size, recvcounts, displs, current_level, mpi.comm, net_timing);

		int num_crits = 0;
		for (int i = 0; i < mpi.comm_size; ++i) {
			num_crits += recvcounts[i];
		}

		float *all_crits;
		int idx;

		if (mpi.rank == 0) {
			crit_path_delay = analyze_timing(net_timing);

			all_crits = new float[num_crits];

			idx = 0;
			for (int pi = 0; pi < mpi.comm_size; ++pi) {
				idx = 0;
				for (auto &net : partition_nets[pi]) {
					zlog_level(delta_log, ROUTER_V3, "Net %d send crits:\n", net->vpr_id);
					for (int k = 1; k <= net->sinks.size(); ++k) {
						all_crits[displs[pi] + idx] = net_timing[net->vpr_id].timing_criticality[k];
						zlog_level(delta_log, ROUTER_V3, "\t%g\n", all_crits[displs[pi] + idx]);
						++idx;
					}
				}
				assert(idx == recvcounts[pi]);
			}
		} else {
			all_crits = nullptr;
		}

		float *crits = new float[recvcounts[mpi.rank]];

		MPI_Scatterv(all_crits, recvcounts, displs, MPI_FLOAT, crits, recvcounts[mpi.rank], MPI_FLOAT, 0, mpi.comm);

		idx = 0;
		for (auto &net : partition_nets[mpi.rank]) {
			zlog_level(delta_log, ROUTER_V3, "Net %d recv crits:\n", net->vpr_id);

			for (int s = 1; s <= net->sinks.size(); ++s) {
				net_timing[net->vpr_id].timing_criticality[s] = crits[idx];
				zlog_level(delta_log, ROUTER_V3, "\t%g\n", net_timing[net->vpr_id].timing_criticality[s]);
				++idx;
			}
		}

		if (mpi.rank == 0) {
			delete [] all_crits;
		}
		delete [] crits;

        analyze_timing_time = clock::now()-analyze_timing_start;

		auto update_cost_start = clock::now();

		if (iter == 0) {
			pres_fac = opts->initial_pres_fac;
			//update_costs_mpi(partitioner.orig_g[0], partitioner.result_pid_by_level[current_level], mpi.rank, congestion, win, pres_fac, 0);
			update_costs(partitioner.orig_g, congestion, pres_fac, 0);
		} else {
			pres_fac *= opts->pres_fac_mult;

			/* Avoid overflow for high iteration counts, even if acc_cost is big */
			pres_fac = std::min(pres_fac, static_cast<float>(HUGE_POSITIVE_FLOAT / 1e5));

			//update_costs_mpi(partitioner.orig_g[0], partitioner.result_pid_by_level[current_level], mpi.rank, congestion, win, pres_fac, opts->acc_fac);
			update_costs(partitioner.orig_g, congestion, pres_fac, opts->acc_fac);
		}

		update_cost_time = clock::now()-update_cost_start;

		total_actual_route_time += actual_route_time;
        total_route_time += route_time;
		total_wait_time += wait_time;
        total_analyze_timing_time += analyze_timing_time;
        total_update_cost_time += update_cost_time;

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

			static const char *name_type[] = { "SOURCE", "SINK", "IPIN", "OPIN",
				"CHANX", "CHANY", "INTRA_CLUSTER_EDGE" };
            zlog_info(delta_log, "Num overused nodes: %lu/%d (%.2f)\n", num_overused_nodes, num_vertices(partitioner.orig_g), num_overused_nodes*100.0/num_vertices(partitioner.orig_g));
			for (int i = 0; i < overused_nodes_by_type.size(); ++i) {
				zlog_info(delta_log, "\t%s: %d (%g)\n", name_type[i], overused_nodes_by_type[i], overused_nodes_by_type[i]*100.0/num_overused_nodes);
			}

			int *all_overused_nodes_by_type = new int[mpi.comm_size*NUM_RR_TYPES];
			int *overused_nodes_by_type_send = new int[NUM_RR_TYPES];
			for (int i = 0; i < overused_nodes_by_type.size(); ++i) {
				overused_nodes_by_type_send[i] = overused_nodes_by_type[i];
			}

			MPI_Gather(overused_nodes_by_type_send, NUM_RR_TYPES, MPI_INT, all_overused_nodes_by_type, NUM_RR_TYPES, MPI_INT, 0, mpi.comm);

			unsigned long *all_num_overused_nodes = new unsigned long[mpi.comm_size];
			MPI_Gather(&num_overused_nodes, 1, MPI_UNSIGNED_LONG, all_num_overused_nodes, 1, MPI_UNSIGNED_LONG, 0, mpi.comm);

			if (mpi.rank == 0) {
				printf("Num overused nodes: ");
				for (int i = 0; i < mpi.comm_size; ++i) {
					printf("%lu/%d (%.2f) ", all_num_overused_nodes[i], num_vertices(partitioner.orig_g), all_num_overused_nodes[i]*100.0/num_vertices(partitioner.orig_g));
				}
				printf("\n");

				for (int i = 0; i < overused_nodes_by_type.size(); ++i) {
					printf("\t%s: ", name_type[i]);
					for (int j = 0; j < mpi.comm_size; ++j) {
						printf("%d (%g) ", all_overused_nodes_by_type[j*NUM_RR_TYPES+i], all_overused_nodes_by_type[j*NUM_RR_TYPES+i]*100.0/all_num_overused_nodes[j]);
					}
					printf("\n");
				}
			} 

			delete [] all_overused_nodes_by_type;
			delete [] overused_nodes_by_type_send;
			delete [] all_num_overused_nodes;

			int not_decreasing = (num_overused_nodes >= prev_num_overused_nodes && iter > 10) ? 1 : 0;
			//int not_decreasing = 1;
			//int not_decreasing = current_level+1 < partitioner.result_pid_by_level.size(); [> testing <]
			//int not_decreasing = current_level+1 <= std::log2(initial_comm_size); [> testing <]
			//int reduced_not_decreasing;
			//MPI_Allreduce(&not_decreasing, &reduced_not_decreasing, 1, MPI_INT, MPI_LOR, mpi.comm);
			//int reduced_not_decreasing = iter > 1;

			prev_num_overused_nodes = num_overused_nodes;

			//zlog_level(delta_log, ROUTER_V1, "not_decreasing: %d reduced_not_decreasing: %d\n", not_decreasing, reduced_not_decreasing);

			if (not_decreasing && initial_comm_size > 1 && current_level < (int)std::log2(initial_comm_size)) {
				/* need to send route tree over */
				for (int i = 0; i < mpi.comm_size; ++i) {
					zlog_level(delta_log, ROUTER_V2, "Freeing comm for %d\n", i);
					int error = MPI_Comm_free(&mpi.ibcast_comm[i]);
					assert(error == MPI_SUCCESS);
				}

				//vector<set<net_t *>> old_partitions(mpi.comm_size);
				//for (int i = 0; i < partitions.size(); ++i) {
					//for (const auto &net : partitions[i]) {
						//old_partitions[i].insert(net);
					//}
				//}
				//
				auto combine_start = clock::now();

				vector<int> displs(mpi.comm_size);
				vector<int> recvcounts(mpi.comm_size);

				assert(partition_nets.size() == mpi.comm_size);

				for (int i = 0; i < mpi.comm_size; ++i) {
					recvcounts[i] = partition_nets[i].size();
				}

				int displ = 0;
				for (int i = 0; i < mpi.comm_size; ++i) {
					displs[i] = displ;
					displ += recvcounts[i];
				}
				
				vector<int> num_congested_nodes(partition_nets[mpi.rank].size());
				for (int i = 0; i < partition_nets[mpi.rank].size(); ++i) {
					net_t *net = partition_nets[mpi.rank][i];

					num_congested_nodes[i] = route_tree_mark_congested_nodes_to_be_ripped(route_trees[net->local_id], partitioner.orig_g, congestion);

					if (num_congested_nodes[i] > 0) {
						zlog_level(delta_log, ROUTER_V3, "Net %d is congested\n", net->local_id);
					}
				}

				vector<int> all_num_congested_nodes(nets.size());

				int error = MPI_Allgatherv(num_congested_nodes.data(), recvcounts[mpi.rank], MPI_INT, all_num_congested_nodes.data(), recvcounts.data(), displs.data(), MPI_INT, mpi.comm);
				assert(error == MPI_SUCCESS);

				int cur_num_non_congested_nets = non_congested_nets_ptr.size();
				nets_ptr.clear();
				for (int i = 0; i < mpi.comm_size; ++i) {
					for (int j = 0; j < partition_nets[i].size(); ++j) {
						assert(j < recvcounts[i]);
						if (all_num_congested_nodes[displs[i]+j] > 0) {
							nets_ptr.push_back(partition_nets[i][j]);

							zlog_level(delta_log, ROUTER_V3, "Recv net %d is congested\n", partition_nets[i][j]->local_id);
						} else {
							assert(all_num_congested_nodes[displs[i]+j] == 0);
							assert(find(begin(non_congested_nets_ptr), end(non_congested_nets_ptr), partition_nets[i][j]) == end(non_congested_nets_ptr));

							non_congested_nets_ptr.push_back(partition_nets[i][j]);
						}
					}
				}

				for (int i = cur_num_non_congested_nets; i < non_congested_nets_ptr.size(); ++i) {
					move_route_tree(non_congested_nets_ptr[i], net_partition_id[non_congested_nets_ptr[i]->local_id], 0, partitioner.orig_g, route_trees, net_route_trees, net_timing, mpi);
				}

				repartition(nets_ptr, opts->load_balanced, mpi.comm_size/2, current_level+1, initial_comm_size, net_partition_id, partition_nets, partitioner.orig_g, route_trees, net_route_trees, net_timing, mpi);

				assert(mpi.comm_size % 2 == 0);
				mpi.comm_size /= 2;

				++current_level;

				idling = mpi.rank >= mpi.comm_size;

				MPI_Comm new_comm;
				MPI_Comm_split(mpi.comm, idling ? 1 : 0, mpi.rank, &new_comm);
				mpi.comm = new_comm;
				MPI_Comm_rank(mpi.comm, &mpi.rank);
				int new_comm_size;
				MPI_Comm_size(mpi.comm, &new_comm_size);
				assert(mpi.comm_size == new_comm_size);

				combine_time = clock::now()-combine_start;

				if (!idling) {
					//assert(current_level < partitioner.result_pid_by_level.size());
					printf("[%d] Transitioned to level %d at iteration %d\n", mpi.rank, current_level, iter);
					zlog_level(delta_log, ROUTER_V1, "Transitioned to level %d at iteration %d\n", current_level, iter);
					zlog_level(delta_log, ROUTER_V1, "New pid %d for initial pid %d\n", mpi.rank, initial_rank);

					prev_num_overused_nodes = std::numeric_limits<unsigned long>::max();

					for (int i = 0; i < num_vertices(partitioner.orig_g); ++i) {
						congestion[i].recalc_occ = 0; 
					}

					if (mpi.rank == 0) {
						for (const auto &net : non_congested_nets_ptr) {
							zlog_level(delta_log, ROUTER_V3, "Checking non-congested net vpr id %d\n", net->vpr_id);

							check_route_tree(route_trees[net->local_id], *net, partitioner.orig_g);
							recalculate_occ(route_trees[net->local_id], partitioner.orig_g, congestion);
						}
					}

					for (const auto &net : partition_nets[mpi.rank]) {
						zlog_level(delta_log, ROUTER_V3, "Checking net vpr id %d\n", net->vpr_id);

						check_route_tree(route_trees[net->local_id], *net, partitioner.orig_g);
						recalculate_occ(route_trees[net->local_id], partitioner.orig_g, congestion);
					}

					sync_recalc_occ(congestion, num_vertices(partitioner.orig_g),  mpi.rank, mpi.comm_size, mpi.comm);

					bool valid = true;
					for (int i = 0; i < num_vertices(partitioner.orig_g); ++i) {
						sprintf_rr_node_impl(i, buffer);
						if (congestion[i].recalc_occ != congestion[i].occ) {
							printf("Node %s occ mismatch, recalc: %d original: %d\n", buffer, congestion[i].recalc_occ, congestion[i].occ);
							valid = false;
						}
					}
					assert(valid);

					for (int i = 0; i < mpi.comm_size; ++i) {
						zlog_level(delta_log, ROUTER_V2, "Duplicating comm for %d\n", i);
						int error = MPI_Comm_dup(mpi.comm, &mpi.ibcast_comm[i]);
						assert(error == MPI_SUCCESS);
					}

					for (int i = mpi.comm_size; i < mpi.comm_size*2; ++i) {
						zlog_level(delta_log, ROUTER_V2, "Freeing req and data %d\n", i);
						assert(mpi.pending_send_req[i] == MPI_REQUEST_NULL);
						mpi.free_send_req_index.push(i);
						assert(mpi.pending_send_data_ref_count[i] == 0);
						mpi.free_send_data_index.push(i);
					}
				}
			} 

			//if ((float)total_num_interpartition_sinks/total_num_sinks_to_route > 0.1) {
				//++current_level;
				//assert(current_level < partitioner.result_pid_by_level.size());
			//}
        } /* feasible_routing */

        iter_time += clock::now()-iter_start;

        total_iter_time += iter_time;
		total_combine_time += combine_time;

		if (prev_rank == 0) {
			printf("Iteration time: %g s.\n", duration_cast<nanoseconds>(iter_time).count() / 1e9);

			vector<float> all_route_time(prev_comm_size, 0);

			printf("\tRoute time: %g ", duration_cast<nanoseconds>(route_time).count() / 1e9);
			for (int i = 1; i < prev_comm_size; ++i) {
				float f_route_time;
				MPI_Recv(&f_route_time, 1, MPI_FLOAT, i, i, prev_comm, MPI_STATUS_IGNORE);
				//unsigned long u_route_time;
				//MPI_Recv(&u_route_time, 1, MPI_UNSIGNED_LONG, i, i, prev_comm, MPI_STATUS_IGNORE);
				printf("%g ", f_route_time);
				all_route_time[i] = f_route_time;
			}
			printf("\n");

			printf("\t\tActual route time: %g (%g) ", duration_cast<nanoseconds>(actual_route_time).count() / 1e9, duration_cast<nanoseconds>(actual_route_time).count() * 100.0 / duration_cast<nanoseconds>(route_time).count());
			for (int i = 1; i < prev_comm_size; ++i) {
				float f_actual_route_time;
				MPI_Recv(&f_actual_route_time, 1, MPI_FLOAT, i, i, prev_comm, MPI_STATUS_IGNORE);
				//unsigned long u_actual_route_time;
				//MPI_Recv(&u_actual_route_time, 1, MPI_UNSIGNED_LONG, i, i, prev_comm, MPI_STATUS_IGNORE);
				printf("%g (%g) ", f_actual_route_time, f_actual_route_time * 100 / all_route_time[i]);
			}
			printf("\n");

			float f_sync_time = duration_cast<nanoseconds>(mpi_perf.total_sync_time).count() / 1e9;
			printf("\t\t\tSync time: %g (%g) ", f_sync_time, f_sync_time * 100.0 / (duration_cast<nanoseconds>(route_time).count() / 1e9));
			for (int i = 1; i < prev_comm_size; ++i) {
				MPI_Recv(&f_sync_time, 1, MPI_FLOAT, i, i, prev_comm, MPI_STATUS_IGNORE);
				printf("%g (%g) ", f_sync_time, f_sync_time * 100 / all_route_time[i]);
			}
			printf("\n");

			//float f_recv_time = duration_cast<nanoseconds>(mpi_perf.total_recv_time).count() / 1e9;
			//printf("\t\t\t\tRecv time: %g (%g) ", f_recv_time, f_recv_time * 100.0 / (duration_cast<nanoseconds>(route_time).count() / 1e9));
			//for (int i = 1; i < prev_comm_size; ++i) {
				//MPI_Recv(&f_recv_time, 1, MPI_FLOAT, i, i, prev_comm, MPI_STATUS_IGNORE);
				//printf("%g (%g) ", f_recv_time, f_recv_time * 100 / all_route_time[i]);
			//}
			//printf("\n");

			//float f_iprobe_time = duration_cast<nanoseconds>(mpi_perf.total_iprobe_time).count() / 1e9;
			//printf("\t\t\t\tIprobe time: %g (%g) ", f_iprobe_time, f_iprobe_time * 100.0 / (duration_cast<nanoseconds>(route_time).count() / 1e9));
			//for (int i = 1; i < prev_comm_size; ++i) {
				//MPI_Recv(&f_iprobe_time, 1, MPI_FLOAT, i, i, prev_comm, MPI_STATUS_IGNORE);
				//printf("%g (%g) ", f_iprobe_time, f_iprobe_time * 100 / all_route_time[i]);
			//}
			//printf("\n");

			//float f_update_time = duration_cast<nanoseconds>(mpi_perf.total_update_time).count() / 1e9;
			//printf("\t\t\t\tUpdate time: %g (%g) ", f_update_time, f_update_time * 100.0 / (duration_cast<nanoseconds>(route_time).count() / 1e9));
			//for (int i = 1; i < prev_comm_size; ++i) {
				//MPI_Recv(&f_update_time, 1, MPI_FLOAT, i, i, prev_comm, MPI_STATUS_IGNORE);
				//printf("%g (%g) ", f_update_time, f_update_time * 100 / all_route_time[i]);
			//}
			//printf("\n");

			float f_broadcast_time = duration_cast<nanoseconds>(mpi_perf.total_broadcast_time).count() / 1e9;
			printf("\t\t\tBroadcast time: %g (%g) ", f_broadcast_time, f_broadcast_time * 100.0 / (duration_cast<nanoseconds>(route_time).count() / 1e9));
			for (int i = 1; i < prev_comm_size; ++i) {
				MPI_Recv(&f_broadcast_time, 1, MPI_FLOAT, i, i, prev_comm, MPI_STATUS_IGNORE);
				printf("%g (%g) ", f_broadcast_time, f_broadcast_time * 100 / all_route_time[i]);
			}
			printf("\n");

			float f_testsome_time = duration_cast<nanoseconds>(mpi_perf.total_send_testsome_time).count() / 1e9;
			printf("\t\t\tSend testsome time: %g (%g) ", f_testsome_time, f_testsome_time * 100.0 / (duration_cast<nanoseconds>(route_time).count() / 1e9));
			for (int i = 1; i < prev_comm_size; ++i) {
				MPI_Recv(&f_testsome_time, 1, MPI_FLOAT, i, i, prev_comm, MPI_STATUS_IGNORE);
				printf("%g (%g) ", f_testsome_time, f_testsome_time * 100 / all_route_time[i]);
			}
			printf("\n");

			float f_last_progress_time = duration_cast<nanoseconds>(last_progress_time).count() / 1e9;
			printf("\t\tLast progress time: %g (%g) ", f_last_progress_time, f_last_progress_time * 100.0 / (duration_cast<nanoseconds>(route_time).count() / 1e9));
			for (int i = 1; i < prev_comm_size; ++i) {
				MPI_Recv(&f_last_progress_time, 1, MPI_FLOAT, i, i, prev_comm, MPI_STATUS_IGNORE);
				printf("%g (%g) ", f_last_progress_time, f_last_progress_time * 100 / all_route_time[i]);
			}
			printf("\n");

			float f_last_sync_time = duration_cast<nanoseconds>(last_sync_time).count() / 1e9;
			printf("\t\tLast sync time: %g (%g) ", f_last_sync_time, f_last_sync_time * 100.0 / (duration_cast<nanoseconds>(route_time).count() / 1e9));
			for (int i = 1; i < prev_comm_size; ++i) {
				MPI_Recv(&f_last_sync_time, 1, MPI_FLOAT, i, i, prev_comm, MPI_STATUS_IGNORE);
				printf("%g (%g) ", f_last_sync_time, f_last_sync_time * 100 / all_route_time[i]);
			}
			printf("\n");

			printf("num last syncs: %d ", num_last_syncs);
			for (int i = 1; i < prev_comm_size; ++i) {
				int tmp_num_last_syncs;
				MPI_Recv(&tmp_num_last_syncs, 1, MPI_INT, i, i, prev_comm, MPI_STATUS_IGNORE);
				printf("%d ", tmp_num_last_syncs);
			}
			printf("\n");

			printf("\t\tWait time: %g (%g) ", duration_cast<nanoseconds>(wait_time).count() / 1e9, duration_cast<nanoseconds>(wait_time).count() * 100.0 / duration_cast<nanoseconds>(route_time).count());
			for (int i = 1; i < prev_comm_size; ++i) {
				float f_wait_time;
				MPI_Recv(&f_wait_time, 1, MPI_FLOAT, i, i, prev_comm, MPI_STATUS_IGNORE);
				printf("%g (%g) ", f_wait_time, f_wait_time * 100 / all_route_time[i]);
			}
			printf("\n");

			printf("\tUpdate cost time: %g s.\n", duration_cast<nanoseconds>(update_cost_time).count() / 1e9);
			printf("\tAnalyze timing time: %g s.\n", duration_cast<nanoseconds>(analyze_timing_time).count() / 1e9);
			printf("\tCombine time: %g s.\n", duration_cast<nanoseconds>(combine_time).count() / 1e9);
			printf("Critical path: %g ns\n", crit_path_delay);
		} else {
			float f_route_time = duration_cast<nanoseconds>(route_time).count() / 1e9;
			MPI_Send(&f_route_time, 1, MPI_FLOAT, 0, prev_rank, prev_comm);

			//unsigned long u_route_time = route_time.count();
			//MPI_Send(&u_route_time, 1, MPI_UNSIGNED_LONG, 0, prev_rank, prev_comm);

			float f_actual_route_time = duration_cast<nanoseconds>(actual_route_time).count() / 1e9;
			MPI_Send(&f_actual_route_time, 1, MPI_FLOAT, 0, prev_rank, prev_comm);

			//unsigned long u_actual_route_time = actual_route_time.count();
			//MPI_Send(&u_actual_route_time, 1, MPI_UNSIGNED_LONG, 0, prev_rank, prev_comm);

			float f_sync_time = duration_cast<nanoseconds>(mpi_perf.total_sync_time).count() / 1e9;
			MPI_Send(&f_sync_time, 1, MPI_FLOAT, 0, prev_rank, prev_comm);

			//float f_recv_time = duration_cast<nanoseconds>(mpi_perf.total_recv_time).count() / 1e9;
			//MPI_Send(&f_recv_time, 1, MPI_FLOAT, 0, prev_rank, prev_comm);

			//float f_iprobe_time = duration_cast<nanoseconds>(mpi_perf.total_iprobe_time).count() / 1e9;
			//MPI_Send(&f_iprobe_time, 1, MPI_FLOAT, 0, prev_rank, prev_comm);

			//float f_update_time = duration_cast<nanoseconds>(mpi_perf.total_update_time).count() / 1e9;
			//MPI_Send(&f_update_time, 1, MPI_FLOAT, 0, prev_rank, prev_comm);

			float f_broadcast_time = duration_cast<nanoseconds>(mpi_perf.total_broadcast_time).count() / 1e9;
			MPI_Send(&f_broadcast_time, 1, MPI_FLOAT, 0, prev_rank, prev_comm);

			float f_testsome_time = duration_cast<nanoseconds>(mpi_perf.total_send_testsome_time).count() / 1e9;
			MPI_Send(&f_testsome_time, 1, MPI_FLOAT, 0, prev_rank, prev_comm);

			float f_last_progress_time = duration_cast<nanoseconds>(last_progress_time).count() / 1e9;
			MPI_Send(&f_last_progress_time, 1, MPI_FLOAT, 0, prev_rank, prev_comm);

			float f_last_sync_time = duration_cast<nanoseconds>(last_sync_time).count() / 1e9;
			MPI_Send(&f_last_sync_time, 1, MPI_FLOAT, 0, prev_rank, prev_comm);

			MPI_Send(&num_last_syncs, 1, MPI_INT, 0, prev_rank, prev_comm);

			float f_wait_time = duration_cast<nanoseconds>(wait_time).count() / 1e9;
			MPI_Send(&f_wait_time, 1, MPI_FLOAT, 0, prev_rank, prev_comm);
		}

		MPI_Barrier(prev_comm);

		prev_comm = mpi.comm;
		prev_comm_size = mpi.comm_size;
		prev_rank = mpi.rank;

        //clock::time_point greedy_earliest_end_time = *std::min_element(begin(greedy_end_time), begin(greedy_end_time)+mpi.comm_size);

        //printf("greedy wait time: ");
        //for (int i = 0; i < mpi.comm_size; ++i) {
			//printf("%g (%g), ", duration_cast<nanoseconds>(greedy_end_time[i]-greedy_earliest_end_time).count() / 1e9, 100.0*(greedy_end_time[i]-greedy_earliest_end_time)/greedy_route_time);
        //}
		//printf("\n");
    }
	
	MPI_Barrier(MPI_COMM_WORLD);
    //sprintf(buffer, "%s_run_%d_rr_graph_occ.txt", s_circuit_name, run);
    //dump_rr_graph_occ(congestion, num_vertices(g), buffer);
	if (initial_rank == 0) {
		if (routed) {
			printf("Routed in %d iterations. Total iteration time: %g\n", iter, duration_cast<nanoseconds>(total_iter_time).count() / 1e9);
		} else {
			printf("Failed to route in %d iterations. Total iteration time: %g\n", opts->max_router_iterations, duration_cast<nanoseconds>(total_iter_time).count() / 1e9);
		}

		vector<float> all_total_route_time(initial_comm_size, 0);

		printf("\tTotal route time: %g ", duration_cast<nanoseconds>(total_route_time).count() / 1e9);
		for (int i = 1; i < initial_comm_size; ++i) {
			float f_total_route_time;
			MPI_Recv(&f_total_route_time, 1, MPI_FLOAT, i, i, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			all_total_route_time[i] = f_total_route_time;
			printf("%g ", f_total_route_time);
		}
		printf("\n");
		printf("\t\tTotal actual route time: %g (%g) ", duration_cast<nanoseconds>(total_actual_route_time).count() / 1e9, duration_cast<nanoseconds>(total_actual_route_time).count() * 100.0 / duration_cast<nanoseconds>(total_route_time).count());
		for (int i = 1; i < initial_comm_size; ++i) {
			float f_total_actual_route_time;
			MPI_Recv(&f_total_actual_route_time, 1, MPI_FLOAT, i, i, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			printf("%g (%g) ", f_total_actual_route_time, f_total_actual_route_time * 100.0 / all_total_route_time[i]);
		}
		printf("\n");
		printf("\t\t\tTotal sync time: %g (%g) ", duration_cast<nanoseconds>(total_sync_time).count() / 1e9, duration_cast<nanoseconds>(total_sync_time).count() * 100.0 / duration_cast<nanoseconds>(total_route_time).count());
		for (int i = 1; i < initial_comm_size; ++i) {
			float f_total_sync_time;
			MPI_Recv(&f_total_sync_time, 1, MPI_FLOAT, i, i, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			printf("%g (%g) ", f_total_sync_time, f_total_sync_time * 100.0 / all_total_route_time[i]);
		}
		printf("\n");
		//printf("\t\t\t\tTotal recv time: %g (%g) ", duration_cast<nanoseconds>(total_recv_time).count() / 1e9, duration_cast<nanoseconds>(total_recv_time).count() * 100.0 / duration_cast<nanoseconds>(total_route_time).count());
		//for (int i = 1; i < initial_comm_size; ++i) {
			//float f_total_recv_time;
			//MPI_Recv(&f_total_recv_time, 1, MPI_FLOAT, i, i, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			//printf("%g (%g) ", f_total_recv_time, f_total_recv_time * 100.0 / all_total_route_time[i]);
		//}
		//printf("\n");
		//printf("\t\t\t\tTotal probe time: %g (%g) ", duration_cast<nanoseconds>(total_iprobe_time).count() / 1e9, duration_cast<nanoseconds>(total_iprobe_time).count() * 100.0 / duration_cast<nanoseconds>(total_route_time).count());
		//for (int i = 1; i < initial_comm_size; ++i) {
			//float f_total_iprobe_time;
			//MPI_Recv(&f_total_iprobe_time, 1, MPI_FLOAT, i, i, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			//printf("%g (%g) ", f_total_iprobe_time, f_total_iprobe_time * 100.0 / all_total_route_time[i]);
		//}
		//printf("\n");
		//printf("\t\t\t\tTotal update time: %g (%g) ", duration_cast<nanoseconds>(total_update_time).count() / 1e9, duration_cast<nanoseconds>(total_update_time).count() * 100.0 / duration_cast<nanoseconds>(total_route_time).count());
		//for (int i = 1; i < initial_comm_size; ++i) {
			//float f_total_update_time;
			//MPI_Recv(&f_total_update_time, 1, MPI_FLOAT, i, i, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			//printf("%g (%g) ", f_total_update_time, f_total_update_time * 100.0 / all_total_route_time[i]);
		//}
		//printf("\n");
		printf("Total calls %d\n", total_calls);
		printf("\t\t\tTotal broadcast time: %g (%g) ", duration_cast<nanoseconds>(total_broadcast_time).count() / 1e9, duration_cast<nanoseconds>(total_broadcast_time).count() * 100.0 / duration_cast<nanoseconds>(total_route_time).count());
		for (int i = 1; i < initial_comm_size; ++i) {
			float f_total_broadcast_time;
			MPI_Recv(&f_total_broadcast_time, 1, MPI_FLOAT, i, i, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			printf("%g (%g) ", f_total_broadcast_time, f_total_broadcast_time * 100.0 / all_total_route_time[i]);
		}
		printf("\n");
		printf("\t\t\tTotal send testsome time: %g (%g) ", duration_cast<nanoseconds>(total_send_testsome_time).count() / 1e9, duration_cast<nanoseconds>(total_send_testsome_time).count() * 100.0 / duration_cast<nanoseconds>(total_route_time).count());
		for (int i = 1; i < initial_comm_size; ++i) {
			float f_total_send_testsome_time;
			MPI_Recv(&f_total_send_testsome_time, 1, MPI_FLOAT, i, i, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			printf("%g (%g) ", f_total_send_testsome_time, f_total_send_testsome_time * 100.0 / all_total_route_time[i]);
		}
		printf("\n");
		printf("\t\tTotal last progress time: %g (%g) ", duration_cast<nanoseconds>(total_last_progress_time).count() / 1e9, duration_cast<nanoseconds>(total_last_progress_time).count() * 100.0 / duration_cast<nanoseconds>(total_route_time).count());
		for (int i = 1; i < initial_comm_size; ++i) {
			float f_total_last_progress_time;
			MPI_Recv(&f_total_last_progress_time, 1, MPI_FLOAT, i, i, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			printf("%g (%g) ", f_total_last_progress_time, f_total_last_progress_time * 100.0 / all_total_route_time[i]);
		}
		printf("\n");
		printf("\t\tTotal last sync time: %g (%g) ", duration_cast<nanoseconds>(total_last_sync_time).count() / 1e9, duration_cast<nanoseconds>(total_last_sync_time).count() * 100.0 / duration_cast<nanoseconds>(total_route_time).count());
		for (int i = 1; i < initial_comm_size; ++i) {
			float f_total_last_sync_time;
			MPI_Recv(&f_total_last_sync_time, 1, MPI_FLOAT, i, i, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			printf("%g (%g) ", f_total_last_sync_time, f_total_last_sync_time * 100.0 / all_total_route_time[i]);
		}
		printf("\n");
		printf("\t\tTotal wait time: %g (%g) ", duration_cast<nanoseconds>(total_wait_time).count() / 1e9, duration_cast<nanoseconds>(total_wait_time).count() * 100.0 / duration_cast<nanoseconds>(total_route_time).count());
		for (int i = 1; i < initial_comm_size; ++i) {
			float f_total_wait_time;
			MPI_Recv(&f_total_wait_time, 1, MPI_FLOAT, i, i, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			printf("%g (%g) ", f_total_wait_time, f_total_wait_time * 100.0 / all_total_route_time[i]);
		}
		printf("\n");
		printf("\tTotal update cost time: %g s.\n", duration_cast<nanoseconds>(total_update_cost_time).count() / 1e9);
		printf("\tTotal analyze timing time: %g s.\n", duration_cast<nanoseconds>(total_analyze_timing_time).count() / 1e9);
		printf("\tTotal combine time: %g s.\n", duration_cast<nanoseconds>(total_combine_time).count() / 1e9);

		if (routed) {
			printf("Final critical path: %g ns\n", crit_path_delay);
		}
	} else {
		float f_total_route_time = duration_cast<nanoseconds>(total_route_time).count() / 1e9;
		MPI_Send(&f_total_route_time, 1, MPI_FLOAT, 0, initial_rank, MPI_COMM_WORLD);

		float f_total_actual_route_time = duration_cast<nanoseconds>(total_actual_route_time).count() / 1e9;
		MPI_Send(&f_total_actual_route_time, 1, MPI_FLOAT, 0, initial_rank, MPI_COMM_WORLD);

		float f_total_sync_time = duration_cast<nanoseconds>(total_sync_time).count() / 1e9;
		MPI_Send(&f_total_sync_time, 1, MPI_FLOAT, 0, initial_rank, MPI_COMM_WORLD);

		//float f_total_recv_time = duration_cast<nanoseconds>(total_recv_time).count() / 1e9;
		//MPI_Send(&f_total_recv_time, 1, MPI_FLOAT, 0, initial_rank, MPI_COMM_WORLD);

		//float f_total_iprobe_time = duration_cast<nanoseconds>(total_iprobe_time).count() / 1e9;
		//MPI_Send(&f_total_iprobe_time, 1, MPI_FLOAT, 0, initial_rank, MPI_COMM_WORLD);

		//float f_total_update_time = duration_cast<nanoseconds>(total_update_time).count() / 1e9;
		//MPI_Send(&f_total_update_time, 1, MPI_FLOAT, 0, initial_rank, MPI_COMM_WORLD);

		float f_total_broadcast_time = duration_cast<nanoseconds>(total_broadcast_time).count() / 1e9;
		MPI_Send(&f_total_broadcast_time, 1, MPI_FLOAT, 0, initial_rank, MPI_COMM_WORLD);

		float f_total_send_testsome_time = duration_cast<nanoseconds>(total_send_testsome_time).count() / 1e9;
		MPI_Send(&f_total_send_testsome_time, 1, MPI_FLOAT, 0, initial_rank, MPI_COMM_WORLD);

		float f_total_last_progress_time = duration_cast<nanoseconds>(total_last_progress_time).count() / 1e9;
		MPI_Send(&f_total_last_progress_time, 1, MPI_FLOAT, 0, initial_rank, MPI_COMM_WORLD);

		float f_total_last_sync_time = duration_cast<nanoseconds>(total_last_sync_time).count() / 1e9;
		MPI_Send(&f_total_last_sync_time, 1, MPI_FLOAT, 0, initial_rank, MPI_COMM_WORLD);

		float f_total_wait_time = duration_cast<nanoseconds>(total_wait_time).count() / 1e9;
		MPI_Send(&f_total_wait_time, 1, MPI_FLOAT, 0, initial_rank, MPI_COMM_WORLD);
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
