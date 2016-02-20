#include "pch.h"
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

int zlog_custom_output(zlog_msg_t *msg);
int zlog_sched_custom_output(zlog_msg_t *msg);

void test_polygon()
{
	polygon p{
		{ {0, 0}, {0, 4}, {4, 4}, {4, 2}, {2, 2}, {2, 0}, {0,0} }
	};

	polygon p2{
		{ {0, 0}, {0, 4}, {4, 4}, {4, 2}, {2, 2}, {2, 0}, {0,0} }
	};
	assert(bg::intersects(p, p2));
}

static void *worker_thread_4(void *args)
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


		net_t *net = current_virtual_net->net;

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

		route_tree_rip_up_marked(wargs->route_trees[net->local_id], wargs->g, nullptr, wargs->params.pres_fac, false, nullptr);

		wargs->perf.total_rip_up_time += clock::now()-rip_up_start;

		auto route_start = clock::now();

		//route_net_2(wargs->g, net->vpr_id, current_virtual_net->source, current_virtual_net->current_sinks, wargs->params, wargs->state, wargs->route_trees[net->local_id], wargs->net_timing[net->vpr_id], false, &wargs->perf, nullptr);

		/*for (const auto &sink : current_virtual_net->sinks) {*/
			/*net->sink_routed[sink->id] = true;*/
		/*}*/
		assert(!current_virtual_net->routed);
		current_virtual_net->routed = true;

		wargs->perf.total_route_time += clock::now()-route_start;

		if (locked) {
			wargs->net_locks[net->local_id].unlock();
		}

		auto push_start = clock::now();

		wargs->routed_virtual_nets->push(current_virtual_net);
		zlog_level(delta_log, ROUTER_V3, "Done routing net %d virtual_net %d in %lld\n", net->vpr_id, current_virtual_net->id,
				std::chrono::duration_cast<std::chrono::nanoseconds>(clock::now()-route_start).count());

		wargs->perf.total_push_time += clock::now()-push_start;

		/*assert(!sem_post(wargs->produce_sem));*/
		/*zlog_info(scheduler_log, "Done routing in %lld ns\n", std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::high_resolution_clock::now()-route_start_time).count());*/
	}
}

//void new_worker()
//{
//}

//class foo {
	//public:
		//void operator() const
		//{
			//const vector<virtual_net_t *> &virtual_nets;
			//virtual	
		//}
//}
//

void wait_for_virtual_net_and_update_scheduler_state(tbb::concurrent_bounded_queue<virtual_net_t *> &routed_virtual_nets, vector<virtual_net_t *> &in_flight_virtual_nets, sched_perf_t *perf)
{
	using clock = std::chrono::high_resolution_clock;

	zlog_level(scheduler_log, ROUTER_V3, "Going to wait for net to finish routing\n");
	/* we need to wait for existing in flight nets to complete */
	auto wait_start = clock::now();

	virtual_net_t *routed_virtual_net;
	routed_virtual_nets.pop(routed_virtual_net);

	if (perf) {
		perf->total_wait_time += clock::now()-wait_start;
	}

	auto update_start = clock::now();

	net_t *net = routed_virtual_net->net;

	//bool locked;
	//if ((*sargs->net_locks)[net->local_id].try_lock()) {
		//locked = true;
	//} else {
		//assert(false);
	//}

	assert(routed_virtual_net->routed);

	zlog_level(scheduler_log, ROUTER_V3, "Removing routed net %d virtual net %d bb %s\n", routed_virtual_net->net->vpr_id, routed_virtual_net->id, static_cast<const std::ostringstream &>(std::ostringstream() << bg::dsv(routed_virtual_net->saved_scheduler_bounding_box)).str().c_str());

	/*assert(pending_rtree.remove(make_pair(routed_virtual_net->saved_scheduler_bounding_box, routed_virtual_net)));*/

	auto iter = find(begin(in_flight_virtual_nets), end(in_flight_virtual_nets), routed_virtual_net);
	assert(iter != end(in_flight_virtual_nets));
	in_flight_virtual_nets.erase(iter);
	/*std::remove(begin(in_flight_virtual_nets), end(in_flight_virtual_nets), routed_virtual_net);*/

	//if (locked) {
		//(*sargs->net_locks)[net->local_id].unlock();
	//}
	/* do we need to remark all the unrouted virtual net's route tree again? */
	/* need to remove the unrouted virtual nets from the pending rtree and reinsert because the bounding boxes have been update */
	/*for (const auto &net : *sargs->nets) {*/
	/*if (!net.ripped_up) {*/
	/*route_tree_mark_nodes_to_be_ripped((*sargs->route_trees)[net.local_id], *sargs->g, 3);*/
	/*update_sinks(net->virtual_nets, (*sargs->route_trees)[net->local_id], *sargs->g, sargs->params->astar_fac);*/
	/*}*/
	/*}*/
}

int dispatch_virtual_nets(vector<virtual_net_t *> &virtual_nets, vector<virtual_net_t *> &in_flight_virtual_nets, tbb::concurrent_bounded_queue<virtual_net_t *> &pending_virtual_nets, sched_perf_t *perf)
{
	using clock = std::chrono::high_resolution_clock;

	extern unsigned long num_leaf_node_pred_calls;
	extern unsigned long num_internal_node_pred_calls;
	/*zlog_level(scheduler_log, ROUTER_V3, "Going to dispatching virtual net. Pending rtree size: %lu\n", pending_rtree.size());*/

#ifdef __linux__
	__itt_frame_begin_v3(pD, NULL);
#endif
	auto dispatch_start = clock::now();
	num_leaf_node_pred_calls = 0;
	num_internal_node_pred_calls = 0;
	/*vector<box> boxes;*/
	int num_dispatched_nets = 0;

	/* try to parallelize this part */
	/* divide and conquer */
	int type = 1;
	if (type == 0) {
		tbb::combinable<vector<virtual_net_t *>> local_dispatched;

		vector<virtual_net_t *> unrouted_virtual_nets;
		unrouted_virtual_nets.reserve(virtual_nets.size());

		for (int i = 0; i < virtual_nets.size(); ++i) {
			if (!virtual_nets[i]->dispatched) {
				unrouted_virtual_nets.push_back(virtual_nets[i]);
			}
		}

		tbb::parallel_for(tbb::blocked_range<size_t>(0, unrouted_virtual_nets.size(), 1024),
				[&local_dispatched, &unrouted_virtual_nets, &in_flight_virtual_nets] (const tbb::blocked_range<size_t> &range) -> void {
				for (int i = range.begin(); i != range.end(); ++i) {
				bool overlap_in_flight = std::any_of(begin(in_flight_virtual_nets), end(in_flight_virtual_nets), [&unrouted_virtual_nets, &i] (virtual_net_t *in_flight) -> bool {
						return bg::intersects(unrouted_virtual_nets[i]->scheduler_bounding_box, in_flight->scheduler_bounding_box);
						});

				bool overlap_dispatched = std::any_of(begin(local_dispatched.local()), end(local_dispatched.local()), [&unrouted_virtual_nets, &i] (virtual_net_t *local_dispatched) -> bool {
						return bg::intersects(unrouted_virtual_nets[i]->scheduler_bounding_box, local_dispatched->scheduler_bounding_box);
						});

				if (!overlap_in_flight && !overlap_dispatched) {
				local_dispatched.local().push_back(unrouted_virtual_nets[i]);
				}
				}
				});

		/*const auto &dispatched = local_dispatched.combine([] (const vector<virtual_net_t *> &a, const vector<virtual_net_t *> &b) -> vector<virtual_net_t *> {*/
		/*if (a.size() > b.size()) {*/
		/*return a;*/
		/*} else {*/
		/*return b;*/
		/*}*/
		/*});*/

		vector<virtual_net_t *> dispatched;
		local_dispatched.combine_each([&dispatched] (const vector<virtual_net_t *> &local_dispatched) -> void {
				for (const auto &a : local_dispatched) {
				bool overlap = std::any_of(begin(dispatched), end(dispatched), [&a] (virtual_net_t *d) -> bool {
						return bg::intersects(d->scheduler_bounding_box, a->scheduler_bounding_box);
						});
				if (!overlap) {
				dispatched.push_back(a);
				}
				}
				}
				);

		for (const auto &d : dispatched) {
			pending_virtual_nets.push(d);

			in_flight_virtual_nets.push_back(d);

			d->dispatched = true;
		}

		num_dispatched_nets = dispatched.size();
	} else if (type == 1) {
		for (int i = 0; i < virtual_nets.size(); ++i) {
			bool overlap_in_flight = std::any_of(begin(in_flight_virtual_nets), end(in_flight_virtual_nets), [&virtual_nets, &i] (virtual_net_t *in_flight) -> bool {
					return !bg::disjoint(virtual_nets[i]->scheduler_bounding_box, in_flight->scheduler_bounding_box);
					});

			if (!virtual_nets[i]->dispatched && !overlap_in_flight) {
				pending_virtual_nets.push(virtual_nets[i]);

				in_flight_virtual_nets.push_back(virtual_nets[i]);

				virtual_nets[i]->dispatched = true;

				++num_dispatched_nets;
			}
		}
	}

	/*for (auto iter = pending_rtree.qbegin(bgi::disjoint(&in_flight_virtual_nets)); iter != pending_rtree.qend(); ++iter) {*/
	/*[> remove the iter from pending_rtree <]*/
	/*[> add the iter to in_flight_virtual_nets <]*/
	/*[> loop again <]*/
	/*[>auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::high_resolution_clock::now() - dispatch_start);<]*/
	/*[>zlog_info(scheduler_log, "Dispatch took %lld ns\n", elapsed.count());<]*/

	/*#ifdef __linux__*/
	/*[>__itt_frame_begin_v3(dispatch_domain, NULL);<]*/
	/*#endif*/
	/*bool not_in_flight = std::find_if(begin(in_flight_virtual_nets), end(in_flight_virtual_nets), [&iter] (virtual_net_t *in_flight) -> bool { return in_flight->sinks[0]->net == iter->second->sinks[0]->net; }) == end(in_flight_virtual_nets);*/

	/*if (not_in_flight) {*/
	/*sargs->pending_virtual_nets->push(iter->second);*/

	/*[>assert(!sem_post(sargs->consume_sem));<]*/

	/*in_flight_virtual_nets.push_back(iter->second);*/

	/*[>printf("Dispatching net %d virtual net %d\n", iter->second->sinks[0]->net->vpr_id, iter->second->id);<]*/
	/*[>printf("Num internal %d Num leaf %d\n", num_internal_node_pred_calls, num_leaf_node_pred_calls);<]*/
	/*++num_dispatched_nets;*/

	/*[>dispatched_virtual_nets.push_back(iter->second);<]*/

	/*[>iter = pending_rtree.qbegin(bgi::disjoint(&in_flight_virtual_nets));<]*/
	/*}*/

	/*[>zlog_level(scheduler_log, ROUTER_V3, "Dispatching net %d virtual net %d\n", iter->second->sinks[0]->net->vpr_id, iter->second->id);<]*/

	/*#ifdef __linux__*/
	/*[>__itt_frame_end_v3(dispatch_domain, NULL);<]*/
	/*#endif*/

	/*[>auto before_post = std::chrono::high_resolution_clock::now();<]*/
	/*[>auto post_time = std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::high_resolution_clock::now()-before_post);<]*/

	/*[>zlog_info(scheduler_log, "Posted at %lld ns in %lld\n", std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::high_resolution_clock::now()-start_time).count(), post_time.count());<]*/

	/*[>dispatch_start = std::chrono::high_resolution_clock::now();<]*/
	/*}*/
	/*printf("Num virtual nets: %d Num pred calls: %d\n", num_virtual_nets, num_pred_calls);*/
	if (perf) {
		perf->total_dispatch_time += std::chrono::high_resolution_clock::now()-dispatch_start;
		perf->num_leaf_node_pred_calls += num_leaf_node_pred_calls;
		perf->num_internal_node_pred_calls += num_internal_node_pred_calls;
	}

#ifdef __linux__
	__itt_frame_end_v3(pD, NULL);
#endif

	/*for (const auto &virtual_net : dispatched_virtual_nets) {*/
	/*zlog_level(scheduler_log, ROUTER_V3, "Removing dispatched net %d virtual net %d bb %s\n", virtual_net->sinks[0]->net->vpr_id, virtual_net->id, static_cast<const std::ostringstream &>(std::ostringstream() << bg::dsv(virtual_net->saved_scheduler_bounding_box)).str().c_str());*/
	/*assert(pending_rtree.remove(make_pair(virtual_net->saved_scheduler_bounding_box, virtual_net)));*/
	/*}*/

	/* verification, to remove later */
	//for (int i = 0; i < in_flight_virtual_nets.size(); ++i) {
	//for (int j = i+1; j < in_flight_virtual_nets.size(); ++j) {
	//assert(!bg::intersects(in_flight_virtual_nets[i]->scheduler_bounding_box, in_flight_virtual_nets[j]->scheduler_bounding_box));
	//}
	//}

#ifdef __linux__
	/*__itt_frame_begin_v3(update_domain, NULL);*/
#endif

	//for (int i = 0; i < num_dispatched_nets; ++i) {
	//}
	//
	return num_dispatched_nets;
}

static void *scheduler_thread_4(SchedulerArgs *sargs, WorkerArgs **wargs, const route_parameters_t &params, RRGraph &g, route_state_t *state, vector<route_tree_t> &route_trees, t_net_timing *net_timing)
{
	char filename[256];
	sprintf(filename, "/Volumes/DATA/scheduler_%d.log", *sargs->iteration);

	extern FILE *sched_output_log;
	extern unsigned long num_leaf_node_pred_calls;
	extern unsigned long num_internal_node_pred_calls;

	using clock = std::chrono::high_resolution_clock;

	sched_output_log = fopen(filename, "w");

	zlog_level(scheduler_log, ROUTER_V1, "Scheduler thread iter %d\n", *sargs->iteration);

	sargs->perf.num_updates = 0;
	sargs->perf.num_leaf_node_pred_calls = 0;
	sargs->perf.num_internal_node_pred_calls = 0;
	sargs->perf.total_build_time = clock::duration::zero();
	sargs->perf.total_dispatch_time = clock::duration::zero();
	sargs->perf.total_rtree_update_time = clock::duration::zero();
	sargs->perf.total_wait_time = clock::duration::zero();

	int total_num_virtual_nets = 0;
	for (auto &net : *sargs->nets) {
		route_tree_mark_nodes_to_be_ripped((*sargs->route_trees)[net.local_id], *sargs->g, nullptr, 10000);
		for (const auto &virtual_net : net.virtual_nets) {
			update_virtual_net_current_sinks(*virtual_net, (*sargs->route_trees)[net.local_id]);
			virtual_net->routed = virtual_net->current_sinks.empty();
			if (!virtual_net->routed) {
				++total_num_virtual_nets;
				virtual_net->dispatched = false;
			}
		}
		std::sort(begin(net.virtual_nets), end(net.virtual_nets), [] (const virtual_net_t *a, const virtual_net_t *b) -> bool {
				return bg::distance(a->centroid, point(a->source->x, a->source->y)) < bg::distance(b->centroid, point(b->source->x, b->source->y));
				});
	}

	tbb::atomic<int> total_num_routed_virtual_nets = 0;
	int round = 0;
	while (total_num_routed_virtual_nets < total_num_virtual_nets) {
		sched_perf_t sched_perf;
		sched_perf.total_build_time = clock::duration::zero();
		sched_perf.total_dispatch_time = clock::duration::zero();
		sched_perf.total_rtree_update_time = clock::duration::zero();
		sched_perf.total_wait_time = clock::duration::zero();

		auto round_start = clock::now();

		auto build_start = clock::now();

		tbb::spin_mutex bulk_lock;
		tbb::atomic<unsigned int> total_bb_area = 0;
		vector<virtual_net_t *> unrouted_virtual_nets;
		unrouted_virtual_nets.reserve(sargs->nets->size());
		for (const auto &net : *sargs->nets) {
			if (round < net.virtual_nets.size() && !net.virtual_nets[round]->routed) {
				unrouted_virtual_nets.push_back(net.virtual_nets[round]);
			} else {
				zlog_level(scheduler_log, ROUTER_V3, "NOT adding net %d virtual net %d because it is routed\n", net.virtual_nets[round]->net->vpr_id, net.virtual_nets[round]->id);
			}
		}

		/* update bounding boxes */
		tbb::parallel_for(tbb::blocked_range<size_t>(0, unrouted_virtual_nets.size(), 1024), [&sargs, &unrouted_virtual_nets, &round, &total_bb_area, &bulk_lock] (const tbb::blocked_range<size_t> &range) -> void {
				for (int i = range.begin(); i != range.end(); ++i) {
					virtual_net_t *virtual_net = unrouted_virtual_nets[i];
					const net_t &net = *virtual_net->net;

					/* the sequence of this 2 calls are important */
					update_virtual_net_bounding_box(*virtual_net, (*sargs->route_trees)[net.local_id], *sargs->g, sargs->params->astar_fac, nullptr);
					update_virtual_net_scheduler_bounding_box(*virtual_net, round == 0 ? (*sargs->route_trees)[net.local_id].scheduler_bounding_box : bg::make_inverse<box>());

					zlog_level(delta_log, ROUTER_V3, "\n");

					zlog_level(scheduler_log, ROUTER_V3, "Adding net %d virtual net %d bb %s\n", virtual_net->net->vpr_id, virtual_net->id, static_cast<const std::ostringstream &>(std::ostringstream() << bg::dsv(virtual_net->scheduler_bounding_box)).str().c_str());

					total_bb_area += bg::area(virtual_net->scheduler_bounding_box);

					virtual_net->saved_scheduler_bounding_box = virtual_net->scheduler_bounding_box;
					/*pending_rtree.insert(make_pair(virtual_net->scheduler_bounding_box, virtual_net));*/
				}
		});

		vector<virtual_net_t *> overlapping_nets;
		vector<pair<box, virtual_net_t *>> nets_to_partition;
		for (const auto &virtual_net : unrouted_virtual_nets) {
			extern int nx;
			extern int ny;
			if (bg::area(virtual_net->scheduler_bounding_box) > 0.05*nx*ny) {
				overlapping_nets.push_back(virtual_net);
			} else {
				nets_to_partition.push_back(make_pair(virtual_net->scheduler_bounding_box, virtual_net));
			}
		}

		//printf("area start:\n");
		//for (int i = 0; i < unrouted_virtual_nets.size(); ++i) {
			//extern int nx;
			//extern int ny;
			//printf("%f\n", bg::area(unrouted_virtual_nets[i]->scheduler_bounding_box)/((nx+1)*(ny+1))*100);
		//}
		//printf("area end:\n");

		sargs->perf.total_build_time = clock::now() - build_start;
		sched_perf.total_build_time = sargs->perf.total_build_time;

		vector<vector<int>> partitions;
		vector<vector<int>> overlaps;
		vector<bool> has_interpartition_overlap;

		auto partitioning_start = clock::now();

		assert(!nets_to_partition.empty());
		partition_nets(nets_to_partition, sargs->opts->num_threads, 1, overlaps, partitions, has_interpartition_overlap);
		//partition_nets_by_clustering(nets_to_partition, sargs->opts->num_threads, partitions, has_interpartition_overlap);

		sargs->perf.total_partitioning_time = clock::now() - partitioning_start;
		sched_perf.total_partitioning_time = sargs->perf.total_partitioning_time;

		printf("Num large nets: %lu\n", overlapping_nets.size());

		//FILE *rect_part = fopen("/Volumes/DATA/clusters/rect_part.txt", "w");
		//for (int i = 0; i < partitions.size(); ++i) {
			//for (int j = 0; j < partitions[i].size(); ++j) {
				//virtual_net_t *vnet = nets_to_partition[partitions[i][j]];
				//int x = vnet->scheduler_bounding_box.min_corner().get<0>();
				//int y = vnet->scheduler_bounding_box.min_corner().get<1>();
				//int width = vnet->scheduler_bounding_box.max_corner().get<0>()-vnet->scheduler_bounding_box.min_corner().get<0>();
				//int height = vnet->scheduler_bounding_box.max_corner().get<1>()-vnet->scheduler_bounding_box.min_corner().get<1>();
				//fprintf(rect_part, "%d %d %d %d %d\n", x, y, width, height, i); 
			//}
		//}
		//fclose(rect_part);

		tbb::atomic<bool> has_idle_threads = false;
		tbb::spin_mutex pending_lock;
		vector<virtual_net_t *> pending_virtual_nets;
		vector<virtual_net_t *> in_flight_virtual_nets;
		tbb::spin_mutex overlapping_lock;

		vector<perf_t> perf(partitions.size());
		for (int i = 0; i < partitions.size(); ++i) {
			perf[i].total_route_time = clock::duration::zero();
			perf[i].total_rip_up_time = clock::duration::zero();
		}

		tbb::parallel_for(tbb::blocked_range<size_t>(0, partitions.size(), 1),
				[&]
				(const tbb::blocked_range<size_t> &range) -> void {
				assert(range.end()-range.begin() == 1);

				int tid = range.begin();
				//printf("tid: %lu\n", tid);

				const vector<int> &virtual_net_indices = partitions[tid];
				for (const auto &virtual_net_index : virtual_net_indices) {
					if (has_interpartition_overlap[virtual_net_index]) { 
						overlapping_lock.lock();
						overlapping_nets.push_back(nets_to_partition[virtual_net_index].second);
						overlapping_lock.unlock();
					} else {
						printf("tid %d vnet %d independent\n", tid, virtual_net_index);
					/*TODO: implement work balancing later */
						//bool dispatched = false;
						//for (int i = 0; i < idling.size() && !dispatched; ++i) {
							//if (idling[i]) {
							//}
						//}
						virtual_net_t *virtual_net = nets_to_partition[virtual_net_index].second;
						const net_t *net = virtual_net->net;

						auto rip_up_start = clock::now();

						route_tree_rip_up_marked(route_trees[net->local_id], g, nullptr, params.pres_fac, false, nullptr);

						//wargs[tid]->perf.total_rip_up_time += clock::now()-rip_up_start;
						perf[tid].total_rip_up_time += clock::now()-rip_up_start;

						auto route_start = clock::now();

						//route_net_2(g, net->vpr_id, virtual_net->source, virtual_net->current_sinks, params, state, route_trees[net->local_id], net_timing[net->vpr_id], false, &(wargs[tid]->perf), nullptr);

						assert(!virtual_net->routed);
						virtual_net->routed = true;

						perf[tid].total_route_time += clock::now()-route_start;
						//wargs[tid]->perf.total_route_time += clock::now()-route_start;

						++total_num_routed_virtual_nets;
					}
				}

				/* steal job from other threads */
				//pending_lock.lock();
				//pending_lock.unlock();

				//route_tree_rip_up_marked(route_trees[net->local_id], g, params.pres_fac);

				//route_net_2(g, net->vpr_id, virtual_net->source, virtual_net->current_sinks, params, state, route_trees[net->local_id], net_timing[net->vpr_id], &perf);
			}
		);

		int total_num_dispatched_nets = 0;
		int num_dispatch_rounds = 0;
		while (total_num_dispatched_nets < overlapping_nets.size()) {
			total_num_dispatched_nets += dispatch_virtual_nets(overlapping_nets, in_flight_virtual_nets, *sargs->pending_virtual_nets, nullptr);
			++num_dispatch_rounds;
			do {
				wait_for_virtual_net_and_update_scheduler_state(*sargs->routed_virtual_nets, in_flight_virtual_nets, nullptr);
				++total_num_routed_virtual_nets;
			} while (sargs->routed_virtual_nets->size() > 0);
		}

		assert(sargs->pending_virtual_nets->empty());
		assert(sargs->routed_virtual_nets->empty());

		auto round_time = clock::now()-round_start;

		extern int nx;
		extern int ny;

		printf("Round %d Time: %g s Num virtual nets: %lu Num overlapping nets: %lu Num dispatch rounds: %d Average concurrency: %g Total bb area: %lu Average bb area: %g (%g)\n", round, std::chrono::duration_cast<std::chrono::nanoseconds>(round_time).count() / 1e9, unrouted_virtual_nets.size(), overlapping_nets.size(), num_dispatch_rounds, (float)overlapping_nets.size()/num_dispatch_rounds, total_bb_area, (float)total_bb_area/overlapping_nets.size(), (float)total_bb_area/overlapping_nets.size()/((nx+1)*(ny+1))*100);
		for (int i = 0; i < partitions.size(); ++i) {
			printf("Partition %d size %lu\n", i, partitions[i].size());
		}
		printf("Round %d time:\n", round);
		printf("Total build time: %g s\n", std::chrono::duration_cast<std::chrono::nanoseconds>(sched_perf.total_build_time).count() / 1e9);
		printf("Partitioning time: %g s\n", std::chrono::duration_cast<std::chrono::nanoseconds>(sched_perf.total_partitioning_time).count() / 1e9);
		for (int i = 0; i < sargs->opts->num_threads; ++i) {
			printf("Thread %d total rip up time: %g s\n", i, std::chrono::duration_cast<std::chrono::nanoseconds>(perf[i].total_rip_up_time).count() / 1e9);
			printf("Thread %d total route time: %g s\n", i, std::chrono::duration_cast<std::chrono::nanoseconds>(perf[i].total_route_time).count() / 1e9);
		}


#ifdef __linux__ 
		/*__itt_frame_end_v3(update_domain, NULL);*/
#endif
		++round;
	}

	if (total_num_routed_virtual_nets != total_num_virtual_nets) {
		zlog_error(delta_log, "Num routed (%d) != Total (%d)\n", total_num_routed_virtual_nets, total_num_virtual_nets);
		assert(false);
	}
	/*assert(sargs->perf.num_updates == num_virtual_nets);*/

	return nullptr;
}

static void *scheduler_thread_3(void *args)
{
	SchedulerArgs *sargs = (SchedulerArgs *)args;

	char filename[256];
	sprintf(filename, "/Volumes/DATA/scheduler_%d.log", *sargs->iteration);

	extern FILE *sched_output_log;
	extern unsigned long num_leaf_node_pred_calls;
	extern unsigned long num_internal_node_pred_calls;

	using clock = std::chrono::high_resolution_clock;

	sched_output_log = fopen(filename, "w");

	zlog_level(scheduler_log, ROUTER_V1, "Scheduler thread iter %d\n", *sargs->iteration);

	sargs->perf.num_updates = 0;
	sargs->perf.num_leaf_node_pred_calls = 0;
	sargs->perf.num_internal_node_pred_calls = 0;
	sargs->perf.total_build_time = clock::duration::zero();
	sargs->perf.total_dispatch_time = clock::duration::zero();
	sargs->perf.total_rtree_update_time = clock::duration::zero();
	sargs->perf.total_wait_time = clock::duration::zero();

	using value = pair<box, virtual_net_t *>;
	struct equal {
		bool operator()(const value &a, const value &b) const
		{
			return bg::equals(a.first, b.first) && a.second == b.second;
		}
	};

#ifdef __linux__
	/*__itt_frame_begin_v3(pD, NULL);*/
#endif

	int total_num_virtual_nets = 0;
	for (const auto &net : *sargs->nets) {
		route_tree_mark_nodes_to_be_ripped((*sargs->route_trees)[net.local_id], *sargs->g, nullptr, 10000);
		for (const auto &virtual_net : net.virtual_nets) {
			update_virtual_net_current_sinks(*virtual_net, (*sargs->route_trees)[net.local_id]);
			virtual_net->routed = virtual_net->current_sinks.empty();
			if (!virtual_net->routed) {
				++total_num_virtual_nets;
				virtual_net->dispatched = false;
			}
		}
	}

	int total_num_routed_virtual_nets = 0;
	int round = 0;
	while (total_num_routed_virtual_nets < total_num_virtual_nets) {
		auto round_start = clock::now();

		auto build_start = clock::now();

		vector<value> virtual_nets;
		tbb::spin_mutex bulk_lock;
		tbb::atomic<unsigned int> total_bb_area = 0;
		vector<virtual_net_t *> temp;
		temp.reserve(sargs->nets->size());
		for (const auto &net : *sargs->nets) {
			if (round < net.virtual_nets.size() && !net.virtual_nets[round]->routed) {
				temp.push_back(net.virtual_nets[round]);
			} else {
				zlog_level(scheduler_log, ROUTER_V3, "NOT adding net %d virtual net %d because it is routed\n", net.virtual_nets[round]->net->vpr_id, net.virtual_nets[round]->id);
			}
		}

		tbb::parallel_for(tbb::blocked_range<size_t>(0, temp.size(), 32), [&sargs, &temp, &round, &virtual_nets, &total_bb_area, &bulk_lock] (const tbb::blocked_range<size_t> &range) -> void {
				for (int i = range.begin(); i != range.end(); ++i) {
					virtual_net_t *virtual_net = temp[i];
					const net_t &net = *virtual_net->net;

					/* the sequence of this 2 calls are important */
					update_virtual_net_bounding_box(*virtual_net, (*sargs->route_trees)[net.local_id], *sargs->g, sargs->params->astar_fac, nullptr);
					update_virtual_net_scheduler_bounding_box(*virtual_net, round == 0 ? (*sargs->route_trees)[net.local_id].scheduler_bounding_box : bg::make_inverse<box>());

					zlog_level(delta_log, ROUTER_V3, "\n");

					zlog_level(scheduler_log, ROUTER_V3, "Adding net %d virtual net %d bb %s\n", virtual_net->net->vpr_id, virtual_net->id, static_cast<const std::ostringstream &>(std::ostringstream() << bg::dsv(virtual_net->scheduler_bounding_box)).str().c_str());

					total_bb_area += bg::area(virtual_net->scheduler_bounding_box);

					virtual_net->saved_scheduler_bounding_box = virtual_net->scheduler_bounding_box;
					/*pending_rtree.insert(make_pair(virtual_net->scheduler_bounding_box, virtual_net));*/
					bulk_lock.lock();
					virtual_nets.push_back(make_pair(virtual_net->scheduler_bounding_box, virtual_net));
					/*com.local().push_back(make_pair(virtual_net->scheduler_bounding_box, virtual_net));*/

					bulk_lock.unlock();
				}
		});


		std::sort(begin(virtual_nets), end(virtual_nets), [] (const value &a, const value &b) -> bool {
				return make_pair(a.second->centroid.get<0>(), a.second->centroid.get<1>()) <
				make_pair(b.second->centroid.get<0>(), b.second->centroid.get<1>());
				});
		/*vector<value> virtual_nets = com.combine([] (const vector<value> &a, const vector<value> &b) -> value { vector<value> con; con.insert(*/

		/*bg::index::rtree<value, bgi::rstar<64>, bgi::indexable<value>, equal> pending_rtree(virtual_nets);*/
		vector<virtual_net_t *> in_flight_virtual_nets;

		sargs->perf.total_build_time += clock::now()-build_start;

		int num_dispatch_rounds = 0;
		int num_routed_virtual_nets = 0;
		int num_virtual_nets = virtual_nets.size();

		/*vector<bool> routed(virtual_nets.size(), false);*/

		while (num_routed_virtual_nets < num_virtual_nets) {
			/*zlog_level(scheduler_log, ROUTER_V3, "Going to dispatching virtual net. Pending rtree size: %lu\n", pending_rtree.size());*/

#ifdef __linux__
	__itt_frame_begin_v3(pD, NULL);
#endif
			auto dispatch_start = clock::now();
			num_leaf_node_pred_calls = 0;
			num_internal_node_pred_calls = 0;
			/*vector<box> boxes;*/
			int num_dispatched_nets = 0;

			/* try to parallelize this part */
			/* divide and conquer */
			int type = 0;
			if (type == 0) {
				tbb::combinable<vector<virtual_net_t *>> local_dispatched;

				vector<virtual_net_t *> unrouted_virtual_nets;
				unrouted_virtual_nets.reserve(virtual_nets.size());

				for (int i = 0; i < virtual_nets.size(); ++i) {
					if (!virtual_nets[i].second->dispatched) {
						unrouted_virtual_nets.push_back(virtual_nets[i].second);
					}
				}

				tbb::parallel_for(tbb::blocked_range<size_t>(0, unrouted_virtual_nets.size(), 1024),
						[&local_dispatched, &unrouted_virtual_nets, &in_flight_virtual_nets] (const tbb::blocked_range<size_t> &range) -> void {
						for (int i = range.begin(); i != range.end(); ++i) {
						bool overlap_in_flight = std::any_of(begin(in_flight_virtual_nets), end(in_flight_virtual_nets), [&unrouted_virtual_nets, &i] (virtual_net_t *in_flight) -> bool {
								return bg::intersects(unrouted_virtual_nets[i]->scheduler_bounding_box, in_flight->scheduler_bounding_box);
								});

						bool overlap_dispatched = std::any_of(begin(local_dispatched.local()), end(local_dispatched.local()), [&unrouted_virtual_nets, &i] (virtual_net_t *local_dispatched) -> bool {
								return bg::intersects(unrouted_virtual_nets[i]->scheduler_bounding_box, local_dispatched->scheduler_bounding_box);
								});

						if (!overlap_in_flight && !overlap_dispatched) {
						local_dispatched.local().push_back(unrouted_virtual_nets[i]);
						}
						}
						});

				/*const auto &dispatched = local_dispatched.combine([] (const vector<virtual_net_t *> &a, const vector<virtual_net_t *> &b) -> vector<virtual_net_t *> {*/
				/*if (a.size() > b.size()) {*/
				/*return a;*/
				/*} else {*/
				/*return b;*/
				/*}*/
				/*});*/

				vector<virtual_net_t *> dispatched;
				local_dispatched.combine_each([&dispatched] (const vector<virtual_net_t *> &local_dispatched) -> void {
						for (const auto &a : local_dispatched) {
						bool overlap = std::any_of(begin(dispatched), end(dispatched), [&a] (virtual_net_t *d) -> bool {
								return bg::intersects(d->scheduler_bounding_box, a->scheduler_bounding_box);
								});
						if (!overlap) {
						dispatched.push_back(a);
						}
						}
						}
						);

				for (const auto &d : dispatched) {
					sargs->pending_virtual_nets->push(d);

					in_flight_virtual_nets.push_back(d);

					d->dispatched = true;
				}

				int num_dispatched_nets = dispatched.size();
			} else if (type == 1) {
				for (int i = 0; i < virtual_nets.size(); ++i) {
					bool overlap_in_flight = std::any_of(begin(in_flight_virtual_nets), end(in_flight_virtual_nets), [&virtual_nets, &i] (virtual_net_t *in_flight) -> bool {
							return !bg::disjoint(virtual_nets[i].first, in_flight->scheduler_bounding_box);
							});

					if (!virtual_nets[i].second->dispatched && !overlap_in_flight) {
						sargs->pending_virtual_nets->push(virtual_nets[i].second);

						in_flight_virtual_nets.push_back(virtual_nets[i].second);

						virtual_nets[i].second->dispatched = true;
					}
				}
				num_dispatched_nets = in_flight_virtual_nets.size();
			} else {
				int num_samples = sqrt(virtual_nets.size());
				for (int i = 0; i < num_samples; ++i) {
				}
			}

			/*for (auto iter = pending_rtree.qbegin(bgi::disjoint(&in_flight_virtual_nets)); iter != pending_rtree.qend(); ++iter) {*/
				/*[> remove the iter from pending_rtree <]*/
				/*[> add the iter to in_flight_virtual_nets <]*/
				/*[> loop again <]*/
				/*[>auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::high_resolution_clock::now() - dispatch_start);<]*/
				/*[>zlog_info(scheduler_log, "Dispatch took %lld ns\n", elapsed.count());<]*/

/*#ifdef __linux__*/
				/*[>__itt_frame_begin_v3(dispatch_domain, NULL);<]*/
/*#endif*/
				/*bool not_in_flight = std::find_if(begin(in_flight_virtual_nets), end(in_flight_virtual_nets), [&iter] (virtual_net_t *in_flight) -> bool { return in_flight->sinks[0]->net == iter->second->sinks[0]->net; }) == end(in_flight_virtual_nets);*/

				/*if (not_in_flight) {*/
					/*sargs->pending_virtual_nets->push(iter->second);*/

					/*[>assert(!sem_post(sargs->consume_sem));<]*/

					/*in_flight_virtual_nets.push_back(iter->second);*/

					/*[>printf("Dispatching net %d virtual net %d\n", iter->second->sinks[0]->net->vpr_id, iter->second->id);<]*/
					/*[>printf("Num internal %d Num leaf %d\n", num_internal_node_pred_calls, num_leaf_node_pred_calls);<]*/
					/*++num_dispatched_nets;*/

					/*[>dispatched_virtual_nets.push_back(iter->second);<]*/

					/*[>iter = pending_rtree.qbegin(bgi::disjoint(&in_flight_virtual_nets));<]*/
				/*}*/

				/*[>zlog_level(scheduler_log, ROUTER_V3, "Dispatching net %d virtual net %d\n", iter->second->sinks[0]->net->vpr_id, iter->second->id);<]*/

/*#ifdef __linux__*/
				/*[>__itt_frame_end_v3(dispatch_domain, NULL);<]*/
/*#endif*/

				/*[>auto before_post = std::chrono::high_resolution_clock::now();<]*/
				/*[>auto post_time = std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::high_resolution_clock::now()-before_post);<]*/

				/*[>zlog_info(scheduler_log, "Posted at %lld ns in %lld\n", std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::high_resolution_clock::now()-start_time).count(), post_time.count());<]*/

				/*[>dispatch_start = std::chrono::high_resolution_clock::now();<]*/
			/*}*/
			/*printf("Num virtual nets: %d Num pred calls: %d\n", num_virtual_nets, num_pred_calls);*/
			sargs->perf.total_dispatch_time += std::chrono::high_resolution_clock::now()-dispatch_start;
			sargs->perf.num_leaf_node_pred_calls += num_leaf_node_pred_calls;
			sargs->perf.num_internal_node_pred_calls += num_internal_node_pred_calls;

			++num_dispatch_rounds;
#ifdef __linux__
	__itt_frame_end_v3(pD, NULL);
#endif

			/*for (const auto &virtual_net : dispatched_virtual_nets) {*/
			/*zlog_level(scheduler_log, ROUTER_V3, "Removing dispatched net %d virtual net %d bb %s\n", virtual_net->sinks[0]->net->vpr_id, virtual_net->id, static_cast<const std::ostringstream &>(std::ostringstream() << bg::dsv(virtual_net->saved_scheduler_bounding_box)).str().c_str());*/
			/*assert(pending_rtree.remove(make_pair(virtual_net->saved_scheduler_bounding_box, virtual_net)));*/
			/*}*/

			/* verification, to remove later */
			//for (int i = 0; i < in_flight_virtual_nets.size(); ++i) {
				//for (int j = i+1; j < in_flight_virtual_nets.size(); ++j) {
					//assert(!bg::intersects(in_flight_virtual_nets[i]->scheduler_bounding_box, in_flight_virtual_nets[j]->scheduler_bounding_box));
				//}
			//}

#ifdef __linux__
			/*__itt_frame_begin_v3(update_domain, NULL);*/
#endif

	//for (int i = 0; i < num_dispatched_nets; ++i) {
			do {
				zlog_level(scheduler_log, ROUTER_V3, "Going to wait for net to finish routing\n");
				/* we need to wait for existing in flight nets to complete */
				auto wait_start = clock::now();

				virtual_net_t *routed_vnet;
				sargs->routed_virtual_nets->pop(routed_vnet);
				/*while (sem_wait(sargs->produce_sem) == -1 && errno == EINTR);*/

				sargs->perf.total_wait_time += clock::now()-wait_start;
				/* remove the completed net from in_flight_virtual_net */

				auto update_start = clock::now();

				/*assert(sargs->routed_virtual_nets->try_pop(routed_vnet));*/

				net_t *net = routed_vnet->net;

				bool locked;
				if ((*sargs->net_locks)[net->local_id].try_lock()) {
					locked = true;
				} else {
					assert(false);
				}

				assert(routed_vnet->routed);

				zlog_level(scheduler_log, ROUTER_V3, "Removing routed net %d virtual net %d bb %s\n", routed_vnet->net->vpr_id, routed_vnet->id, static_cast<const std::ostringstream &>(std::ostringstream() << bg::dsv(routed_vnet->saved_scheduler_bounding_box)).str().c_str());

				/*assert(pending_rtree.remove(make_pair(routed_vnet->saved_scheduler_bounding_box, routed_vnet)));*/

				auto iter = find(begin(in_flight_virtual_nets), end(in_flight_virtual_nets), routed_vnet);
				assert(iter != end(in_flight_virtual_nets));
				in_flight_virtual_nets.erase(iter);
				/*std::remove(begin(in_flight_virtual_nets), end(in_flight_virtual_nets), routed_vnet);*/

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

				++total_num_routed_virtual_nets;
				++num_routed_virtual_nets;

				sargs->perf.total_rtree_update_time += clock::now()-update_start;

				++sargs->perf.num_updates;
			} while (sargs->routed_virtual_nets->size() > 0);
	//}
		}

		assert(sargs->pending_virtual_nets->empty());
		assert(sargs->routed_virtual_nets->empty());

		auto round_time = clock::now()-round_start;

		printf("Round %d Time: %g s Num virtual nets: %lu Num dispatch rounds: %d Average concurrency: %g Total bb area: %lu\n", round, std::chrono::duration_cast<std::chrono::nanoseconds>(round_time).count() / 1e9, virtual_nets.size(), num_dispatch_rounds, (float)virtual_nets.size()/num_dispatch_rounds, total_bb_area);

#ifdef __linux__ 
		/*__itt_frame_end_v3(update_domain, NULL);*/
#endif
		++round;
	}

	if (total_num_routed_virtual_nets != total_num_virtual_nets) {
		zlog_error(delta_log, "Num routed (%d) != Total (%d)\n", total_num_routed_virtual_nets, total_num_virtual_nets);
		assert(false);
	}
	/*assert(sargs->perf.num_updates == num_virtual_nets);*/

	return nullptr;
}

void test_partition();

bool phased_greedy_route(t_router_opts *opts)
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

	test_polygon();
	test_partition();

	/*test_rtree();*/

	/*test_route_tree();*/

	/*test_geometry();*/

	/*test_clustering(opts->num_threads, opts->bb_factor);*/
	/*return;*/

	/*test_perf();*/
	/*exit(-1);*/

	/*test_quadrant();*/

	/*auto time = clock::now();*/

	/*printf("test\n");*/

	/*auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(clock::now()-time);*/
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
		pthread_create(&threads[i], NULL, worker_thread_4, args[i]);
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

	//for_all_vertices(g, [] (RRNode &v) -> void {
			//v.properties.acc_cost = 1;
			//v.properties.pres_cost = 1;
			//v.properties.occ = 0;
			//});

	char buffer[256];

	if (current_output_log) {
		fclose(current_output_log);
	}

	std::mt19937 mt(time(NULL));

	bool routed = false;
	clock::duration total_route_time = clock::duration::zero();

	for (int iter = 0; iter < opts->max_router_iterations && !routed; ++iter) {
		sprintf(buffer, "/Volumes/DATA/iter_%d.txt", iter);
		current_output_log = fopen(buffer, "w");

		zlog_info(delta_log, "Routing iteration: %d\n", iter);

		for (int i = 0; i < opts->num_threads; ++i) {
			args[i]->iteration = &iter;

			args[i]->perf.num_heap_pushes = 0;
			args[i]->perf.num_heap_pops = 0;
			args[i]->perf.num_neighbor_visits = 0;
			args[i]->perf.total_wait_time = clock::duration::zero();
			args[i]->perf.total_rip_up_time = clock::duration::zero();
			args[i]->perf.total_route_time = clock::duration::zero();
			args[i]->perf.total_update_time = clock::duration::zero();
			args[i]->perf.total_push_time = clock::duration::zero();
			args[i]->perf.total_centroid_time = clock::duration::zero();
			args[i]->perf.total_get_nearest_time = clock::duration::zero();
			args[i]->perf.total_verificaion_time = clock::duration::zero();
			args[i]->perf.total_expansion_time = clock::duration::zero();
			args[i]->perf.total_scheduler_box_time = clock::duration::zero();
		}
		sargs.iteration = &iter;

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

		auto route_start = clock::now();

		/*std::shuffle(begin(virtual_nets), end(virtual_nets), mt);*/
		
		scheduler_thread_4(&sargs, args, params, g, state, route_trees, net_timing);

		auto route_time = clock::now()-route_start;
		total_route_time += route_time;
		/*zlog_info(delta_log, "Scheduling took %s\n", format(elapsed).c_str());*/
		zlog_info(delta_log, "Routing took %g\n", std::chrono::duration_cast<std::chrono::nanoseconds>(route_time).count() / 1e9);
		printf("Routing took %g\n", std::chrono::duration_cast<std::chrono::nanoseconds>(route_time).count() / 1e9);
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
			get_overused_nodes(route_trees[net.local_id], get_vertex(route_trees[net.local_id].graph, route_trees[net.local_id].root_rt_node_id), g, nullptr, overused_rr_node);
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

		if (feasible_routing(g, nullptr)) {
			zlog_info(delta_log, "Routed in %d iterations. Total route time: %g\n", iter+1, std::chrono::duration_cast<std::chrono::nanoseconds>(total_route_time).count() / 1e9);
			printf("Routed in %d iterations. Total route time: %g\n", iter+1, std::chrono::duration_cast<std::chrono::nanoseconds>(total_route_time).count() / 1e9);
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

			/* swap the pointers */
			std::swap(current_traces_ptr, prev_traces_ptr);
			for (const auto &trace : *current_traces_ptr) {
				assert(trace_empty(trace));
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

bool phased_greedy_route_old(t_router_opts *opts)
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

	test_partition();

	/*test_rtree();*/

	/*test_route_tree();*/

	/*test_geometry();*/

	/*test_clustering(opts->num_threads, opts->bb_factor);*/
	/*return;*/

	/*test_perf();*/
	/*exit(-1);*/

	/*test_quadrant();*/

	/*auto time = clock::now();*/

	/*printf("test\n");*/

	/*auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(clock::now()-time);*/
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
		pthread_create(&threads[i], NULL, worker_thread_4, args[i]);
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

	//for_all_vertices(g, [] (RRNode &v) -> void {
			//v.properties.acc_cost = 1;
			//v.properties.pres_cost = 1;
			//v.properties.occ = 0;
			//});

	char buffer[256];

	if (current_output_log) {
		fclose(current_output_log);
	}

	std::mt19937 mt(time(NULL));

	bool routed = false;
	clock::duration total_route_time = clock::duration::zero();

	for (int iter = 0; iter < opts->max_router_iterations && !routed; ++iter) {
		sprintf(buffer, "/Volumes/DATA/iter_%d.txt", iter);
		current_output_log = fopen(buffer, "w");

		zlog_info(delta_log, "Routing iteration: %d\n", iter);

		for (int i = 0; i < opts->num_threads; ++i) {
			args[i]->iteration = &iter;

			args[i]->perf.num_heap_pushes = 0;
			args[i]->perf.num_heap_pops = 0;
			args[i]->perf.num_neighbor_visits = 0;
			args[i]->perf.total_wait_time = clock::duration::zero();
			args[i]->perf.total_rip_up_time = clock::duration::zero();
			args[i]->perf.total_route_time = clock::duration::zero();
			args[i]->perf.total_update_time = clock::duration::zero();
			args[i]->perf.total_push_time = clock::duration::zero();
			args[i]->perf.total_centroid_time = clock::duration::zero();
			args[i]->perf.total_get_nearest_time = clock::duration::zero();
			args[i]->perf.total_verificaion_time = clock::duration::zero();
			args[i]->perf.total_expansion_time = clock::duration::zero();
			args[i]->perf.total_scheduler_box_time = clock::duration::zero();
		}
		sargs.iteration = &iter;

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

		auto route_start = clock::now();

		/*std::shuffle(begin(virtual_nets), end(virtual_nets), mt);*/
		scheduler_thread_3(&sargs);

		auto route_time = clock::now()-route_start;
		total_route_time += route_time;
		/*zlog_info(delta_log, "Scheduling took %s\n", format(elapsed).c_str());*/
		zlog_info(delta_log, "Routing took %g\n", std::chrono::duration_cast<std::chrono::nanoseconds>(route_time).count() / 1e9);
		printf("Routing took %g\n", std::chrono::duration_cast<std::chrono::nanoseconds>(route_time).count() / 1e9);
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
			get_overused_nodes(route_trees[net.local_id], get_vertex(route_trees[net.local_id].graph, route_trees[net.local_id].root_rt_node_id), g, nullptr, overused_rr_node);
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

		if (feasible_routing(g, nullptr)) {
			zlog_info(delta_log, "Routed in %d iterations. Total route time: %g\n", iter+1, std::chrono::duration_cast<std::chrono::nanoseconds>(total_route_time).count() / 1e9);
			printf("Routed in %d iterations. Total route time: %g\n", iter+1, std::chrono::duration_cast<std::chrono::nanoseconds>(total_route_time).count() / 1e9);
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

			/* swap the pointers */
			std::swap(current_traces_ptr, prev_traces_ptr);
			for (const auto &trace : *current_traces_ptr) {
				assert(trace_empty(trace));
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
