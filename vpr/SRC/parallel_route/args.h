#ifndef ARGS_H
#define ARGS_H

#include <semaphore.h>
#include <tbb/tbb.h>
#include "route.h"
#include "route_tree.h"

struct WorkerArgs {
	int tid;
	int num_threads;
	int *iteration;

	pthread_mutex_t *start_cond_lock;
	pthread_cond_t *start_cond;
	my_pthread_barrier_t *barrier;

	std::mutex &lock;
	QuadTree<virtual_net_t *> &in_flight_qt;
	/*vector<virtual_net_t> &virtual_nets;*/
	std::vector<std::vector<virtual_net_t>> *virtual_nets_by_net;
	int &num_routed_virtual_nets;

	std::vector<std::mutex> &net_locks;

	RRGraph &g;
	route_parameters_t &params;
	route_state_t *state;
	std::vector<route_tree_t> &route_trees;
	std::vector<trace_t> **prev_traces_ptr;
	std::vector<trace_t> **current_traces_ptr;
	t_net_timing *net_timing;

/*#if __linux__*/
	sem_t *consume_sem;
	sem_t *produce_sem;
/*#elif __APPLE__*/
	/*dispatched*/
/*#else*/
/*#error "Unsupported OS!"*/
/*#endif*/

	tbb::concurrent_bounded_queue<virtual_net_t *> *pending_virtual_nets;
	tbb::concurrent_bounded_queue<virtual_net_t *> *routed_virtual_nets;

	perf_t perf;

	WorkerArgs(
			int tid,
			int num_threads,

			pthread_mutex_t *start_cond_lock,
			pthread_cond_t *start_cond,
			my_pthread_barrier_t *barrier,

			std::mutex &lock,
			QuadTree<virtual_net_t *> &in_flight_qt,
			/*vector<virtual_net_t> &virtual_nets,*/
			int &num_routed_virtual_nets,

			std::vector<std::mutex> &net_locks,

			RRGraph &g,
			route_parameters_t &params,
			route_state_t *state,
			std::vector<route_tree_t> &route_trees,
			std::vector<trace_t> **prev_traces_ptr,
			std::vector<trace_t> **current_traces_ptr,
			t_net_timing *net_timing,

			sem_t *consume_sem,
			sem_t *produce_sem,

			tbb::concurrent_bounded_queue<virtual_net_t *> *pending_virtual_nets,
			tbb::concurrent_bounded_queue<virtual_net_t *> *routed_virtual_nets
				) :
					tid(tid),
					num_threads(num_threads),

					start_cond_lock(start_cond_lock),
					start_cond(start_cond),
					barrier(barrier),

					lock(lock),
					in_flight_qt(in_flight_qt),
					/*virtual_nets(virtual_nets),*/
					num_routed_virtual_nets(num_routed_virtual_nets),

					net_locks(net_locks),

					g(g),
					params(params),
					state(state),
					route_trees(route_trees),
					prev_traces_ptr(prev_traces_ptr),
					current_traces_ptr(current_traces_ptr),
					net_timing(net_timing),

					consume_sem(consume_sem),
					produce_sem(produce_sem),

					pending_virtual_nets(pending_virtual_nets),
					routed_virtual_nets(routed_virtual_nets)
					{
					}
};


struct SchedulerArgs {
	int *iteration;
	std::vector<virtual_net_t *> current_virtual_nets;
	std::vector<virtual_net_t> *virtual_nets;
	std::vector<net_t> *nets;
	std::vector<std::vector<virtual_net_t>> *virtual_nets_by_net;

	sem_t *produce_sem;	
	sem_t *consume_sem;	
	tbb::concurrent_bounded_queue<virtual_net_t *> *routed_virtual_nets;
	tbb::concurrent_bounded_queue<virtual_net_t *> *pending_virtual_nets;
	QuadTree<virtual_net_t *> *in_flight_qt;
	RRGraph *g;
	route_parameters_t *params;
	std::vector<route_tree_t> *route_trees;
	int num_nets;
	t_router_opts *opts;

	sched_perf_t perf;

	std::vector<std::mutex> *net_locks;
};

#endif
