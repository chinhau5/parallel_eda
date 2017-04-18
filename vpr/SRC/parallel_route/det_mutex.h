#ifndef DET_MUTEX_H

#define DET_MUTEX_H

#include <atomic>
#include <chrono>
#include <vector>
#include <tbb/spin_mutex.h>
#include <x86intrin.h>
#include "log.h"
#include "clock.h"
#include "tp.h"

#define PROFILE_WAIT_LEVEL 2
#define SET_CONTEXT
//#define TRACE_WAIT
//#define PMC
#define TSC
#define INC_THRESHOLD 1000000

//using mtimer = std::chrono::high_resolution_clock;
using mtimer = myclock;

struct alignas(64) exec_state_t {
	//std::atomic<unsigned long> *logical_clock;
	//std::atomic<int> *inet;
	int tid;
	std::atomic<unsigned long> logical_clock;
	std::atomic<int> inet;
	int lock_count;
	int context;
	int region;
	int level;
	unsigned long prev_lock_clock;
	std::atomic<unsigned long> previous_logical_clock;
	long long total_counter;
};

struct copyable_state_t {
	int tid;
	unsigned long logical_clock;
	int inet;
	int lock_count;
	int context;
	int region;
	int level;
};

struct det_mutex_stats_t {
	int index;

	unsigned long num_incs;

	unsigned long released_logical_clock;

	mtimer::duration wait_time;
	int num_rounds;

	unsigned long logical_clock;
	int context;
	int inet;
	int level;

	//std::vector<std::pair<unsigned long, int>> max_clock_diff;
	std::vector<copyable_state_t> others;
	
	unsigned long max;
	int max_thread;
	int max_lock_count;

	unsigned long max_logical_clock;
	int max_context;
	int max_inet;
	int max_level;

	mtimer::time_point start;
	mtimer::duration length;
};

const std::vector<det_mutex_stats_t> &get_wait_stats(int level, int tid);
void init_wait_stats(int num_levels, int num_threads);
void clear_wait_stats();

#if defined(TSC)
using inc_time_t = unsigned long long;
#else
using inc_time_t = mtimer::time_point;
#endif

const std::vector<std::pair<inc_time_t, int>> &get_inc_time(int tid);
void add_inc_time(int tid, const inc_time_t &time, int context);
void init_inc_time(int num_threads);
void clear_inc_time();

void inline set_context(exec_state_t *state, int context)
{
	state->context = context;
}

int inline get_context(exec_state_t *state)
{
	return state->context;
}

#define mem_order std::memory_order_seq_cst

void inline inc_logical_clock_no_profile(exec_state_t *state, int count)
{
	//state->logical_clock->fetch_add(count, std::memory_order_relaxed);
	if (state->logical_clock.load() - state->prev_lock_clock < INC_THRESHOLD) {
#if defined(TRACE_WAIT)
		tracepoint(hello_world, my_first_tracepoint, state->tid, state->context, 0, state->lock_count);
#endif
	
		state->logical_clock.fetch_add(count, mem_order);
	}

	zlog_level(delta_log, ROUTER_V3, "Increasing logical clock to %lu from %d\n", state->logical_clock.load(), get_context(state));
}

void inline inc_logical_clock(exec_state_t *state, int count)
{
	inc_logical_clock_no_profile(state, count);

#if PROFILE_WAIT_LEVEL >= 3
	int context = get_context(state);
#if defined(TSC)
	//unsigned int aux = 0;
	//unsigned long long tsc = __rdtscp(&aux);
	//add_inc_time(tsc, state->tid, context);
	add_inc_time(state->tid, 0, context);
#else
	add_inc_time(state->tid, mtimer::now(), context);
#endif
#endif
}

void inline set_logical_clock(exec_state_t *state, int count)
{
	//state->logical_clock->fetch_add(count, std::memory_order_relaxed);
#if defined(TRACE_WAIT)
	tracepoint(hello_world, my_first_tracepoint, state->tid, state->context, 0, state->lock_count);
#endif

	state->logical_clock.store(count, mem_order);
}

unsigned long inline get_logical_clock(exec_state_t *state)
{
	return state->logical_clock.load(mem_order);
}

void inline set_inet(exec_state_t *state, int inet)
{
	state->inet.store(inet, mem_order);
}

int inline get_inet(exec_state_t *state)
{
	return state->inet.load(mem_order);
}

struct det_mutex_t {
	tbb::spin_mutex *lock;
	exec_state_t *e_state;
	int num_threads;
	unsigned long released_logical_clock;
	int released_thread;
	int released_no_wait;
	int no_wait;
	//tbb::spin_mutex lock;
	//int last_tid;
	//vec_clock_t last_released;
};

void det_mutex_reset(det_mutex_t &m);

void det_mutex_init(det_mutex_t &m, exec_state_t *e_state, int num_threads);

void det_mutex_destroy(det_mutex_t &m);

void det_mutex_lock(det_mutex_t &m, int tid, int num_nodes);
void det_mutex_lock_no_wait(det_mutex_t &m, int tid);

void det_mutex_unlock(det_mutex_t &m, int tid);
void det_mutex_unlock_no_wait(det_mutex_t &m, int tid);

#endif
