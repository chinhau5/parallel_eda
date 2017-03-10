#ifndef DET_MUTEX_H

#define DET_MUTEX_H

#include <atomic>
#include <chrono>
#include <vector>
#include <tbb/spin_mutex.h>
#include "clock.h"

//using mtimer = std::chrono::high_resolution_clock;
using mtimer = myclock;

struct exec_state_t {
	//std::atomic<unsigned long> *logical_clock;
	//std::atomic<int> *inet;
	int tid;
	alignas(64) std::atomic<unsigned long> logical_clock;
	alignas(64) std::atomic<int> inet;
	alignas(64) int context;
};

struct det_mutex_stats_t {
	std::vector<std::pair<unsigned long, int>> max_clock_diff;
	mtimer::duration wait_time;
};

const std::vector<det_mutex_stats_t> &get_mutex_stats(int tid);
void init_wait_stats(int num_threads);
void clear_wait_stats();

const std::vector<std::pair<mtimer::time_point, int>> &get_inc_time(int tid);
void add_inc_time(const mtimer::time_point &time, int tid, int context);
void init_inc_time(int num_threads);
void clear_inc_time();

void inline inc_logical_clock(exec_state_t *state, int count, int context)
{
	//state->logical_clock->fetch_add(count, std::memory_order_relaxed);
	state->logical_clock.fetch_add(count, std::memory_order_relaxed);
	add_inc_time(mtimer::now(), state->tid, context);
	state->context = context;
}

void inline set_logical_clock(exec_state_t *state, int count)
{
	//state->logical_clock->fetch_add(count, std::memory_order_relaxed);
	state->logical_clock.store(count, std::memory_order_relaxed);
}

unsigned long inline get_logical_clock(exec_state_t *state)
{
	return state->logical_clock.load(std::memory_order_relaxed);
}

void inline set_inet(exec_state_t *state, int inet)
{
	state->inet.store(inet, std::memory_order_relaxed);
}

int inline get_inet(exec_state_t *state)
{
	return state->inet.load(std::memory_order_relaxed);
}

struct det_mutex_t {
	tbb::spin_mutex *lock;
	exec_state_t *e_state;
	int num_threads;
	unsigned long released_logical_clock;
	//tbb::spin_mutex lock;
	//int last_tid;
	//vec_clock_t last_released;
};

void det_mutex_init(det_mutex_t &m, exec_state_t *e_state, int num_threads);

void det_mutex_destroy(det_mutex_t &m);

void det_mutex_lock(det_mutex_t &m, int tid);

void det_mutex_unlock(det_mutex_t &m, int tid);

#endif
