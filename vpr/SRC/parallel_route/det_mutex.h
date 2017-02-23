#ifndef DET_MUTEX_H

#define DET_MUTEX_H

#include <atomic>
#include <tbb/spin_mutex.h>

struct exec_state_t {
	std::atomic<unsigned long> *logical_clock;
	std::atomic<int> *inet;
};

struct det_mutex_t {
	tbb::spin_mutex *lock;
	exec_state_t *e_state;
	int num_threads;
	//tbb::spin_mutex lock;
	//int last_tid;
	//vec_clock_t last_released;
};

void det_mutex_init(det_mutex_t &m, exec_state_t *e_state, int num_threads);

void det_mutex_lock(det_mutex_t &m, int tid);

void det_mutex_unlock(det_mutex_t &m, int tid);

#endif
