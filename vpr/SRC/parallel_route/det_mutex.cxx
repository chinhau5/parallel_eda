#include "det_mutex.h"
#include "log.h"
#include <assert.h>

void det_mutex_init(det_mutex_t &m, exec_state_t *e_state, int num_threads)
{
	m.lock = new tbb::spin_mutex();
	m.e_state = e_state;
	m.num_threads = num_threads;
}

static void det_mutex_wait_for_turn(det_mutex_t &m, int tid)
{
	int cur = *m.e_state[tid].logical_clock;
	int cur_inet = *m.e_state[tid].inet;
	bool other_less;
	zlog_level(delta_log, ROUTER_V3, "Thread %d acquiring lock, cur clock %d\n", tid, cur);

	do {
		other_less = false;
		for (int i = 0; i < m.num_threads && !other_less; ++i) {
			if (i == tid) {
				continue;
			}
			int other = *m.e_state[i].logical_clock;
			int other_inet = *m.e_state[i].inet;
			//if (other_inet == cur_inet) {
				//printf("%d %d\n", other_inet, cur_inet);
				//assert(false);
			//}
			//assert(other_inet != cur_inet || (other_inet == cur_inet &&other_inet == -1));
			assert(other_inet != cur_inet);

			other_less = other < cur || (other == cur && other_inet < cur_inet);
		}
	} while (other_less);
	assert(m.lock->try_lock());
}

void det_mutex_lock(det_mutex_t &m, int tid)
{
	det_mutex_wait_for_turn(m, tid);
	//m.lock.lock();
}

void det_mutex_unlock(det_mutex_t &m, int tid)
{
	m.lock->unlock();
	++(*m.e_state[tid].logical_clock);
	zlog_level(delta_log, ROUTER_V3, "Thread %d released lock, new clock %d\n", tid, m.e_state[tid].logical_clock->load());
}
