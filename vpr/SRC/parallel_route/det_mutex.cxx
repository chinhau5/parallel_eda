#include "det_mutex.h"
#include "log.h"

void det_mutex_init(det_mutex_t &m, exec_state_t *e_state, int num_threads)
{
	m.e_state = e_state;
	m.num_threads = num_threads;
}

void det_mutex_wait_for_turn(det_mutex_t &m, int tid)
{
	int cur = *m.e_state[tid].logical_clock;
	bool other_less;
	zlog_level(delta_log, ROUTER_V3, "Thread %d acquiring lock, cur clock %d\n", tid, cur);

	do {
		other_less = false;
		for (int i = 0; i < m.num_threads && !other_less; ++i) {
			if (i == tid) {
				continue;
			}
			other_less = *m.e_state[i].logical_clock < cur || (*m.e_state[i].logical_clock == cur && i <= tid);
		}
	} while (other_less);
}

void det_mutex_lock(det_mutex_t &m, int tid)
{
	det_mutex_wait_for_turn(m, tid);
	//m.lock.lock();
}

void det_mutex_unlock(det_mutex_t &m, int tid)
{
	++(*m.e_state[tid].logical_clock);
	zlog_level(delta_log, ROUTER_V3, "Thread %d released lock, new clock %d\n", tid, m.e_state[tid].logical_clock->load());
	//m.lock.unlock();
}
