#include "det_mutex.h"
#include "log.h"
#include <assert.h>
#include <tbb/cache_aligned_allocator.h>

void start_logical_clock();
void stop_logical_clock();

void start_logical_clock(int es);
void stop_logical_clock(int es);

static std::vector<std::vector<det_mutex_stats_t>, tbb::cache_aligned_allocator<std::vector<det_mutex_stats_t>>> stats;

static std::vector<std::vector<std::pair<inc_time_t, int>>, tbb::cache_aligned_allocator<std::vector<std::pair<inc_time_t, int>>>> inc_time;

const std::vector<det_mutex_stats_t> &get_mutex_stats(int tid)
{
	return stats[tid];
}

void init_wait_stats(int num_threads)
{
	stats.resize(num_threads);
	for (auto &s : stats) {
		s.reserve(65536);
	}
}

void clear_wait_stats()
{
	for (auto &s : stats) {
		s.clear();
	}
}

const std::vector<std::pair<inc_time_t, int>> &get_inc_time(int tid)
{
	return inc_time[tid];
}

void add_inc_time(int tid, const inc_time_t &time, int context)
{
	inc_time[tid].emplace_back(std::make_pair(time, context));
}

void init_inc_time(int num_threads)
{
	inc_time.resize(num_threads);
	for (auto &i : inc_time) {
		i.reserve(4194304);
	}
}

void clear_inc_time()
{
	for (auto &i : inc_time) {
		i.clear();
	}
}

void det_mutex_init(det_mutex_t &m, exec_state_t *e_state, int num_threads)
{
	m.lock = new tbb::spin_mutex();
	m.e_state = e_state;
	m.num_threads = num_threads;
	m.released_logical_clock = 0;
	m.no_wait = -1;
}

void det_mutex_destroy(det_mutex_t &m)
{
	delete m.lock;
}

static void det_mutex_wait_for_turn_fixed(det_mutex_t &m, int tid, det_mutex_stats_t &stat)
{

}

static void det_mutex_wait_for_turn(det_mutex_t &m, int tid)
{
	unsigned long cur = get_logical_clock(&m.e_state[tid]);
	int cur_inet = get_inet(&m.e_state[tid]);
	bool other_less;

	do {
		other_less = false;
		for (int i = 0; i < m.num_threads && !other_less; ++i) {
			if (i == tid) {
				continue;
			}
			unsigned long other = get_logical_clock(&m.e_state[i]);
			int other_inet = get_inet(&m.e_state[i]);
			//if (other_inet == cur_inet) {
				//printf("%d %d\n", other_inet, cur_inet);
				//assert(false);
			//}
			//assert(other_inet != cur_inet || (other_inet == cur_inet &&other_inet == -1));
			assert(other_inet != cur_inet);

			other_less = other < cur || (other == cur && other_inet < cur_inet);
		}
	} while (other_less);
}

void det_mutex_lock(det_mutex_t &m, int tid)
{
#if PROFILE_WAIT_LEVEL >= 1
	stats[tid].emplace_back();
	auto &stat = stats[tid].back();
	stat.index = m.e_state[tid].lock_count;
	++(m.e_state[tid].lock_count);
	stat.inet = m.e_state[tid].inet;
	stat.previous_context = m.e_state[tid].context;
#endif

#if PROFILE_WAIT_LEVEL >= 2
	unsigned long cur = get_logical_clock(&m.e_state[tid]);

	unsigned long min_clock = cur;
	int min_thread = tid;
	int min_inet = m.e_state[tid].inet;
	int min_context = m.e_state[tid].context;
	int min_lock_count = m.e_state[tid].lock_count;

	for (int i = 0; i < m.num_threads; ++i) {
		if (i == tid) {
			continue;
		}
		unsigned long other = get_logical_clock(&m.e_state[i]);
		if (other < cur) {
			if (other < min_clock) {
				min_clock = other;
				min_thread = i;
				min_inet = m.e_state[i].inet;
				min_context = m.e_state[i].context;
				min_lock_count = m.e_state[i].lock_count;
			}
		}
	}

	//if (cur < min_clock) {
		//printf("%lu %lu\n", cur, min_clock);
		//assert(false);
	//}
	assert(cur >= min_clock);
	//assert(min_thread != tid);

	//stat.max_clock_diff.push_back(std::make_pair(cur-min_clock, min_context));
	stat.logical_clock = cur;
	stat.max = cur-min_clock;
	stat.max_logical_clock = min_clock;
	stat.max_thread = min_thread;
	stat.max_inet = min_inet;
	stat.max_context = min_context;
	stat.max_lock_count = min_lock_count;
#endif

#if PROFILE_WAIT_LEVEL >= 1
	auto wait_start = mtimer::now();
#endif

#if defined(INSTRUMENT)
	stop_logical_clock();
#elif defined(PMC)
	stop_logical_clock(tid);
#endif

	zlog_level(delta_log, ROUTER_V3, "[%d] Thread %d acquiring lock, region %d released %lu cur clock %lu [e_state %lX]\n", stat.index, tid, m.e_state[tid].region, m.released_logical_clock, get_logical_clock(&m.e_state[tid]), m.e_state);

#if PROFILE_WAIT_LEVEL >= 2
	zlog_level(delta_log, ROUTER_V3, "\tMax diff %lu Thread %d Current lock index %d Context %d\n", stat.max, stat.max_thread, stat.max_lock_count, stat.max_context);

#endif

#if PROFILE_WAIT_LEVEL >= 1
	stat.num_rounds = 0;
#endif
#if PROFILE_WAIT_LEVEL >= 1
	set_context(&m.e_state[tid], 50);
#endif
	bool locked = false;
	bool contested = false;
	//bool managed_to_lock = false;
	while (!locked) {
		det_mutex_wait_for_turn(m, tid);
#if PROFILE_WAIT_LEVEL >= 1
		++stat.num_rounds;
#endif
		if (m.lock->try_lock()) {
			assert(m.no_wait == -1);

			unsigned long cur = get_logical_clock(&m.e_state[tid]);

			if (cur > m.released_logical_clock) {
				zlog_level(delta_log, ROUTER_V3, "\tManaged to lock\n");

				if (contested) {
					assert(cur == m.released_logical_clock + 1);
				} 
				
				locked = true;

				m.no_wait = 0;
			} else {
				zlog_level(delta_log, ROUTER_V3, "\tLocked in logical time, %lu <= %lu\n", cur, m.released_logical_clock);

				//managed_to_lock = true;

				m.lock->unlock();

				inc_logical_clock(&m.e_state[tid], 1);
			}
		} else {
			zlog_level(delta_log, ROUTER_V3, "\tLocked in physical time\n");

			/* the assertion below fails because we only assign no_wait
			 * when lock is successful in both logical and physical time (see above) */
			//assert(m.no_wait == 0 || m.no_wait == 1);

			/* only increase the clock when waiting for a normal lock to prevent 
			 * race condition */
			if (m.no_wait == 0) {
				inc_logical_clock(&m.e_state[tid], 1);
				if (!contested) {
					contested = true;
				}
			}

			/* there is a possibility that we might lose our turn after unlocking while it's locked in
			 * logical time, so the assert below fails */
			//assert(!managed_to_lock);
			
			//m.e_state[tid].logical_clock->fetch_add(1, std::memory_order_relaxed);
		}
	}

#if PROFILE_WAIT_LEVEL >= 1
	set_context(&m.e_state[tid], 51);
#endif
	inc_logical_clock(&m.e_state[tid], 1);
	//m.e_state[tid].logical_clock->fetch_add(1, std::memory_order_relaxed);
#if defined(INSTRUMENT)
	start_logical_clock();
#elif defined(PMC)
	start_logical_clock(tid);
#endif

#if PROFILE_WAIT_LEVEL >= 1
	stat.wait_time = mtimer::now()-wait_start;
#endif
}

void det_mutex_lock_no_wait(det_mutex_t &m, int tid)
{
#if PROFILE_WAIT_LEVEL >= 1
	stats[tid].emplace_back();
	auto &stat = stats[tid].back();
	stat.index = m.e_state[tid].lock_count;
	++(m.e_state[tid].lock_count);
	stat.inet = m.e_state[tid].inet;
	stat.previous_context = m.e_state[tid].context;
#endif

#if PROFILE_WAIT_LEVEL >= 2
	unsigned long cur = get_logical_clock(&m.e_state[tid]);

	unsigned long min_clock = cur;
	int min_thread = tid;
	int min_inet = m.e_state[tid].inet;
	int min_context = m.e_state[tid].context;
	int min_lock_count = m.e_state[tid].lock_count;

	for (int i = 0; i < m.num_threads; ++i) {
		if (i == tid) {
			continue;
		}
		unsigned long other = get_logical_clock(&m.e_state[i]);
		if (other < cur) {
			if (other < min_clock) {
				min_clock = other;
				min_thread = i;
				min_inet = m.e_state[i].inet;
				min_context = m.e_state[i].context;
				min_lock_count = m.e_state[i].lock_count;
			}
		}
	}

	//if (cur < min_clock) {
		//printf("%lu %lu\n", cur, min_clock);
		//assert(false);
	//}
	assert(cur >= min_clock);
	//assert(min_thread != tid);

	//stat.max_clock_diff.push_back(std::make_pair(cur-min_clock, min_context));
	stat.logical_clock = cur;
	stat.max = cur-min_clock;
	stat.max_logical_clock = min_clock;
	stat.max_thread = min_thread;
	stat.max_inet = min_inet;
	stat.max_context = min_context;
	stat.max_lock_count = min_lock_count;
#endif

#if PROFILE_WAIT_LEVEL >= 1
	auto wait_start = mtimer::now();
#endif

#if defined(INSTRUMENT)
	stop_logical_clock();
#elif defined(PMC)
	stop_logical_clock(tid);
#endif

	zlog_level(delta_log, ROUTER_V3, "[%d] Thread %d acquiring lock, region %d released %lu cur clock %lu no wait [e_state %lX]\n", stat.index, tid, m.e_state[tid].region, m.released_logical_clock, get_logical_clock(&m.e_state[tid]), m.e_state);

#if PROFILE_WAIT_LEVEL >= 2
	zlog_level(delta_log, ROUTER_V3, "\tMax diff %lu Thread %d Current lock index %d Context %d\n", stat.max, stat.max_thread, stat.max_lock_count, stat.max_context);

#endif

#if PROFILE_WAIT_LEVEL >= 1
	stat.num_rounds = 0;
#endif
#if PROFILE_WAIT_LEVEL >= 1
	set_context(&m.e_state[tid], 52);
#endif
#if PROFILE_WAIT_LEVEL >= 1
	++stat.num_rounds;
#endif
	m.lock->lock();

	m.no_wait = 1;

	inc_logical_clock(&m.e_state[tid], 1);
	//set_logical_clock(&m.e_state[tid], m.released_logical_clock+1);
	//m.e_state[tid].logical_clock->fetch_add(1, std::memory_order_relaxed);

	//m.e_state[tid].logical_clock->fetch_add(1, std::memory_order_relaxed);
#if defined(INSTRUMENT)
	start_logical_clock();
#elif defined(PMC)
	start_logical_clock(tid);
#endif

#if PROFILE_WAIT_LEVEL >= 1
	stat.wait_time = mtimer::now()-wait_start;
#endif
}

void det_mutex_unlock(det_mutex_t &m, int tid)
{
#if defined(INSTRUMENT)
	stop_logical_clock();
#elif defined(PMC)
	stop_logical_clock(tid);
#endif

	unsigned long cur = get_logical_clock(&m.e_state[tid]);

	if (cur <= m.released_logical_clock) {
		zlog_error(delta_log, "%lu <= %lu\n", cur, m.released_logical_clock);
		assert(false);
	}

	m.released_logical_clock = cur;
	m.no_wait = -1;

	m.lock->unlock();

#if PROFILE_WAIT_LEVEL >= 1
	set_context(&m.e_state[tid], 53);
#endif
	inc_logical_clock(&m.e_state[tid], 1);
	//m.e_state[tid].logical_clock->fetch_add(1, std::memory_order_relaxed);
#if defined(INSTRUMENT)
	start_logical_clock();
#elif defined(PMC)
	start_logical_clock(tid);
#endif

	zlog_level(delta_log, ROUTER_V3, "Thread %d released lock, region %d new clock %lu [e_state %lX]\n", tid, m.e_state[tid].region, get_logical_clock(&m.e_state[tid]), m.e_state);
}

void det_mutex_unlock_no_wait(det_mutex_t &m, int tid)
{
#if defined(INSTRUMENT)
	stop_logical_clock();
#elif defined(PMC)
	stop_logical_clock(tid);
#endif

	m.no_wait = -1;

	m.lock->unlock();

#if PROFILE_WAIT_LEVEL >= 1
	set_context(&m.e_state[tid], 54);
#endif
	inc_logical_clock(&m.e_state[tid], 1);
	//m.e_state[tid].logical_clock->fetch_add(1, std::memory_order_relaxed);
	
#if defined(INSTRUMENT)
	start_logical_clock();
#elif defined(PMC)
	start_logical_clock(tid);
#endif

	zlog_level(delta_log, ROUTER_V3, "Thread %d released lock, region %d new clock %lu no wait [e_state %lX]\n", tid, m.e_state[tid].region, get_logical_clock(&m.e_state[tid]), m.e_state);
}
