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

void add_inc_time(const inc_time_t &time, int tid, int context)
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
	stat.previous_context = m.e_state[tid].context;
	stat.inet = m.e_state[tid].inet;
#endif

#if PROFILE_WAIT_LEVEL >= 2
	unsigned long cur = get_logical_clock(&m.e_state[tid]);

	unsigned long min_clock = cur;
	int min_context = m.e_state[tid].context;
	int min_lock_count = m.e_state[tid].lock_count;
	int min_thread = tid;
	int min_inet = m.e_state[tid].inet;

	for (int i = 0; i < m.num_threads; ++i) {
		if (i == tid) {
			continue;
		}
		unsigned long other = get_logical_clock(&m.e_state[i]);
		if (other < cur) {
			if (other < min_clock) {
				min_clock = other;
				min_context = m.e_state[i].context;
				min_lock_count = m.e_state[i].lock_count;
				min_thread = i;
				min_inet = m.e_state[i].inet;
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
	stat.max_context = min_context;
	stat.max_lock_count = min_lock_count;
	stat.max_inet = min_inet;
#endif

#if PROFILE_WAIT_LEVEL >= 1
	auto wait_start = mtimer::now();
#endif

#if defined(INSTRUMENT)
	stop_logical_clock();
#elif defined(PMC)
	stop_logical_clock(tid);
#endif

	zlog_level(delta_log, ROUTER_V3, "Thread %d acquiring lock, cur clock %d [e_state %lX]\n", tid, get_logical_clock(&m.e_state[tid]), m.e_state);

#if PROFILE_WAIT_LEVEL >= 2
	zlog_level(delta_log, ROUTER_V3, "\tIndex %d Max diff %lu Thread %d Current lock index %d Context %d\n", stat.index, stat.max, stat.max_thread, stat.max_lock_count, stat.max_context);

#endif

#if PROFILE_WAIT_LEVEL >= 1
	stat.num_rounds = 0;
#endif
	while (true) {
		det_mutex_wait_for_turn(m, tid);
#if PROFILE_WAIT_LEVEL >= 1
		++stat.num_rounds;
#endif
		if (m.lock->try_lock()) {
			if (m.released_logical_clock >= get_logical_clock(&m.e_state[tid])) {
				m.lock->unlock();
			} else {
				break;
			}
		}
		inc_logical_clock(&m.e_state[tid], 1, 50);
		//m.e_state[tid].logical_clock->fetch_add(1, std::memory_order_relaxed);
	}

	inc_logical_clock(&m.e_state[tid], 1, 51);
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

	m.released_logical_clock = get_logical_clock(&m.e_state[tid]);

	m.lock->unlock();

	inc_logical_clock(&m.e_state[tid], 1, 52);
	//m.e_state[tid].logical_clock->fetch_add(1, std::memory_order_relaxed);
#if defined(INSTRUMENT)
	start_logical_clock();
#elif defined(PMC)
	start_logical_clock(tid);
#endif

	zlog_level(delta_log, ROUTER_V3, "Thread %d released lock, new clock %d [e_state %lX]\n", tid, get_logical_clock(&m.e_state[tid]), m.e_state);
}
