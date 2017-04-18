#include "det_mutex.h"
#include "log.h"
#include <assert.h>
#include <tbb/cache_aligned_allocator.h>

void start_logical_clock();
void stop_logical_clock();

void start_logical_clock(int es);
void stop_logical_clock(int es);

static std::vector<std::vector<std::vector<det_mutex_stats_t>>> stats;

static std::vector<std::vector<std::pair<inc_time_t, int>>, tbb::cache_aligned_allocator<std::vector<std::pair<inc_time_t, int>>>> inc_time;

const std::vector<det_mutex_stats_t> &get_wait_stats(int level, int tid)
{
	assert(level < stats.size());
	assert(tid < stats[level].size());

	return stats[level][tid];
}

void init_wait_stats(int num_levels, int num_threads)
{
	stats.resize(num_levels);

	for (int level = 0; level < num_levels; ++level) {
		stats[level].resize(num_threads);

		for (int thread = 0; thread < num_threads; ++thread) {
			stats[level][thread].reserve(65536);
		}
	}
}

void clear_wait_stats()
{
	for (int level = 0; level < stats.size(); ++level) {
		for (int thread = 0; thread < stats[level].size(); ++thread) {
			stats[level][thread].clear();
		}
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

void det_mutex_reset(det_mutex_t &m)
{
	m.released_logical_clock = 0;

	m.released_thread = -1;
	m.released_no_wait = -1;
}

void det_mutex_init(det_mutex_t &m, exec_state_t *e_state, int num_threads)
{
	m.lock = new tbb::spin_mutex();
	m.e_state = e_state;
	m.num_threads = num_threads;

	det_mutex_reset(m);

	m.no_wait = -1;
}

void det_mutex_destroy(det_mutex_t &m)
{
	delete m.lock;
}

static void det_mutex_wait_for_turn_fixed(det_mutex_t &m, int tid, det_mutex_stats_t &stat)
{

}

static void det_mutex_wait_for_turn(det_mutex_t &m, int tid, int num_nodes)
{
	unsigned long cur = get_logical_clock(&m.e_state[tid]);
	int cur_inet = get_inet(&m.e_state[tid]);
	bool other_less;

	do {
		other_less = false;

		int threads_per_node = m.num_threads/num_nodes;
		int base = tid - (tid % threads_per_node);

		for (int i = base; i < base+threads_per_node && !other_less; ++i) {
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
			if (other_inet == cur_inet) {
				zlog_level(delta_log, ROUTER_V3, "thread %d inet %d thread %d inet %d\n", tid, cur_inet, i, other_inet);
				assert(false);
			}
			//assert(other_inet != cur_inet);

			other_less = other < cur || (other == cur && other_inet < cur_inet);
		}
	} while (other_less);
}

void det_mutex_lock(det_mutex_t &m, int tid, int num_nodes)
{
#if PROFILE_WAIT_LEVEL >= 1
	int level = m.e_state[tid].level;

	stats[level][tid].emplace_back();
	auto &stat = stats[level][tid].back();

	stat.index = m.e_state[tid].lock_count;
	stat.inet = get_inet(&m.e_state[tid]);
	stat.context = m.e_state[tid].context;
	stat.logical_clock = get_logical_clock(&m.e_state[tid]);
	stat.level = level;
	stat.released_logical_clock = m.released_logical_clock;
	stat.num_incs = stat.logical_clock-m.e_state[tid].prev_lock_clock;
#endif

#if PROFILE_WAIT_LEVEL >= 2
	unsigned long min_clock = stat.logical_clock;
	int min_inet = stat.inet;
	int min_context = stat.context;
	int min_thread = tid;
	int min_lock_count = stat.index;
	int min_level = stat.level;

	int threads_per_node = m.num_threads/num_nodes;
	int base = tid - (tid % threads_per_node);

	for (int i = base; i < base+threads_per_node; ++i) {
		if (i == tid) {
			continue;
		}

		//stat.others.emplace_back();
		//stat.others.back().tid = i;
		//stat.others.back().logical_clock = other;
		//stat.others.back().inet = get_inet(&m.e_state[i]);
		//stat.others.back().lock_count = i;
		//stat.others.back().tid = i;
		//min_lock_count = m.e_state[i].lock_count;
		//min_level = m.e_state[i].level;

		unsigned long other = get_logical_clock(&m.e_state[i]);
		if (other < stat.logical_clock) {
			if (other < min_clock) {
				min_clock = other;
				min_thread = i;
				min_inet = get_inet(&m.e_state[i]);
				min_context = m.e_state[i].context;
				min_lock_count = m.e_state[i].lock_count;
				min_level = m.e_state[i].level;
			}
		}

	}

	//if (cur < min_clock) {
		//printf("%lu %lu\n", cur, min_clock);
		//assert(false);
	//}
	assert(stat.logical_clock >= min_clock);
	//assert(min_thread != tid);

	//stat.max_clock_diff.push_back(std::make_pair(cur-min_clock, min_context));
	stat.max = stat.logical_clock-min_clock;
	stat.max_logical_clock = min_clock;
	stat.max_thread = min_thread;
	stat.max_inet = min_inet;
	stat.max_context = min_context;
	stat.max_lock_count = min_lock_count;
	stat.max_level = min_level;
#endif

#if PROFILE_WAIT_LEVEL >= 1
	auto wait_start = mtimer::now();
#endif

#if defined(TRACE_WAIT)
	tracepoint(hello_world, my_first_tracepoint, tid, m.e_state[tid].context, 1, m.e_state[tid].lock_count);
#endif

	m.e_state[tid].prev_lock_clock = get_logical_clock(&m.e_state[tid]);

#if defined(INSTRUMENT)
	stop_logical_clock();
#elif defined(PMC)
	stop_logical_clock(tid);
#endif

	zlog_level(delta_log, ROUTER_V3, "[%d] Thread %d acquiring lock, region %d cur clock %lu [e_state %lX]\n", stat.index, tid, m.e_state[tid].region, get_logical_clock(&m.e_state[tid]), m.e_state);

#if PROFILE_WAIT_LEVEL >= 2
	zlog_level(delta_log, ROUTER_V3, "\tMax diff %lu Thread %d Current lock index %d Context %d\n", stat.max, stat.max_thread, stat.max_lock_count, stat.max_context);

#endif

#if PROFILE_WAIT_LEVEL >= 1
	stat.num_rounds = 0;
#endif
#if defined(SET_CONTEXT)
	set_context(&m.e_state[tid], 50);
#endif
	bool locked = false;
	bool contested = false;
	//bool managed_to_lock = false;
	while (!locked) {
		det_mutex_wait_for_turn(m, tid, num_nodes);
#if PROFILE_WAIT_LEVEL >= 1
		++stat.num_rounds;
#endif
		if (m.lock->try_lock()) {
			assert(m.no_wait == -1);

			unsigned long cur = get_logical_clock(&m.e_state[tid]);

			if (cur > m.released_logical_clock) {
				zlog_level(delta_log, ROUTER_V3, "\tManaged to lock, released_thread %d released_logical_clock %lu released_no_wait %d\n", m.released_thread, m.released_logical_clock, m.released_no_wait);

				if (contested) {
					assert(cur == m.released_logical_clock + 1);
				} 
				
				locked = true;

				m.no_wait = 0;
			} else {
				zlog_level(delta_log, ROUTER_V3, "\tLocked in logical time, %lu <= %lu\n", cur, m.released_logical_clock);

				//managed_to_lock = true;

				m.lock->unlock();

				//inc_logical_clock(&m.e_state[tid], 1);
				set_logical_clock(&m.e_state[tid], m.released_logical_clock+1);
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

#if defined(TRACE_WAIT)
	tracepoint(hello_world, my_first_tracepoint, tid, m.e_state[tid].context, 2, m.e_state[tid].lock_count);
#endif

#if defined(SET_CONTEXT)
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
	stat.start = mtimer::now();
	stat.wait_time = stat.start-wait_start;
#endif

	++(m.e_state[tid].lock_count);
}

void det_mutex_lock_no_wait(det_mutex_t &m, int tid)
{
#if PROFILE_WAIT_LEVEL >= 1
	int level = m.e_state[tid].level;

	stats[level][tid].emplace_back();
	auto &stat = stats[level][tid].back();

	stat.index = m.e_state[tid].lock_count;
	stat.inet = get_inet(&m.e_state[tid]);
	stat.context = m.e_state[tid].context;
	stat.logical_clock = get_logical_clock(&m.e_state[tid]);
	stat.level = level;
	stat.released_logical_clock = 0;
	stat.num_incs = stat.logical_clock-m.e_state[tid].prev_lock_clock;

	stat.max = 0;
	stat.max_logical_clock = 0;
	stat.max_thread = -1;
	stat.max_inet = -1;
	stat.max_context = -2;
	stat.max_lock_count = -1;
	stat.max_level = -1;
#endif

#if PROFILE_WAIT_LEVEL >= 1
	auto wait_start = mtimer::now();
#endif

#if defined(TRACE_WAIT)
	tracepoint(hello_world, my_first_tracepoint, tid, m.e_state[tid].context, 3, m.e_state[tid].lock_count);
#endif

	m.e_state[tid].prev_lock_clock = get_logical_clock(&m.e_state[tid]);

#if defined(INSTRUMENT)
	stop_logical_clock();
#elif defined(PMC)
	stop_logical_clock(tid);
#endif

	zlog_level(delta_log, ROUTER_V3, "[%d] Thread %d acquiring lock, region %d cur clock %lu no wait [e_state %lX]\n", stat.index, tid, m.e_state[tid].region, get_logical_clock(&m.e_state[tid]), m.e_state);

#if defined(SET_CONTEXT)
	set_context(&m.e_state[tid], 52);
#endif
	m.lock->lock();

#if defined(TRACE_WAIT)
	tracepoint(hello_world, my_first_tracepoint, tid, m.e_state[tid].context, 4, m.e_state[tid].lock_count);
#endif

	zlog_level(delta_log, ROUTER_V3, "\tManaged to lock, released_thread %d released_logical_clock %lu released_no_wait %d\n", m.released_thread, m.released_logical_clock, m.released_no_wait);

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
	stat.num_rounds = 1;
	stat.start = mtimer::now();
	stat.wait_time = stat.start-wait_start;
#endif

	++(m.e_state[tid].lock_count);
}

void det_mutex_unlock(det_mutex_t &m, int tid)
{
#if defined(INSTRUMENT)
	stop_logical_clock();
#elif defined(PMC)
	stop_logical_clock(tid);
#endif

#if PROFILE_WAIT_LEVEL >= 1
	auto &stat = stats[m.e_state[tid].level][tid].back();
	stat.length = mtimer::now()-stat.start;
#endif

	unsigned long cur = get_logical_clock(&m.e_state[tid]);

	if (cur <= m.released_logical_clock) {
		zlog_error(delta_log, "%lu <= %lu\n", cur, m.released_logical_clock);
		assert(false);
	}

	m.released_logical_clock = cur;
	m.released_thread = tid;
	m.released_no_wait = 0;

	m.no_wait = -1;

	m.lock->unlock();

	zlog_level(delta_log, ROUTER_V3, "Thread %d released lock, region %d clock %lu [e_state %lX]\n", tid, m.e_state[tid].region, get_logical_clock(&m.e_state[tid]), m.e_state);

#if defined(SET_CONTEXT)
	set_context(&m.e_state[tid], 53);
#endif
	inc_logical_clock(&m.e_state[tid], 1);
	//m.e_state[tid].logical_clock->fetch_add(1, std::memory_order_relaxed);
#if defined(INSTRUMENT)
	start_logical_clock();
#elif defined(PMC)
	start_logical_clock(tid);
#endif
}

void det_mutex_unlock_no_wait(det_mutex_t &m, int tid)
{
#if defined(INSTRUMENT)
	stop_logical_clock();
#elif defined(PMC)
	stop_logical_clock(tid);
#endif

#if PROFILE_WAIT_LEVEL >= 1
	auto &stat = stats[m.e_state[tid].level][tid].back();
	stat.length = mtimer::now()-stat.start;
#endif

	m.released_thread = tid;
	m.released_no_wait = 1;

	m.no_wait = -1;

	m.lock->unlock();

	zlog_level(delta_log, ROUTER_V3, "Thread %d released lock, region %d clock %lu no wait [e_state %lX]\n", tid, m.e_state[tid].region, get_logical_clock(&m.e_state[tid]), m.e_state);

#if defined(SET_CONTEXT)
	set_context(&m.e_state[tid], 54);
#endif
	inc_logical_clock(&m.e_state[tid], 1);
	//m.e_state[tid].logical_clock->fetch_add(1, std::memory_order_relaxed);
	
#if defined(INSTRUMENT)
	start_logical_clock();
#elif defined(PMC)
	start_logical_clock(tid);
#endif
}
