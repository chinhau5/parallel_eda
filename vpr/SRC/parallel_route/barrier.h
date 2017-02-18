#ifndef PTHREAD_BARRIER_H_
#define PTHREAD_BARRIER_H_

#include <pthread.h>
#include <errno.h>
#include <atomic>

typedef int my_pthread_barrierattr_t;
typedef struct
{
    pthread_mutex_t mutex;
    pthread_cond_t cond;
	pthread_t first_waiting_thread;
    int count;
    int tripCount;
} my_pthread_barrier_t;

int my_pthread_barrier_init(my_pthread_barrier_t *barrier, const my_pthread_barrierattr_t *attr, unsigned int count);

int my_pthread_barrier_destroy(my_pthread_barrier_t *barrier);

int my_pthread_barrier_wait(my_pthread_barrier_t *barrier, int tid);

class SpinningBarrier {
	private:
		int _count;
		std::atomic<int> _bar; // Counter of threads, faced barrier.
		std::atomic<int> _passed; // Number of barriers, passed by all threads.

	public:
		SpinningBarrier(int count)
			: _count(count), _bar(0), _passed(0)
		{
		}

		void wait()
		{
			int passed_old = _passed.load(std::memory_order_relaxed);

			if (_bar.fetch_add(1) == (_count - 1)) {
				// The last thread, faced barrier.
				_bar.store(0, std::memory_order_relaxed);
				// Synchronize and store in one operation.
				_passed.store(passed_old + 1, std::memory_order_relaxed);
			} else {
				// Not the last thread. Wait others.
				while (_passed.load(std::memory_order_relaxed) == passed_old) {
				}
				// Need to synchronize cache with other threads, passed barrier.
				//std::atomic_thread_fence(std::memory_order_acquire);
			}
		}

		template<typename Func>
		void wait(const Func &f)
		{
			int passed_old = _passed.load(std::memory_order_relaxed);

			if (_bar.fetch_add(1) == (_count - 1)) {
				// The last thread, faced barrier.
				_bar.store(0, std::memory_order_relaxed);
				// Synchronize and store in one operation.
				_passed.store(passed_old + 1, std::memory_order_relaxed);
			} else {
				// Not the last thread. Wait others.
				while (_passed.load(std::memory_order_relaxed) == passed_old) {
					f();
				}
				// Need to synchronize cache with other threads, passed barrier.
				//std::atomic_thread_fence(std::memory_order_acquire);
			}
		}
};

#endif // PTHREAD_BARRIER_H_
