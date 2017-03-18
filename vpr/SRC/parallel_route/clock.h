#ifndef CCLOCK_H
#define CCLOCK_H

//using namespace std;
#include <assert.h>

struct myclock {
	typedef std::chrono::nanoseconds duration;
	typedef duration::rep rep;
	typedef duration::period period;
	typedef std::chrono::time_point<myclock, duration> time_point;

	static time_point now() {
		timespec time;
#ifdef __linux__
		assert(!clock_gettime(CLOCK_MONOTONIC, &time));
#else
		assert(false);
#endif
		return time_point(std::chrono::nanoseconds((time.tv_sec * (int)1e9) + time.tv_nsec));
	}
};

#endif
