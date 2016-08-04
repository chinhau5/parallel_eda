#ifndef CCLOCK_H
#define CCLOCK_H

using namespace std;

struct myclock {
	typedef chrono::nanoseconds duration;
	typedef duration::rep rep;
	typedef duration::period period;
	typedef chrono::time_point<myclock, duration> time_point;

	static time_point now() {
		timespec time;
#ifdef __linux__
		assert(!clock_gettime(CLOCK_MONOTONIC, &time));
#endif
		return time_point(chrono::nanoseconds((time.tv_sec * (int)1e9) + time.tv_nsec));
	}
};

#endif
