#include <stdio.h>
#include <zlog.h>
#include "barrier.h"

int my_pthread_barrier_init(my_pthread_barrier_t *barrier, const my_pthread_barrierattr_t *attr, unsigned int count)
{
    if(count == 0)
    {
        errno = EINVAL;
        return -1;
    }
    if(pthread_mutex_init(&barrier->mutex, 0) < 0)
    {
        return -1;
    }
    if(pthread_cond_init(&barrier->cond, 0) < 0)
    {
        pthread_mutex_destroy(&barrier->mutex);
        return -1;
    }
    barrier->tripCount = count;
    barrier->count = 0;

    return 0;
}

int my_pthread_barrier_destroy(my_pthread_barrier_t *barrier)
{
    pthread_cond_destroy(&barrier->cond);
    pthread_mutex_destroy(&barrier->mutex);
    return 0;
}

int my_pthread_barrier_wait(my_pthread_barrier_t *barrier, int tid)
{
    pthread_mutex_lock(&barrier->mutex);
	if (barrier->count == 0) {
		/* first thread to start waiting on barrier */
		barrier->first_waiting_thread = pthread_self();
	}
    ++(barrier->count);

	int ret;
    if(barrier->count >= barrier->tripCount)
    {
        barrier->count = 0;
		dzlog_debug("%d Waking all threads\n", tid);
        pthread_cond_broadcast(&barrier->cond);
        pthread_mutex_unlock(&barrier->mutex);
		ret = 0;
    }
    else
    {
		dzlog_debug("%d Waiting with count: %d\n", tid, barrier->count);
        pthread_cond_wait(&barrier->cond, &(barrier->mutex));
		if (barrier->first_waiting_thread == pthread_self()) {
			ret = 1;
		} else {
			ret = 0;
		}
        pthread_mutex_unlock(&barrier->mutex);
    }
	return ret;
}
