#ifndef PTHREAD_BARRIER_H_
#define PTHREAD_BARRIER_H_

#include <pthread.h>
#include <errno.h>

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

#endif // PTHREAD_BARRIER_H_
