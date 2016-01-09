#ifndef PCH_H
#define PCH_H

#include <stdio.h>
#include <assert.h>

#include <fcntl.h>
#include <unistd.h>
#include <semaphore.h>

#include <iostream>
#include <string>
#include <sstream>
#include <random>
#include <memory>
#include <ctime>
#include <chrono>
#include <mutex>

//#define TBB_USE_DEBUG 1
#include <tbb/task_scheduler_init.h>
#include <tbb/concurrent_queue.h>
#include <tbb/tbb.h>

#ifdef __linux__
#include <sys/syscall.h>
#include <ittnotify.h>
#endif

#include <boost/numeric/interval.hpp>
#include <boost/geometry.hpp>

#endif
