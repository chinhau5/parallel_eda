#ifndef PCH_H
#define PCH_H

#include <stdio.h>
#include <stdlib.h>
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
#include <cmath>
#include <map>
#include <atomic>

//#define TBB_USE_DEBUG 1
#include <tbb/task_scheduler_init.h>
#include <tbb/concurrent_queue.h>
#include <tbb/tbb.h>

#ifdef __linux__
#include <sys/syscall.h>
//#include <ittnotify.h>
#endif

#include <boost/numeric/interval.hpp>
#include <boost/geometry.hpp>
#include <boost/iterator/counting_iterator.hpp>
#include <boost/range/counting_range.hpp>
#include <boost/range/adaptor/transformed.hpp>
#include <boost/range/adaptor/sliced.hpp>
#include <boost/range/algorithm/copy.hpp>
#include <boost/thread/barrier.hpp>
#include <boost/property_map/vector_property_map.hpp>
#include <boost/property_map/shared_array_property_map.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/sequential_vertex_coloring.hpp>
#include <boost/graph/smallest_last_ordering.hpp>
#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics/stats.hpp>
#include <boost/accumulators/statistics/mean.hpp>
#include <boost/accumulators/statistics/min.hpp>
#include <boost/accumulators/statistics/max.hpp>
#include <boost/accumulators/statistics/moment.hpp>
#include <boost/accumulators/statistics/variance.hpp>
#include <boost/heap/binomial_heap.hpp>

#endif
