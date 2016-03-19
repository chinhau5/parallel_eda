#ifndef LOG_H
#define LOG_H

#include <zlog.h>

enum {
	ROUTER_V1 = ZLOG_LEVEL_DEBUG+3,
	ROUTER_V2 = ZLOG_LEVEL_DEBUG+2, 
	ROUTER_V3 = ZLOG_LEVEL_DEBUG+1
};

extern zlog_category_t *sort_log;
extern zlog_category_t *delta_log;
extern zlog_category_t *rr_log;
extern zlog_category_t *net_log;
extern zlog_category_t *schedule_log;
extern zlog_category_t *scheduler_log;
extern zlog_category_t *independent_log;
extern zlog_category_t *static_log;
extern zlog_category_t *dynamic_log;
extern zlog_category_t *missing_edge_log;
extern zlog_category_t *ss_log;

#define zlog_level(cat, level, ...) \
	zlog(cat, __FILE__, sizeof(__FILE__)-1, __func__, sizeof(__func__)-1, __LINE__, \
	level, __VA_ARGS__)

#define zlog_level(cat, level, ...)
#define zlog_debug(cat, ...)
#define zlog_info(cat, ...)
#define zlog_warn(cat, ...)

#endif
