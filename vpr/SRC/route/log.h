#ifndef LOG_H
#define LOG_H

#include <zlog.h>

enum {
	ROUTER_V1 = ZLOG_LEVEL_DEBUG+3,
	ROUTER_V2 = ZLOG_LEVEL_DEBUG+2, 
	ROUTER_V3 = ZLOG_LEVEL_DEBUG+1
};

#define zlog_level(cat, level, ...) \
	zlog(cat, __FILE__, sizeof(__FILE__)-1, __func__, sizeof(__func__)-1, __LINE__, \
	level, __VA_ARGS__)

//#define zlog_level(cat, level, ...)

#endif
