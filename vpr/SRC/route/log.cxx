#include "log.h"

zlog_category_t *sort_log;
zlog_category_t *delta_log;
zlog_category_t *rr_log;
zlog_category_t *net_log;
zlog_category_t *schedule_log;
zlog_category_t *scheduler_log;
zlog_category_t *independent_log;
zlog_category_t *static_log;
zlog_category_t *dynamic_log;
zlog_category_t *missing_edge_log;
zlog_category_t *ss_log;

void init_logging()
{
	delta_log = zlog_get_category("delta");
	rr_log = zlog_get_category("rr");
	net_log = zlog_get_category("net");
	schedule_log = zlog_get_category("schedule");
	scheduler_log = zlog_get_category("scheduler");
	independent_log = zlog_get_category("independent");
	sort_log = zlog_get_category("sort");
	dynamic_log = zlog_get_category("dynamic");
	static_log = zlog_get_category("static");
	missing_edge_log = zlog_get_category("missing_edge");
	ss_log = zlog_get_category("second_stage");
}

