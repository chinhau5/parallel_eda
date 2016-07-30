#include "pch.h"
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

using namespace std;

static map<string, FILE *> log_files;

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

#define LOG_PATH_PREFIX "/tmp/"

int concurrent_log_impl(zlog_msg_t *msg)
{
	auto iter = log_files.find(string(msg->path));
	if (iter == end(log_files)) {
		char filename[256];
		sprintf(filename, "%s%s", LOG_PATH_PREFIX, msg->path);

		FILE *file = fopen(filename, "w");
		if (!file) {
			perror(nullptr);
			assert(false);
		}

		auto res = log_files.insert(make_pair(string(msg->path), file));
		assert(res.second);

		iter = res.first;
	}
	fprintf(iter->second, "%s", msg->buf);
	fflush(iter->second);

	return 0;
}

