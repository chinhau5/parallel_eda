#include <assert.h>
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

void concurrent_log_impl(zlog_msg_t *msg, std::vector<std::vector<FILE *>> &log_files, int iter, int tid)
{
	assert(iter >= 0 && iter < log_files.size());
	assert(tid >= 0 && tid < log_files[iter].size());
	FILE *file = log_files[iter][tid];
	if (!file) {
		char filename[256];
		sprintf(filename, "%s%s", LOG_PATH_PREFIX, msg->path);

		file = fopen(filename, "w");
		if (!file) {
			perror(nullptr);
			assert(false);
		}

		log_files[iter][tid] = file;
	}
	fprintf(file, "%s", msg->buf);
	fflush(file);
}

