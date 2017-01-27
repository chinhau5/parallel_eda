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
static vector<vector<FILE *>> log_files_2;
static std::mutex log_files_lock;

void init_logging(int num_threads)
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

	log_files_2.resize(100, vector<FILE *>(num_threads, nullptr));
	//std::fill(begin(log_files_2), end(log_files_2), nullptr);
}

int concurrent_log_impl_2(zlog_msg_t *msg)
{
	char *citer = zlog_get_mdc("iter");
	assert(citer);

	char *ctid = zlog_get_mdc("tid");
	assert(ctid);

	int iter = atoi(citer);
    int tid	= atoi(ctid);
	assert(tid >= 0 && tid < log_files_2.size());

	FILE *file;
	if (!log_files_2[iter][tid]) {
		file = fopen(msg->path, "w");
		if (!file) {
			perror(nullptr);
			assert(false);
		}
		log_files_2[iter][tid] = file;
	} else {
		file = log_files_2[iter][tid];
	}

	fprintf(file, "%s", msg->buf);
	fflush(file);

	return 0;
}

int concurrent_log_impl(zlog_msg_t *msg)
{
	log_files_lock.lock();
	auto iter = log_files.find(string(msg->path));
	if (iter == end(log_files)) {
		char filename[256];
		sprintf(filename, "%s", msg->path);

		FILE *file = fopen(filename, "w");
		if (!file) {
			perror(nullptr);
			assert(false);
		}

		auto res = log_files.insert(make_pair(string(msg->path), file));
		assert(res.second);

		iter = res.first;
	}
	log_files_lock.unlock();
	fprintf(iter->second, "%s", msg->buf);
	fflush(iter->second);

	return 0;
}

