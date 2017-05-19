#include "pch.h"
#include "log.h"
#include "delta_stepping_common.h" 

using namespace std;

void bucket_insert(Buckets &buckets, float delta, vector<int> &in_bucket, int v, float distance)
{
	int b_i = floor(distance/delta);

	zlog_level(delta_log, ROUTER_V3, "\t\t\tInserting vertex %d dist: %g into bucket %d\n", v, distance, b_i);

	assert(b_i >= 0);

	if (b_i >= buckets.size()) {
		buckets.resize(b_i+1, nullptr);
	}
	if (!buckets[b_i]) {
		buckets[b_i] = new vector<int>();
	}
	buckets[b_i]->push_back(v);

	assert(in_bucket[v] == -1);
	in_bucket[v] = b_i;
}

bool bucket_remove(Buckets &buckets, float delta, vector<int> &in_bucket, const vector<bool> &vertex_deleted, int v, float *distance)
{
	if (distance[v] == std::numeric_limits<float>::max() || vertex_deleted[v]) {
		assert(in_bucket[v] == -1);
		return false;
	}

	int b_i = floor(distance[v]/delta);

	zlog_level(delta_log, ROUTER_V3, "\t\t\tRemoving vertex %d from bucket %d\n", v, b_i);

	assert(b_i >= 0 && b_i < buckets.size());

	auto *bucket = buckets[b_i];
	auto iter = find(bucket->begin(), bucket->end(), v);
	assert(iter != bucket->end());
	bucket->erase(iter);

	assert(in_bucket[v] >= 0);
	in_bucket[v] = -1;

	return true;
}

int bucket_get_non_empty(const Buckets &buckets)
{
	int non_empty = -1;

	for (int i = 0; i < buckets.size() && non_empty < 0; ++i) {
		if (buckets[i] && !buckets[i]->empty()) {
			non_empty = i;
		}
	}

	return non_empty;
}

void bucket_insert(delta_stepping_t &ds, int v, float distance)
{
	int b_i = floor(distance/ds.delta);

	DEBUG_PRINTS("\t\t\tInserting vertex %d dist: %g into bucket %d\n", v, distance, b_i);

	assert(b_i >= 0);

	if (b_i >= ds.buckets.size()) {
		ds.buckets.resize(b_i+1, nullptr);
	}
	if (!ds.buckets[b_i]) {
		ds.buckets[b_i] = new vector<int>();
	}
	ds.buckets[b_i]->push_back(v);

	assert(ds.in_bucket[v] == -1);
	ds.in_bucket[v] = b_i;
}

bool bucket_remove(delta_stepping_t &ds, int v, const float *distance)
{
	if (ds.in_bucket[v] == -1) {
		return false;
	}

	int b_i = floor(distance[v]/ds.delta);

	DEBUG_PRINTS("\t\t\tRemoving vertex %d from bucket %d\n", v, b_i);

	assert(b_i >= 0 && b_i < ds.buckets.size());

	auto *bucket = ds.buckets[b_i];
	auto iter = find(bucket->begin(), bucket->end(), v);
	assert(iter != bucket->end());
	bucket->erase(iter);

	assert(ds.in_bucket[v] >= 0);
	ds.in_bucket[v] = -1;

	return true;
}
