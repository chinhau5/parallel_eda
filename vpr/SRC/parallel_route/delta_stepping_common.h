#ifndef DELTA_STEPPING_COMMON_H
#define DELTA_STEPPING_COMMON_H

using Buckets = std::vector<std::vector<int> *>;

struct delta_stepping_t {
	Buckets buckets;
	float delta;
	std::vector<int> in_bucket;
	std::vector<bool> vertex_deleted;
};

//#define DEBUG_PRINTS(msg, ...) printf((msg), __VA_ARGS__)
//#define DEBUG_PRINT(msg) printf((msg))
#define DEBUG_PRINTS(msg, ...) 
#define DEBUG_PRINT(msg) 

void bucket_insert(Buckets &buckets, float delta, std::vector<int> &in_bucket, int v, float distance);

bool bucket_remove(Buckets &buckets, float delta, std::vector<int> &in_bucket, const std::vector<bool> &vertex_deleted, int v, const float *distance);

int bucket_get_non_empty(const Buckets &buckets);

void bucket_insert(delta_stepping_t &ds, int v, float distance);

bool bucket_remove(delta_stepping_t &ds, int v, const float *distance);

#endif
