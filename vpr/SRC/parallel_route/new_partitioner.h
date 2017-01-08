#ifndef NEW_PARTITIONER_H
#define NEW_PARTITIONER_H

#include "geometry.h"

template<typename Item>
struct partition_t {
	std::vector<Item> items;
	box bounding_box;
};

void split(const std::vector<std::pair<point, int>> &sitems,
		float mid, bool xsort,
		std::vector<std::pair<point, int>> &left, std::vector<std::pair<point, int>> &right);

float get_median_mid(const std::vector<std::pair<point, int>> &xsorted_items, const std::vector<std::pair<point, int>> &ysorted_items, bool xcut);

void make_partition(const std::vector<std::pair<point, int>> &items, std::vector<partition_t<int>> &partitions);

void partition_internal(const std::vector<std::pair<point, int>> &xsorted_items, const std::vector<std::pair<point, int>> &ysorted_items, bool xcut, bool one_sided, float bb_area_threshold, int level, std::vector<partition_t<int>> &partitions);

template<typename Item, typename ItemToPoint>
void partition(const std::vector<Item> &items, const ItemToPoint &to_point, float bb_area_threshold, std::vector<partition_t<int>> &partitions)
{
	std::vector<std::pair<point, int>> xsorted_items;
	std::vector<std::pair<point, int>> ysorted_items;

	box initial_bb = bg::make_inverse<box>();

	for (int i = 0; i < items.size(); ++i) {
		const auto &point = to_point(items[i]);

		xsorted_items.push_back(std::make_pair(point, i));
		ysorted_items.push_back(std::make_pair(point, i));

		bg::expand(initial_bb, point);
	}

	//std::sort(begin(xsorted_items), end(xsorted_items));
	//std::sort(begin(ysorted_items), end(ysorted_items));
	std::sort(begin(xsorted_items), end(xsorted_items), [&] (const std::pair<point, int> &a, const std::pair<point, int> &b) -> bool { return std::make_pair(bg::get<0>(a.first), a.second) < std::make_pair(bg::get<0>(b.first), b.second); });
	std::sort(begin(ysorted_items), end(ysorted_items), [&] (const std::pair<point, int> &a, const std::pair<point, int> &b) -> bool { return std::make_pair(bg::get<1>(a.first), a.second) < std::make_pair(bg::get<1>(b.first), b.second); });
	//
	assert(partitions.empty());

	int width = bg::get<bg::max_corner, 0>(initial_bb) - bg::get<bg::min_corner, 0>(initial_bb);
	int height = bg::get<bg::max_corner, 1>(initial_bb) - bg::get<bg::min_corner, 1>(initial_bb);
	bool xcut = width > height;
	//PRINT("width %d height %d\n", width, height);
	//if (xcut) {
		//PRINTS("using xcut initially\n");
	//} else {
		//PRINTS("not using xcut initially\n");
	//}

	partition_internal(xsorted_items, ysorted_items, xcut, false, bb_area_threshold, 0, partitions);
}


#endif
