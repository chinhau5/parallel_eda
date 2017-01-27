#include "new_partitioner.h"

void split(const std::vector<std::pair<point, int>> &sitems,
		float mid, bool xsort,
		std::vector<std::pair<point, int>> &left, std::vector<std::pair<point, int>> &right)
{
	for (const auto &item : sitems) {
		if (xsort) {
			if (bg::get<0>(item.first) < mid) {
				left.push_back(item);
			} else {
				assert(bg::get<0>(item.first) >= mid);
				right.push_back(item);
			}
		} else {
			if (bg::get<1>(item.first) < mid) {
				left.push_back(item);
			} else {
				assert(bg::get<1>(item.first) >= mid);
				right.push_back(item);
			}
		}
	}
}

float get_median_mid(const std::vector<std::pair<point, int>> &xsorted_items, const std::vector<std::pair<point, int>> &ysorted_items, bool xcut)
{
	float mid;
	if ((xsorted_items.size() % 2) == 0) {
		if (xcut) {
			mid = (float)(bg::get<0>(xsorted_items[xsorted_items.size()/2 - 1].first) + bg::get<0>(xsorted_items[xsorted_items.size()/2].first))/2;
		} else {
			mid = (float)(bg::get<1>(ysorted_items[ysorted_items.size()/2 - 1].first) + bg::get<1>(ysorted_items[ysorted_items.size()/2].first))/2;
		}
	} else {
		if (xcut) {
			mid = bg::get<0>(xsorted_items[xsorted_items.size()/2].first);
		} else {
			mid = bg::get<1>(ysorted_items[ysorted_items.size()/2].first);
		}
	}
	return mid;
}

void make_partition(const std::vector<std::pair<point, int>> &items, std::vector<partition_t<int>> &partitions)
{
	assert(!items.empty());

	partition_t<int> part;
	part.bounding_box = bg::make_inverse<box>();
	for (const auto &item : items) {
		bg::expand(part.bounding_box, item.first);
		part.items.push_back(item.second);
	}
	partitions.emplace_back(part);
}

void partition_internal(const std::vector<std::pair<point, int>> &xsorted_items, const std::vector<std::pair<point, int>> &ysorted_items, bool xcut, bool one_sided, float bb_area_threshold, int level, std::vector<partition_t<int>> &partitions)
{
	assert(xsorted_items.size() == ysorted_items.size());
	assert(!xsorted_items.empty());

	for (int i = 0; i < xsorted_items.size(); ++i) {
		if (i > 0) {
			printf("x[%d] = %d x[%d] = %d\n", i-1, bg::get<0>(xsorted_items[i-1].first), i, bg::get<0>(xsorted_items[i].first));
			assert(bg::get<0>(xsorted_items[i-1].first) <= bg::get<0>(xsorted_items[i].first));
		}
	}
	for (int i = 0; i < ysorted_items.size(); ++i) {
		if (i > 0) {
			printf("y[%d] = %d y[%d] = %d\n", i-1, bg::get<1>(ysorted_items[i-1].first), i, bg::get<1>(ysorted_items[i].first));
			assert(bg::get<1>(ysorted_items[i-1].first) <= bg::get<1>(ysorted_items[i].first));
		}
	}

	box bb = bg::make_inverse<box>();
	for (const auto &item : xsorted_items) {
		bg::expand(bb, item.first);
	}

	//assert((bg::get<bg::min_corner, 0>(bb) > 0));
	//assert((bg::get<bg::min_corner, 1>(bb) > 0));

	bg::set<bg::min_corner, 0>(bb, bg::get<bg::min_corner, 0>(bb)-1);
	bg::set<bg::min_corner, 1>(bb, bg::get<bg::min_corner, 1>(bb)-1);
	assert(bg::is_valid(bb));

	double area = bg::area(bb);
	assert(area > 0);

	/* base case */
	if (area <= bb_area_threshold) {
		make_partition(xsorted_items, partitions);

		return;
	}

	float mid = get_median_mid(xsorted_items, ysorted_items, xcut);

	//PRINT("level %d xcut %d mid %g\n", level, xcut ? 1 : 0, mid);

	std::vector<std::pair<point, int>> xsorted_left;
	std::vector<std::pair<point, int>> xsorted_right;

	split(xsorted_items, 
			mid, xcut,
			xsorted_left, xsorted_right);

	std::vector<std::pair<point, int>> ysorted_left;
	std::vector<std::pair<point, int>> ysorted_right;

	split(ysorted_items,
			mid, xcut,
			ysorted_left, ysorted_right);

	assert(xsorted_left.size() == ysorted_left.size());  
	assert(xsorted_right.size() == ysorted_right.size());  

	if (one_sided) {
		if (!xsorted_left.empty() && xsorted_right.empty()) {
			/* second time one sided, stop recursing */
			//PRINTS("\t1double one sided stop recursing\n");
			make_partition(xsorted_left, partitions);
		} else if (xsorted_left.empty() && !xsorted_right.empty()) {
			/* second time one sided, stop recursing */
			//PRINTS("\t2double one sided stop recursing\n");
			make_partition(xsorted_right, partitions);
		} else {
			assert(!xsorted_left.empty() && !xsorted_right.empty());

			partition_internal(xsorted_left, ysorted_left, !xcut, false, bb_area_threshold, level+1, partitions);
			partition_internal(xsorted_right, ysorted_right, !xcut, false, bb_area_threshold, level+1, partitions);
		}
	} else {
		if (!xsorted_left.empty() && xsorted_right.empty()) {
			/* first time one sided, continue recursing */
			//PRINTS("\t1first one sided continue recursing\n");
			partition_internal(xsorted_left, ysorted_left, !xcut, true, bb_area_threshold, level+1, partitions);
		} else if (xsorted_left.empty() && !xsorted_right.empty()) {
			/* first time one sided, continue recursing */
			//PRINTS("\t2first one sided continue recursing\n");
			partition_internal(xsorted_right, ysorted_right, !xcut, true, bb_area_threshold, level+1, partitions);
		} else {
			assert(!xsorted_left.empty() && !xsorted_right.empty());

			partition_internal(xsorted_left, ysorted_left, !xcut, false, bb_area_threshold, level+1, partitions);
			partition_internal(xsorted_right, ysorted_right, !xcut, false, bb_area_threshold, level+1, partitions);
		}
	}
}
