#include "pch.h"
#include "geometry.h"
#include "route.h"

using namespace std;

void partition_fpga(const box &current, int level, bool horizontal, vector<box> &all_boxes)
{
	//printf("current %d->%d,%d->%d\n", bg::get<bg::min_corner, 0>(current), bg::get<bg::max_corner, 0>(current),
			//bg::get<bg::min_corner, 1>(current), bg::get<bg::max_corner, 1>(current));

	if (level == 0) {
		//printf("done\n");
		all_boxes.push_back(current);
	} else {
		assert(level > 0);
		
		box split_0;
		box split_1;

		if (horizontal) {
			bg::set<bg::min_corner, 0>(split_0, bg::get<bg::min_corner, 0>(current));
			bg::set<bg::max_corner, 0>(split_0, bg::get<bg::max_corner, 0>(current));
			bg::set<bg::min_corner, 1>(split_0, bg::get<bg::min_corner, 1>(current));
			bg::set<bg::max_corner, 1>(split_0, (bg::get<bg::min_corner, 1>(current)+bg::get<bg::max_corner, 1>(current))/2);

			bg::set<bg::min_corner, 0>(split_1, bg::get<bg::min_corner, 0>(current));
			bg::set<bg::max_corner, 0>(split_1, bg::get<bg::max_corner, 0>(current));
			bg::set<bg::min_corner, 1>(split_1, bg::get<bg::max_corner, 1>(split_0)+1);
			bg::set<bg::max_corner, 1>(split_1, bg::get<bg::max_corner, 1>(current));

			//printf("hsplit 0 %d->%d,%d->%d 1 %d->%d,%d->%d\n",
					//bg::get<bg::min_corner, 0>(split_0), bg::get<bg::max_corner, 0>(split_0),
					//bg::get<bg::min_corner, 1>(split_0), bg::get<bg::max_corner, 1>(split_0),
					//bg::get<bg::min_corner, 0>(split_1), bg::get<bg::max_corner, 0>(split_1),
					//bg::get<bg::min_corner, 1>(split_1), bg::get<bg::max_corner, 1>(split_1));

			assert((bg::get<bg::min_corner, 1>(split_0) < bg::get<bg::max_corner, 1>(split_0)));
			assert((bg::get<bg::min_corner, 1>(split_1) < bg::get<bg::max_corner, 1>(split_1)));
		} else {
			bg::set<bg::min_corner, 1>(split_0, bg::get<bg::min_corner, 1>(current));
			bg::set<bg::max_corner, 1>(split_0, bg::get<bg::max_corner, 1>(current));
			bg::set<bg::min_corner, 0>(split_0, bg::get<bg::min_corner, 0>(current));
			bg::set<bg::max_corner, 0>(split_0, (bg::get<bg::min_corner, 0>(current)+bg::get<bg::max_corner, 0>(current))/2);

			bg::set<bg::min_corner, 1>(split_1, bg::get<bg::min_corner, 1>(current));
			bg::set<bg::max_corner, 1>(split_1, bg::get<bg::max_corner, 1>(current));
			bg::set<bg::min_corner, 0>(split_1, bg::get<bg::max_corner, 0>(split_0)+1);
			bg::set<bg::max_corner, 0>(split_1, bg::get<bg::max_corner, 0>(current));

			//printf("vsplit 0 %d->%d,%d->%d 1 %d->%d,%d->%d\n",
					//bg::get<bg::min_corner, 0>(split_0), bg::get<bg::max_corner, 0>(split_0),
					//bg::get<bg::min_corner, 1>(split_0), bg::get<bg::max_corner, 1>(split_0),
					//bg::get<bg::min_corner, 0>(split_1), bg::get<bg::max_corner, 0>(split_1),
					//bg::get<bg::min_corner, 1>(split_1), bg::get<bg::max_corner, 1>(split_1));

			assert((bg::get<bg::min_corner, 0>(split_0) < bg::get<bg::max_corner, 0>(split_0)));
			assert((bg::get<bg::min_corner, 0>(split_1) < bg::get<bg::max_corner, 0>(split_1)));
		}

		partition_fpga(split_0, level-1, !horizontal, all_boxes); 
		partition_fpga(split_1, level-1, !horizontal, all_boxes); 
	}
}

int get_num_interpartition_nets(const vector<net_t> &nets, int num_partitions)
{
	vector<box> all_boxes;
	extern int nx, ny;
	box fpga_box;
	bg::set<bg::min_corner, 0>(fpga_box, 0);
	bg::set<bg::max_corner, 0>(fpga_box, nx+1);
	bg::set<bg::min_corner, 1>(fpga_box, 0);
	bg::set<bg::max_corner, 1>(fpga_box, ny+1);

	//printf("nx: %d ny: %d\n", nx, ny);

	partition_fpga(fpga_box, std::log2(num_partitions), true, all_boxes);

	//for (const auto &box : all_boxes) {
		//printf("box l %d r %d b %d t %d\n",
				//bg::get<bg::min_corner, 0>(box),
				//bg::get<bg::max_corner, 0>(box),
				//bg::get<bg::min_corner, 1>(box),
				//bg::get<bg::max_corner, 1>(box));
	//}

	int num_interpartition_nets = 0;
	for (const auto &net : nets) {
		int num_intersection = 0;
		//box net_bb;

		//bg::set<bg::min_corner, 0>(net_bb, net.bounding_box.xmin);
		//bg::set<bg::max_corner, 0>(net_bb, net.bounding_box.xmax);
		//bg::set<bg::min_corner, 1>(net_bb, net.bounding_box.ymin);
		//bg::set<bg::max_corner, 1>(net_bb, net.bounding_box.ymax);

		for (const auto &box : all_boxes) {
			if (bg::intersects(net.bounding_box, box)) {
				++num_intersection;
			}
		}
		if (num_intersection > 1) {
			++num_interpartition_nets;
		}
	}

	return num_interpartition_nets;
}
