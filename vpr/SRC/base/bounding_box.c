#include <assert.h>
#include "vpr_types.h"
#include "parallel_route_timing.h"

extern struct s_bb *route_bb;

int get_bounding_box_area(const struct s_bb *bb)
{
	int area = (bb->xmax - bb->xmin + 1) * (bb->ymax - bb->ymin + 1);
	assert(area >= 0);
	return area;
}

bool bounding_box_overlap(const struct s_bb &bb_a, const struct s_bb &bb_b)
{
	Interval<int> a_hor(bb_a.xmin, bb_a.xmax);
	Interval<int> b_hor(bb_b.xmin, bb_b.xmax);

	Interval<int> a_vert(bb_a.ymin, bb_a.ymax);
	Interval<int> b_vert(bb_b.ymin, bb_b.ymax);

	return a_hor.intersects(b_hor) && a_vert.intersects(b_vert);
}

int get_bounding_box_overlap_area(const struct s_bb &bb_a, const struct s_bb &bb_b)
{
	Interval<int> a_hor(bb_a.xmin, bb_a.xmax);
	Interval<int> b_hor(bb_b.xmin, bb_b.xmax);

	Interval<int> a_vert(bb_a.ymin, bb_a.ymax);
	Interval<int> b_vert(bb_b.ymin, bb_b.ymax);

	return (a_hor.intersect_size(b_hor) + 1) * (a_vert.intersect_size(b_vert) + 1);
}

bool bounding_box_overlap(int net_a, int net_b)
{
	Interval<int> a_hor(route_bb[net_a].xmin, route_bb[net_a].xmax);
	Interval<int> b_hor(route_bb[net_b].xmin, route_bb[net_b].xmax);

	Interval<int> a_vert(route_bb[net_a].ymin, route_bb[net_a].ymax);
	Interval<int> b_vert(route_bb[net_b].ymin, route_bb[net_b].ymax);

	return a_hor.intersects(b_hor) && a_vert.intersects(b_vert);
}
