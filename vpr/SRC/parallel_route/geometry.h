#ifndef GEOMETRY_H
#define GEOMETRY_H

#include <boost/numeric/interval.hpp>
#include <boost/geometry.hpp>
#include <boost/geometry/geometries/point.hpp>
#include <boost/geometry/geometries/box.hpp>
#include <boost/geometry/index/rtree.hpp>

namespace bg = boost::geometry;
namespace bgi = boost::geometry::index;

using point = typename bg::model::point<int, 2, bg::cs::cartesian>; 
using box = typename bg::model::box<point>;
using multi_point = typename bg::model::multi_point<point>;
using segment = typename bg::model::segment<point>;
using polygon = typename bg::model::polygon<point>;

template<typename T>
struct point_t {
	T x;
	T y;
};

typedef struct bounding_box_t {
	int xmin;
	int xmax;
	int ymin;
	int ymax;
	bool operator==(const bounding_box_t &other) const {
		return other.xmin == xmin && 
			other.xmax == xmax && 
			other.ymin == ymin && 
			other.ymax == ymax;
	}
} bounding_box_t;

template<typename PointA, typename PointB>
bounding_box_t get_bounding_box(const PointA &a, const PointB &b, int bb_factor)
{
	bounding_box_t bb;

	bb.xmin = std::min(a.x, b.x);
	bb.xmax = std::max(a.x, b.x);
	bb.ymin = std::min(a.y, b.y);
	bb.ymax = std::max(a.y, b.y);

	bb.xmin -= 1;
	bb.ymin -= 1;

	extern int nx;
	extern int ny;

	bb.xmin = std::max(bb.xmin - bb_factor, 0);
	bb.xmax = std::min(bb.xmax + bb_factor, nx + 1);
	bb.ymin = std::max(bb.ymin - bb_factor, 0);
	bb.ymax = std::min(bb.ymax + bb_factor, ny + 1);

	return bb;
}

template<typename BoundingBox>
int get_bounding_box_area(const BoundingBox &bb)
{
	assert(bb.xmax >= bb.xmin && bb.ymax >= bb.ymin);
	int area = (bb.xmax - bb.xmin + 1) * (bb.ymax - bb.ymin + 1);
	assert(area >= 0);
	return area;
}

template<typename BoundingBox>
bool box_overlap(const BoundingBox &box_a, const BoundingBox &box_b)
{
	using namespace boost::numeric;
	interval<int> a_hor(box_a.xmin, box_a.xmax);
	interval<int> b_hor(box_b.xmin, box_b.xmax);

	interval<int> a_vert(box_a.ymin, box_a.ymax);
	interval<int> b_vert(box_b.ymin, box_b.ymax);

	return overlap(a_hor, b_hor) && overlap(a_vert, b_vert);
}

void expand_and_clip(box &b, const std::pair<int, int> &xexpand, const std::pair<int, int> &yexpand, const std::pair<int, int> &max);

void make_l_segment(const point &source, const point &sink, bool xfirst, const std::pair<int, int> &xexpand, const std::pair<int, int> &yexpand, const std::pair<int, int> &max, std::vector<box> &boxes);

#endif
