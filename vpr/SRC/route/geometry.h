#ifndef GEOMETRY_H
#define GEOMETRY_H

#include <boost/geometry.hpp>
#include <boost/geometry/geometries/point.hpp>
#include <boost/geometry/geometries/box.hpp>
#include <boost/geometry/index/rtree.hpp>

namespace bg = boost::geometry;
namespace bgi = boost::geometry::index;

using point = bg::model::point<int, 2, bg::cs::cartesian>; 
using box = bg::model::box<point>;
using multi_point = bg::model::multi_point<point>;
using segment = bg::model::segment<point>;

template<typename T>
struct point_t {
	T x;
	T y;
};

#endif
