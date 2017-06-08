#include "geometry.h"

using namespace std;

void box_to_polygon(const box &b, polygon &p)
{
	assert(p.outer().empty());
	assert(p.inners().empty());

	bg::append(p, b.min_corner());
	bg::append(p, point(bg::get<0>(b.min_corner()), bg::get<1>(b.max_corner())));
	bg::append(p, b.max_corner());
	bg::append(p, point(bg::get<0>(b.max_corner()), bg::get<1>(b.min_corner())));
	bg::append(p, b.min_corner());

	assert(bg::is_valid(p));
	assert(bg::area(b) == bg::area(p));
}

void clip(box &bb, const pair<int, int> &max)
{
	box clipped;
	bg::intersection(bg::make<box>(0, 0, max.first, max.second), bb, clipped);
	bb = clipped;
}

void expand(box &b, const pair<int, int> &xexpand, const pair<int, int> &yexpand)
{
	bg::set<bg::min_corner, 0>(b, bg::get<bg::min_corner, 0>(b) - xexpand.first);
	bg::set<bg::max_corner, 0>(b, bg::get<bg::max_corner, 0>(b) + xexpand.second);
	bg::set<bg::min_corner, 1>(b, bg::get<bg::min_corner, 1>(b) - yexpand.first);
	bg::set<bg::max_corner, 1>(b, bg::get<bg::max_corner, 1>(b) + yexpand.second);
}

void make_l_segment(const point &source, const point &sink, bool xfirst, const pair<int, int> &xexpand, const pair<int, int> &yexpand, const pair<int, int> &max, bool c, vector<box> &boxes)
{
	box hb;
	box vb;

	if (xfirst) {
		bg::assign_values(hb, bg::get<0>(source), bg::get<1>(source), bg::get<0>(sink), bg::get<1>(source));
		bg::assign_values(vb, bg::get<0>(sink), bg::get<1>(source), bg::get<0>(sink), bg::get<1>(sink));
	} else {
		bg::assign_values(hb, bg::get<0>(source), bg::get<1>(source), bg::get<0>(source), bg::get<1>(sink));
		bg::assign_values(vb, bg::get<0>(source), bg::get<1>(sink), bg::get<0>(sink), bg::get<1>(sink));
	}

	/* horizontal segment */
	//cout << "hseg before correct " << bg::dsv(hb) << endl;
	bg::correct(hb);
	//cout << "hseg after correct " << bg::dsv(hb) << endl;
	expand(hb, xexpand, yexpand);
	if (c) clip(hb, max);
	//cout << "hseg after expand " << bg::dsv(hb) << endl << endl;
	boxes.push_back(hb);

	//cout << "vseg before correct " << bg::dsv(vb) << endl;
	bg::correct(vb);
	//cout << "vseg after correct " << bg::dsv(vb) << endl;
	expand(vb, xexpand, yexpand);
	if (c) clip(vb, max);
	//cout << "vseg after expand " << bg::dsv(vb) << endl << endl;
	boxes.push_back(vb);
}
