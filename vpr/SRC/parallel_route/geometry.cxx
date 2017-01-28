#include "geometry.h"

using namespace std;

void expand_and_clip(box &b, const pair<int, int> &xexpand, const pair<int, int> &yexpand, const pair<int, int> &max)
{
	bg::set<bg::min_corner, 0>(b, std::max(0, bg::get<bg::min_corner, 0>(b) - xexpand.first));
	bg::set<bg::max_corner, 0>(b, std::min(max.first, bg::get<bg::max_corner, 0>(b) + xexpand.second));
	bg::set<bg::min_corner, 1>(b, std::max(0, bg::get<bg::min_corner, 1>(b) - yexpand.first));
	bg::set<bg::max_corner, 1>(b, std::min(max.second, bg::get<bg::max_corner, 1>(b) + yexpand.second));
}

void make_l_segment(const point &source, const point &sink, bool xfirst, const pair<int, int> &xexpand, const pair<int, int> &yexpand, const pair<int, int> &max, vector<box> &boxes)
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
	expand_and_clip(hb, xexpand, yexpand, max);
	//cout << "hseg after expand " << bg::dsv(hb) << endl << endl;
	boxes.push_back(hb);

	//cout << "vseg before correct " << bg::dsv(vb) << endl;
	bg::correct(vb);
	//cout << "vseg after correct " << bg::dsv(vb) << endl;
	expand_and_clip(vb, xexpand, yexpand, max);
	//cout << "vseg after expand " << bg::dsv(vb) << endl << endl;
	boxes.push_back(vb);
}
