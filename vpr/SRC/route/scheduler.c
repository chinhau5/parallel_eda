#include <boost/numeric/interval.hpp>
#include <boost/timer/timer.hpp>
#include <boost/geometry.hpp>
#include <boost/geometry/geometries/point.hpp>
#include <boost/geometry/geometries/box.hpp>
#include <boost/geometry/index/rtree.hpp>
#include "zlog.h"
#include "vpr_types.h"
#include "route.h"

namespace bg = boost::geometry;
namespace bgi = boost::geometry::index;

typedef bg::model::point<int, 2, bg::cs::cartesian> point;
typedef bg::model::box<point> box;
typedef std::pair<box, const net_t *> value;

/* TODO: check whether nets are global before routing */
using namespace boost::timer;

void verify_rtree_overlap(vector<net_t *> &nets)
{
	vector<vector<const net_t *>> overlap(nets.size());
	vector<vector<const net_t *>> non_overlap(nets.size());

	for (int i = 0; i < nets.size(); ++i) {
		for (int j = i + 1;  j < nets.size(); ++j) {
			if (box_overlap(nets[i]->current_bounding_box, nets[j]->current_bounding_box)) {
				overlap[nets[i]->current_local_id].push_back(nets[j]);
				overlap[nets[j]->current_local_id].push_back(nets[i]);
			} else {
				non_overlap[nets[i]->current_local_id].push_back(nets[j]);
				non_overlap[nets[j]->current_local_id].push_back(nets[i]);
			}
		}
	}

	for (int i = 0; i < nets.size(); ++i) {
		sort(overlap[nets[i]->current_local_id].begin(), overlap[nets[i]->current_local_id].end());
		sort(non_overlap[nets[i]->current_local_id].begin(), non_overlap[nets[i]->current_local_id].end());

		/*assert(i == nets[i]->current_local_id);*/

		sort(nets[i]->overlapping_nets_vec.begin(), nets[i]->overlapping_nets_vec.end());
		sort(nets[i]->non_overlapping_nets_vec.begin(), nets[i]->non_overlapping_nets_vec.end());

		assert(overlap[nets[i]->current_local_id] == nets[i]->overlapping_nets_vec);
		assert(non_overlap[nets[i]->current_local_id] == nets[i]->non_overlapping_nets_vec);
	}
}

void load_overlapping_nets_rtree(vector<net_t *> &nets)
{
	extern zlog_category_t *delta_log;

	cpu_timer timer;
	timer.start();

	vector<value> bulk;
	for (const auto &net : nets) {
		// create a box
		box b(point(net->current_bounding_box.xmin, net->current_bounding_box.ymin), point(net->current_bounding_box.xmax, net->current_bounding_box.ymax));
		// insert new value
		bulk.push_back(std::make_pair(b, net));
	}
	bgi::rtree< value, bgi::rstar<16> > t(bulk);

	timer.stop();
	cpu_times elapsed = timer.elapsed();
	/*zlog_info(delta_log, "Scheduling took %s\n", format(elapsed).c_str());*/
	zlog_info(delta_log, "R-tree creation took %g\n", elapsed.wall / 1e9);

	struct intersect_inserter {
		vector<bool> &overlapping_nets;
		intersect_inserter(vector<bool> &overlapping_nets) : overlapping_nets(overlapping_nets) {}

		intersect_inserter &operator=(const value &val)
		{
			assert(val.second->current_local_id < overlapping_nets.size());
			overlapping_nets[val.second->current_local_id] = true;
			return *this;
		}
		intersect_inserter &operator*() {
			return *this;
		}
		intersect_inserter &operator++() {
			return *this;
		}
		intersect_inserter &operator++(int i) {
			return *this;
		}
	};

	timer.start();
	for (auto &net : nets) {
		box query_box(point(net->current_bounding_box.xmin, net->current_bounding_box.ymin), point(net->current_bounding_box.xmax, net->current_bounding_box.ymax));
		std::vector<value> result_s;

		vector<bool> overlapping_nets(nets.size(), false);
		t.query(bgi::intersects(query_box), intersect_inserter(overlapping_nets));

		net->overlapping_nets_vec.clear();
		net->non_overlapping_nets_vec.clear();
		for (auto &other_net : nets) {
			//assert(v.second != &net);
			if (overlapping_nets[other_net->current_local_id]) {
				if (&other_net != &net) {
					net->overlapping_nets_vec.push_back(other_net);
				}
			} else {
				assert(&other_net != &net);
				net->non_overlapping_nets_vec.push_back(other_net);
			}
		}

		/*result_s.clear();*/
		/*timer.start();*/
		/*t.query(!bgi::intersects(query_box), std::back_inserter(result_s));*/
		/*timer.stop();*/
		/*elapsed = timer.elapsed();*/
		/*zlog_info(delta_log, "non intersect query took %g\n", elapsed.wall / 1e9);*/
		/*for (const auto &v : result_s) {*/
			/*assert(v.second != &net);*/
			/*net->non_overlapping_nets_vec.push_back(v.second);*/
		/*}*/

		assert(net->overlapping_nets_vec.size() + net->non_overlapping_nets_vec.size() == nets.size()-1);
	}
	timer.stop();
	elapsed = timer.elapsed();
	/*zlog_info(delta_log, "Scheduling took %s\n", format(elapsed).c_str());*/
	zlog_info(delta_log, "R-tree query took %g\n", elapsed.wall / 1e9);

	/* verify */
	verify_rtree_overlap(nets);
}

void one_d_independent()
{
}



