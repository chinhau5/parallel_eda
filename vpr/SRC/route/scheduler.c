#include <boost/numeric/interval.hpp>
#include <boost/timer/timer.hpp>
#include <boost/geometry.hpp>
#include <boost/geometry/geometries/point.hpp>
#include <boost/geometry/geometries/box.hpp>
#include <boost/geometry/index/rtree.hpp>
#include "zlog.h"
#include "vpr_types.h"
#include "route.h"
#include "scheduler.h"

namespace bg = boost::geometry;
namespace bgi = boost::geometry::index;

typedef bg::model::point<int, 2, bg::cs::cartesian> point;
typedef bg::model::box<point> box;
typedef std::pair<box, const net_t *> value;

/* TODO: check whether nets are global before routing */
using namespace boost::timer;

void sprintf_rr_node(int rr_node, char *buffer);

void verify_rtree_overlap(vector<net_t *> &nets)
{
	vector<vector<const net_t *>> overlap(nets.size());
	vector<vector<const net_t *>> non_overlap(nets.size());

	for (int i = 0; i < nets.size(); ++i) {
		for (int j = i + 1;  j < nets.size(); ++j) {
			if (box_overlap(nets[i]->sinks[nets[i]->current_sink_index].current_bounding_box, nets[j]->sinks[nets[j]->current_sink_index].current_bounding_box)) {
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
		box b(point(net->sinks[net->current_sink_index].current_bounding_box.xmin, net->sinks[net->current_sink_index].current_bounding_box.ymin), point(net->sinks[net->current_sink_index].current_bounding_box.xmax, net->sinks[net->current_sink_index].current_bounding_box.ymax));
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
		box query_box(point(net->sinks[net->current_sink_index].current_bounding_box.xmin, net->sinks[net->current_sink_index].current_bounding_box.ymin), point(net->sinks[net->current_sink_index].current_bounding_box.xmax, net->sinks[net->current_sink_index].current_bounding_box.ymax));
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

void verify_ind(const vector<const net_t *> &nets)
{
	for (int i = 0; i < nets.size(); ++i) {
		for (int j = i + 1;  j < nets.size(); ++j) {
			assert(!box_overlap(nets[i]->sinks[nets[i]->current_sink_index].current_bounding_box, nets[j]->sinks[nets[j]->current_sink_index].current_bounding_box));
		}
	}
}

void schedule_nets_ind(vector<net_t *> &nets, vector<vector<const net_t *>> &net_scheduled_at_time)
{
	cpu_timer timer;
	timer.start();

	vector<pair<const bounding_box_t *, net_t *>> bbs(nets.size());

	for (int i = 0; i < nets.size(); ++i) {
		bbs[i].first = &(nets[i]->sinks[nets[i]->current_sink_index].current_bounding_box);
		bbs[i].second = nets[i];
	}

	vector<bool> scheduled(nets.size(), false);
	net_scheduled_at_time.reserve(nets.size()); //worst case
	int num_scheduled_nets = 0;
	int time = 0;
	while (num_scheduled_nets < nets.size()) {
		vector<pair<const bounding_box_t *, net_t *> *> bbs_ptr;
		for (int i = 0; i < nets.size(); ++i) {
			if (!scheduled[nets[i]->current_local_id]) {
				bbs_ptr.push_back(&bbs[i]);
			}
		}
		vector<pair<const bounding_box_t *, net_t *> *> chosen;
		max_independent_rectangles(bbs_ptr, chosen);

		net_scheduled_at_time.push_back(vector<const net_t *>());
		for (const auto &c : chosen) {
			scheduled[c->second->current_local_id] = true;
			net_scheduled_at_time[time].push_back(c->second);
		}
		num_scheduled_nets += chosen.size();
		++time;
	}
	timer.stop();
	cpu_times elapsed = timer.elapsed();
	extern zlog_category_t *delta_log;
	zlog_info(delta_log, "Scheduling took %g\n", elapsed.wall / 1e9);
	assert(all_of(scheduled.begin(), scheduled.end(), [] (const bool &item) -> bool { return item; }));

	for (const auto &at_time : net_scheduled_at_time) {
		verify_ind(at_time);
	}
}

void create_simple_nets(vector<net_t> &nets) 
{
	/*for (const auto &net) {*/
	/*}*/
}

void schedule_nets_bounding_box(vector<net_t *> &nets, vector<pair<sink_t *, net_t *>> &scheduled_sinks)
{
	cpu_timer timer;
	timer.start();

	vector<pair<sink_t *, net_t *>> sinks;
	for (auto &net : nets) {
		int num_sinks = 0;
		for (auto &sink : net->sinks) {
			assert(sink.id < net->sink_routed.size() && sink.id >= 0);
			if (!net->sink_routed[sink.id]) {
				sinks.push_back(make_pair(&sink, net));
				++num_sinks;
			}
		}
		assert(num_sinks > 0);
	}

	std::sort(sinks.begin(), sinks.end(), [] (const pair<sink_t *, net_t *> &a, const pair<sink_t *, net_t *> &b) -> bool {
			return get_bounding_box_area(a.first->current_bounding_box) > get_bounding_box_area(b.first->current_bounding_box);
			});


	extern zlog_category_t *schedule_log;
	zlog_debug(schedule_log, "Largest %d sinks:\n", std::min(20, (int)sinks.size()));
	for (int i = 0; i < std::min(20, (int)sinks.size()); ++i) {
		char source[256];
		char sink[256];
		sprintf_rr_node(sinks[i].first->source.rr_node, source);
		sprintf_rr_node(sinks[i].first->rr_node, sink);
		zlog_debug(schedule_log, "Net %d Source: %s Sink: %s BB area: %d\n", sinks[i].second->vpr_id, source, sink, get_bounding_box_area(sinks[i].first->current_bounding_box));
	}

	/*printf("Num overlapping_nets:\n");*/
	/*for (const auto &net : nets) {*/
		/*printf("Net %d: %d\n", net->current_local_id, net->overlapping_nets_vec.size());*/
	/*}*/

	/*vector<bool> scheduled_sinks(sinks.size(), false);*/

	for (int i = 0; i < sinks.size(); ++i) {
		bool overlap_scheduled_sinks = any_of(scheduled_sinks.begin(), scheduled_sinks.end(), [&i, &sinks] (const pair<sink_t *, net_t *> &scheduled_sink) -> bool {
				return box_overlap(sinks[i].first->current_bounding_box, scheduled_sink.first->current_bounding_box);
				});

		if (!overlap_scheduled_sinks) {
			scheduled_sinks.push_back(sinks[i]);
			net_t *net;
			sink_t *sink;
			std::tie(sink, net)	= sinks[i];
			net->current_sink = sink;
		}
	}
	timer.stop();
	cpu_times elapsed = timer.elapsed();
	extern zlog_category_t *delta_log;
	zlog_info(delta_log, "Scheduling took %g\n", elapsed.wall / 1e9);
}
