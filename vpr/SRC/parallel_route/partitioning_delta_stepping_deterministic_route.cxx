#include "pch.h"
#include <thread>
#include <condition_variable>
#include <sys/stat.h>

#include "vpr_types.h"

#include "log.h"
#include "graph.h"
#include "route.h"
#include "misr.h"
#include "utility.h"
#include "init.h"
#include "delta_stepping.h"
#include "route_tree.h"
#include "congestion.h"
#include "router.h"
#include "geometry.h"
#include "clock.h"
#include "new_partitioner.h"
#include "barrier.h"
#include "blockingconcurrentqueue.h"
#include "dijkstra.h"

using namespace std;

//#include <mlpack/methods/kmeans/kmeans.hpp>
//#include <mlpack/methods/kmeans/refined_start.hpp>

float get_timing_driven_expected_cost(const rr_node_property_t &current, const rr_node_property_t &target, float criticality_fac, float R_upstream);
void update_sink_criticalities(net_t &net, const t_net_timing &net_timing, const route_parameters_t &params);

using timer = std::chrono::high_resolution_clock;
using std::chrono::duration_cast;
using std::chrono::nanoseconds;

typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::directedS> OverlapGraph;
typedef boost::graph_traits<OverlapGraph>::vertex_descriptor vertex_descriptor;
typedef boost::graph_traits<OverlapGraph>::vertices_size_type vertices_size_type;
typedef boost::property_map<OverlapGraph, boost::vertex_index_t>::const_type vertex_index_map;

typedef struct new_virtual_net_t {
	int global_index;
	int index;
	net_t *net;
	vector<sink_t *> sinks;
	box bounding_box;
	bool l_segment;
	vector<box> bounding_boxes;
	point last_point;
	/* set when build the overlap graph */
	int parent;
	/* set by the scheduler */
	//int num_incoming_edges;
} new_virtual_net_t;

using rtree_value = pair<box, net_t *>;

struct net_to_rtree_item {
	rtree_value operator()(net_t *net) const
	{
		return make_pair(box(point(net->bounding_box.xmin, net->bounding_box.ymin), point(net->bounding_box.xmax, net->bounding_box.ymax)), net);
	}
};

using virtual_rtree_value = pair<box, new_virtual_net_t *>;

struct virtual_net_to_rtree_item {
	virtual_rtree_value operator()(new_virtual_net_t *vnet) const
	{
		return make_pair(vnet->bounding_box, vnet);
	}
};

using point_rtree_value = pair<point, sink_t *>;

struct sink_to_point_rtree_item {
	point_rtree_value operator()(sink_t &sink) const
	{
		return make_pair(point(sink.x, sink.y), &sink);
	}
};

struct rtree_value_equal {
	template<typename GeometryPtrPair>
	bool operator()(const GeometryPtrPair &a, const GeometryPtrPair &b) const
	{
		return bg::equals(a.first, b.first) && a.second == b.second;
	}
};

//template<typename Net>
//struct partition_t {
	//vector<Net *> nets;
	//partition_t *left;
	//partition_t *right;
//};

//void create_partitions_internal(const vector<net_t *> nets, vector<> &sorted_ver_edges, vector<> &sorted_hor_edges, bool horizontal_split)
//{
//}

//void create_partitions()
//{
//}
//

static char dirname[256];

void split_bb_7(net_t &net, const RRGraph &g, float bb_area_threshold, vector<new_virtual_net_t> &virtual_nets)
{
	vector<partition_t<int>> partitions;
	partition(net.sinks, [] (const sink_t &sink) -> point { return point(sink.x, sink.y); }, bb_area_threshold, partitions);
	//partition_t<int> p;
	//bg::assign_inverse(p.bounding_box);
	//for (int i = 0; i < net.sinks.size(); ++i) {
		//p.items.push_back(i);
		//bg::expand(p.bounding_box, point(net.sinks[i].x, net.sinks[i].y));
	//}
	//partitions.push_back(p);

	bgi::rtree<pair<point, int>, bgi::rstar<16>> tree;
	tree.insert(make_pair(point(net.source.x, net.source.y), -1));

	printf("current net %d num sinks %d\n", net.local_id, net.sinks.size());

	vector<bool> added(partitions.size(), false);
	vector<bool> inserted(net.sinks.size(), false);
	int vnet_index = 0;
	for (int i = 0; i < partitions.size(); ++i) {
		double min_distance = std::numeric_limits<double>::max();
		int nearest_source = -2;
		int nearest_partition = -1;

		for (int j = 0; j < partitions.size(); ++j) {
			if (added[j]) {
				continue;
			}

			vector<pair<point, int>> nearest;
			tree.query(bgi::nearest(partitions[j].bounding_box, 1), std::back_inserter(nearest));
			assert(nearest.size() == 1);

			double d = bg::comparable_distance(nearest[0].first, partitions[j].bounding_box);
			assert(d >= 0);

			printf("\tpartition %d distance %g\n", j, d);

			if (d < min_distance) {
				min_distance = d;
				nearest_source = nearest[0].second;
				nearest_partition = j;
			}
		}

		assert(min_distance != std::numeric_limits<double>::max());
		assert(nearest_source != -2);
		assert(nearest_partition != -1);

		assert(!added[nearest_partition]);
		added[nearest_partition] = true;

		printf("\tnearest partition %d\n", nearest_partition);

		point nearest_source_p;
		box source_box;
		if (nearest_source == -1) {
			/* source */
			bg::assign_values(nearest_source_p, net.source.x, net.source.y);
			assert(bg::get<0>(nearest_source_p) == net.source.x && bg::get<1>(nearest_source_p) == net.source.y);

			source_box = bg::make<box>(net.source.x, net.source.y, net.source.x, net.source.y);
			//const auto &source_props = get_vertex_props(g, net.source.rr_node);
			//source_box = bg::make<box>(source_props.xlow, source_props.ylow, source_props.xhigh, source_props.yhigh);

			//cout << "\t\texpanding " << bg::dsv(vnet.bounding_box) << " with source box " << bg::dsv(source_box) << endl;

			//bg::assign(vnet.last_point, nearest_source_p);
			cout << "\tuse net source " << bg::dsv(nearest_source_p) << " as source, box " << bg::dsv(source_box) << endl;
		} else {
			assert(nearest_source >= 0);

			bg::assign_values(nearest_source_p, net.sinks[nearest_source].x, net.sinks[nearest_source].y);
			assert(bg::get<0>(nearest_source_p) == net.sinks[nearest_source].x && bg::get<1>(nearest_source_p) == net.sinks[nearest_source].y);

			extern int nx, ny;
			int xmin = std::max(bg::get<0>(nearest_source_p)-1, 0);
			int xmax = std::min(bg::get<0>(nearest_source_p)+1, nx+2);

			int ymin = std::max(bg::get<1>(nearest_source_p)-1, 0);
			int ymax = std::min(bg::get<1>(nearest_source_p)+1, ny+2);

			source_box = bg::make<box>(xmin, ymin, xmax, ymax);

			cout << "\tuse sink " << bg::dsv(nearest_source_p) << " as source, box " << bg::dsv(source_box) << endl;
			
			//source_box = bg::make_inverse<box>();
			//bg::expand(source_box, p);

			//cout << "\t\texpanding " << bg::dsv(vnet.bounding_box) << " with sink box " << bg::dsv(source_box) << " of sink " << bg::dsv(p) << endl;

			//bg::assign(vnet.last_point, nearest_source_p);
		}
		//bg::assign(vnet.bounding_box, source_box);

		int nearest_sink = -1;
		min_distance = std::numeric_limits<double>::max();

		box sinks_bb = bg::make_inverse<box>();
		for (const auto &item : partitions[nearest_partition].items) {
			auto &sink = net.sinks[item];

			point p(sink.x, sink.y);

			tree.insert(make_pair(p, item));

			assert(!inserted[item]);
			inserted[item] = true;

			//const auto &sink_props = get_vertex_props(g, sink.rr_node);
			//bg::expand(sinks_bb, bg::make<box>(sink_props.xlow, sink_props.ylow, sink_props.xhigh, sink_props.yhigh));
			bg::expand(sinks_bb, p);

			double d = bg::comparable_distance(source_box, p);
			assert(d >= 0);
			if (d < min_distance) {
				min_distance = d;
				nearest_sink = item;
			}
		}

		assert(nearest_sink != -1);

		assert(bg::equals(sinks_bb, partitions[nearest_partition].bounding_box));

		point nearest_sink_p(net.sinks[nearest_sink].x, net.sinks[nearest_sink].y);

		cout << "\tneaerest sink " << nearest_sink << " " << bg::dsv(nearest_sink_p) << endl;

		box smallest_bb = bg::make_inverse<box>();
		bg::expand(smallest_bb, source_box);
		bg::expand(smallest_bb, nearest_sink_p);

		double smallest_bb_area = bg::area(smallest_bb);

		cout << "\tsmallest bb " << bg::dsv(smallest_bb) << " area " << smallest_bb_area << endl;

		bool has_l_segment = smallest_bb_area > bb_area_threshold;
		if (has_l_segment) {
			/* need to add one more L segment virtual net */

			new_virtual_net_t lvnet;
			lvnet.index = vnet_index++;
			lvnet.net = &net;
			lvnet.l_segment = true;

			lvnet.sinks.push_back(&net.sinks[nearest_sink]);

			extern int nx, ny;
			make_l_segment(nearest_source_p, nearest_sink_p, false, make_pair(1, 0), make_pair(1, 0), make_pair(nx+2, ny+2), lvnet.bounding_boxes);
			assert(lvnet.bounding_boxes.size() == 2);

			cout << "\thas l segment from " << bg::dsv(nearest_source_p) << " to " << bg::dsv(nearest_sink_p) << endl;
			cout << "\tl segment boxes " << bg::dsv(lvnet.bounding_boxes[0]) << " " << bg::dsv(lvnet.bounding_boxes[1]) << endl;

			box intersect_box;
			bg::intersection(lvnet.bounding_boxes[0], lvnet.bounding_boxes[1], intersect_box);
			cout << "\tintersect box " << bg::dsv(intersect_box) << endl;
			cout << "\treduced area from " << smallest_bb_area << " to " << bg::area(lvnet.bounding_boxes[0])+bg::area(lvnet.bounding_boxes[1])-bg::area(intersect_box) << endl;

			virtual_nets.emplace_back(lvnet);
		}

		/* TODO: should recursively create l segments instead of just for the first one */
		if (!has_l_segment || (has_l_segment && partitions[nearest_partition].items.size() > 1)) {
			new_virtual_net_t vnet;
			vnet.index = vnet_index++;
			vnet.net = &net;
			vnet.l_segment = false;

			//bg::assign_inverse(sinks_bb);
			for (const auto &item : partitions[nearest_partition].items) {
				if (has_l_segment && item == nearest_sink) {
					continue;
				}

				auto &sink = net.sinks[item];
				vnet.sinks.push_back(&sink);
				//bg::expand(sinks_bb, point(sink.x, sink.y));
			}
			assert(!vnet.sinks.empty());

			assert((has_l_segment && vnet.sinks.size() == partitions[nearest_partition].items.size()-1)
					|| (!has_l_segment && vnet.sinks.size() == partitions[nearest_partition].items.size()));

			//bg::expand(sinks_bb, source_box);
			//vnet.bounding_boxes.push_back(sinks_bb);
			if (!has_l_segment) {
				bg::expand(sinks_bb, source_box);
			}
			extern int nx, ny;
			expand_and_clip(sinks_bb, make_pair(1, 1), make_pair(1, 1), make_pair(nx+2, ny+2));

			if (has_l_segment) {
				printf("\tconnecting to rest of sinks with bb area %g\n", bg::area(sinks_bb));
			}

			vnet.bounding_boxes.push_back(sinks_bb);
			assert(vnet.bounding_boxes.size() == 1);

			virtual_nets.emplace_back(vnet);
		}
	}

	assert(tree.size() == net.sinks.size()+1);
}

void split_bb_6(net_t &net, const RRGraph &g, float bb_area_threshold, vector<new_virtual_net_t> &virtual_nets)
{
	vector<partition_t<int>> partitions;
	partition(net.sinks, [] (const sink_t &sink) -> point { return point(sink.x, sink.y); }, bb_area_threshold, partitions);
	//partition_t<int> p;
	//bg::assign_inverse(p.bounding_box);
	//for (int i = 0; i < net.sinks.size(); ++i) {
		//p.items.push_back(i);
		//bg::expand(p.bounding_box, point(net.sinks[i].x, net.sinks[i].y));
	//}
	//partitions.push_back(p);

	bgi::rtree<pair<point, int>, bgi::rstar<16>> tree;
	tree.insert(make_pair(point(net.source.x, net.source.y), -1));

	printf("current net %d\n", net.local_id);

	vector<bool> added(partitions.size(), false);
	vector<bool> inserted(net.sinks.size(), false);
	for (int i = 0; i < partitions.size(); ++i) {
		double min_distance = std::numeric_limits<double>::max();
		int nearest_point = -2;
		int nearest_partition = -1;

		printf("\tcurrent partition %d\n", i);

		for (int j = 0; j < partitions.size(); ++j) {
			if (added[j]) {
				continue;
			}

			vector<pair<point, int>> nearest;
			tree.query(bgi::nearest(partitions[j].bounding_box, 1), std::back_inserter(nearest));
			assert(nearest.size() == 1);

			double d = bg::comparable_distance(nearest[0].first, partitions[j].bounding_box);
			if (d < min_distance) {
				min_distance = d;
				nearest_point = nearest[0].second;
				nearest_partition = j;
			}
		}

		assert(min_distance != std::numeric_limits<double>::max());
		assert(nearest_point != -2);
		assert(nearest_partition != -1);

		assert(!added[nearest_partition]);
		added[nearest_partition] = true;

		new_virtual_net_t vnet;
		vnet.index = i;
		vnet.net = &net;

		box sinks_bb = bg::make_inverse<box>();
		for (const auto &item : partitions[nearest_partition].items) {
			auto &sink = net.sinks[item];
			vnet.sinks.push_back(&sink);

			assert(!inserted[item]);

			point p(sink.x, sink.y);

			tree.insert(make_pair(p, item));
			inserted[item] = true;

			//const auto &sink_props = get_vertex_props(g, sink.rr_node);
			//bg::expand(sinks_bb, bg::make<box>(sink_props.xlow, sink_props.ylow, sink_props.xhigh, sink_props.yhigh));
			bg::expand(sinks_bb, p);
		}

		assert(bg::equals(sinks_bb, partitions[nearest_partition].bounding_box));

		box source_box;
		if (nearest_point == -1) {
			/* source */
			source_box = bg::make<box>(net.source.x, net.source.y, net.source.x, net.source.y);
			//const auto &source_props = get_vertex_props(g, net.source.rr_node);
			//source_box = bg::make<box>(source_props.xlow, source_props.ylow, source_props.xhigh, source_props.yhigh);

			//cout << "\t\texpanding " << bg::dsv(vnet.bounding_box) << " with source box " << bg::dsv(source_box) << endl;

			bg::set<0>(vnet.last_point, net.source.x);
			bg::set<1>(vnet.last_point, net.source.y);
		} else {
			assert(nearest_point >= 0);

			point p(net.sinks[nearest_point].x, net.sinks[nearest_point].y);

			extern int nx, ny;
			int xmin = std::max(bg::get<0>(p)-1, 0);
			int xmax = std::min(bg::get<0>(p)+1, nx+2);

			int ymin = std::max(bg::get<1>(p)-1, 0);
			int ymax = std::min(bg::get<1>(p)+1, ny+2);

			source_box = bg::make<box>(xmin, ymin, xmax, ymax);
			
			//source_box = bg::make_inverse<box>();
			//bg::expand(source_box, p);

			//cout << "\t\texpanding " << bg::dsv(vnet.bounding_box) << " with sink box " << bg::dsv(source_box) << " of sink " << bg::dsv(p) << endl;

			bg::assign(vnet.last_point, p);
		}
		bg::expand(sinks_bb, source_box);
		//cout << "\t\tafter expansion with box " << bg::dsv(source_box) << " = " << bg::dsv(vnet.bounding_box) << endl;

		extern int nx, ny;
		expand_and_clip(sinks_bb, make_pair(1, 1), make_pair(1, 1), make_pair(nx+2, ny+2));

		vnet.bounding_boxes.push_back(sinks_bb);

		//cout << "\t\tafter further expansion " << bg::dsv(vnet.bounding_box) << endl;

		//for (auto &sink : vnet.sinks) {
			//sink->current_bounding_box.xmin = bg::get<bg::min_corner, 0>(vnet.bounding_box);
			//sink->current_bounding_box.ymin = bg::get<bg::min_corner, 1>(vnet.bounding_box);
			//sink->current_bounding_box.xmax = bg::get<bg::max_corner, 0>(vnet.bounding_box);
			//sink->current_bounding_box.ymax = bg::get<bg::max_corner, 1>(vnet.bounding_box);
		//}

		virtual_nets.emplace_back(vnet);
	}

	assert(tree.size() == net.sinks.size()+1);
}

void split_bb_5(net_t &net, float bb_area_threshold, vector<new_virtual_net_t> &virtual_nets)
{
	using point_f = typename bg::model::point<float, 2, bg::cs::cartesian>;

	point_f centroid(net.source.x, net.source.y);
	point_f total = centroid;
	int num_points = 1;

	sink_to_point_rtree_item to_point_rtree_item;

	bgi::rtree<point_rtree_value, bgi::rstar<16>, bgi::indexable<point_rtree_value>, rtree_value_equal> tree(net.sinks | boost::adaptors::transformed(to_point_rtree_item));

	//int num_acc = 0;

	//box previous_box;
	int debug = 0;

	printf("net id %d\n", net.local_id);
	printf("source %d %d\n", net.source.x, net.source.y);

	while (!tree.empty()) {
		new_virtual_net_t vnet;
		vnet.index = debug;
		vnet.net = &net;
		vnet.bounding_box = bg::make_inverse<box>();

		bg::assign(vnet.last_point, centroid);

		bg::expand(vnet.bounding_box, centroid);

		printf("centroid: %g %g\n", centroid.get<0>(), centroid.get<1>());

		//if ((num_acc % 2) == 0) {
			//total = centroid;
			//num_points = 1;
			//num_acc = 0;

			////printf("resetting total to %g %g\n", total.get<0>(), total.get<1>());
		//}

		while (bg::area(vnet.bounding_box) < bb_area_threshold && !tree.empty()) {
			vector<point_rtree_value> nearest;

			assert(tree.query(bgi::nearest(centroid, 1), std::back_inserter(nearest)) == 1);

			assert(nearest.size() == 1);
			assert(nearest[0].first.get<0>() == nearest[0].second->x && nearest[0].first.get<1>() == nearest[0].second->y);

			vnet.sinks.push_back(nearest[0].second);
			bg::expand(vnet.bounding_box, nearest[0].first);

			printf("\t%g %g + %d %d = ", total.get<0>(), total.get<1>(), nearest[0].first.get<0>(), nearest[0].first.get<1>());

			bg::add_point(total, nearest[0].first);
			++num_points;

			printf("%g %g\n", total.get<0>(), total.get<1>());

			/* recalc centroid */
			centroid = total;
			bg::divide_value(centroid, num_points); 

			printf("\tnew centroid %g %g\n", centroid.get<0>(), centroid.get<1>());

			assert(tree.remove(nearest[0]) == 1);
		}

		printf("bb %d,%d -> %d,%d\n", bg::get<0>(vnet.bounding_box.min_corner()), bg::get<1>(vnet.bounding_box.min_corner()),
				bg::get<0>(vnet.bounding_box.max_corner()), bg::get<1>(vnet.bounding_box.max_corner()));

		assert(!vnet.sinks.empty());

		//previous_box = vnet.bounding_box;
		
		extern int nx, ny;
		bg::set<0>(vnet.bounding_box.min_corner(), std::max(0, bg::get<0>(vnet.bounding_box.min_corner())-1));
		bg::set<1>(vnet.bounding_box.min_corner(), std::max(0, bg::get<1>(vnet.bounding_box.min_corner())-1));
		bg::set<0>(vnet.bounding_box.max_corner(), std::min(nx+2, bg::get<0>(vnet.bounding_box.max_corner())+1));
		bg::set<1>(vnet.bounding_box.max_corner(), std::min(ny+2, bg::get<1>(vnet.bounding_box.max_corner())+1));

		for (auto &sink : vnet.sinks) {
			sink->current_bounding_box.xmin = bg::get<0>(vnet.bounding_box.min_corner());
			sink->current_bounding_box.ymin = bg::get<1>(vnet.bounding_box.min_corner());
			sink->current_bounding_box.xmax = bg::get<0>(vnet.bounding_box.max_corner());
			sink->current_bounding_box.ymax = bg::get<1>(vnet.bounding_box.max_corner());
		}

		//++num_acc;
		++debug;

		virtual_nets.emplace_back(std::move(vnet));
	}
}

void split_bb_4(net_t &net, float bb_area_threshold, vector<new_virtual_net_t> &virtual_nets)
{
	using point_f = typename bg::model::point<float, 2, bg::cs::cartesian>;

	sink_to_point_rtree_item to_point_rtree_item;

	bgi::rtree<point_rtree_value, bgi::rstar<16>, bgi::indexable<point_rtree_value>, rtree_value_equal> tree(net.sinks | boost::adaptors::transformed(to_point_rtree_item));

	int debug = 0;

	//printf("source %d %d\n", net.source.x, net.source.y);

	while (!tree.empty()) {
		new_virtual_net_t vnet;
		vnet.net = &net;
		vnet.index = debug;
		//vnet.num_incoming_edges = 0;

		box source_box;

		if (debug > 0) {
			vector<point_rtree_value> nearest;
			assert(tree.query(bgi::nearest(virtual_nets.back().bounding_box, 1), std::back_inserter(nearest)) == 1);

			float min_distance = std::numeric_limits<float>::max();
			point nearest_point;
			for (const auto &sink : virtual_nets.back().sinks) {
				point p(sink->x, sink->y);
				float d = bg::comparable_distance(p, nearest[0].first);
				if (d < min_distance) {
					min_distance = d;
					nearest_point = p;
				}
			}

			extern int nx, ny;
			int xmin = std::max(bg::get<0>(nearest_point)-2, 0);
			int xmax = std::min(bg::get<0>(nearest_point)+2, nx+2);

			int ymin = std::max(bg::get<1>(nearest_point)-2, 0);
			int ymax = std::min(bg::get<1>(nearest_point)+2, ny+2);

			//printf("xmin %d xmax %d ymin %d ymax %d\n", xmin, xmax, ymin, ymax);

			source_box = bg::make<box>(xmin, ymin, xmax, ymax);

			vnet.last_point = nearest_point;
		} else {
			source_box = bg::make_inverse<box>();
			bg::expand(source_box, point(net.source.x, net.source.y));

			vnet.last_point = point(net.source.x, net.source.y);
		}

		//printf("source box %d %d, %d %d\n", bg::get<0>(source_box.min_corner()), bg::get<1>(source_box.min_corner()),
				//bg::get<0>(source_box.max_corner()), bg::get<1>(source_box.max_corner()));

		vnet.bounding_box = source_box;

		while (bg::area(vnet.bounding_box) < bb_area_threshold && !tree.empty()) {
			vector<point_rtree_value> nearest;

			assert(tree.query(bgi::nearest(source_box, 16), std::back_inserter(nearest)) > 0);

			std::sort(begin(nearest), end(nearest), [&source_box] (const point_rtree_value &a, const point_rtree_value &b) -> bool {
						return bg::comparable_distance(a.first, source_box) < bg::comparable_distance(b.first, source_box);
						});

			for (int i = 0; i < nearest.size() && bg::area(vnet.bounding_box) < bb_area_threshold; ++i) {
				assert(nearest[i].first.get<0>() == nearest[i].second->x && nearest[i].first.get<1>() == nearest[i].second->y);
				vnet.sinks.push_back(nearest[i].second);
				bg::expand(vnet.bounding_box, nearest[i].first);

				//printf("expanding with point %d %d\n", nearest[i].first.get<0>(), nearest[i].first.get<1>());

				assert(tree.remove(nearest[i]) == 1);
			}
		}

		assert(!vnet.sinks.empty());

		extern int nx, ny;
		bg::set<0>(vnet.bounding_box.min_corner(), std::max(0, bg::get<0>(vnet.bounding_box.min_corner())-1));
		bg::set<1>(vnet.bounding_box.min_corner(), std::max(0, bg::get<1>(vnet.bounding_box.min_corner())-1));
		bg::set<0>(vnet.bounding_box.max_corner(), std::min(nx+2, bg::get<0>(vnet.bounding_box.max_corner())+1));
		bg::set<1>(vnet.bounding_box.max_corner(), std::min(ny+2, bg::get<1>(vnet.bounding_box.max_corner())+1));
		//box expanded = vnet.bounding_box;
		//bg::set<0>(expanded.min_corner(), std::max(0, bg::get<0>(expanded.min_corner())-1));
		//bg::set<1>(expanded.min_corner(), std::max(0, bg::get<1>(expanded.min_corner())-1));
		//bg::set<0>(expanded.max_corner(), std::min(nx+2, bg::get<0>(expanded.max_corner())+1));
		//bg::set<1>(expanded.max_corner(), std::min(ny+2, bg::get<1>(expanded.max_corner())+1));

		for (auto &sink : vnet.sinks) {
			sink->current_bounding_box.xmin = bg::get<0>(vnet.bounding_box.min_corner());
			sink->current_bounding_box.ymin = bg::get<1>(vnet.bounding_box.min_corner());
			sink->current_bounding_box.xmax = bg::get<0>(vnet.bounding_box.max_corner());
			sink->current_bounding_box.ymax = bg::get<1>(vnet.bounding_box.max_corner());
			//sink->current_bounding_box.xmin = bg::get<0>(expanded.min_corner());
			//sink->current_bounding_box.ymin = bg::get<1>(expanded.min_corner());
			//sink->current_bounding_box.xmax = bg::get<0>(expanded.max_corner());
			//sink->current_bounding_box.ymax = bg::get<1>(expanded.max_corner());
		}

		virtual_nets.emplace_back(std::move(vnet));

		++debug;
	}
}

#ifdef OLD_CODE
void split_bb_3(net_t &net, float bb_area_threshold, vector<new_virtual_net_t> &virtual_nets)
{
	using point_f = typename bg::model::point<float, 2, bg::cs::cartesian>;

	sink_to_point_rtree_item to_point_rtree_item;

	bgi::rtree<point_rtree_value, bgi::rstar<16>, bgi::indexable<point_rtree_value>, rtree_value_equal> tree(net.sinks | boost::adaptors::transformed(to_point_rtree_item));

	box source_box = bg::make_inverse<box>();
	bg::expand(source_box, point(net.source.x, net.source.y));

	int debug = 0;

	printf("source %d %d\n", net.source.x, net.source.y);

	while (!tree.empty()) {
		new_virtual_net_t vnet;
		vnet.net = &net;

		vnet.bounding_box = source_box;

		printf("source box %d %d, %d %d\n", bg::get<0>(vnet.bounding_box.min_corner()), bg::get<1>(vnet.bounding_box.min_corner()),
				bg::get<0>(vnet.bounding_box.max_corner()), bg::get<1>(vnet.bounding_box.max_corner()));

		point last_point;

		while (bg::area(vnet.bounding_box) < bb_area_threshold && !tree.empty()) {
			vector<point_rtree_value> nearest;

			assert(tree.query(bgi::nearest(source_box, 16), std::back_inserter(nearest)) > 0);

			//assert(std::is_sorted(begin(nearest), end(nearest), [&source_box] (const point_rtree_value &a, const point_rtree_value &b) -> bool {
						//return bg::comparable_distance(a.first, source_box) > bg::comparable_distance(b.first, source_box);
						//}));

			std::sort(begin(nearest), end(nearest), [&source_box] (const point_rtree_value &a, const point_rtree_value &b) -> bool {
						return bg::comparable_distance(a.first, source_box) < bg::comparable_distance(b.first, source_box);
						});

			for (int i = 0; i < nearest.size() && bg::area(vnet.bounding_box) < bb_area_threshold; ++i) {
				assert(nearest[i].first.get<0>() == nearest[i].second->x && nearest[i].first.get<1>() == nearest[i].second->y);
				vnet.sinks.push_back(nearest[i].second);
				bg::expand(vnet.bounding_box, nearest[i].first);

				printf("expanding with point %d %d\n", nearest[i].first.get<0>(), nearest[i].first.get<1>());

				last_point = nearest[i].first;

				assert(tree.remove(nearest[i]) == 1);
			}
		}

		assert(!vnet.sinks.empty());
		vnet.last_point = last_point;

		/* TODO: clip the values here */
		source_box.min_corner() = last_point;
		bg::subtract_value(source_box.min_corner(), 2);

		source_box.max_corner() = last_point;
		bg::add_value(source_box.max_corner(), 2);

		printf("last point %d %d\n", bg::get<0>(last_point), bg::get<1>(last_point));
		printf("new source box %d %d, %d %d\n", bg::get<0>(source_box.min_corner()), bg::get<1>(source_box.min_corner()),
				bg::get<0>(source_box.max_corner()), bg::get<1>(source_box.max_corner()));

		assert(!vnet.sinks.empty());

		++debug;

		virtual_nets.emplace_back(std::move(vnet));
	}
}

void split_bb_2(net_t &net, float bb_area_threshold, vector<new_virtual_net_t> &virtual_nets)
{
	using point_f = typename bg::model::point<float, 2, bg::cs::cartesian>;

	sink_to_point_rtree_item to_point_rtree_item;

	bgi::rtree<point_rtree_value, bgi::rstar<16>, bgi::indexable<point_rtree_value>, rtree_value_equal> tree(net.sinks | boost::adaptors::transformed(to_point_rtree_item));

	box previous_box;
	int debug = 0;

	printf("source %d %d\n", net.source.x, net.source.y);

	while (!tree.empty()) {
		new_virtual_net_t vnet;
		vnet.net = &net;

		if (debug > 0) {
			vnet.bounding_box = previous_box;
		} else {
			vnet.bounding_box = bg::make_inverse<box>();
			bg::expand(vnet.bounding_box, point(net.source.x, net.source.y));
		}

		printf("box %d %d, %d %d\n", bg::get<0>(vnet.bounding_box.min_corner()), bg::get<1>(vnet.bounding_box.min_corner()),
				bg::get<0>(vnet.bounding_box.max_corner()), bg::get<1>(vnet.bounding_box.max_corner()));

		box source_box = vnet.bounding_box;
		point last_point;

		while (bg::area(vnet.bounding_box) < bb_area_threshold && !tree.empty()) {
			vector<point_rtree_value> nearest;

			assert(tree.query(bgi::nearest(source_box, 16), std::back_inserter(nearest)) > 0);

			//assert(std::is_sorted(begin(nearest), end(nearest), [&source_box] (const point_rtree_value &a, const point_rtree_value &b) -> bool {
						//return bg::comparable_distance(a.first, source_box) > bg::comparable_distance(b.first, source_box);
						//}));

			std::sort(begin(nearest), end(nearest), [&source_box] (const point_rtree_value &a, const point_rtree_value &b) -> bool {
						return bg::comparable_distance(a.first, source_box) < bg::comparable_distance(b.first, source_box);
						});

			for (int i = 0; i < nearest.size() && bg::area(vnet.bounding_box) < bb_area_threshold; ++i) {
				assert(nearest[i].first.get<0>() == nearest[i].second->x && nearest[i].first.get<1>() == nearest[i].second->y);
				vnet.sinks.push_back(nearest[i].second);
				bg::expand(vnet.bounding_box, nearest[i].first);

				printf("expanding with point %d %d\n", nearest[i].first.get<0>(), nearest[i].first.get<1>());

				last_point = nearest[i].first;

				assert(tree.remove(nearest[i]) == 1);
			}
		}

		bg::assign(previous_box.min_corner(), last_point);
		bg::centroid(vnet.bounding_box, previous_box.max_corner());
		bg::correct(previous_box);

		printf("previous box %d %d, %d %d\n", bg::get<0>(previous_box.min_corner()), bg::get<1>(previous_box.min_corner()),
				bg::get<0>(previous_box.max_corner()), bg::get<1>(previous_box.max_corner()));

		assert(!vnet.sinks.empty());

		++debug;

		virtual_nets.emplace_back(std::move(vnet));
	}
}

void split_bb(net_t &net, float bb_area_threshold, vector<new_virtual_net_t> &virtual_nets)
{
	using point_f = typename bg::model::point<float, 2, bg::cs::cartesian>;

	sink_to_point_rtree_item to_point_rtree_item;

	bgi::rtree<point_rtree_value, bgi::rstar<16>, bgi::indexable<point_rtree_value>, rtree_value_equal> tree(net.sinks | boost::adaptors::transformed(to_point_rtree_item));

	box previous_box;
	int debug = 0;

	//printf("source %d %d\n", net.source.x, net.source.y);

	while (!tree.empty()) {
		new_virtual_net_t vnet;
		vnet.net = &net;

		point_f centroid;

		if (debug > 0) {
			bg::centroid(previous_box, centroid);

			vnet.bounding_box = previous_box;
		} else {
			bg::assign_values(centroid, net.source.x, net.source.y);

			vnet.bounding_box = bg::make_inverse<box>();
			bg::expand(vnet.bounding_box, centroid);
		}

		printf("box %d %d, %d %d\n", bg::get<0>(vnet.bounding_box.min_corner()), bg::get<1>(vnet.bounding_box.min_corner()),
				bg::get<0>(vnet.bounding_box.max_corner()), bg::get<1>(vnet.bounding_box.max_corner()));
		printf("centroid: %g %g\n", centroid.get<0>(), centroid.get<1>());

		point last_point;
		point_f total = centroid;
		int num_points = 1;

		while (bg::area(vnet.bounding_box) < bb_area_threshold && !tree.empty()) {
			vector<point_rtree_value> nearest;

			assert(tree.query(bgi::nearest(centroid, 1), std::back_inserter(nearest)) == 1);

			assert(nearest.size() == 1);
			assert(nearest[0].first.get<0>() == nearest[0].second->x && nearest[0].first.get<1>() == nearest[0].second->y);

			vnet.sinks.push_back(nearest[0].second);
			bg::expand(vnet.bounding_box, nearest[0].first);

			printf("%g %g + %d %d = ", total.get<0>(), total.get<1>(), nearest[0].first.get<0>(), nearest[0].first.get<1>());

			bg::add_point(total, nearest[0].first);
			++num_points;

			printf("%g %g\n", total.get<0>(), total.get<1>());

			/* recalc centroid */
			centroid = total;
			bg::divide_value(centroid, num_points); 

			last_point = nearest[0].first;

			printf("new centroid %g %g\n", centroid.get<0>(), centroid.get<1>());

			assert(tree.remove(nearest[0]) == 1);
		}

		bg::assign(previous_box.min_corner(), last_point);
		bg::assign(previous_box.max_corner(), centroid);
		bg::correct(previous_box);

		printf("previous box %d %d, %d %d\n", bg::get<0>(previous_box.min_corner()), bg::get<1>(previous_box.min_corner()),
				bg::get<0>(previous_box.max_corner()), bg::get<1>(previous_box.max_corner()));

		assert(!vnet.sinks.empty());

		++debug;

		virtual_nets.emplace_back(std::move(vnet));
	}
}

void split(net_t &net, float bb_area_threshold, vector<new_virtual_net_t> &virtual_nets)
{
	using point_f = typename bg::model::point<float, 2, bg::cs::cartesian>;

	point_f centroid(net.source.x, net.source.y);
	point_f total;
	int num_points;

	sink_to_point_rtree_item to_point_rtree_item;

	bgi::rtree<point_rtree_value, bgi::rstar<16>, bgi::indexable<point_rtree_value>, rtree_value_equal> tree(net.sinks | boost::adaptors::transformed(to_point_rtree_item));

	int num_acc = 0;

	box previous_box;
	int debug = 0;

	//printf("source %d %d\n", net.source.x, net.source.y);

	while (!tree.empty()) {
		new_virtual_net_t vnet;
		vnet.net = &net;
		vnet.bounding_box = bg::make_inverse<box>();

		bg::expand(vnet.bounding_box, centroid);

		//printf("centroid: %g %g\n", centroid.get<0>(), centroid.get<1>());

		if ((num_acc % 2) == 0) {
			total = centroid;
			num_points = 1;
			num_acc = 0;

			//printf("resetting total to %g %g\n", total.get<0>(), total.get<1>());
		}

		while (bg::area(vnet.bounding_box) < bb_area_threshold && !tree.empty()) {
			vector<point_rtree_value> nearest;

			assert(tree.query(bgi::nearest(centroid, 1), std::back_inserter(nearest)) == 1);

			assert(nearest.size() == 1);
			assert(nearest[0].first.get<0>() == nearest[0].second->x && nearest[0].first.get<1>() == nearest[0].second->y);

			vnet.sinks.push_back(nearest[0].second);
			bg::expand(vnet.bounding_box, nearest[0].first);

			//printf("%g %g + %d %d = ", total.get<0>(), total.get<1>(), nearest[0].first.get<0>(), nearest[0].first.get<1>());

			bg::add_point(total, nearest[0].first);
			++num_points;

			//printf("%g %g\n", total.get<0>(), total.get<1>());

			/* recalc centroid */
			centroid = total;
			bg::divide_value(centroid, num_points); 

			//printf("new centroid %g %g\n", centroid.get<0>(), centroid.get<1>());

			assert(tree.remove(nearest[0]) == 1);
		}

		assert(!vnet.sinks.empty());

		if (debug > 0) {
			assert(bg::intersects(previous_box, vnet.bounding_box));
		}
		previous_box = vnet.bounding_box;

		++num_acc;
		++debug;

		virtual_nets.emplace_back(std::move(vnet));
	}
}
#endif

int create_virtual_nets(vector<net_t> &nets, const RRGraph &g, float threshold_scale, vector<vector<new_virtual_net_t>> &all_virtual_nets, vector<new_virtual_net_t *> &all_virtual_nets_ptr)
{
	extern int nx, ny;
	int fpga_area = (nx+2) * (ny+2);
	float bb_area_threshold = (float)fpga_area / threshold_scale;

	all_virtual_nets.resize(nets.size());
	for (int i = 0; i < nets.size(); ++i) {
		assert(all_virtual_nets[i].empty());
		split_bb_7(nets[i], g, bb_area_threshold, all_virtual_nets[i]);
	}

	assert(all_virtual_nets_ptr.empty());
	for (auto &virtual_nets : all_virtual_nets) {
		for (auto &virtual_net : virtual_nets) {
			all_virtual_nets_ptr.push_back(&virtual_net);
		}
	}

	int global_index = 0;
	for (auto &virtual_net : all_virtual_nets_ptr) {
		virtual_net->global_index = global_index++;
	}

	for (auto &virtual_net : all_virtual_nets_ptr) {
		for (auto &sink : virtual_net->sinks) {
			sink->vnet = virtual_net;
		}
	}

	int num_l_segments = 0;
	for (auto &virtual_net : all_virtual_nets_ptr) {
		if (virtual_net->l_segment) {
			++num_l_segments;
		}
	}
	printf("num l segments %d\n", num_l_segments);

	return global_index;
}

void best_case(vector<net_t> &nets, vector<vector<net_t *>> &phase_nets)
{
	phase_nets.emplace_back();
	for (auto &net : nets) {
		phase_nets.back().push_back(&net);
	}
}

int build_overlap_graph_4(vector<vector<new_virtual_net_t>> &all_virtual_nets, const vector<new_virtual_net_t *> &all_virtual_nets_ptr, vector<vector<int>> &overlap)
{
	overlap.resize(all_virtual_nets_ptr.size());
	for (int i = 0; i < all_virtual_nets_ptr.size(); ++i) {
		assert(overlap[i].empty());
	}

	int num_edges = 0;

	for (int i = 0; i < all_virtual_nets_ptr.size(); ++i) {
		for (int j = i+1; j < all_virtual_nets_ptr.size(); ++j) {
			if (all_virtual_nets_ptr[i]->net != all_virtual_nets_ptr[j]->net) {
				bool intersect = false;

				for (int k = 0; k < all_virtual_nets_ptr[i]->bounding_boxes.size() && !intersect; ++k) {
					for (int l = 0; l < all_virtual_nets_ptr[j]->bounding_boxes.size() && !intersect; ++l) {
						intersect = bg::intersects(all_virtual_nets_ptr[i]->bounding_boxes[k], all_virtual_nets_ptr[j]->bounding_boxes[l]);
					} 
				}

				if (intersect) {
					overlap[i].push_back(j);
					overlap[j].push_back(i);

					num_edges += 2;
				}
			}
		}
	}

	/* TODO: this is very pessimistic. we can remove some edges if we are sure that
	 * the bounding box overlaps with the existing route tree sufficiently */
	//for (auto &virtual_nets : all_virtual_nets) {
		//for (int i = 0; i < virtual_nets.size(); ++i) {
			//for (int j = 0; j < virtual_nets.size(); ++j) {
				//if (i != j) {
					//auto e = edge(virtual_nets[i].v, virtual_nets[j].v, g);
					//if (!e.second) {
						//add_edge(virtual_nets[i].v, virtual_nets[j].v, g);
					//}
				//}
			//}
		//}
	//}

	/* Don't need this here because we are checking the parent based on the virtual net index */
	//for (auto &virtual_nets : all_virtual_nets) {
		//for (int i = virtual_nets.size()-1; i > 0; --i) {
			//auto e = edge(virtual_nets[i].v, virtual_nets[i-1].v, g);
			//assert(!e.second);
			//add_edge(virtual_nets[i].v, virtual_nets[i-1].v, g);
		//}
	//}
	
	printf("starting to find parents\n");

	//for (auto &virtual_nets : all_virtual_nets) {
		//printf("virtual nets for net %d num sinks %d\n", virtual_nets[0].net->local_id, virtual_nets[0].net->sinks.size());

		//for (int i = 0; i < virtual_nets.size(); ++i) {
			//printf("vnet %d\n", i);

			//if (i == 0) {
				//virtual_nets[i].parent = -1;
			//} else {
				//int num_overlaps = 0;
				//double max_overlap_area = std::numeric_limits<double>::lowest();
				//int best = -1;
				//for (int j = 0; j < i; ++j) {
					//box intersection;
					//bool intersects = bg::intersection(virtual_nets[i].bounding_box, virtual_nets[j].bounding_box, intersection);
					//printf("bb0 %d,%d -> %d,%d bb1 %d,%d -> %d,%d\n", bg::get<0>(virtual_nets[i].bounding_box.min_corner()), bg::get<1>(virtual_nets[i].bounding_box.min_corner()),
							//bg::get<0>(virtual_nets[i].bounding_box.max_corner()), bg::get<1>(virtual_nets[i].bounding_box.max_corner()),
							//bg::get<0>(virtual_nets[j].bounding_box.min_corner()), bg::get<1>(virtual_nets[j].bounding_box.min_corner()),
							//bg::get<0>(virtual_nets[j].bounding_box.max_corner()), bg::get<1>(virtual_nets[j].bounding_box.max_corner()));
					//if (intersects) {
						//assert(bg::intersects(virtual_nets[i].bounding_box, virtual_nets[j].bounding_box));
						//++num_overlaps;
						//double area = bg::area(intersection);
						//printf("\toverlap between vnet %d and %d area %g max area %g\n", i, j, area, max_overlap_area);
						//if (area > max_overlap_area) {
							//max_overlap_area = area;
							//best = virtual_nets[j].global_index;
							//printf("\tsetting max area to %g and best vnet to %d\n", area, best);
						//}
					//} else {
						//assert(!bg::intersects(virtual_nets[i].bounding_box, virtual_nets[j].bounding_box));
						//printf("\tno overlap between vnet %d and %d\n", i, j);
					//}
				//}
				//assert(num_overlaps > 0);
				//assert(best != -1);
				//virtual_nets[i].parent = best;
				
				//overlap[best].push_back(virtual_nets[i].global_index);
			//}
		//}
	//}
	for (auto &virtual_nets : all_virtual_nets) {
		for (int i = 0; i < virtual_nets.size(); ++i) {
			if (i == 0) {
				virtual_nets[i].parent = -1;
			} else {
				virtual_nets[i].parent = virtual_nets[i-1].global_index;

				overlap[virtual_nets[i-1].global_index].push_back(virtual_nets[i].global_index);
				++num_edges;
			}
		}
	}
	//int num_orig_directed_edges = 0;
	//for (int i = 0; i < all_virtual_nets_ptr.size(); ++i) {
		//assert(all_virtual_nets_ptr[i]->v == i);
		//assert(all_virtual_nets_ptr[i]->global_index == i);

		//if (all_virtual_nets_ptr[i]->index > 0) {
			//overlap[all_virtual_nets_ptr[i]->parent].push_back(all_virtual_nets_ptr[i]->v);
			//++num_orig_directed_edges;
		//} 
	//}

	return num_edges;
}

int build_overlap_graph_3(vector<vector<new_virtual_net_t>> &all_virtual_nets, const vector<new_virtual_net_t *> &all_virtual_nets_ptr, vector<vector<int>> &overlap)
{
	overlap.resize(all_virtual_nets_ptr.size());
	for (int i = 0; i < all_virtual_nets_ptr.size(); ++i) {
		assert(overlap[i].empty());
	}

	int num_edges = 0;

	for (int i = 0; i < all_virtual_nets_ptr.size(); ++i) {
		for (int j = i+1; j < all_virtual_nets_ptr.size(); ++j) {
			if (all_virtual_nets_ptr[i]->net != all_virtual_nets_ptr[j]->net
					&& bg::intersects(all_virtual_nets_ptr[i]->bounding_box, all_virtual_nets_ptr[j]->bounding_box)) {
				overlap[i].push_back(j);
				overlap[j].push_back(i);

				num_edges += 2;
			}
		}
	}

	/* TODO: this is very pessimistic. we can remove some edges if we are sure that
	 * the bounding box overlaps with the existing route tree sufficiently */
	//for (auto &virtual_nets : all_virtual_nets) {
		//for (int i = 0; i < virtual_nets.size(); ++i) {
			//for (int j = 0; j < virtual_nets.size(); ++j) {
				//if (i != j) {
					//auto e = edge(virtual_nets[i].v, virtual_nets[j].v, g);
					//if (!e.second) {
						//add_edge(virtual_nets[i].v, virtual_nets[j].v, g);
					//}
				//}
			//}
		//}
	//}

	/* Don't need this here because we are checking the parent based on the virtual net index */
	//for (auto &virtual_nets : all_virtual_nets) {
		//for (int i = virtual_nets.size()-1; i > 0; --i) {
			//auto e = edge(virtual_nets[i].v, virtual_nets[i-1].v, g);
			//assert(!e.second);
			//add_edge(virtual_nets[i].v, virtual_nets[i-1].v, g);
		//}
	//}
	
	printf("starting to find parents\n");

	//for (auto &virtual_nets : all_virtual_nets) {
		//printf("virtual nets for net %d num sinks %d\n", virtual_nets[0].net->local_id, virtual_nets[0].net->sinks.size());

		//for (int i = 0; i < virtual_nets.size(); ++i) {
			//printf("vnet %d\n", i);

			//if (i == 0) {
				//virtual_nets[i].parent = -1;
			//} else {
				//int num_overlaps = 0;
				//double max_overlap_area = std::numeric_limits<double>::lowest();
				//int best = -1;
				//for (int j = 0; j < i; ++j) {
					//box intersection;
					//bool intersects = bg::intersection(virtual_nets[i].bounding_box, virtual_nets[j].bounding_box, intersection);
					//printf("bb0 %d,%d -> %d,%d bb1 %d,%d -> %d,%d\n", bg::get<0>(virtual_nets[i].bounding_box.min_corner()), bg::get<1>(virtual_nets[i].bounding_box.min_corner()),
							//bg::get<0>(virtual_nets[i].bounding_box.max_corner()), bg::get<1>(virtual_nets[i].bounding_box.max_corner()),
							//bg::get<0>(virtual_nets[j].bounding_box.min_corner()), bg::get<1>(virtual_nets[j].bounding_box.min_corner()),
							//bg::get<0>(virtual_nets[j].bounding_box.max_corner()), bg::get<1>(virtual_nets[j].bounding_box.max_corner()));
					//if (intersects) {
						//assert(bg::intersects(virtual_nets[i].bounding_box, virtual_nets[j].bounding_box));
						//++num_overlaps;
						//double area = bg::area(intersection);
						//printf("\toverlap between vnet %d and %d area %g max area %g\n", i, j, area, max_overlap_area);
						//if (area > max_overlap_area) {
							//max_overlap_area = area;
							//best = virtual_nets[j].global_index;
							//printf("\tsetting max area to %g and best vnet to %d\n", area, best);
						//}
					//} else {
						//assert(!bg::intersects(virtual_nets[i].bounding_box, virtual_nets[j].bounding_box));
						//printf("\tno overlap between vnet %d and %d\n", i, j);
					//}
				//}
				//assert(num_overlaps > 0);
				//assert(best != -1);
				//virtual_nets[i].parent = best;
				
				//overlap[best].push_back(virtual_nets[i].global_index);
			//}
		//}
	//}
	for (auto &virtual_nets : all_virtual_nets) {
		for (int i = 0; i < virtual_nets.size(); ++i) {
			if (i == 0) {
				virtual_nets[i].parent = -1;
			} else {
				virtual_nets[i].parent = virtual_nets[i-1].global_index;

				overlap[virtual_nets[i-1].global_index].push_back(virtual_nets[i].global_index);
				++num_edges;
			}
		}
	}
	//int num_orig_directed_edges = 0;
	//for (int i = 0; i < all_virtual_nets_ptr.size(); ++i) {
		//assert(all_virtual_nets_ptr[i]->v == i);
		//assert(all_virtual_nets_ptr[i]->global_index == i);

		//if (all_virtual_nets_ptr[i]->index > 0) {
			//overlap[all_virtual_nets_ptr[i]->parent].push_back(all_virtual_nets_ptr[i]->v);
			//++num_orig_directed_edges;
		//} 
	//}

	return num_edges;
}

int build_overlap_graph_2(vector<vector<new_virtual_net_t>> &all_virtual_nets, vector<vector<int>> &overlap, vector<new_virtual_net_t *> &all_virtual_nets_ptr)
{
	for (auto &virtual_nets : all_virtual_nets) {
		for (auto &virtual_net : virtual_nets) {
			all_virtual_nets_ptr.push_back(&virtual_net);
		}
	}

	virtual_net_to_rtree_item to_rtree_item;
	bgi::rtree<virtual_rtree_value, bgi::rstar<64>, bgi::indexable<virtual_rtree_value>, rtree_value_equal> tree(all_virtual_nets_ptr | boost::adaptors::transformed(to_rtree_item));

	overlap.resize(all_virtual_nets_ptr.size());

	int num_edges = 0;

	for (const auto &current : all_virtual_nets_ptr) {
		//vector<virtual_rtree_value> overlapping_nets;
		//tree.query(bgi::intersects(current->bounding_box), std::back_inserter(overlapping_nets));

		assert(overlap[current->global_index].empty());

		int num_equal = 0;
		//for (const auto &overlapping : overlapping_nets) {
			//assert(bg::equals(overlapping.first, overlapping.second->bounding_box));
			//assert(bg::intersects(overlapping.first, current->bounding_box));

		for (auto iter = tree.qbegin(bgi::intersects(current->bounding_box));
				iter != tree.qend();
				++iter) {
			const auto &overlapping = *iter;

			if (overlapping.second != current && overlapping.second->net != current->net) {
				//add_edge(current->v, overlapping.second->v, g);
				overlap[current->global_index].push_back(overlapping.second->global_index);

				++num_edges;
			}

			if (overlapping.second == current) {
				++num_equal;
			}
		}
		assert(num_equal == 1);
	}

	/* TODO: this is very pessimistic. we can remove some edges if we are sure that
	 * the bounding box overlaps with the existing route tree sufficiently */
	//for (auto &virtual_nets : all_virtual_nets) {
		//for (int i = 0; i < virtual_nets.size(); ++i) {
			//for (int j = 0; j < virtual_nets.size(); ++j) {
				//if (i != j) {
					//auto e = edge(virtual_nets[i].v, virtual_nets[j].v, g);
					//if (!e.second) {
						//add_edge(virtual_nets[i].v, virtual_nets[j].v, g);
					//}
				//}
			//}
		//}
	//}

	/* Don't need this here because we are checking the parent based on the virtual net index */
	//for (auto &virtual_nets : all_virtual_nets) {
		//for (int i = virtual_nets.size()-1; i > 0; --i) {
			//auto e = edge(virtual_nets[i].v, virtual_nets[i-1].v, g);
			//assert(!e.second);
			//add_edge(virtual_nets[i].v, virtual_nets[i-1].v, g);
		//}
	//}
	
	for (auto &virtual_nets : all_virtual_nets) {
		for (int i = 0; i < virtual_nets.size(); ++i) {
			if (i == 0) {
				virtual_nets[i].parent = -1;
			} else {
				virtual_nets[i].parent = virtual_nets[i-1].global_index;

				overlap[virtual_nets[i-1].global_index].push_back(virtual_nets[i].global_index);
				++num_edges;
			}
		}
	}

	return num_edges;
}

void build_overlap_graph(vector<vector<new_virtual_net_t>> &all_virtual_nets, OverlapGraph &g, vector<new_virtual_net_t *> &all_virtual_nets_ptr)
{
	for (auto &virtual_nets : all_virtual_nets) {
		for (auto &virtual_net : virtual_nets) {
			add_vertex(g);

			all_virtual_nets_ptr.push_back(&virtual_net);
		}
	}

	virtual_net_to_rtree_item to_rtree_item;
	bgi::rtree<virtual_rtree_value, bgi::rstar<64>, bgi::indexable<virtual_rtree_value>, rtree_value_equal> tree(all_virtual_nets_ptr | boost::adaptors::transformed(to_rtree_item));

	//vector<vector<int>> overlap(num_vertices(g));

	for (const auto &current : all_virtual_nets_ptr) {
		vector<virtual_rtree_value> overlapping_nets;
		tree.query(bgi::intersects(current->bounding_box), std::back_inserter(overlapping_nets));

		int num_equal = 0;
		for (const auto &overlapping_net : overlapping_nets) {
			assert(bg::equals(overlapping_net.first, overlapping_net.second->bounding_box));
			assert(bg::intersects(overlapping_net.first, current->bounding_box));

		//for (auto iter = tree.qbegin(bgi::intersects(current->bounding_box));
				//iter != tree.qend();
				//++iter) {
			//const auto &overlapping_net = *iter;

			if (overlapping_net.second != current && overlapping_net.second->net != current->net) {
				add_edge(current->global_index, overlapping_net.second->global_index, g);
				//overlap[current->v].push_back(overlapping_net.second->v);
			}

			if (overlapping_net.second == current) {
				++num_equal;
			}
		}
		assert(num_equal == 1);
	}

	/* TODO: this is very pessimistic. we can remove some edges if we are sure that
	 * the bounding box overlaps with the existing route tree sufficiently */
	//for (auto &virtual_nets : all_virtual_nets) {
		//for (int i = 0; i < virtual_nets.size(); ++i) {
			//for (int j = 0; j < virtual_nets.size(); ++j) {
				//if (i != j) {
					//auto e = edge(virtual_nets[i].v, virtual_nets[j].v, g);
					//if (!e.second) {
						//add_edge(virtual_nets[i].v, virtual_nets[j].v, g);
					//}
				//}
			//}
		//}
	//}

	/* Don't need this here because we are checking the parent based on the virtual net index */
	//for (auto &virtual_nets : all_virtual_nets) {
		//for (int i = virtual_nets.size()-1; i > 0; --i) {
			//auto e = edge(virtual_nets[i].v, virtual_nets[i-1].v, g);
			//assert(!e.second);
			//add_edge(virtual_nets[i].v, virtual_nets[i-1].v, g);
		//}
	//}
	
	for (auto &virtual_nets : all_virtual_nets) {
		for (int i = 0; i < virtual_nets.size(); ++i) {
			if (i == 0) {
				virtual_nets[i].parent = boost::graph_traits<OverlapGraph>::null_vertex();
			} else {
				virtual_nets[i].parent = virtual_nets[i-1].global_index;
			}
		}
	}
}

#ifdef OLD_CODE
template <class VirtualNetsPtr, class Load, class OrderPA, class ColorMap>
	typename boost::property_traits<ColorMap>::value_type
custom_vertex_coloring_3(const vector<vector<int>> &G, const VirtualNetsPtr &all_virtual_nets_ptr, const Load &load, int num_threads, OrderPA &order, 
		ColorMap &color)
{
	typedef typename boost::property_traits<ColorMap>::value_type size_type;

	size_type max_color = 0;
	const size_type V = G.size();

	// We need to keep track of which colors are used by
	// adjacent vertices. We do this by marking the colors
	// that are used. The mark array contains the mark
	// for each color. The length of mark is the
	// number of vertices since the maximum possible number of colors
	// is the number of vertices.

	//Initialize colors 
	for (int v = 0; v < V; ++v)
		put(color, v, V-1);

	bool has_uncolored = true;
	int num_rounds = 0;

	using LoadValueType = decltype(load(0));
	using HeapElement = pair<LoadValueType, int>;
	vector<std::priority_queue<HeapElement, vector<HeapElement>, std::greater<HeapElement>>> heaps(V);
	vector<LoadValueType> max_load(V, std::numeric_limits<LoadValueType>::min());

	while (has_uncolored) {
		std::vector<size_type> mark(V, 
				std::numeric_limits<size_type>::max BOOST_PREVENT_MACRO_SUBSTITUTION());
		//Determine the color for every vertex one by one
		has_uncolored = false;
		for ( size_type i = 0; i < V; i++) {
			int current = get(order,i);

			//printf("current %d\n", current);

			/* we are done coloring this node */
			if (get(color, current) != V-1) {
				//printf("\tdone coloring\n");
				continue;
			}

			/* don't color this node if its parent is uncolored */
			int parent = all_virtual_nets_ptr[current]->parent;
			if (parent != -1 
					&& all_virtual_nets_ptr[current]->index > all_virtual_nets_ptr[parent]->index
					&& get(color, parent) == V-1) {
				has_uncolored = true;
				//printf("\tparent uncolored\n");
				continue;
			}

			//Mark the colors of vertices adjacent to current.
			//i can be the value for marking since i increases successively
			for (int j = 0; j < G[current].size(); ++j) {
				mark[get(color, G[current][j])] = i; 
			}

			//Next step is to assign the smallest un-marked color
			//to the current vertex.
			size_type j;
			if (parent != -1) {
				j = get(color, parent)+1;
				assert(j <= max_color);
			} else {
				j = 0;
			}

			//Scan through all useable colors, find the smallest possible
			//color that is not used by neighbors.  Note that if mark[j]
			//is equal to i, color j is used by one of the current vertex's
			//neighbors.
			int best_color = -1;
			LoadValueType min_load_delta = std::numeric_limits<LoadValueType>::max();

			LoadValueType l = load(current);
			//printf("\tnew load %d\n", l);

			for (; j < max_color; ++j) {
				if (mark[j] != i) {
					assert(!heaps[j].empty());
					assert(max_load[j] != std::numeric_limits<LoadValueType>::min());

					/* choose the best color here */
					const auto &min_load = heaps[j].top();

					LoadValueType delta = min_load.first + l - max_load[j];

					//printf("\tc %d max load %d min load %d thread %d delta %d\n", j, max_load[j], min_load.first, min_load.second, delta);

					if (delta < min_load_delta) {
						best_color = j;
						min_load_delta = delta;
					}
				}
			}

			if (best_color != -1) {
				j = best_color;

				auto min_load = heaps[j].top();
				heaps[j].pop();

				min_load.first += l;

				heaps[j].push(min_load);

				//printf("\tbest color %d pushing min load %d thread %d to heap. old max load %d ", j, min_load.first, min_load.second, max_load[j]);

				max_load[j] = std::max(max_load[j], min_load.first);

				//printf("new max load %d\n", max_load[j]);
			} else {
				assert(j == max_color);

				assert(heaps[j].empty());

				heaps[j].emplace(l, 0);
				for (int thread = 1; thread < num_threads; ++thread) {
					heaps[j].emplace(0, thread);
				}

				max_load[j] = l;

				//printf("\tnew color %d max load %d\n", j, max_load[j]);

				++max_color;
			}

			//while ( max_color <= j)  //All colors are used up. Add one more color
				//++max_color;

			//At this point, j is the smallest possible color
			put(color, current, j);  //Save the color of vertex current
		}

		++num_rounds;
	}

	printf("num rounds %d\n", num_rounds);

	return max_color;
}

template <class OrderPA, class ColorMap, class VirtualNetsPtr>
	typename boost::property_traits<ColorMap>::value_type
custom_vertex_coloring_2(const vector<vector<int>> &G, OrderPA &order, 
		ColorMap &color, const VirtualNetsPtr &all_virtual_nets_ptr)
{
	typedef typename boost::property_traits<ColorMap>::value_type size_type;

	size_type max_color = 0;
	const size_type V = G.size();

	// We need to keep track of which colors are used by
	// adjacent vertices. We do this by marking the colors
	// that are used. The mark array contains the mark
	// for each color. The length of mark is the
	// number of vertices since the maximum possible number of colors
	// is the number of vertices.

	//Initialize colors 
	for (int v = 0; v < V; ++v)
		put(color, v, V-1);

	bool has_uncolored = true;
	int num_rounds = 0;

	while (has_uncolored) {
		std::vector<size_type> mark(V, 
				std::numeric_limits<size_type>::max BOOST_PREVENT_MACRO_SUBSTITUTION());
		//Determine the color for every vertex one by one
		has_uncolored = false;
		for ( size_type i = 0; i < V; i++) {
			int current = get(order,i);

			/* we are done coloring this node */
			if (get(color, current) != V-1) {
				continue;
			}

			/* don't color this node if its parent is uncolored */
			int parent = all_virtual_nets_ptr[current]->parent;
			if (parent != -1 
					&& all_virtual_nets_ptr[current]->index > all_virtual_nets_ptr[parent]->index
					&& get(color, parent) == V-1) {
				has_uncolored = true;
				continue;
			}

			//Mark the colors of vertices adjacent to current.
			//i can be the value for marking since i increases successively
			for (int j = 0; j < G[current].size(); ++j) {
				mark[get(color, G[current][j])] = i; 
			}

			//Next step is to assign the smallest un-marked color
			//to the current vertex.
			size_type j;
			if (parent != -1) {
				j = get(color, parent)+1;
				//j = get(color, parent);
				assert(j <= max_color);
				/* just to make sure we don't use our parents color */
				//mark[j] = i;
			} else {
				j = 0;
			}

			//Scan through all useable colors, find the smallest possible
			//color that is not used by neighbors.  Note that if mark[j]
			//is equal to i, color j is used by one of the current vertex's
			//neighbors.
			while ( j < max_color && mark[j] == i ) 
				++j;

			assert(j <= max_color);
			if (j == max_color) {
				++max_color;
			}
			//while ( max_color <= j)  //All colors are used up. Add one more color
				//++max_color;

			//At this point, j is the smallest possible color
			put(color, current, j);  //Save the color of vertex current
		}

		++num_rounds;
	}

	printf("num rounds %d\n", num_rounds);

	return max_color;
}

template <class VertexListGraph, class OrderPA, class ColorMap, class VirtualNetsPtr>
	typename boost::property_traits<ColorMap>::value_type
custom_vertex_coloring(const VertexListGraph& G, OrderPA &order, 
		ColorMap &color, const VirtualNetsPtr &all_virtual_nets_ptr)
{
	typedef boost::graph_traits<VertexListGraph> GraphTraits;
	typedef typename GraphTraits::vertex_descriptor Vertex;
	typedef typename boost::property_traits<ColorMap>::value_type size_type;

	size_type max_color = 0;
	const size_type V = num_vertices(G);

	// We need to keep track of which colors are used by
	// adjacent vertices. We do this by marking the colors
	// that are used. The mark array contains the mark
	// for each color. The length of mark is the
	// number of vertices since the maximum possible number of colors
	// is the number of vertices.

	//Initialize colors 
	typename GraphTraits::vertex_iterator v, vend;
	for (boost::tie(v, vend) = vertices(G); v != vend; ++v)
		put(color, *v, V-1);

	bool has_uncolored = true;
	int num_rounds = 0;

	while (has_uncolored) {
		std::vector<size_type> mark(V, 
				std::numeric_limits<size_type>::max BOOST_PREVENT_MACRO_SUBSTITUTION());
		//Determine the color for every vertex one by one
		has_uncolored = false;
		for ( size_type i = 0; i < V; i++) {
			Vertex current = get(order,i);

			/* we are done coloring this node */
			if (get(color, current) != V-1) {
				continue;
			}

			/* don't color this node if its parent is uncolored */
			Vertex parent = all_virtual_nets_ptr[current]->parent;
			if (parent != GraphTraits::null_vertex()
					&& all_virtual_nets_ptr[current]->index > all_virtual_nets_ptr[parent]->index
					&& get(color, parent) == V-1) {
				has_uncolored = true;
				continue;
			}

			//Mark the colors of vertices adjacent to current.
			//i can be the value for marking since i increases successively
			typename GraphTraits::adjacency_iterator v, vend;
			for (boost::tie(v,vend) = adjacent_vertices(current, G); v != vend; ++v) {
				mark[get(color,*v)] = i; 
			}

			//Next step is to assign the smallest un-marked color
			//to the current vertex.
			size_type j;
			if (parent != GraphTraits::null_vertex()) {
				j = get(color, parent)+1;
			} else {
				j = 0;
			}

			//Scan through all useable colors, find the smallest possible
			//color that is not used by neighbors.  Note that if mark[j]
			//is equal to i, color j is used by one of the current vertex's
			//neighbors.
			while ( j < max_color && mark[j] == i ) 
				++j;

			max_color = std::max(j+1, max_color);
			//while ( max_color <= j)  //All colors are used up. Add one more color
				//++max_color;

			//At this point, j is the smallest possible color
			put(color, current, j);  //Save the color of vertex current
		}

		++num_rounds;
	}

	printf("num rounds %d\n", num_rounds);

	return max_color;
}
#endif

template <class Order, class Degree, 
		 class Marker, class BucketSorter>
void 
custom_smallest_last_vertex_ordering(const vector<vector<int>> & G, Order order, 
		Degree degree, Marker marker,
		BucketSorter& degree_buckets) {
	typedef std::size_t size_type;

	const size_type num = G.size();

	for (int v = 0; v < G.size(); ++v) {
		put(marker, v, num);
		put(degree, v, G[v].size());
		degree_buckets.push(v);
	}

	size_type minimum_degree = 0;
	size_type current_order = num - 1;

	while ( 1 ) {
		typedef typename BucketSorter::stack MDStack;
		MDStack minimum_degree_stack = degree_buckets[minimum_degree];
		while (minimum_degree_stack.empty())
			minimum_degree_stack = degree_buckets[++minimum_degree];

		int node = minimum_degree_stack.top();
		put(order, current_order, node);

		if ( current_order == 0 ) //find all vertices
			break;

		minimum_degree_stack.pop();
		put(marker, node, 0); //node has been ordered.

		for (int v = 0; v < G[node].size(); ++v) {
			int vn = G[node][v];
			if ( get(marker, vn) > current_order ) { //vn is unordered vertex
				put(marker, vn, current_order);  //mark the columns adjacent to node

				//delete vn from the bucket sorter         
				degree_buckets.remove(vn);

				//It is possible minimum degree goes down
				//Here we keep tracking it.
				put(degree, vn, get(degree, vn) - 1); 
				minimum_degree = std::min(minimum_degree, get(degree, vn)); 

				//reinsert vn in the bucket sorter with the new degree
				degree_buckets.push(vn);
			}
		}

		current_order--;
	}

	//at this point, order[i] = v_i;
}

template <class Order, class Degree, class Marker>
void 
custom_smallest_last_vertex_ordering(const vector<vector<int>> & G, Order order, 
		Degree degree, Marker marker) {
	typedef std::size_t size_type;

	const size_type num = G.size();

	typedef boost::bucket_sorter<size_type, int, Degree, boost::identity_property_map> BucketSorter;

	BucketSorter degree_bucket_sorter(num, num, degree);

	custom_smallest_last_vertex_ordering(G, order, degree, marker, degree_bucket_sorter);
}

template <class Order>
void 
custom_smallest_last_vertex_ordering(const vector<vector<int>> & G, Order order) {
	custom_smallest_last_vertex_ordering(G, order,
			make_shared_array_property_map(G.size(), std::size_t(0), boost::identity_property_map()),
			make_shared_array_property_map(G.size(), 0, boost::identity_property_map()));
}

void build_and_order(vector<vector<new_virtual_net_t>> &all_virtual_nets, const vector<new_virtual_net_t *> &all_virtual_nets_ptr,
		vector<vector<int>> &overlap, vector<vertex_descriptor> &order)
{
	auto build_start = timer::now();

	int num_edges = build_overlap_graph_4(all_virtual_nets, all_virtual_nets_ptr, overlap);

	auto build_time = timer::now()-build_start;

	printf("overlap graph build time %g\n", duration_cast<nanoseconds>(build_time).count() / 1e9);
	printf("overlap num vertices %d num edges %d\n", overlap.size(), num_edges);

	//auto ep = edges(g);
	//for (auto ei = ep.first; ei != ep.second; ++ei) {
		//printf("edge %d -> %d\n", source(*ei, g), target(*ei, g));
	//}

	order.resize(overlap.size());
	auto order_map = make_iterator_property_map(begin(order), boost::identity_property_map());

#ifdef LARGEST_FIRST
	for (int i = 0; i < overlap.size(); ++i) {
		order[i] = i;
	}

	std::sort(begin(order), end(order), [&] (int a, int b) -> bool {
			assert(a == all_virtual_nets_ptr[a]->v); 
			assert(b == all_virtual_nets_ptr[b]->v); 

			/* BUG: the relation is not transitive */
			//if (all_virtual_nets_ptr[a]->net == all_virtual_nets_ptr[b]->net) {
			//return all_virtual_nets_ptr[a]->index < all_virtual_nets_ptr[b]->index;
			//} else {
			//return out_degree(a, g) > out_degree(b, g);
			//}
			//return all_virtual_nets_ptr[a]->index < all_virtual_nets_ptr[b]->index || (all_virtual_nets_ptr[a]->index == all_virtual_nets_ptr[b]->index && out_degree(a, g) > out_degree(b, g));
			return overlap[a].size() > overlap[b].size();
			});
#else
	/* this gives better result than largest first order */
	custom_smallest_last_vertex_ordering(overlap, order_map); 
#endif

	//for (int i = 0; i < overlap.size(); ++i) {
		//printf("order %d\n", order[i]);
	//}
}

void topological_sort_2(const vector<vector<int>> &g, const vector<int> &order, vector<int> &sorted)
{
	vector<bool> marked(g.size(), false);

	for (const auto &current : order) {

	}

	vector<int> num_incoming_edges(g.size(), 0);
	for (const auto &adj : g) {
		for (const auto &to : adj) {
			++num_incoming_edges[to];
		}
	}
	queue<int> s;
	for (int i = 0; i < num_incoming_edges.size(); ++i) {
		if (num_incoming_edges[i] == 0) {
			s.push(i);
		}
	}
	while (!s.empty()) {
		int item = s.front(); s.pop();
		sorted.push_back(item);
		for (const auto &to : g[item]) {
			if (--num_incoming_edges[to] == 0) {
				s.push(to);
			}
		}
	}
	assert(all_of(begin(num_incoming_edges), end(num_incoming_edges), [] (int val) -> bool { return val == 0; }));
}

void topological_sort(const vector<vector<int>> &g, vector<int> &sorted)
{
	vector<int> num_incoming_edges(g.size(), 0);
	for (const auto &adj : g) {
		for (const auto &to : adj) {
			++num_incoming_edges[to];
		}
	}
	queue<int> s;
	for (int i = 0; i < num_incoming_edges.size(); ++i) {
		if (num_incoming_edges[i] == 0) {
			s.push(i);
		}
	}
	while (!s.empty()) {
		int item = s.front(); s.pop();
		sorted.push_back(item);
		for (const auto &to : g[item]) {
			if (--num_incoming_edges[to] == 0) {
				s.push(to);
			}
		}
	}
	assert(all_of(begin(num_incoming_edges), end(num_incoming_edges), [] (int val) -> bool { return val == 0; }));
}

template<typename LoadFunc>
void schedule_virtual_nets_6(const vector<vector<int>> &g, const vector<new_virtual_net_t *> &all_virtual_nets_ptr, const LoadFunc &load, int num_partitions, vector<vector<int>> &directed_g, vector<int> &num_incoming_edges)
{
	vector<int> sorted_nodes(g.size());
	for (int i = 0; i < g.size(); ++i) {
		sorted_nodes[i] = i;
	}
	std::sort(begin(sorted_nodes), end(sorted_nodes), [&g] (int a, int b) -> bool { return g[a].size() > g[b].size(); });

	struct node_item_t {
		int id;
		float start_time;
		unsigned long num_out_edges;

		//bool operator<(const node_item_t &other) const {
			//return start_time > other.start_time || (start_time == other.start_time && num_out_edges < other.num_out_edges);
		//}
		bool operator<(const node_item_t &other) const {
			return start_time > other.start_time;
		}
	};

	boost::heap::binomial_heap<node_item_t, boost::heap::stable<true>> min_node_heap;
	using node_handle_t = typename boost::heap::binomial_heap<node_item_t, boost::heap::stable<true>>::handle_type;

	vector<node_handle_t> node_handles(g.size());

	vector<float> start_time(g.size());

	for (const auto &n : sorted_nodes) {
		float stime;
		if (all_virtual_nets_ptr[n]->index > 0) {
			stime = std::numeric_limits<float>::max();
		} else {
			assert(all_virtual_nets_ptr[n]->index == 0);
			stime = 0;
		}
		start_time[n] = stime;
		node_handles[n] = min_node_heap.push({ n, stime, g[n].size() });
	}

	struct partition_item_t {
		int id;
		float start_time;

		bool operator<(const partition_item_t &other) const {
			return start_time > other.start_time;
		}
	};
	boost::heap::binomial_heap<partition_item_t, boost::heap::stable<true>> min_partition_heap;

	for (int i = 0; i < num_partitions; ++i) {
		min_partition_heap.push({ i, 0 });
	}

	vector<bool> scheduled(g.size(), false);

	directed_g.resize(g.size());
	for (auto &u : directed_g) {
		u.clear();
	}

	int num_edges = 0;
	for (const auto &u : g) {
		num_edges += u.size();
	}

	int num_true_directed_edges = 0;

	num_incoming_edges.resize(g.size());
	std::fill(begin(num_incoming_edges), end(num_incoming_edges), 0);

	while (!min_node_heap.empty()) {
		node_item_t min_node = min_node_heap.top();
		min_node_heap.pop();

		int u = min_node.id;

		assert(!scheduled[u]);
		scheduled[u] = true;

		start_time[u] = min_node.start_time;

		partition_item_t min_partition = min_partition_heap.top();
		min_partition_heap.pop();

		float uload = load(u);
		assert(uload > 0);
		float new_start = std::max(min_node.start_time, min_partition.start_time) + uload;
		//float new_start = min_node.start_time + uload;

		min_partition_heap.push({ min_partition.id, new_start });

		PRINT("Current %d Parent %d [start time %g load %g] sent to partition %d with start time %g\n", u, all_virtual_nets_ptr[u]->parent, min_node.start_time, uload, min_partition.id, min_partition.start_time);

		PRINT("New start %g = max(item %g, part %g) + load %g\n", new_start, min_node.start_time, min_partition.start_time, uload);

		/* it is possible that we have less amount of threads than the amount of parallelism available */
		//assert(min_node.start_time >= min_partition.start_time);
		
		int parent = all_virtual_nets_ptr[u]->parent;
		if (parent >= 0) {
			assert(scheduled[parent]);
		} else {
			assert(parent == -1);
		}	

		for (const auto &v : g[u]) {
			if (scheduled[v]) {
				assert(start_time[v] < start_time[u]);
				continue;
			}

			assert((*node_handles[v]).id == v);
			assert((*node_handles[v]).start_time == start_time[v]);

			PRINT("\tNeighbor %d Current start %g\n", v, (*node_handles[v]).start_time);

			directed_g[u].push_back(v);
			++num_incoming_edges[v];
			if (all_virtual_nets_ptr[u]->net == all_virtual_nets_ptr[v]->net) {
				assert(all_virtual_nets_ptr[v]->index == all_virtual_nets_ptr[u]->index+1);
				++num_true_directed_edges;
			}

			if ((*node_handles[v]).start_time == std::numeric_limits<float>::max() && all_virtual_nets_ptr[v]->net == all_virtual_nets_ptr[u]->net) {
				assert(all_virtual_nets_ptr[v]->index == all_virtual_nets_ptr[u]->index+1);
				assert(new_start < (*node_handles[v]).start_time);

				PRINT("\t\tUpdating virtual child %d start time from %g to %g\n", v, (*node_handles[v]).start_time, new_start);

				min_node_heap.update(node_handles[v], { v, new_start, g[v].size() });

				start_time[v] = new_start;
			} else if (new_start > (*node_handles[v]).start_time) {
				PRINT("\t\tUpdating neighbor %d start time from %g to %g\n", v, (*node_handles[v]).start_time, new_start);

				min_node_heap.update(node_handles[v], { v, new_start, g[v].size() });

				start_time[v] = new_start;
			}
		}
	}

	assert(std::all_of(begin(scheduled), end(scheduled), [] (const bool &s) -> bool { return s; }));

	int num_directed_edges = 0;
	for (const auto &u : directed_g) {
		num_directed_edges += u.size();
	}

	printf("%d %d %d\n", num_edges, num_true_directed_edges, num_directed_edges);
	assert(((num_edges-num_true_directed_edges) % 2) == 0);
	assert((num_edges-num_true_directed_edges)/2 == num_directed_edges-num_true_directed_edges);
}

template<typename LoadFunc>
void schedule_virtual_nets_5(const vector<vector<int>> &g, const vector<new_virtual_net_t *> &all_virtual_nets_ptr, const LoadFunc &load, int num_partitions, vector<vector<int>> &directed_g, vector<int> &num_incoming_edges)
{
	vector<int> sorted_nodes(g.size());
	for (int i = 0; i < g.size(); ++i) {
		sorted_nodes[i] = i;
	}
	std::sort(begin(sorted_nodes), end(sorted_nodes), [&g] (int a, int b) -> bool { return g[a].size() > g[b].size(); });

	struct node_item_t {
		int id;
		float start_time;
		unsigned long num_out_edges;

		//bool operator<(const node_item_t &other) const {
			//return start_time > other.start_time || (start_time == other.start_time && num_out_edges < other.num_out_edges);
		//}
		bool operator<(const node_item_t &other) const {
			return start_time > other.start_time;
		}
	};

	boost::heap::binomial_heap<node_item_t, boost::heap::stable<true>> min_node_heap;
	using node_handle_t = typename boost::heap::binomial_heap<node_item_t, boost::heap::stable<true>>::handle_type;

	vector<node_handle_t> node_handles(g.size());

	for (const auto &n : sorted_nodes) {
		float start_time;
		if (all_virtual_nets_ptr[n]->index > 0) {
			start_time = std::numeric_limits<float>::max();
		} else {
			assert(all_virtual_nets_ptr[n]->index == 0);
			start_time = 0;
		}
		node_handles[n] = min_node_heap.push({ n, start_time, g[n].size() });
	}

	struct partition_item_t {
		int id;
		float start_time;

		bool operator<(const partition_item_t &other) const {
			return start_time > other.start_time;
		}
	};
	boost::heap::binomial_heap<partition_item_t, boost::heap::stable<true>> min_partition_heap;

	for (int i = 0; i < num_partitions; ++i) {
		min_partition_heap.push({ i, 0 });
	}

	vector<bool> scheduled(g.size(), false);

	directed_g.resize(g.size());
	for (auto &u : directed_g) {
		u.clear();
	}

	int num_edges = 0;
	for (const auto &u : g) {
		num_edges += u.size();
	}

	int num_directed_edges = 0;
	int num_true_directed_edges = 0;

	num_incoming_edges.resize(g.size());
	std::fill(begin(num_incoming_edges), end(num_incoming_edges), 0);

	while (!min_node_heap.empty()) {
		node_item_t min_node = min_node_heap.top();
		min_node_heap.pop();

		int u = min_node.id;

		assert(!scheduled[u]);
		scheduled[u] = true;

		partition_item_t min_partition = min_partition_heap.top();
		min_partition_heap.pop();

		float uload = load(u);
		float new_start = std::max(min_node.start_time, min_partition.start_time) + uload;
		//float new_start = min_node.start_time + uload;

		min_partition_heap.push({ min_partition.id, new_start });

		PRINT("Current %d Parent %d [start time %g load %g] sent to partition %d with start time %g\n", u, all_virtual_nets_ptr[u]->parent, min_node.start_time, uload, min_partition.id, min_partition.start_time);

		PRINT("New start %g = max(item %g, part %g) + load %g\n", new_start, min_node.start_time, min_partition.start_time, uload);

		/* it is possible that we have less amount of threads than the amount of parallelism available */
		//assert(min_node.start_time >= min_partition.start_time);
		
		int parent = all_virtual_nets_ptr[u]->parent;
		if (parent >= 0) {
			assert(scheduled[parent]);
			directed_g[parent].push_back(u);
			++num_incoming_edges[u];
			++num_directed_edges;
			++num_true_directed_edges;
		} else {
			assert(parent == -1);
		}	

		for (const auto &v : g[u]) {
			if (scheduled[v]) {
				directed_g[v].push_back(u);
				++num_incoming_edges[u];
				++num_directed_edges;
				continue;
			}

			assert((*node_handles[v]).id == v);

			PRINT("\tNeighbor %d Current start %g\n", v, (*node_handles[v]).start_time);

			if ((*node_handles[v]).start_time == std::numeric_limits<float>::max() && all_virtual_nets_ptr[v]->net == all_virtual_nets_ptr[u]->net) {
				assert(all_virtual_nets_ptr[v]->index == all_virtual_nets_ptr[u]->index+1);
				assert(new_start < (*node_handles[v]).start_time);

				PRINT("\t\tUpdating virtual child %d start time from %g to %g\n", v, (*node_handles[v]).start_time, new_start);

				min_node_heap.update(node_handles[v], { v, new_start, g[v].size() });
			} else if (new_start > (*node_handles[v]).start_time) {
				PRINT("\t\tUpdating neighbor %d start time from %g to %g\n", v, (*node_handles[v]).start_time, new_start);

				min_node_heap.update(node_handles[v], { v, new_start, g[v].size() });
			}
		}
	}

	assert(std::all_of(begin(scheduled), end(scheduled), [] (const bool &s) -> bool { return s; }));

	assert(((num_edges-num_true_directed_edges) % 2) == 0);
	assert((num_edges-num_true_directed_edges)/2 == num_directed_edges-num_true_directed_edges);
}

#ifdef OLD_CODE
void schedule_virtual_nets_4(vector<vector<int>> &overlap, vector<new_virtual_net_t *> &all_virtual_nets_ptr, vector<vertex_descriptor> &order, int num_threads, const vector<timer::duration> *vnet_route_time, vector<vector<int>> &directed, vector<new_virtual_net_t *> &topo_nets)
{
	boost::vector_property_map<vertices_size_type> color;

	auto coloring_start = timer::now();

	auto order_map = make_iterator_property_map(begin(order), boost::identity_property_map());

	int num_colors;
	if (vnet_route_time) {
		num_colors = custom_vertex_coloring_3(overlap, all_virtual_nets_ptr, [&vnet_route_time] (int v) -> long long { return (*vnet_route_time)[v].count(); }, num_threads, order_map, color);
	} else {
		num_colors = custom_vertex_coloring_3(overlap, all_virtual_nets_ptr, [&all_virtual_nets_ptr] (int v) -> int { return bg::area(all_virtual_nets_ptr[v]->bounding_box); }, num_threads, order_map, color);
	}

	auto coloring_time = timer::now()-coloring_start;

	/* checking validity of coloring */
	int vnet = 0;
	net_t *current_net = all_virtual_nets_ptr[0]->net;

	for (int i = 0; i < all_virtual_nets_ptr.size(); ) {
		set<int> colors;
		int previous_color = -1;
		int previous_index = -1;
		int num_virtual_nets = 0;

		for (; i < all_virtual_nets_ptr.size() && all_virtual_nets_ptr[i]->net == current_net; ++i) {
			const auto &virtual_net = *all_virtual_nets_ptr[i];

			int c = color[virtual_net.global_index];
			colors.insert(c);

			printf("vnet %d index %d v %d color %d\n", vnet, virtual_net.index, virtual_net.global_index, c);

			assert(virtual_net.index > previous_index);
			previous_index = virtual_net.index;

			assert(c > previous_color);
			previous_color = c;

			++num_virtual_nets;
		}

		assert(colors.size() == num_virtual_nets);

		if (i < all_virtual_nets_ptr.size()) {
			assert(all_virtual_nets_ptr[i]->net != current_net);
			assert(all_virtual_nets_ptr[i-1]->net == current_net);
			current_net = all_virtual_nets_ptr[i]->net;
		}

		++vnet;
	}

	vector<vector<new_virtual_net_t *>> phase_nets(num_colors);

	for (const auto &ptr : all_virtual_nets_ptr) {
		phase_nets[get(color, ptr->global_index)].push_back(ptr);
	}

	for (const auto &phase : phase_nets) {
		verify_ind(phase, [] (new_virtual_net_t *net) -> box { return net->bounding_box; });
	}

	/* coloring stats */
	vector<int> num_v(num_colors, 0);
	for (const auto &ptr : all_virtual_nets_ptr) {
		++num_v[get(color, ptr->global_index)];
	}

	using namespace boost::accumulators;
	accumulator_set<int, stats<tag::mean, tag::max, tag::min, tag::variance > > acc;

	for (int i = 0; i < num_v.size(); ++i) {
		printf("color %d %d\n", i, num_v[i]);
		acc(num_v[i]);
	}

	printf("num colors %d min %d max %d mean %g std dev %g\n", num_colors, boost::accumulators::min(acc), boost::accumulators::max(acc), mean(acc), std::sqrt(variance(acc)));
	printf("coloring time %g\n", duration_cast<nanoseconds>(coloring_time).count() / 1e9);

	/* topo sort */
	directed.resize(overlap.size());

	for (int i = 0; i < overlap.size(); ++i) {
		assert(directed[i].empty());
		assert(all_virtual_nets_ptr[i]->global_index == i);

		int p = all_virtual_nets_ptr[i]->parent;
		if (p != -1) {
			directed[p].push_back(i);
			//++all_virtual_nets_ptr[i]->num_incoming_edges;
		}

		for (int j = 0; j < overlap[i].size(); ++j) {
			int to = overlap[i][j];
			assert(color[i] != color[to]);
			if (color[to] > color[i]) {
				directed[i].push_back(to);

				//++all_virtual_nets_ptr[to]->num_incoming_edges;
			}
		}
	}

	vector<int> sorted;
	topological_sort(directed, sorted);

	for (const auto &s : sorted) {
		topo_nets.push_back(all_virtual_nets_ptr[s]);
	}
}

void schedule_virtual_nets_3(vector<vector<int>> &overlap, vector<new_virtual_net_t *> &all_virtual_nets_ptr, vector<vertex_descriptor> &order, int num_threads, const vector<timer::duration> *vnet_route_time, vector<vector<new_virtual_net_t *>> &phase_nets)
{
	boost::vector_property_map<vertices_size_type> color;

	auto coloring_start = timer::now();

	auto order_map = make_iterator_property_map(begin(order), boost::identity_property_map());

	int num_colors;
	if (vnet_route_time) {
		num_colors = custom_vertex_coloring_3(overlap, all_virtual_nets_ptr, [&vnet_route_time] (int v) -> long long { return (*vnet_route_time)[v].count(); }, num_threads, order_map, color);
	} else {
		num_colors = custom_vertex_coloring_3(overlap, all_virtual_nets_ptr, [&all_virtual_nets_ptr] (int v) -> int { return bg::area(all_virtual_nets_ptr[v]->bounding_box); }, num_threads, order_map, color);
	}

	auto coloring_time = timer::now()-coloring_start;

	int vnet = 0;
	net_t *current_net = all_virtual_nets_ptr[0]->net;

	for (int i = 0; i < all_virtual_nets_ptr.size(); ) {
		set<int> colors;
		int previous_color = -1;
		int previous_index = -1;
		int num_virtual_nets = 0;

		for (; i < all_virtual_nets_ptr.size() && all_virtual_nets_ptr[i]->net == current_net; ++i) {
			const auto &virtual_net = *all_virtual_nets_ptr[i];

			int c = color[virtual_net.global_index];
			colors.insert(c);

			printf("vnet %d index %d v %d color %d\n", vnet, virtual_net.index, virtual_net.global_index, c);

			assert(virtual_net.index > previous_index);
			previous_index = virtual_net.index;

			assert(c > previous_color);
			previous_color = c;

			++num_virtual_nets;
		}

		assert(colors.size() == num_virtual_nets);

		if (i < all_virtual_nets_ptr.size()) {
			assert(all_virtual_nets_ptr[i]->net != current_net);
			assert(all_virtual_nets_ptr[i-1]->net == current_net);
			current_net = all_virtual_nets_ptr[i]->net;
		}

		++vnet;
	}

	vector<int> num_v(num_colors, 0);
	for (const auto &ptr : all_virtual_nets_ptr) {
		++num_v[get(color, ptr->global_index)];
	}

	using namespace boost::accumulators;
	accumulator_set<int, stats<tag::mean, tag::max, tag::min, tag::variance > > acc;

	for (int i = 0; i < num_v.size(); ++i) {
		printf("color %d %d\n", i, num_v[i]);
		acc(num_v[i]);
	}

	printf("num colors %d min %d max %d mean %g std dev %g\n", num_colors, boost::accumulators::min(acc), boost::accumulators::max(acc), mean(acc), std::sqrt(variance(acc)));
	printf("coloring time %g\n", duration_cast<nanoseconds>(coloring_time).count() / 1e9);
	
	phase_nets.resize(num_colors);

	for (const auto &ptr : all_virtual_nets_ptr) {
		phase_nets[get(color, ptr->global_index)].push_back(ptr);
	}

	for (const auto &phase : phase_nets) {
		verify_ind(phase, [] (new_virtual_net_t *net) -> box { return net->bounding_box; });
	}
}

void schedule_virtual_nets_2(vector<vector<new_virtual_net_t>> &all_virtual_nets, vector<vector<new_virtual_net_t *>> &phase_nets)
{
	OverlapGraph g;
	vector<new_virtual_net_t *> all_virtual_nets_ptr;

	auto build_start = timer::now();

	build_overlap_graph(all_virtual_nets, g, all_virtual_nets_ptr);

	auto build_time = timer::now()-build_start;

	printf("overlap graph build time %g\n", duration_cast<nanoseconds>(build_time).count() / 1e9);
	printf("overlap num vertices %d num edges %d\n", num_vertices(g), num_edges(g));

	//auto ep = edges(g);
	//for (auto ei = ep.first; ei != ep.second; ++ei) {
		//printf("edge %d -> %d\n", source(*ei, g), target(*ei, g));
	//}

	vector<vertex_descriptor> order(num_vertices(g));
	auto order_map = make_iterator_property_map(begin(order), get(boost::vertex_index, g));

#ifdef LARGEST_FIRST
	for (int i = 0; i < num_vertices(g); ++i) {
		order[i] = i;
	}

	std::sort(begin(order), end(order), [&] (int a, int b) -> bool {
			//assert(a == all_virtual_nets_ptr[a]->v); 
			//assert(b == all_virtual_nets_ptr[b]->v); 

			/* BUG: the relation is not transitive */
			//if (all_virtual_nets_ptr[a]->net == all_virtual_nets_ptr[b]->net) {
			//return all_virtual_nets_ptr[a]->index < all_virtual_nets_ptr[b]->index;
			//} else {
			//return out_degree(a, g) > out_degree(b, g);
			//}
			//return all_virtual_nets_ptr[a]->index < all_virtual_nets_ptr[b]->index || (all_virtual_nets_ptr[a]->index == all_virtual_nets_ptr[b]->index && out_degree(a, g) > out_degree(b, g));
			return out_degree(a, g) > out_degree(b, g);
			});
#else
	/* this gives better result than largest first order */
	smallest_last_vertex_ordering(g, order_map); 
#endif

	for (int i = 0; i < num_vertices(g); ++i) {
		printf("order %d\n", order[i]);
	}

	boost::vector_property_map<vertices_size_type> color;

	auto coloring_start = timer::now();

	int num_colors = custom_vertex_coloring(g, order_map, color, all_virtual_nets_ptr);

	auto coloring_time = timer::now()-coloring_start;

	int vnet = 0;
	for (const auto &virtual_nets : all_virtual_nets) {
		set<int> colors;
		int previous_color = -1;
		int previous_index = -1;
		for (const auto &virtual_net : virtual_nets) {
			int c = color[virtual_net.global_index];
			colors.insert(c);

			//printf("vnet %d index %d v %d color %d\n", vnet, virtual_net.index, virtual_net.global_index, c);

			assert(virtual_net.index > previous_index);
			previous_index = virtual_net.index;

			assert(c > previous_color);
			previous_color = c;
		}
		assert(colors.size() == virtual_nets.size());
		++vnet;
	}

	vector<int> num_v(num_colors, 0);
	for (const auto &ptr : all_virtual_nets_ptr) {
		++num_v[get(color, ptr->global_index)];
	}

	using namespace boost::accumulators;
	accumulator_set<int, stats<tag::mean, tag::max, tag::min, tag::variance > > acc;

	for (int i = 0; i < num_v.size(); ++i) {
		printf("color %d %d\n", i, num_v[i]);
		acc(num_v[i]);
	}

	printf("num colors %d min %d max %d mean %g std dev %g\n", num_colors, boost::accumulators::min(acc), boost::accumulators::max(acc), mean(acc), std::sqrt(variance(acc)));
	printf("coloring time %g\n", duration_cast<nanoseconds>(coloring_time).count() / 1e9);
	
	phase_nets.resize(num_colors);
	for (const auto &ptr : all_virtual_nets_ptr) {
		phase_nets[get(color, ptr->global_index)].push_back(ptr);
	}

	for (const auto &phase : phase_nets) {
		verify_ind(phase, [] (new_virtual_net_t *net) -> box { return net->bounding_box; });
	}
}

void schedule_virtual_nets(vector<vector<new_virtual_net_t>> &all_virtual_nets, vector<vector<new_virtual_net_t *>> &phase_nets)
{
	OverlapGraph g;
	vector<new_virtual_net_t *> all_virtual_nets_ptr;
	build_overlap_graph(all_virtual_nets, g, all_virtual_nets_ptr);

	auto ep = edges(g);
	for (auto ei = ep.first; ei != ep.second; ++ei) {
		printf("edge %d -> %d\n", source(*ei, g), target(*ei, g));
	}

	vector<vertex_descriptor> order(num_vertices(g));

	for (int i = 0; i < num_vertices(g); ++i) {
		order[i] = i;
	}

	std::sort(begin(order), end(order), [&] (int a, int b) -> bool {
			//assert(a == all_virtual_nets_ptr[a]->v); 
			//assert(b == all_virtual_nets_ptr[b]->v); 

			/* BUG: the relation is not transitive */
			//if (all_virtual_nets_ptr[a]->net == all_virtual_nets_ptr[b]->net) {
			//return all_virtual_nets_ptr[a]->index < all_virtual_nets_ptr[b]->index;
			//} else {
			//return out_degree(a, g) > out_degree(b, g);
			//}
			return all_virtual_nets_ptr[a]->index < all_virtual_nets_ptr[b]->index || (all_virtual_nets_ptr[a]->index == all_virtual_nets_ptr[b]->index && out_degree(a, g) > out_degree(b, g));
			});

	auto order_map = make_iterator_property_map(begin(order), get(boost::vertex_index, g));

	boost::vector_property_map<vertices_size_type> color;

	int n = boost::sequential_vertex_coloring(g, order_map, color);
	printf("num colors %d\n", n);
	
	int vnet = 0;
	for (const auto &virtual_nets : all_virtual_nets) {
		set<int> colors;
		int previous_color = -1;
		int previous_index = -1;
		for (const auto &virtual_net : virtual_nets) {
			int c = color[virtual_net.global_index];
			colors.insert(c);

			printf("vnet %d index %d v %d color %d\n", vnet, virtual_net.index, virtual_net.global_index, c);

			assert(virtual_net.index > previous_index);
			previous_index = virtual_net.index;

			assert(c > previous_color);
			previous_color = c;
		}
		assert(colors.size() == virtual_nets.size());
		++vnet;
	}
}
#endif

void test_rtree(vector<net_t> &nets, vector<vector<net_t *>> &phase_nets)
{
	vector<net_t *> sorted_nets;
	for (auto &net : nets) {
		sorted_nets.push_back(&net);
	}
	std::sort(begin(sorted_nets), end(sorted_nets), [] (const net_t *a, const net_t *b) -> bool {
			return a->sinks.size() > b->sinks.size();
			});

	auto rtree_start = timer::now();
	//bgi::rtree<rtree_value, bgi::rstar<16>> tree(nets | boost::adaptors::transformed([] (const net_t &net) -> pair<box, int> {
				//return make_pair(box(point(net.bounding_box.xmin, net.bounding_box.ymin), point(net.bounding_box.xmax, net.bounding_box.ymax)), net.local_id);
				//}
				//));

	net_to_rtree_item to_rtree_item;
	bgi::rtree<rtree_value, bgi::rstar<16>, bgi::indexable<rtree_value>, rtree_value_equal> tree(sorted_nets | boost::adaptors::transformed(to_rtree_item));

	auto rtree_time = timer::now()-rtree_start;

	printf("rtree build time %g\n", duration_cast<nanoseconds>(rtree_time).count() / 1e9);

	auto to_box = [] (const net_t &net) -> box { return box(point(net.bounding_box.xmin, net.bounding_box.ymin), point(net.bounding_box.xmax, net.bounding_box.ymax)); };

	//auto lol = tree.qbegin(bgi::satisfies([] (const rtree_value &val) -> bool { return true; }));
	vector<bool> net_scheduled(sorted_nets.size(), false);
	int inet = 0;
	int num_rounds = 0;
	int max_con = 0;
	auto sched_start = timer::now();
	while (!tree.empty()) {
		net_t *net = sorted_nets[inet];

		if (!net_scheduled[net->local_id]) {
			vector<rtree_value> disjoint_nets;
			vector<rtree_value> dispatched;

			dispatched.push_back(make_pair(to_box(*net), net));

			if (tree.query(bgi::disjoint(to_box(*net)), std::back_inserter(disjoint_nets)) > 0) {
				/* we could do a brute force search here instead of using MISR if sizeof disjoin_nets is small enough to consider all
				 * possible combinations */
				/* run graph coloring on disjoint_nets */
				/* for each item in disjoint_nets, find out their disjoint nets and do a set union or set covering */
				std::sort(begin(disjoint_nets), end(disjoint_nets), [] (const rtree_value &a, const rtree_value &b) -> bool {
						return bg::area(a.first) > bg::area(b.first);
						});

				for (int i = 0; i < disjoint_nets.size(); ++i) {
					const auto &cur_dis = disjoint_nets[i];
					bool intersect = std::any_of(begin(dispatched), end(dispatched), [&cur_dis] (const rtree_value &n) -> bool {
							return bg::intersects(cur_dis.first, n.first);
							});
					if (!intersect) {
						dispatched.push_back(cur_dis);
					}
				}


				//int max_area = -1;
				//int max = -1;
				//for (int i = 0; i < disjoint_nets.size(); ++i) {
					//int area = bg::area(to_box(disjoint_nets[i]));
					//if (area > max_area) {
						//max_area = area;
						//max = i;
					//}
				//}

				//assert(max >= 0);
				//tree.remove(make_pair(to_box(disjoint_nets[max]), disjoint_nets[max].local_id));
			}

			assert(tree.remove(begin(dispatched), end(dispatched)) == dispatched.size());

			//printf("Num dispatched: %d\n", dispatched.size());
			verify_ind(dispatched, [] (const rtree_value &a) -> box { return a.first; });

			for (const auto &d : dispatched) {
				assert(!net_scheduled[d.second->local_id]);
				net_scheduled[d.second->local_id] = true;
			}

			phase_nets.emplace_back();
			for (const auto &d : dispatched) {
				phase_nets.back().push_back(d.second);
			}

			++num_rounds;
			max_con = std::max(max_con, (int)dispatched.size());
		}

		++inet;
	}
	auto sched_time = timer::now() - sched_start;

	printf("sched_time: %g max_con: %d ave_con: %g\n", duration_cast<nanoseconds>(sched_time).count()/1e9, max_con, sorted_nets.size()*1.0/num_rounds);

	assert(std::all_of(begin(net_scheduled), end(net_scheduled), [] (bool a) -> bool { return a; }));
}

void test_misr(const vector<net_t> &nets)
{
	set<int> all_nets;

	auto to_box = [] (const net_t &net) -> box { return box(point(net.bounding_box.xmin, net.bounding_box.ymin), point(net.bounding_box.xmax, net.bounding_box.ymax)); };

	int subiter = 0;
	int max_con = 0;

	vector<net_t> nets_copy = nets;

	while (!nets_copy.empty()) {
		vector<net_t *> chosen;

		max_independent_rectangles(nets_copy, to_box, [] (const net_t &net) -> pair<int, int> { return make_pair(net.bounding_box.ymin, net.bounding_box.ymax); }, chosen);

		verify_ind(chosen, [&to_box] (const net_t *net) -> box { return to_box(*net); });

		extern int nx, ny;

		FILE *file;
		if (subiter < 10) {
			char buffer[256];
			sprintf(buffer, "%d_boxes.txt", subiter);
			file = fopen(buffer, "w");
			fprintf(file, "0 0 %d %d 0\n", nx+2, ny+2);
		}

		int total_area = 0;
		for (const auto &c : chosen) {
			assert(all_nets.find(c->vpr_id) == all_nets.end());
			all_nets.insert(c->vpr_id);

			total_area += bg::area(to_box(*c));

			if (subiter < 10) {
				fprintf(file, "%d %d %d %d 0\n", c->bounding_box.xmin, c->bounding_box.ymin, c->bounding_box.xmax-c->bounding_box.xmin, c->bounding_box.ymax-c->bounding_box.ymin);
			}
		}

		vector<net_t *> new_chosen;
		for (auto &net : nets_copy) {
			const auto &box = to_box(net);
			bool dis = std::all_of(begin(chosen), end(chosen), [&to_box, &box] (const net_t *c) -> bool {
						return bg::disjoint(to_box(*c), box);
					}) && std::all_of(begin(new_chosen), end(new_chosen), [&to_box, &box] (const net_t *c) -> bool {
						return bg::disjoint(to_box(*c), box);
					});
			if (dis) {
				if (subiter < 10) {
					fprintf(file, "%d %d %d %d 1\n", net.bounding_box.xmin, net.bounding_box.ymin, net.bounding_box.xmax-net.bounding_box.xmin, net.bounding_box.ymax-net.bounding_box.ymin);
				}
				new_chosen.push_back(&net);
			}
		}

		if (subiter < 10) {
			fclose(file);
		}

		++subiter;

		printf("chosen size: %d %d %g\n", chosen.size(), total_area, 100.0*total_area/((nx+2)*(ny+2)));
		max_con = std::max((int)chosen.size(), max_con);

		nets_copy.erase(std::remove_if(begin(nets_copy), end(nets_copy), [&chosen] (const net_t &other) -> bool {
					return std::any_of(begin(chosen), end(chosen), [&other] (const net_t *c) -> bool {
							return c->vpr_id == other.vpr_id;
							});
					}), end(nets_copy));
	}

	printf("max_con: %d ave_con: %g\n", max_con, nets.size()*1.0/subiter);
}

typedef struct extra_route_state_t {
	float upstream_R;
	float delay;
} extra_route_state_t;

typedef struct dijkstra_stats_t {
	int num_heap_pops;
	int num_heap_pushes;
	int num_neighbor_visits;
} dijkstra_stats_t;

bool outside_bb(const rr_node_property_t &node, const vector<box> &boxes)
{
	bool outside = true;

	for (int i = 0; i < boxes.size() && outside; ++i) {
		if (node.type == CHANX || node.type == CHANY) {
			int x, y;
			if (node.inc_direction) {
				x = node.xlow;
				y = node.ylow;
			} else {
				x = node.xhigh;
				y = node.yhigh;
			}

			if (x >= bg::get<bg::min_corner, 0>(boxes[i])
					&& x <= bg::get<bg::max_corner, 0>(boxes[i])
					&& y >= bg::get<bg::min_corner, 1>(boxes[i])
					&& y <= bg::get<bg::max_corner, 1>(boxes[i])) {
				outside = false;
			}
		} else {
			if (node.xhigh >= bg::get<bg::min_corner, 0>(boxes[i])
					&& node.xlow <= bg::get<bg::max_corner, 0>(boxes[i])
					&& node.yhigh >= bg::get<bg::min_corner, 1>(boxes[i])
					&& node.ylow <= bg::get<bg::max_corner, 1>(boxes[i])) {
				outside = false;
			}
		}
	}

	return outside;
}

void rand_delay()
{
	timespec delay;
	delay.tv_sec = 0;
	int r = rand() % 100;
	delay.tv_nsec = r * 10;
	assert(!nanosleep(&delay, nullptr));
}

class DeltaSteppingRouter {
	private:
		RRGraph &_g;

		float _astar_fac;

		float *_known_distance;
		float *_distance;
		RREdge *_prev_edge;
		extra_route_state_t *_state;

		congestion_t *_congestion;
		const float &_pres_fac;

		vector<int> _modified_nodes;
		vector<bool> _modified_node_added;

		const sink_t *_current_sink;
		RRNode _existing_opin;

		route_tree_t *_current_rt;

		dijkstra_stats_t _stats;

		static tbb::atomic<int> _num_instances;
		int _instance;

	private:
		void popped_node(int v)
		{
			char buffer[256];
			sprintf_rr_node(v, buffer);
			zlog_level(delta_log, ROUTER_V3, "%s\n", buffer);
			++_stats.num_heap_pops;

			rand_delay();
		}

		void relax_node(int v, const RREdge &e)
		{
			//RouteTreeNode rt_node = route_tree_get_rt_node(*_current_rt, v);

			//if (rt_node != RouteTree::null_vertex()) {
				//const auto &rt_node_p = get_vertex_props(_current_rt->graph, rt_node);
				//if (!valid(_prev_edge[v])) {
					//assert(!valid(rt_node_p.rt_edge_to_parent));
				//} else {
					//int rr = get_vertex_props(_current_rt->graph, get_source(_current_rt->graph, rt_node_p.rt_edge_to_parent)).rr_node;
					//assert(rr == get_source(_g, _prev_edge[v]));
				//}
			//}

			const auto &v_p = get_vertex_props(_g, v);

			if (valid(e)) {
				assert(v == get_target(_g, e));

				int u = get_source(_g, e);

				float u_delay;
				float u_upstream_R;
				if (_state[u].upstream_R != std::numeric_limits<float>::max()) {
					assert(_state[u].delay != std::numeric_limits<float>::max());
					u_upstream_R = _state[u].upstream_R;
					u_delay = _state[u].delay;
				} else {
					auto rt_node = route_tree_get_rt_node(*_current_rt, u);
					assert(rt_node != RouteTree::null_vertex());
					const auto &u_rt_node_p = get_vertex_props(_current_rt->graph, rt_node);
					u_upstream_R = u_rt_node_p.upstream_R;
					u_delay = u_rt_node_p.delay;
				}

				const auto &e_p = get_edge_props(_g, e);

				extern struct s_switch_inf *switch_inf;
				const struct s_switch_inf *sw = &switch_inf[e_p.switch_index];

				_state[v].upstream_R = sw->R + v_p.R;
				if (!sw->buffered)  {
					_state[v].upstream_R += u_upstream_R;
				} 

				float delay;
				if (sw->buffered) {
					delay = sw->Tdel + v_p.C * (sw->R + 0.5f * v_p.R);
				} else {
					delay = sw->Tdel + v_p.C * (u_upstream_R + sw->R + 0.5f * v_p.R);
				}
				_state[v].delay = u_delay + delay;
			} else {
				_state[v].upstream_R = v_p.R;
				_state[v].delay = v_p.C * 0.5f * v_p.R;
			}

			if (!_modified_node_added[v]) {
				_modified_nodes.push_back(v);
				_modified_node_added[v] = true;
			}
		}

		bool expand_node(int node)
		{
			const auto &prop = get_vertex_props(_g, node);

			char buffer[256];
			sprintf_rr_node(node, buffer);

			zlog_level(delta_log, ROUTER_V3, "\tChecking whether to expand %s ", buffer);

			//const auto &sink_prop = get_vertex_props(_g, _current_sink->rr_node);

			// this has to come before the bounding box checking 
			if (prop.type == IPIN) { 
				assert(num_out_edges(_g, node) == 1);
				int target = get_target(_g, *begin(get_out_edges(_g, node)));
				assert(get_vertex_props(_g, target).type == SINK);

				/* this condition is sufficient because the bounding box would have covered
				 * the SINK */
				if (target != _current_sink->rr_node) {
					zlog_level(delta_log, ROUTER_V3, "not target IPIN\n");
					return false;
				}
			}

			if (_existing_opin != RRGraph::null_vertex() && prop.type == OPIN && node != _existing_opin) {
				zlog_level(delta_log, ROUTER_V3, "not expanding other OPIN\n");
				return false;
			}

			//if (prop.xlow >= _current_sink->current_bounding_box.xmin
					//&& prop.xhigh <= _current_sink->current_bounding_box.xmax
					//&& prop.ylow >= _current_sink->current_bounding_box.ymin
					//&& prop.yhigh <= _current_sink->current_bounding_box.ymax) {

			//if (prop.xlow < _current_sink->current_bounding_box.xmin
					//|| prop.xhigh > _current_sink->current_bounding_box.xmax
					//|| prop.ylow < _current_sink->current_bounding_box.ymin
					//|| prop.yhigh > _current_sink->current_bounding_box.ymax) {
				//zlog_level(delta_log, ROUTER_V3, "outside of bounding box\n");
				//return false;
			//}
			
			//if (prop.xhigh < _current_sink->current_bounding_box.xmin
					//|| prop.xlow > _current_sink->current_bounding_box.xmax
					//|| prop.yhigh < _current_sink->current_bounding_box.ymin
					//|| prop.ylow > _current_sink->current_bounding_box.ymax) {
				//zlog_level(delta_log, ROUTER_V3, "outside of bounding box\n");
				//return false;
			//}


			if (outside_bb(prop, _current_sink->vnet->bounding_boxes)) {
				zlog_level(delta_log, ROUTER_V3, "outside of bounding box\n");
				return false;
			}

			zlog_level(delta_log, ROUTER_V3, " OK\n");

			++_stats.num_neighbor_visits;

			return true;
		}

		void push_node(int node)
		{
			++_stats.num_heap_pushes;
		}

		pair<float, float> get_edge_weight(const RREdge &e)
		{
			int u = get_source(_g, e);
			int v = get_target(_g, e);

			const auto &v_p = get_vertex_props(_g, v);
			const auto &e_p = get_edge_props(_g, e);

			extern struct s_switch_inf *switch_inf;
			const struct s_switch_inf *sw = &switch_inf[e_p.switch_index];

			assert(_state[u].upstream_R != std::numeric_limits<float>::max());
			float delay;
			if (sw->buffered) {
				delay = sw->Tdel + v_p.C * (sw->R + 0.5f * v_p.R);
			} else {
				delay = sw->Tdel + v_p.C * (_state[u].upstream_R + sw->R + 0.5f * v_p.R);
			}
			extern t_rr_indexed_data *rr_indexed_data;
			float congestion_cost = rr_indexed_data[v_p.cost_index].base_cost * _congestion[v].acc_cost * _congestion[v].pres_cost;
			float known_cost = _current_sink->criticality_fac * delay + (1 - _current_sink->criticality_fac) * congestion_cost;

			float upstream_R = sw->R + v_p.R;
			if (!sw->buffered) {
				upstream_R += _state[u].upstream_R;
			}
			float expected_cost = get_timing_driven_expected_cost(v_p, get_vertex_props(_g, _current_sink->rr_node), _current_sink->criticality_fac, upstream_R);

			zlog_level(delta_log, ROUTER_V3, "\t%d -> %d delay %g congestion %g crit_fac %g expected %g expected_hex %X known %g predicted %g\n", 
					u, v, delay, congestion_cost, _current_sink->criticality_fac, expected_cost, *(unsigned int *)&expected_cost, known_cost, known_cost + _astar_fac * expected_cost);
			//zlog_level(delta_log, ROUTER_V3, "\t[u: upstream %g] [edge: d %g R %g] [v: R %g C %g]\n",
					//_state[u].upstream_R, sw->Tdel, sw->R, v_p.R, v_p.C);

			return make_pair(known_cost, known_cost + _astar_fac * expected_cost);
		}

		RREdge get_previous_edge(int rr_node_id, const route_tree_t &rt)
		{
			RREdge previous_edge;

			if (!valid(_prev_edge[rr_node_id])) {
				previous_edge = RRGraph::null_edge();
			} else {
				RouteTreeNode rt_node = route_tree_get_rt_node(rt, rr_node_id);

				char buffer[256];
				sprintf_rr_node(rr_node_id, buffer);

				if (rt_node != RouteTree::null_vertex() && valid(get_vertex_props(rt.graph, rt_node).rt_edge_to_parent)) {

					/* already reach an existing route tree node but the current state suggests there's a possible path of 
					 * lower cost */
					/*if (rt_node->rr_edge_to_parent != -1) {*/
					/*char parent[256];*/
					/*sprintf_rr_node(get_source(g, rt_node->rr_edge_to_parent), parent);*/
					/*zlog_error(delta_log, "Error: Existing route tree node %s has non-null rr_edge_to_parent that connects to %s\n", buffer, parent);*/
					/*assert(false);*/
					/*}*/

					char s_state[256];
					char s_rt[256];
					sprintf_rr_node(get_source(_g, _prev_edge[rr_node_id]), s_state);
					sprintf_rr_node(get_vertex_props(rt.graph, get_source(rt.graph, get_vertex_props(rt.graph, rt_node).rt_edge_to_parent)).rr_node, s_rt);
					zlog_warn(delta_log, "Warning: Existing route tree node %s does not have a matching route state. (state.prev_edge: %s rt_node.rr_edge_to_parent: %s) because we have found a shorter path to that node\n", buffer, s_state, s_rt);

					previous_edge = RRGraph::null_edge();
				} else {
					previous_edge = _prev_edge[rr_node_id];
				}
			} 

			return previous_edge;
		}

		void backtrack(int sink_node, route_tree_t &rt, vector<RRNode> &added_rr_nodes)
		{
			char buffer[256];
			sprintf_rr_node(sink_node, buffer);
			zlog_level(delta_log, ROUTER_V3, "Adding %s to route tree\n", buffer);

			RouteTreeNode rt_node = route_tree_add_rr_node(rt, sink_node, _g);
			assert(rt_node != RouteTree::null_vertex());
			assert(_state[sink_node].upstream_R != std::numeric_limits<float>::max() && _state[sink_node].delay != std::numeric_limits<float>::max());
			route_tree_set_node_properties(rt, rt_node, false, _state[sink_node].upstream_R, _state[sink_node].delay);

			const auto &sink_node_p = get_vertex_props(_g, sink_node);
			//assert(sink_node_p.xlow >= _current_sink->current_bounding_box.xmin
					//&& sink_node_p.xhigh <= _current_sink->current_bounding_box.xmax
					//&& sink_node_p.ylow >= _current_sink->current_bounding_box.ymin
					//&& sink_node_p.yhigh <= _current_sink->current_bounding_box.ymax);
			//assert(sink_node_p.xhigh >= _current_sink->current_bounding_box.xmin
					//&& sink_node_p.xlow <= _current_sink->current_bounding_box.xmax
					//&& sink_node_p.yhigh >= _current_sink->current_bounding_box.ymin
					//&& sink_node_p.ylow <= _current_sink->current_bounding_box.ymax);
			assert(!outside_bb(sink_node_p, _current_sink->vnet->bounding_boxes));

			update_one_cost_internal(sink_node, _g, _congestion, 1, _pres_fac);
			added_rr_nodes.push_back(sink_node);

			zlog_level(delta_log, ROUTER_V3, "\n");

			RREdge edge = get_previous_edge(sink_node, rt);

			RRNode child_rr_node = sink_node;

			while (valid(edge)) {
				RRNode parent_rr_node = get_source(_g, edge);

				const auto &parent_rr_node_p = get_vertex_props(_g, parent_rr_node);

				sprintf_rr_node(parent_rr_node, buffer);
				zlog_level(delta_log, ROUTER_V3, "Adding %s to route tree\n", buffer);

				RouteTreeNode parent_rt_node = route_tree_add_rr_node(rt, parent_rr_node, _g);
				if (parent_rt_node != RouteTree::null_vertex()) {
					assert(_state[parent_rr_node].upstream_R != std::numeric_limits<float>::max() && _state[parent_rr_node].delay != std::numeric_limits<float>::max());
					route_tree_set_node_properties(rt, parent_rt_node, parent_rr_node_p.type != IPIN && parent_rr_node_p.type != SINK, _state[parent_rr_node].upstream_R, _state[parent_rr_node].delay);
				} 

				assert(child_rr_node == get_target(_g, edge));

				char buffer2[256];
				sprintf_rr_node(child_rr_node, buffer2);

				zlog_level(delta_log, ROUTER_V3, "Adding edge %s -> %s\n", buffer, buffer2);

				const RouteTreeEdge &rt_edge = route_tree_add_edge_between_rr_node(rt, parent_rr_node, child_rr_node); 
				auto &rt_edge_props = get_edge_props(rt.graph, rt_edge);
				rt_edge_props.rr_edge = edge;

				if (parent_rr_node_p.type == OPIN && _existing_opin == RRGraph::null_vertex()) {
					_existing_opin = parent_rr_node;
				}

				if ((parent_rt_node != RouteTree::null_vertex()) || (get_vertex_props(_g, parent_rr_node).type == SOURCE)) {
					const auto &rr_node_p = get_vertex_props(_g, parent_rr_node);
					//assert(rr_node_p.xlow >= _current_sink->current_bounding_box.xmin
							//&& rr_node_p.xhigh <= _current_sink->current_bounding_box.xmax
							//&& rr_node_p.ylow >= _current_sink->current_bounding_box.ymin
							//&& rr_node_p.yhigh <= _current_sink->current_bounding_box.ymax);
					//assert(rr_node_p.xhigh >= _current_sink->current_bounding_box.xmin
							//&& rr_node_p.xlow <= _current_sink->current_bounding_box.xmax
							//&& rr_node_p.yhigh >= _current_sink->current_bounding_box.ymin
							//&& rr_node_p.ylow <= _current_sink->current_bounding_box.ymax);
					assert(!outside_bb(rr_node_p, _current_sink->vnet->bounding_boxes));

					update_one_cost_internal(parent_rr_node, _g, _congestion, 1, _pres_fac);
					added_rr_nodes.push_back(parent_rr_node);
				}

				zlog_level(delta_log, ROUTER_V3, "\n");

				child_rr_node = parent_rr_node;
				edge = get_previous_edge(parent_rr_node, rt);
			}
		}

		template<typename Graph, typename Edge, typename EdgeWeightFunc, typename Callbacks>
		friend void delta_stepping(const Graph &g, const vector<existing_source_t<Edge>> &sources, int sink, float delta, float *known_distance, float *distance, Edge *prev_edge, const EdgeWeightFunc &edge_weight, Callbacks &callbacks);

		template<typename Graph, typename Edge, typename EdgeWeightFunc, typename ExpandCheckFunc, typename Callbacks>
		friend void dijkstra(const Graph &g, const vector<existing_source_t<Edge>> &sources, int sink, float *known_distance, float *distance, Edge *prev_edge, const ExpandCheckFunc &expand_node, const EdgeWeightFunc &edge_weight, Callbacks &callbacks);

		template<typename Edge, typename Callbacks>
		friend void relax(Buckets &buckets, float delta, vector<bool> &in_bucket, const vector<bool> &vertex_deleted,
					float *known_distance, float *distance, Edge *predecessor,
					int v, float new_known_distance, float new_distance, const Edge &edge,
					Callbacks &callbacks);

	public:
		DeltaSteppingRouter(RRGraph &g, congestion_t *congestion, const float &pres_fac)
			: _g(g), _congestion(congestion), _pres_fac(pres_fac), _modified_node_added(num_vertices(g), false)
		{
			_known_distance = new float[num_vertices(g)];
			_distance = new float[num_vertices(g)];
			_prev_edge = new RREdge[num_vertices(g)];
			_state = new extra_route_state_t[num_vertices(g)];

			for (int i = 0; i < num_vertices(g); ++i) {
				_known_distance[i] = std::numeric_limits<float>::max();
				_distance[i] = std::numeric_limits<float>::max();
				_prev_edge[i] = RRGraph::null_edge();
				_state[i].delay = std::numeric_limits<float>::max();
				_state[i].upstream_R = std::numeric_limits<float>::max();
			}
			
			reset_stats();

			_instance = _num_instances++;
		}

		int get_instance() const
		{
			return _instance;
		}

		int get_num_instances() const
		{
			return _num_instances;
		}

		void reset_stats()
		{
			_stats.num_heap_pops = 0;
			_stats.num_heap_pushes = 0;
			_stats.num_neighbor_visits = 0;
		}

		const dijkstra_stats_t &get_stats() const
		{
			return _stats;
		}

		void route(const source_t *source, const vector<const sink_t *> &sinks, float delta, float astar_fac, route_tree_t &rt, vector<RRNode> &added_rr_nodes, t_net_timing &net_timing)
		{
			char buffer[256];

			_astar_fac = astar_fac;

			/* TODO: we should not reset the exisiting OPIN here */
			_existing_opin = RRGraph::null_vertex();

			if (source) {
				const auto &source_p = get_vertex_props(_g, source->rr_node);
				RouteTreeNode root_rt_node = route_tree_add_rr_node(rt, source->rr_node, _g);
				route_tree_set_node_properties(rt, root_rt_node, true, source_p.R, 0.5f * source_p.R * source_p.C);
				route_tree_add_root(rt, source->rr_node);
			} else {
				assert(rt.root_rt_nodes.size() == 1);
			}

			vector<const sink_t *> sorted_sinks = sinks;

			std::sort(begin(sorted_sinks), end(sorted_sinks), [] (const sink_t *a, const sink_t *b) -> bool {
					return a->criticality_fac > b->criticality_fac;
					});

			_current_rt = &rt;

			int isink = 0;
			for (const auto &sink : sorted_sinks) {
				_current_sink = sink;

				sprintf_rr_node(sink->rr_node, buffer);
				//zlog_level(delta_log, ROUTER_V3, "Current sink: %s BB %d->%d, %d->%d\n", buffer, sink->current_bounding_box.xmin, sink->current_bounding_box.xmax, sink->current_bounding_box.ymin, sink->current_bounding_box.ymax);
				zlog_level(delta_log, ROUTER_V3, "Current sink: %s BB ", buffer);
				for (int i = 0; i < sink->vnet->bounding_boxes.size(); ++i) {
				   	zlog_level(delta_log, ROUTER_V3, "%d->%d,%d->%d ", bg::get<bg::min_corner, 0>(sink->vnet->bounding_boxes[i]), bg::get<bg::max_corner, 0>(sink->vnet->bounding_boxes[i]), bg::get<bg::min_corner, 1>(sink->vnet->bounding_boxes[i]), bg::get<bg::max_corner, 1>(sink->vnet->bounding_boxes[i]));
				}
				zlog_level(delta_log, ROUTER_V3, "\n");

				vector<existing_source_t<RREdge>> sources;

				for (const auto &rt_node : route_tree_get_nodes(rt)) {
					const auto &rt_node_p = get_vertex_props(rt.graph, rt_node);
					RRNode node = rt_node_p.rr_node;
					const auto &node_p = get_vertex_props(_g, node);

					sprintf_rr_node(node, buffer);

					if (rt_node_p.reexpand) {
						//if (node_p.xhigh < sink->current_bounding_box.xmin
								//|| node_p.xlow > sink->current_bounding_box.xmax
								//|| node_p.yhigh < sink->current_bounding_box.ymin
								//|| node_p.ylow > sink->current_bounding_box.ymax) {
						if (outside_bb(node_p, _current_sink->vnet->bounding_boxes)) {
							zlog_level(delta_log, ROUTER_V3, "Existing %s out of bounding box\n", buffer);
						} else {
							float kd = sink->criticality_fac * rt_node_p.delay; 
							float d = kd + _astar_fac * get_timing_driven_expected_cost(node_p, get_vertex_props(_g, sink->rr_node), sink->criticality_fac, rt_node_p.upstream_R); 

							RREdge prev = RRGraph::null_edge();
							if (valid(rt_node_p.rt_edge_to_parent)) {
								prev = get_edge_props(rt.graph, rt_node_p.rt_edge_to_parent).rr_edge;
							}

							zlog_level(delta_log, ROUTER_V3, "Adding %s back to heap [delay %g upstream_R %g] [kd=%g d=%g prev=%d]\n", buffer, rt_node_p.delay, rt_node_p.upstream_R, kd, d, get_source(_g, prev));

							sources.push_back({ node, kd, d, prev });
						}
					} else {
						zlog_level(delta_log, ROUTER_V3, "Not reexpanding %s\n", buffer);
					}
				}

				//delta_stepping(_g, sources, sink->rr_node, delta, _known_distance, _distance, _prev_edge, [this] (const RREdge &e) -> pair<float, float> { return get_edge_weight(e); }, *this);
				dijkstra(_g, sources, sink->rr_node, _known_distance, _distance, _prev_edge, [this] (int rr_node) -> bool { return expand_node(rr_node); }, [this] (const RREdge &e) -> pair<float, float> { return get_edge_weight(e); }, *this);

				backtrack(sink->rr_node, rt, added_rr_nodes);

				assert(get_vertex_props(rt.graph, route_tree_get_rt_node(rt, sink->rr_node)).delay == _state[sink->rr_node].delay);
				net_timing.delay[sink->id+1] = _state[sink->rr_node].delay;

				for (const auto &n : _modified_nodes) {
					_known_distance[n] = std::numeric_limits<float>::max();
					_distance[n] = std::numeric_limits<float>::max();
					_prev_edge[n] = RRGraph::null_edge();
					_state[n].delay = std::numeric_limits<float>::max();
					_state[n].upstream_R = std::numeric_limits<float>::max();

					assert(_modified_node_added[n]);
					_modified_node_added[n] = false;
				}

				_modified_nodes.clear();

				//for (const auto &n : get_vertices(_g)) {
					//assert(_known_distance[n] == std::numeric_limits<float>::max() &&
							//_distance[n] == std::numeric_limits<float>::max() &&
							//_prev_edge[n] == RRGraph::null_edge() &&
							//_state[n].delay == std::numeric_limits<float>::max() &&
							//_state[n].upstream_R == std::numeric_limits<float>::max());
				//}

				++isink;
			}
		}
};

tbb::atomic<int> DeltaSteppingRouter::_num_instances = 0;

//bool partitioning_delta_stepping_deterministic_route_test(t_router_opts *opts)
//{
	//init_logging();
    //zlog_set_record("custom_output", concurrent_log_impl);

	//zlog_put_mdc("iter", "0");
	//zlog_put_mdc("tid", "0");

	//RRGraph g;
	//init_graph(g);

	//init_sprintf_rr_node(&g);

	//vector<net_t> nets;
	//vector<net_t> global_nets;
	//init_nets(nets, global_nets, opts->bb_factor, opts->large_bb);

	//vector<vector<net_t *>> phase_nets;
	//test_rtree(nets, phase_nets);
	//test_misr(nets);
	//return false;

	//congestion_t *congestion = new congestion_t[num_vertices(g)];
    //for (int i = 0; i < num_vertices(g); ++i) {
        //congestion[i].acc_cost = 1;
        //congestion[i].pres_cost = 1;
        //congestion[i].occ = 0;
    //}

	//route_tree_t *route_trees = new route_tree_t[nets.size()];
	//for (int i = 0; i < nets.size(); ++i) {
		//route_tree_init(route_trees[i]);
	//}

    //t_net_timing *net_timing = new t_net_timing[nets.size()+global_nets.size()];
    //init_net_timing(nets, global_nets, net_timing);

    //route_parameters_t params;
    //params.criticality_exp = opts->criticality_exp;
    //params.astar_fac = opts->astar_fac;
    //params.max_criticality = opts->max_criticality;
    //params.bend_cost = opts->bend_cost;

	////std::sort(begin(nets), end(nets), [] (const net_t &a, const net_t &b) -> bool {
			////return a.sinks.size() > b.sinks.size();
			////});

	//float pres_fac = opts->first_iter_pres_fac;

	//float delta;
	//for (const auto &n : get_vertices(g)) {
		//const auto &prop = get_vertex_props(g, n);
		//if (prop.type == CHANX) {
			//extern t_rr_indexed_data *rr_indexed_data;
			//delta = rr_indexed_data[prop.cost_index].base_cost;
			//break;
		//}
	//}

	//DeltaSteppingRouter router(g, congestion, pres_fac);

	//bool routed = false;
	//for (int iter = 0; iter < opts->max_router_iterations && !routed; ++iter) {
		//int inet = 0;
		//for (auto &net : nets) {
			//update_sink_criticalities(net, net_timing[net.vpr_id], params);

			//route_tree_mark_all_nodes_to_be_ripped(route_trees[net.local_id], g);

			//route_tree_rip_up_marked(route_trees[net.local_id], g, congestion, pres_fac);

			//vector<const sink_t *> sinks;
			//for (int i = 0; i < net.sinks.size(); ++i) {
				//sinks.push_back(&net.sinks[i]);
			//}

			//zlog_level(delta_log, ROUTER_V3, "Routing net %d\n", net.local_id);
			//router.route(&net.source, sinks, net.local_id, delta, opts->astar_fac, route_trees[net.local_id], net_timing[net.vpr_id]); 

			//++inet;
		//}

		//[> checking <]
		//for (int i = 0; i < num_vertices(g); ++i) {
			//congestion[i].recalc_occ = 0; 
		//}

		//for (const auto &net : nets) {
			//check_route_tree(route_trees[net.local_id], net, g);
			//recalculate_occ(route_trees[net.local_id], g, congestion);
		//}

		//unsigned long num_overused_nodes = 0;
		//vector<int> overused_nodes_by_type(NUM_RR_TYPES, 0);

		//if (feasible_routing(g, congestion)) {
			////dump_route(*current_traces_ptr, "route.txt");
			//routed = true;
		//} else {
			//for (int i = 0; i < num_vertices(g); ++i) {
				//if (congestion[i].occ > get_vertex_props(g, i).capacity) {
					//++num_overused_nodes;
					//const auto &v_p = get_vertex_props(g, i);
					//++overused_nodes_by_type[v_p.type];
				//}
			//}

			////auto update_cost_start = timer::now();

			//if (iter == 0) {
				//pres_fac = opts->initial_pres_fac;
				//update_costs(g, congestion, pres_fac, 0);
			//} else {
				//pres_fac *= opts->pres_fac_mult;

				//[> Avoid overflow for high iteration counts, even if acc_cost is big <]
				//pres_fac = std::min(pres_fac, static_cast<float>(HUGE_POSITIVE_FLOAT / 1e5));

				//update_costs(g, congestion, pres_fac, opts->acc_fac);
			//}

			////update_cost_time = timer::now()-update_cost_start;
		//}

		////auto analyze_timing_start = timer::now();

		//float crit_path_delay = analyze_timing(net_timing);

		////analyze_timing_time = timer::now()-analyze_timing_start;
		
		//printf("Overused: %lu Crit path delay: %g\n", num_overused_nodes, crit_path_delay);
	//}

	//return false;
//}

class Barrier
{
private:
    std::mutex _mutex;
    std::condition_variable _cv;
    std::size_t _count;
	std::size_t _initial_count;

public:
    explicit Barrier(std::size_t count)
		: _count{count}, _initial_count(count)
	{
	}

    void wait()
    {
        std::unique_lock<std::mutex> lock{_mutex};
        if (--_count == 0) {
			lock.unlock();
            _cv.notify_all();
        } else {
            _cv.wait(lock, [this] { return _count == 0; });
        }
    }

	void reset(std::size_t count)
	{
        std::unique_lock<std::mutex> lock{_mutex};
		_count = count;
	}

	void reset()
	{
		reset(_initial_count);
	}
};

using heap_item_t = pair<int, new_virtual_net_t *>;
using heap_t = boost::heap::binomial_heap<heap_item_t, boost::heap::stable<true>, boost::heap::compare<std::greater<heap_item_t>>>;

template<typename T, std::size_t Align>
struct alignas(Align) aligned_atomic {
	std::atomic<T> val;
};

typedef struct worker_sync_t {
	/*alignas(64)*/ tbb::atomic<int> stop_routing;
	/*alignas(64)*/ tbb::atomic<int> net_index;
#ifdef PTHREAD_BARRIER
	/*alignas(64)*/ pthread_barrier_t barrier;
#else
	/*alignas(64)*/ boost::barrier *barrier;
	//alignas(64) SpinningBarrier *barrier;
#endif
	//std::atomic<int> nums_net_to_route;
	tbb::spin_mutex lock;
	tbb::concurrent_queue<new_virtual_net_t *> global_pending_nets;
	//std::atomic<bool> has_global;
	aligned_atomic<int, 64> *current_num_incoming_edges;
	moodycamel::BlockingConcurrentQueue<new_virtual_net_t *> fast_global_pending_nets;
	heap_t num_incoming_edges;
	vector<heap_t::handle_type> handles;
} worker_sync_t;

typedef struct global_route_state_t {
	congestion_t *congestion;
	route_tree_t *route_trees;
	//route_tree_t *back_route_trees;
	vector<RRNode> *added_rr_nodes;	
	vector<RRNode> *back_added_rr_nodes;	
	t_net_timing *net_timing;
} global_route_state_t;

template<typename Net>
struct /*alignas(64)*/ router_t {
	vector<vector<Net>> nets;
	vector<Net> topo_nets;
	vector<int> num_incoming_edges;
	vector<vector<int>> directed;
	worker_sync_t sync;
	RRGraph g;
	global_route_state_t state;
	route_parameters_t params;
	float pres_fac;
	float delta;
	int num_threads;
	int num_virtual_nets;
	vector<dijkstra_stats_t> stats;
	vector<timer::duration> route_time;
	vector<timer::duration> wait_time;
	vector<timer::duration> pure_route_time;
	vector<timer::duration> update_time;
	vector<timer::duration> last_barrier_wait_time;
	vector<vector<timer::duration>> phase_time;
	vector<vector<timer::duration>> pure_phase_time;
	vector<timer::duration> net_route_time;
	vector<int> net_router;
	vector<dijkstra_stats_t> net_stats;
	tbb::atomic<int> iter;
};

//class Worker {
	//private:
		//const vector<vector<net_t *>> &_nets;
		//worker_sync_t _sync;
		//shared_route &_r;
		//float &_pres_fac;

		//DeltaSteppingRouter _router;

	//public:
		//Worker(tbb::atomic<bool> &stop_routing, const vector<vector<net_t *>> &nets, tbb::atomic<int> &net_index, pthread_barrier_t &barrier, shared_route &r, float &pres_fac)
			//: _stop_routing(stop_routing), _nets(nets), _net_index(net_index), _barrier(barrier), _r(r), _pres_fac(pres_fac), _router(r.g, r.congestion, pres_fac)
		//{
		//}

		//void main(const route_parameters_t &params, float delta)
		//{
		//}

		//void operator()(const route_parameters_t &params, float delta)
		//{
			//while (!_stop_routing) {
				//main(params, delta);
			//}
		//}
//};
//

void do_virtual_work_topo_4(router_t<new_virtual_net_t *> *router, DeltaSteppingRouter &net_router)
{
#ifdef PTHREAD_BARRIER
	pthread_barrier_wait(&router->sync.barrier);
#else
	router->sync.barrier->wait();
#endif

	assert(net_router.get_num_instances() == router->num_threads);

	char buffer[256];

	sprintf(buffer, "%d", router->iter);
	zlog_put_mdc("iter", buffer);

	vector<timer::duration> net_route_time(router->topo_nets.size(), timer::duration::max());
	vector<dijkstra_stats_t> net_stats(router->topo_nets.size(),
			{ std::numeric_limits<int>::max(), std::numeric_limits<int>::max(), std::numeric_limits<int>::max() });
	vector<int> routed_vnets;

	auto pure_route_time = timer::duration::zero();
	auto wait_time = timer::duration::zero();
	auto update_time = timer::duration::zero();

	auto route_start = timer::now();

	int inet;

	queue<new_virtual_net_t *> local_pending_nets;

	while ((inet =router->sync.net_index++) < router->topo_nets.size()) {
		new_virtual_net_t *vnet = nullptr;

		if (!local_pending_nets.empty()) {
			vnet = local_pending_nets.front();
			local_pending_nets.pop();
		}

		if (!vnet) {
			/* get from global */
			auto wait_start = timer::now();
			/* wait */
			//bool expected = true;
			//while (!router->sync.has_global.compare_exchange_weak(expected, false)) {
			//}
			//router->sync.fast_global_pending_nets.wait_dequeue(vnet);
			while (!router->sync.global_pending_nets.try_pop(vnet)) {
			}
			//router->sync.global_pending_nets.pop(vnet);

			wait_time += timer::now()-wait_start;
		}

		assert(vnet);

		//router->sync.lock.lock();

		//assert(!router->sync.global_pending_nets.empty());
		//vnet = router->sync.global_pending_nets.front();
		//router->sync.global_pending_nets.pop();

		//router->sync.lock.unlock();

		auto &net = *vnet->net;

		auto real_net_route_start = timer::now();

		zlog_level(delta_log, ROUTER_V3, "Routing net %d vnet global %d vnet %d num sinks %d l segment %d\n", net.local_id, vnet->global_index, vnet->index, vnet->sinks.size(), vnet->l_segment ? 1 : 0);

		if (router->iter > 0) {
			//for (const auto &rt_node : route_tree_get_nodes(router->state.back_route_trees[net.local_id])) {
			//RRNode rr_node = get_vertex_props(router->state.back_route_trees[net.local_id].graph, rt_node).rr_node;
			//const auto &rr_node_p = get_vertex_props(router->g, rr_node);

			//if (rr_node_p.xhigh >= _current_sink->current_bounding_box.xmin
			//&& rr_node_p.xlow <= _current_sink->current_bounding_box.xmax
			//&& rr_node_p.yhigh >= _current_sink->current_bounding_box.ymin
			//&& rr_node_p.ylow <= _current_sink->current_bounding_box.ymax) {
			//update_one_cost_internal(rr_node, router->g, router->state.congestion, -1, router->pres_fac);
			//}
			//}

			for (const auto &rr_node : router->state.back_added_rr_nodes[vnet->global_index]) {
				const auto &rr_node_p = get_vertex_props(router->g, rr_node);
				//assert(rr_node_p.xlow >= vnet->sinks[0]->current_bounding_box.xmin
						//&& rr_node_p.xhigh <= vnet->sinks[0]->current_bounding_box.xmax
						//&& rr_node_p.ylow >= vnet->sinks[0]->current_bounding_box.ymin
						//&& rr_node_p.yhigh <= vnet->sinks[0]->current_bounding_box.ymax);
				//assert(rr_node_p.xhigh >= vnet->sinks[0]->current_bounding_box.xmin
						//&& rr_node_p.xlow <= vnet->sinks[0]->current_bounding_box.xmax
						//&& rr_node_p.yhigh >= vnet->sinks[0]->current_bounding_box.ymin
						//&& rr_node_p.ylow <= vnet->sinks[0]->current_bounding_box.ymax);
				assert(!outside_bb(rr_node_p, vnet->bounding_boxes));
				update_one_cost_internal(rr_node, router->g, router->state.congestion, -1, router->pres_fac);
			}
		}

		vector<const sink_t *> sinks;
		for (int i = 0; i < vnet->sinks.size(); ++i) {
			sinks.push_back(vnet->sinks[i]);
		}

		source_t *source;
		if (vnet->index == 0) {
			source = &net.source;
		} else {
			source = nullptr;
		}

		net_router.reset_stats();

		net_router.route(source, sinks, router->delta, router->params.astar_fac, router->state.route_trees[net.local_id], router->state.added_rr_nodes[vnet->global_index], router->state.net_timing[net.vpr_id]); 

		net_route_time[vnet->global_index] = timer::now()-real_net_route_start;
		pure_route_time += net_route_time[vnet->global_index];

		net_stats[vnet->global_index] = net_router.get_stats();

		auto update_start = timer::now();
		for (const auto &to : router->directed[vnet->global_index]) {
			int val = --router->sync.current_num_incoming_edges[to].val;
			assert(val >= 0);
			if (val == 0) {
				//router->sync.lock.lock();
				if (local_pending_nets.empty()) {
					local_pending_nets.push(router->topo_nets[to]);
				} else {
					//router->sync.fast_global_pending_nets.enqueue(router->topo_nets[to]);
					router->sync.global_pending_nets.push(router->topo_nets[to]);
				}
				//router->sync.lock.unlock();
				//router->sync.has_global = true;
			}
		}
		update_time += timer::now() - update_start;

		routed_vnets.push_back(vnet->global_index);
	}

	assert(local_pending_nets.size() == 0 || local_pending_nets.size() == 1);
	while (!local_pending_nets.empty()) {
		//router->sync.fast_global_pending_nets.enqueue(local_pending_nets.front());
		router->sync.global_pending_nets.push(local_pending_nets.front());
		local_pending_nets.pop();
	}

	auto last_barrier_wait_start = timer::now();
#ifdef PTHREAD_BARRIER 
	pthread_barrier_wait(&router->sync.barrier);
#else
	router->sync.barrier->wait();
#endif
	router->last_barrier_wait_time[net_router.get_instance()] = timer::now() - last_barrier_wait_start;

	router->route_time[net_router.get_instance()] = timer::now()-route_start;
	router->wait_time[net_router.get_instance()] = wait_time;
	router->pure_route_time[net_router.get_instance()] = pure_route_time; 
	router->update_time[net_router.get_instance()] = update_time; 

	router->stats[net_router.get_instance()].num_heap_pops = 0;
	router->stats[net_router.get_instance()].num_heap_pushes = 0;
	router->stats[net_router.get_instance()].num_neighbor_visits = 0;

	for (const auto &r : routed_vnets) {
		router->net_route_time[r] = net_route_time[r];
		router->net_router[r] = net_router.get_instance();
		router->net_stats[r] = net_stats[r];

		router->stats[net_router.get_instance()].num_heap_pops += net_stats[r].num_heap_pops;
		router->stats[net_router.get_instance()].num_heap_pushes += net_stats[r].num_heap_pushes;
		router->stats[net_router.get_instance()].num_neighbor_visits += net_stats[r].num_neighbor_visits;
	}

#ifdef PTHREAD_BARRIER 
	pthread_barrier_wait(&router->sync.barrier);
#else
	router->sync.barrier->wait();
#endif
}

void do_virtual_work_topo_3(router_t<new_virtual_net_t *> *router, DeltaSteppingRouter &net_router)
{
#ifdef PTHREAD_BARRIER
	pthread_barrier_wait(&router->sync.barrier);
#else
	router->sync.barrier->wait();
#endif

	assert(net_router.get_num_instances() == router->num_threads);

	char buffer[256];

	sprintf(buffer, "%d", router->iter);
	zlog_put_mdc("iter", buffer);

	vector<timer::duration> net_route_time(router->topo_nets.size(), timer::duration::max());
	vector<dijkstra_stats_t> net_stats(router->topo_nets.size(),
			{ std::numeric_limits<int>::max(), std::numeric_limits<int>::max(), std::numeric_limits<int>::max() });
	vector<int> routed_vnets;

#ifdef PTHREAD_BARRIER
	pthread_barrier_wait(&router->sync.barrier);
#else
	router->sync.barrier->wait();
#endif

	auto pure_route_time = timer::duration::zero();
	auto wait_time = timer::duration::zero();
	auto update_time = timer::duration::zero();

	auto route_start = timer::now();

	int inet;

	queue<new_virtual_net_t *> local_pending_nets;

	while ((inet =router->sync.net_index++) < router->topo_nets.size()) {
		new_virtual_net_t *vnet = nullptr;

		if (!local_pending_nets.empty()) {
			vnet = local_pending_nets.front();
			local_pending_nets.pop();
		}

		if (!vnet) {
			/* get from global */
			auto wait_start = timer::now();
			/* wait */
			//bool expected = true;
			//while (!router->sync.has_global.compare_exchange_weak(expected, false)) {
			//}
			while (!router->sync.global_pending_nets.try_pop(vnet)) {
			}
			wait_time += timer::now()-wait_start;
		}

		assert(vnet);

		//router->sync.lock.lock();

		//assert(!router->sync.global_pending_nets.empty());
		//vnet = router->sync.global_pending_nets.front();
		//router->sync.global_pending_nets.pop();

		//router->sync.lock.unlock();

		auto &net = *vnet->net;

		auto real_net_route_start = timer::now();

		zlog_level(delta_log, ROUTER_V3, "Routing net %d vnet %d\n", net.local_id, vnet->index);

		if (router->iter > 0) {
			//for (const auto &rt_node : route_tree_get_nodes(router->state.back_route_trees[net.local_id])) {
			//RRNode rr_node = get_vertex_props(router->state.back_route_trees[net.local_id].graph, rt_node).rr_node;
			//const auto &rr_node_p = get_vertex_props(router->g, rr_node);

			//if (rr_node_p.xhigh >= _current_sink->current_bounding_box.xmin
			//&& rr_node_p.xlow <= _current_sink->current_bounding_box.xmax
			//&& rr_node_p.yhigh >= _current_sink->current_bounding_box.ymin
			//&& rr_node_p.ylow <= _current_sink->current_bounding_box.ymax) {
			//update_one_cost_internal(rr_node, router->g, router->state.congestion, -1, router->pres_fac);
			//}
			//}

			for (const auto &rr_node : router->state.back_added_rr_nodes[vnet->global_index]) {
				const auto &rr_node_p = get_vertex_props(router->g, rr_node);
				//assert(rr_node_p.xlow >= vnet->sinks[0]->current_bounding_box.xmin
						//&& rr_node_p.xhigh <= vnet->sinks[0]->current_bounding_box.xmax
						//&& rr_node_p.ylow >= vnet->sinks[0]->current_bounding_box.ymin
						//&& rr_node_p.yhigh <= vnet->sinks[0]->current_bounding_box.ymax);
				assert(rr_node_p.xhigh >= vnet->sinks[0]->current_bounding_box.xmin
						&& rr_node_p.xlow <= vnet->sinks[0]->current_bounding_box.xmax
						&& rr_node_p.yhigh >= vnet->sinks[0]->current_bounding_box.ymin
						&& rr_node_p.ylow <= vnet->sinks[0]->current_bounding_box.ymax);
				update_one_cost_internal(rr_node, router->g, router->state.congestion, -1, router->pres_fac);
			}
		}

		vector<const sink_t *> sinks;
		for (int i = 0; i < vnet->sinks.size(); ++i) {
			sinks.push_back(vnet->sinks[i]);
		}

		source_t *source;
		if (vnet->index == 0) {
			source = &net.source;
		} else {
			source = nullptr;
		}

		net_router.reset_stats();

		net_router.route(source, sinks, router->delta, router->params.astar_fac, router->state.route_trees[net.local_id], router->state.added_rr_nodes[vnet->global_index], router->state.net_timing[net.vpr_id]); 

		//rand_delay();

		net_route_time[vnet->global_index] = timer::now()-real_net_route_start;
		pure_route_time += net_route_time[vnet->global_index];

		net_stats[vnet->global_index] = net_router.get_stats();

		auto update_start = timer::now();
		for (const auto &to : router->directed[vnet->global_index]) {
			int val = --router->sync.current_num_incoming_edges[to].val;
			assert(val >= 0);
			if (val == 0) {
				//router->sync.lock.lock();
				if (local_pending_nets.empty()) {
					local_pending_nets.push(router->topo_nets[to]);
				} else {
					router->sync.global_pending_nets.push(router->topo_nets[to]);
				}
				//router->sync.lock.unlock();
				//router->sync.has_global = true;
			}
		}
		update_time += timer::now() - update_start;

		routed_vnets.push_back(vnet->global_index);
	}

	while (!local_pending_nets.empty()) {
		router->sync.global_pending_nets.push(local_pending_nets.front());
		local_pending_nets.pop();
	}

	auto last_barrier_wait_start = timer::now();
#ifdef PTHREAD_BARRIER 
	pthread_barrier_wait(&router->sync.barrier);
#else
	router->sync.barrier->wait();
#endif
	router->last_barrier_wait_time[net_router.get_instance()] = timer::now() - last_barrier_wait_start;

	router->route_time[net_router.get_instance()] = timer::now()-route_start;
	router->wait_time[net_router.get_instance()] = wait_time;
	router->pure_route_time[net_router.get_instance()] = pure_route_time; 
	router->update_time[net_router.get_instance()] = update_time; 

	router->stats[net_router.get_instance()].num_heap_pops = 0;
	router->stats[net_router.get_instance()].num_heap_pushes = 0;
	router->stats[net_router.get_instance()].num_neighbor_visits = 0;

	for (const auto &r : routed_vnets) {
		router->net_route_time[r] = net_route_time[r];
		router->net_router[r] = net_router.get_instance();
		router->net_stats[r] = net_stats[r];

		router->stats[net_router.get_instance()].num_heap_pops += net_stats[r].num_heap_pops;
		router->stats[net_router.get_instance()].num_heap_pushes += net_stats[r].num_heap_pushes;
		router->stats[net_router.get_instance()].num_neighbor_visits += net_stats[r].num_neighbor_visits;
	}

#ifdef PTHREAD_BARRIER 
	pthread_barrier_wait(&router->sync.barrier);
#else
	router->sync.barrier->wait();
#endif
}

void do_virtual_work_topo_2(router_t<new_virtual_net_t *> *router, DeltaSteppingRouter &net_router)
{
#ifdef PTHREAD_BARRIER
	pthread_barrier_wait(&router->sync.barrier);
#else
	router->sync.barrier->wait();
#endif

	assert(net_router.get_num_instances() == router->num_threads);

	char buffer[256];

	sprintf(buffer, "%d", router->iter);
	zlog_put_mdc("iter", buffer);

	vector<timer::duration> net_route_time(router->topo_nets.size(), timer::duration::max());
	vector<int> routed_vnets;

	net_router.reset_stats();

#ifdef PTHREAD_BARRIER
	pthread_barrier_wait(&router->sync.barrier);
#else
	router->sync.barrier->wait();
#endif

	auto pure_route_time = timer::duration::zero();
	auto wait_time = timer::duration::zero();
	auto update_time = timer::duration::zero();

	auto route_start = timer::now();

	int inet;

	while (router->sync.net_index++ < router->topo_nets.size()) {
		/* wait */
		new_virtual_net_t *vnet = nullptr;

		auto wait_start = timer::now();
		while (!vnet) {
			//tbb::spin_mutex::scoped_lock lock(router->sync.lock);

			router->sync.lock.lock();
			if (router->sync.num_incoming_edges.top().first == 0) {
				auto item = router->sync.num_incoming_edges.top();
				router->sync.num_incoming_edges.pop();

				assert(item.first == 0);
				vnet = item.second;
			}
			router->sync.lock.unlock();
		}
		wait_time += timer::now()-wait_start;

		auto &net = *vnet->net;

		auto real_net_route_start = timer::now();

		update_sink_criticalities(net, router->state.net_timing[net.vpr_id], router->params);

		zlog_level(delta_log, ROUTER_V3, "Routing net %d vnet %d\n", net.local_id, vnet->index);

		if (router->iter > 0) {
			//for (const auto &rt_node : route_tree_get_nodes(router->state.back_route_trees[net.local_id])) {
			//RRNode rr_node = get_vertex_props(router->state.back_route_trees[net.local_id].graph, rt_node).rr_node;
			//const auto &rr_node_p = get_vertex_props(router->g, rr_node);

			//if (rr_node_p.xhigh >= _current_sink->current_bounding_box.xmin
			//&& rr_node_p.xlow <= _current_sink->current_bounding_box.xmax
			//&& rr_node_p.yhigh >= _current_sink->current_bounding_box.ymin
			//&& rr_node_p.ylow <= _current_sink->current_bounding_box.ymax) {
			//update_one_cost_internal(rr_node, router->g, router->state.congestion, -1, router->pres_fac);
			//}
			//}

			for (const auto &rr_node : router->state.back_added_rr_nodes[vnet->global_index]) {
				update_one_cost_internal(rr_node, router->g, router->state.congestion, -1, router->pres_fac);
			}
		}

		vector<const sink_t *> sinks;
		for (int i = 0; i < vnet->sinks.size(); ++i) {
			sinks.push_back(vnet->sinks[i]);
		}

		source_t *source;
		if (vnet->index == 0) {
			source = &net.source;
		} else {
			source = nullptr;
		}

		net_router.route(source, sinks, router->delta, router->params.astar_fac, router->state.route_trees[net.local_id], router->state.added_rr_nodes[vnet->global_index], router->state.net_timing[net.vpr_id]); 

		net_route_time[vnet->global_index] = timer::now()-real_net_route_start;
		pure_route_time += net_route_time[vnet->global_index];

		{
			router->sync.lock.lock();
			auto update_start = timer::now();
			for (const auto &to : router->directed[vnet->global_index]) {
				//tbb::spin_mutex::scoped_lock lock(router->sync.lock);
				assert((*router->sync.handles[to]).first > 0);
				assert((*router->sync.handles[to]).second->global_index == to);
				assert(!router->sync.num_incoming_edges.empty());

				router->sync.num_incoming_edges.update(router->sync.handles[to], { (*router->sync.handles[to]).first-1, router->topo_nets[to] });
			}
			update_time += timer::now() - update_start;
			router->sync.lock.unlock();
		}

		routed_vnets.push_back(vnet->global_index);
	}

	auto last_barrier_wait_start = timer::now();
#ifdef PTHREAD_BARRIER 
	pthread_barrier_wait(&router->sync.barrier);
#else
	router->sync.barrier->wait();
#endif
	router->last_barrier_wait_time[net_router.get_instance()] = timer::now() - last_barrier_wait_start;

	assert(router->sync.num_incoming_edges.empty());

	router->route_time[net_router.get_instance()] = timer::now()-route_start;
	router->wait_time[net_router.get_instance()] = wait_time;
	router->pure_route_time[net_router.get_instance()] = pure_route_time; 
	router->update_time[net_router.get_instance()] = update_time; 

	for (const auto &r : routed_vnets) {
		router->net_route_time[r] = net_route_time[r];
		router->net_router[r] = net_router.get_instance();
	}

	router->stats[net_router.get_instance()] = net_router.get_stats();

#ifdef PTHREAD_BARRIER 
	pthread_barrier_wait(&router->sync.barrier);
#else
	router->sync.barrier->wait();
#endif
}

void do_virtual_work_topo(router_t<new_virtual_net_t *> *router, DeltaSteppingRouter &net_router)
{
#ifdef PTHREAD_BARRIER
	pthread_barrier_wait(&router->sync.barrier);
#else
	router->sync.barrier->wait();
#endif

	assert(net_router.get_num_instances() == router->num_threads);

	char buffer[256];

	sprintf(buffer, "%d", router->iter);
	zlog_put_mdc("iter", buffer);

	vector<timer::duration> net_route_time(router->topo_nets.size(), timer::duration::max());
	vector<int> routed_vnets;

	net_router.reset_stats();

	router->sync.net_index = 0;

#ifdef PTHREAD_BARRIER
	pthread_barrier_wait(&router->sync.barrier);
#else
	router->sync.barrier->wait();
#endif

	auto pure_route_time = timer::duration::zero();
	auto route_start = timer::now();

	int inet;

	while ((inet = router->sync.net_index++) < router->topo_nets.size()) {
		auto &vnet = *router->topo_nets[inet];
		auto &net = *vnet.net;

		/* wait */
		//while (router->sync.current_num_incoming_edges[vnet.global_index].load(std::memory_order_acquire) > 0) {
			//;
		//}

		//router->sync.current_num_incoming_edges[vnet.global_index] = vnet.num_incoming_edges;

		auto real_net_route_start = timer::now();

		update_sink_criticalities(net, router->state.net_timing[net.vpr_id], router->params);

		zlog_level(delta_log, ROUTER_V3, "Routing net %d vnet %d\n", net.local_id, vnet.index);

		if (router->iter > 0) {
			//for (const auto &rt_node : route_tree_get_nodes(router->state.back_route_trees[net.local_id])) {
			//RRNode rr_node = get_vertex_props(router->state.back_route_trees[net.local_id].graph, rt_node).rr_node;
			//const auto &rr_node_p = get_vertex_props(router->g, rr_node);

			//if (rr_node_p.xhigh >= _current_sink->current_bounding_box.xmin
			//&& rr_node_p.xlow <= _current_sink->current_bounding_box.xmax
			//&& rr_node_p.yhigh >= _current_sink->current_bounding_box.ymin
			//&& rr_node_p.ylow <= _current_sink->current_bounding_box.ymax) {
			//update_one_cost_internal(rr_node, router->g, router->state.congestion, -1, router->pres_fac);
			//}
			//}

			for (const auto &rr_node : router->state.back_added_rr_nodes[vnet.global_index]) {
				update_one_cost_internal(rr_node, router->g, router->state.congestion, -1, router->pres_fac);
			}
		}

		vector<const sink_t *> sinks;
		for (int i = 0; i < vnet.sinks.size(); ++i) {
			sinks.push_back(vnet.sinks[i]);
		}

		source_t *source;
		if (vnet.index == 0) {
			source = &net.source;
		} else {
			source = nullptr;
		}

		net_router.route(source, sinks, router->delta, router->params.astar_fac, router->state.route_trees[net.local_id], router->state.added_rr_nodes[vnet.global_index], router->state.net_timing[net.vpr_id]); 

		for (const auto &to : router->directed[vnet.global_index]) {
			//--router->sync.current_num_incoming_edges[to];
		}

		net_route_time[vnet.global_index] = timer::now()-real_net_route_start;
		pure_route_time += net_route_time[vnet.global_index];
		routed_vnets.push_back(vnet.global_index);
	}

#ifdef PTHREAD_BARRIER 
	pthread_barrier_wait(&router->sync.barrier);
#else
	router->sync.barrier->wait();
#endif

	router->route_time[net_router.get_instance()] = timer::now()-route_start;
	router->pure_route_time[net_router.get_instance()] = pure_route_time; 

	for (const auto &r : routed_vnets) {
		router->net_route_time[r] = net_route_time[r];
		router->net_router[r] = net_router.get_instance();
	}

	router->stats[net_router.get_instance()] = net_router.get_stats();

#ifdef PTHREAD_BARRIER 
	pthread_barrier_wait(&router->sync.barrier);
#else
	router->sync.barrier->wait();
#endif
}

void do_virtual_work(router_t<new_virtual_net_t *> *router, DeltaSteppingRouter &net_router)
{
#ifdef PTHREAD_BARRIER
	pthread_barrier_wait(&router->sync.barrier);
#else
	router->sync.barrier->wait();
#endif

	assert(net_router.get_num_instances() == router->num_threads);

	char buffer[256];

	sprintf(buffer, "%d", router->iter);
	zlog_put_mdc("iter", buffer);

	net_router.reset_stats();

	vector<timer::duration> phase_time(router->nets.size(), timer::duration::max());
	vector<timer::duration> pure_phase_time(router->nets.size(), timer::duration::max());
	vector<timer::duration> net_route_time(router->num_virtual_nets, timer::duration::max());
	vector<int> routed_vnets;

	auto pure_route_time = timer::duration::zero();

	auto route_start = timer::now();

	for (int phase = 0; phase < router->nets.size(); ++phase) {
		int inet;

		router->sync.net_index = 0;

#ifdef PTHREAD_BARRIER
		pthread_barrier_wait(&router->sync.barrier);
#else
		router->sync.barrier->wait();
#endif

		zlog_level(delta_log, ROUTER_V2, "Phase %d\n", phase);

		auto pure_phase_start = timer::now();

		while ((inet = router->sync.net_index++) < router->nets[phase].size()) {
			auto net_route_start = timer::now();

			auto &vnet = *router->nets[phase][inet];
			auto &net = *vnet.net;

			update_sink_criticalities(net, router->state.net_timing[net.vpr_id], router->params);

			zlog_level(delta_log, ROUTER_V3, "Routing net %d vnet %d\n", net.local_id, vnet.index);

			if (router->iter > 0) {
				//for (const auto &rt_node : route_tree_get_nodes(router->state.back_route_trees[net.local_id])) {
					//RRNode rr_node = get_vertex_props(router->state.back_route_trees[net.local_id].graph, rt_node).rr_node;
					//const auto &rr_node_p = get_vertex_props(router->g, rr_node);

					//if (rr_node_p.xhigh >= _current_sink->current_bounding_box.xmin
							//&& rr_node_p.xlow <= _current_sink->current_bounding_box.xmax
							//&& rr_node_p.yhigh >= _current_sink->current_bounding_box.ymin
							//&& rr_node_p.ylow <= _current_sink->current_bounding_box.ymax) {
						//update_one_cost_internal(rr_node, router->g, router->state.congestion, -1, router->pres_fac);
					//}
				//}
				
				for (const auto &rr_node : router->state.back_added_rr_nodes[vnet.global_index]) {
					update_one_cost_internal(rr_node, router->g, router->state.congestion, -1, router->pres_fac);
				}
			}

			vector<const sink_t *> sinks;
			for (int i = 0; i < vnet.sinks.size(); ++i) {
				sinks.push_back(vnet.sinks[i]);
			}

			source_t *source;
			if (vnet.index == 0) {
				source = &net.source;
			} else {
				source = nullptr;
			}

			net_router.route(source, sinks, router->delta, router->params.astar_fac, router->state.route_trees[net.local_id], router->state.added_rr_nodes[vnet.global_index], router->state.net_timing[net.vpr_id]); 

			net_route_time[vnet.global_index] = timer::now()-net_route_start;
			routed_vnets.push_back(vnet.global_index);
		}

		pure_phase_time[phase] = timer::now()-pure_phase_start; 
		pure_route_time += pure_phase_time[phase];

#ifdef PTHREAD_BARRIER 
		pthread_barrier_wait(&router->sync.barrier);
#else
		router->sync.barrier->wait();
#endif

		phase_time[phase] = timer::now()-pure_phase_start;
	}

	router->route_time[net_router.get_instance()] = timer::now()-route_start;
	router->pure_route_time[net_router.get_instance()] = pure_route_time; 

	router->phase_time[net_router.get_instance()] = std::move(phase_time); 
	router->pure_phase_time[net_router.get_instance()] = std::move(pure_phase_time); 

	for (const auto &r : routed_vnets) {
		router->net_route_time[r] = net_route_time[r];
		router->net_router[r] = net_router.get_instance();
	}

	router->stats[net_router.get_instance()] = net_router.get_stats();

#ifdef PTHREAD_BARRIER 
	pthread_barrier_wait(&router->sync.barrier);
#else
	router->sync.barrier->wait();
#endif
}

void do_work(router_t<net_t *> *router, DeltaSteppingRouter &net_router)
{
	char buffer[256];

	sprintf(buffer, "%d", router->iter);
	zlog_put_mdc("iter", buffer);

#ifdef PTHREAD_BARRIER
	pthread_barrier_wait(&router->sync.barrier);
#else
	router->sync.barrier->wait();
#endif

	assert(net_router.get_num_instances() == router->num_threads);

	net_router.reset_stats();

	auto pure_route_time = timer::duration::zero();

	auto route_start = timer::now();

	for (int phase = 0; phase < router->nets.size(); ++phase) {
		int inet;

		router->sync.net_index = 0;

#ifdef PTHREAD_BARRIER
		pthread_barrier_wait(&router->sync.barrier);
#else
		router->sync.barrier->wait();
#endif

		//auto pure_route_start = timer::now();

		while ((inet = router->sync.net_index++) < router->nets[phase].size()) {
			auto &net = *router->nets[phase][inet];

			update_sink_criticalities(net, router->state.net_timing[net.vpr_id], router->params);

			route_tree_mark_all_nodes_to_be_ripped(router->state.route_trees[net.local_id], router->g);

			route_tree_rip_up_marked(router->state.route_trees[net.local_id], router->g, router->state.congestion, router->pres_fac);

			vector<const sink_t *> sinks;
			for (int i = 0; i < net.sinks.size(); ++i) {
				sinks.push_back(&net.sinks[i]);
			}

			zlog_level(delta_log, ROUTER_V3, "Routing net %d\n", net.local_id);
			net_router.route(&net.source, sinks, router->delta, router->params.astar_fac, router->state.route_trees[net.local_id], router->state.added_rr_nodes[net.local_id], router->state.net_timing[net.vpr_id]); 
		}

		//pure_route_time += timer::now()-pure_route_start;

#ifdef PTHREAD_BARRIER
		pthread_barrier_wait(&router->sync.barrier);
#else
		router->sync.barrier->wait();
#endif
	}

	router->route_time[net_router.get_instance()] = timer::now()-route_start;
	router->pure_route_time[net_router.get_instance()] = pure_route_time; 

	router->stats[net_router.get_instance()] = net_router.get_stats();
}

void *virtual_topo_worker_thread(void *args)
{
	router_t<new_virtual_net_t *> *router = static_cast<router_t<new_virtual_net_t *> *>(args);

#ifdef PTHREAD_BARRIER
	pthread_barrier_wait(&router->sync.barrier);
#else
	router->sync.barrier->wait();
#endif

	DeltaSteppingRouter net_router(router->g, router->state.congestion, router->pres_fac);

	char buffer[256];
	sprintf(buffer, "%d", net_router.get_instance());
	zlog_put_mdc("tid", buffer);
	zlog_put_mdc("log_dir", dirname);

	while (!router->sync.stop_routing) {
		do_virtual_work_topo_4(router, net_router);
	}

	return nullptr;
}

void *virtual_worker_thread(void *args)
{
	router_t<new_virtual_net_t *> *router = (router_t<new_virtual_net_t *> *)args;

#ifdef PTHREAD_BARRIER
	pthread_barrier_wait(&router->sync.barrier);
#else
	router->sync.barrier->wait();
#endif

	DeltaSteppingRouter net_router(router->g, router->state.congestion, router->pres_fac);

	char buffer[256];
	sprintf(buffer, "%d", net_router.get_instance());
	zlog_put_mdc("tid", buffer);

	while (!router->sync.stop_routing) {
		do_virtual_work(router, net_router);
	}

	return nullptr;
}

void *worker_thread(void *args)
{
	router_t<net_t *> *router = (router_t<net_t *> *)args;

	DeltaSteppingRouter net_router(router->g, router->state.congestion, router->pres_fac);

	char buffer[256];
	sprintf(buffer, "%d", net_router.get_instance());
	zlog_put_mdc("tid", buffer);

	while (!router->sync.stop_routing) {
		do_work(router, net_router);
	}

	return nullptr;
}

void write_nets(const vector<net_t> &nets)
{
	auto nets_copy = nets;

	std::sort(begin(nets_copy), end(nets_copy), [] (const net_t &a, const net_t &b) -> bool {
			return a.sinks.size() > b.sinks.size();
			});

	for (int i = 0; i < 5; ++i) {
		char buffer[256];

		sprintf(buffer, "net_sinks/net_%d.txt", i);
		FILE *file = fopen(buffer, "w");

		for (const auto &sink : nets_copy[i].sinks) {
			fprintf(file, "%d %d\n", sink.x, sink.y);
		}

		fclose(file);
	}
}

template<typename Net, typename ToBox, typename ToNumSinks>
void write_net_stats(const vector<Net> &nets, const char *filename, const ToBox &to_box, const ToNumSinks &to_num_sinks)
{
	FILE *file = fopen(filename, "w");
	for (const auto &net : nets) {
		const box &bb = to_box(net);
		fprintf(file, "%g %d\n", bg::area(bb), to_num_sinks(net));
	}
	fclose(file);
}

void write_routes(const char *fname, const route_tree_t *route_trees, int num_nets)
{
	FILE *routes = fopen(fname, "w");

	for (int i = 0; i < num_nets; ++i) {
		fprintf(routes, "Net %d\n", i);
		for (const auto &node : route_tree_get_nodes(route_trees[i])) {
			const auto &props = get_vertex_props(route_trees[i].graph, node);
			fprintf(routes, "%d\n", props.rr_node);
		}
	}

	fclose(routes);
}

void write_net_dijkstra_stats(const char *fname, const vector<dijkstra_stats_t> &net_stats)
{
	FILE *stats = fopen(fname, "w");
	int i = 0;
	for (const auto &s : net_stats) {
		fprintf(stats, "Net %d\n", i);
		fprintf(stats, "%d\n", s.num_heap_pops);
		fprintf(stats, "%d\n", s.num_neighbor_visits);
		fprintf(stats, "%d\n", s.num_heap_pushes);

		++i;
	}
	
	fclose(stats);
}

void write_congestion_state(const char *fname, const congestion_t *congestion, int num_rr_nodes)
{
	FILE *state = fopen(fname, "w");

	for (int i = 0; i < num_rr_nodes; ++i) {
		fprintf(state, "%d %d\n", i, congestion[i].occ);
	}
	
	fclose(state);
}

void write_virtual_nets(const char *dirname, const vector<vector<new_virtual_net_t>> &all_virtual_nets)
{
	vector<const vector<new_virtual_net_t> *> sorted_all_virtual_nets;
	for (const auto &virtual_nets : all_virtual_nets) {
		sorted_all_virtual_nets.push_back(&virtual_nets);
	}
	std::sort(begin(sorted_all_virtual_nets), end(sorted_all_virtual_nets), [] (const vector<new_virtual_net_t> *a, const vector<new_virtual_net_t> *b) -> bool {
			return a->at(0).net->sinks.size() > b->at(0).net->sinks.size();
			});

	const auto &virtual_nets = *sorted_all_virtual_nets[0];
	for (int i = 0; i < virtual_nets.size(); ++i) {
		char buffer[256];

		sprintf(buffer, "%s/virtual_nets_0_%d_bb.txt", dirname, i);
		FILE *vnet_bb = fopen(buffer, "w");
		const auto &box = virtual_nets[i].bounding_box;
		extern int nx, ny;
		fprintf(vnet_bb, "0 0 %d %d 0\n", nx+2, ny+2);
		fprintf(vnet_bb, "%d %d %d %d 0\n",
				box.min_corner().get<0>(), box.min_corner().get<1>(),
				box.max_corner().get<0>()-box.min_corner().get<0>(),
				box.max_corner().get<1>()-box.min_corner().get<1>());
		//if (i > 0) {
			//const auto &prev_box = virtual_nets[i-1].bounding_box;
			//fprintf(vnet_bb, "%d %d %d %d 1\n",
					//prev_box.min_corner().get<0>(), prev_box.min_corner().get<1>(),
					//prev_box.max_corner().get<0>()-prev_box.min_corner().get<0>(),
					//prev_box.max_corner().get<1>()-prev_box.min_corner().get<1>());
		//}
		fclose(vnet_bb);

		sprintf(buffer, "%s/virtual_nets_0_%d_p.txt", dirname, i);
		FILE *vnet_p = fopen(buffer, "w");

		fprintf(vnet_p, "%d %d %d 0\n", virtual_nets[i].net->source.x, virtual_nets[i].net->source.y, i);
		for (int j = 0; j < virtual_nets[i].sinks.size(); ++j) {
			const auto &sink = virtual_nets[i].sinks[j];
			//if (j == virtual_nets[i].sinks.size()-1) {
				//assert(virtual_nets[i].last_point.get<0>() == sink->x &&virtual_nets[i].last_point.get<1>() == sink->y);
				//fprintf(vnet_p, "%d %d %d 1\n", sink->x, sink->y, i);
			//} else {
				//fprintf(vnet_p, "%d %d %d 0\n", sink->x, sink->y, i);
			//}
				fprintf(vnet_p, "%d %d %d 0\n", sink->x, sink->y, i);
		}

		fprintf(vnet_p, "%d %d %d 1\n", bg::get<0>(virtual_nets[i].last_point), bg::get<1>(virtual_nets[i].last_point), i);

		fclose(vnet_p);
	}
}

int get_num_interpartition_nets(const vector<net_t> &nets, int num_partitions);

bool partitioning_delta_stepping_deterministic_route_virtual(t_router_opts *opts)
{
	char buffer[256];

	int dir_index;
	if (!opts->log_dir) {
		bool dir_exists;
		dir_index = 0;
		do {
			extern char *s_circuit_name;
			sprintf(dirname, "%s_stats_%d", s_circuit_name, dir_index);
			struct stat s;
			dir_exists = stat(dirname, &s) == 0 ? true : false;
			if (dir_exists) {
				++dir_index;
			}
		} while (dir_exists);
	} else {
		sprintf(dirname, opts->log_dir);

		struct stat s;
		if (stat(dirname, &s) == 0) {
			printf("log dir %s already exists\n", dirname);
			dir_index = -1;
		} else {
			dir_index = 0;
		}
	}

	if (dir_index >= 0) {
		if (mkdir(dirname, 0700) == -1) {
			printf("failed to create dir %s\n", dirname);
			exit(-1);
		}
	}

	init_logging(opts->num_threads);
    zlog_set_record("custom_output", concurrent_log_impl_2);
	sprintf(buffer, "%d", 0);
	zlog_put_mdc("iter", buffer);
	zlog_put_mdc("tid", buffer);
	zlog_put_mdc("log_dir", dirname);

	router_t<new_virtual_net_t *> router;

	//assert(((unsigned long)&router.sync.stop_routing & 63) == 0);
	//assert(((unsigned long)&router.sync.net_index & 63) == 0);
#ifdef __linux__
	//assert(((unsigned long)&router.sync.barrier & 63) == 0);
#endif

	router.num_threads = opts->num_threads;

	init_graph(router.g);

	init_sprintf_rr_node(&router.g);

	vector<net_t> nets;
	vector<net_t> global_nets;
	init_nets(nets, global_nets, opts->bb_factor, opts->large_bb);

	printf("Num_interpatition_nets [%d nets]: ", nets.size());
	for (int i = 2; i <= 32; i *= 2) {
		printf("%d ", get_num_interpartition_nets(nets, i));
	}
	printf("\n");
	exit(0);

	//vector<net_t> nets_copy = nets;
	//std::sort(begin(nets_copy), end(nets_copy), [] (const net_t &a, const net_t &b) -> bool {
			//return a.sinks.size() > b.sinks.size();
			//});
	//vector<net_t> only;
	//boost::copy(nets_copy | boost::adaptors::sliced(0, 5000), std::back_inserter(only));
	//boost::copy(nets_copy, std::back_inserter(only));

	vector<vector<new_virtual_net_t>> all_virtual_nets;
	vector<new_virtual_net_t *> all_virtual_nets_ptr;
	int num_virtual_nets = create_virtual_nets(nets, router.g, opts->bb_area_threshold_scale, all_virtual_nets, all_virtual_nets_ptr);

	write_virtual_nets(dirname, all_virtual_nets);

	printf("start dump\n");
	for (int i = 0; i < num_virtual_nets; ++i) {
		printf("%d %d %d\n", i, all_virtual_nets_ptr[i]->net->local_id, all_virtual_nets_ptr[i]->index);
	}
	printf("end dump\n");
	//best_case(nets, router.nets);
	
	//vector<vector<int>> overlap2;
	//vector<new_virtual_net_t *> all_virtual_nets_ptr2;

	//auto build_start = timer::now();

	//build_overlap_graph_3(all_virtual_nets, all_virtual_nets_ptr2, overlap2);

	//auto build_time = timer::now()-build_start;

	//printf("overlap graph build 3 time %g\n", duration_cast<nanoseconds>(build_time).count() / 1e9);

	//vector<vector<int>> overlap1;
	//vector<new_virtual_net_t *> all_virtual_nets_ptr1;
	//build_start = timer::now();
	//build_overlap_graph_2(all_virtual_nets, overlap1, all_virtual_nets_ptr1);
	//build_time = timer::now()-build_start;

	//printf("overlap graph build 2 time %g\n", duration_cast<nanoseconds>(build_time).count() / 1e9);

	//for (int i = 0; i < overlap1.size(); ++i) {
		//vector<int> o1_temp = overlap1[i];
		//vector<int> o2_temp = overlap2[i];

		//std::sort(begin(o1_temp), end(o1_temp));
		//std::sort(begin(o2_temp), end(o2_temp));

		//assert(o1_temp == o2_temp);
	//}

	vector<vector<int>> overlap;
	vector<vertex_descriptor> order;
	build_and_order(all_virtual_nets, all_virtual_nets_ptr, overlap, order);

	//int num_orig_directed_edges = 0;
	//for (int i = 0; i < all_virtual_nets_ptr.size(); ++i) {
		//assert(all_virtual_nets_ptr[i]->v == i);
		//assert(all_virtual_nets_ptr[i]->global_index == i);

		//if (all_virtual_nets_ptr[i]->index > 0) {
			//overlap[all_virtual_nets_ptr[i]->parent].push_back(all_virtual_nets_ptr[i]->v);
			//++num_orig_directed_edges;
		//} 
	//}

	vector<int> all_num_incoming_edges(overlap.size(), 0);
	for (int i = 0; i < overlap.size(); ++i) {
		assert(overlap[i].size() > 0);
		for (int j = 0; j < overlap[i].size(); ++j) {
			++all_num_incoming_edges[overlap[i][j]];
		}
	}
	for (int i = 0; i < overlap.size(); ++i) {
		assert(all_num_incoming_edges[i] > 0);
	}
	vector<vector<int>> directed;
	vector<int> num_incoming_edges;
	//schedule_virtual_nets_5(overlap, all_virtual_nets_ptr, [&all_virtual_nets_ptr] (int v) -> float { return all_virtual_nets_ptr[v]->sinks.size(); }, opts->num_threads, directed, num_incoming_edges);

	schedule_virtual_nets_6(overlap, all_virtual_nets_ptr, [&all_virtual_nets_ptr] (int v) -> float { return all_virtual_nets_ptr[v]->sinks.size(); }, opts->num_threads, router.directed, router.num_incoming_edges);

	//assert(directed.size() == router.directed.size());
	//for (int i = 0; i < directed.size(); ++i) {
		//auto sorted_directed = directed[i];
		//auto sorted_router_directed = router.directed[i];
		//std::sort(begin(sorted_directed), end(sorted_directed));
		//std::sort(begin(sorted_router_directed), end(sorted_router_directed));
		//assert(sorted_directed == sorted_router_directed);
	//}

	//assert(num_incoming_edges.size() == router.num_incoming_edges.size());
	//for (int i = 0; i < num_incoming_edges.size(); ++i) {
		//assert(num_incoming_edges[i] == router.num_incoming_edges[i]);
	//}

	//int num_edges = 0;
	//for (int i = 0; i < overlap.size(); ++i) {
		//num_edges += overlap[i].size();
	//}
	
	//int num_directed_edges = 0;
	//for (int i = 0; i < router.directed.size(); ++i) {
		//num_directed_edges += router.directed[i].size();
	//}

	//assert(((num_edges-num_orig_directed_edges) % 2) == 0);
	//assert(num_directed_edges-num_orig_directed_edges == (num_edges-num_orig_directed_edges) / 2);

	router.topo_nets = all_virtual_nets_ptr;
	//schedule_virtual_nets_3(overlap, all_virtual_nets_ptr, order, opts->num_threads, nullptr, router.nets);
	//schedule_virtual_nets_4(overlap, all_virtual_nets_ptr, order, opts->num_threads, nullptr, router.directed, router.topo_nets);
	//test_misr(nets);
	//write_nets(nets);
	router.num_virtual_nets = num_virtual_nets;

	//sprintf(buffer, "%s/virtual_nets_bb_%g.txt", dirname, opts->bb_area_threshold_scale);
	//write_net_stats(all_virtual_nets_ptr, buffer, [] (const new_virtual_net_t *vnet) -> box { return vnet->bounding_box; }, [] (const new_virtual_net_t *vnet) -> int { return vnet->net->sinks.size(); });

	sprintf(buffer, "%s/nets_bb.txt", dirname);
	write_net_stats(nets, buffer, [] (const net_t &net) -> box { return box(point(net.bounding_box.xmin, net.bounding_box.ymin), point(net.bounding_box.xmax, net.bounding_box.ymax)); }, [] (const net_t &net) -> int { return net.sinks.size(); });

	congestion_t *congestion_aligned;
//#ifdef __linux__
	//assert(posix_memalign((void **)&congestion_aligned, 64, sizeof(congestion_t)*num_vertices(router.g)) == 0);
	//assert(((unsigned long)congestion_aligned & 63) == 0);
//#else
	//congestion_aligned = (congestion_t *)malloc(sizeof(congestion_t)*num_vertices(router.g));
//#endif 
	//router.state.congestion = new(congestion_aligned) congestion_t[num_vertices(router.g)];
	router.state.congestion = new congestion_t[num_vertices(router.g)];
    for (int i = 0; i < num_vertices(router.g); ++i) {
        router.state.congestion[i].acc_cost = 1;
        router.state.congestion[i].pres_cost = 1;
        router.state.congestion[i].occ = 0;
    }

	router.state.route_trees = new route_tree_t[nets.size()];
	//router.state.back_route_trees = new route_tree_t[nets.size()];
	for (int i = 0; i < nets.size(); ++i) {
		route_tree_init(router.state.route_trees[i]);
		//route_tree_init(router.state.back_route_trees[i]);
	}

	router.state.added_rr_nodes = new vector<RRNode>[num_virtual_nets];
	router.state.back_added_rr_nodes = new vector<RRNode>[num_virtual_nets];

    router.state.net_timing = new t_net_timing[nets.size()+global_nets.size()];
    init_net_timing(nets, global_nets, router.state.net_timing);

    router.params.criticality_exp = opts->criticality_exp;
    router.params.astar_fac = opts->astar_fac;
    router.params.max_criticality = opts->max_criticality;
    router.params.bend_cost = opts->bend_cost;

	//std::sort(begin(nets), end(nets), [] (const net_t &a, const net_t &b) -> bool {
			//return a.sinks.size() > b.sinks.size();
			//});

	router.pres_fac = opts->first_iter_pres_fac;

	for (const auto &n : get_vertices(router.g)) {
		const auto &prop = get_vertex_props(router.g, n);
		if (prop.type == CHANX) {
			extern t_rr_indexed_data *rr_indexed_data;
			router.delta = rr_indexed_data[prop.cost_index].base_cost;
			break;
		}
	}

#ifdef PTHREAD_BARRIER
	assert(pthread_barrier_init(&router.sync.barrier, nullptr, opts->num_threads) == 0);
#else
	router.sync.barrier = new boost::barrier(opts->num_threads);
	//router.sync.barrier = new SpinningBarrier(opts->num_threads);
#endif
	router.sync.stop_routing = 0;

	void *aligned_current_num_incoming_edges;
	assert(posix_memalign(&aligned_current_num_incoming_edges, 64, sizeof(aligned_atomic<int, 64>)*num_virtual_nets) == 0);
	assert(((unsigned long)aligned_current_num_incoming_edges & 63) == 0);
	router.sync.current_num_incoming_edges = new(aligned_current_num_incoming_edges) aligned_atomic<int, 64>[num_virtual_nets];
	//router.sync.current_num_incoming_edges = new std::atomic<int>[num_virtual_nets];

	router.stats.resize(opts->num_threads);
	router.route_time.resize(opts->num_threads);
	router.wait_time.resize(opts->num_threads);
	router.pure_route_time.resize(opts->num_threads);
	router.update_time.resize(opts->num_threads);
	router.last_barrier_wait_time.resize(opts->num_threads);
	router.pure_phase_time.resize(opts->num_threads);
	router.phase_time.resize(opts->num_threads);

	router.net_route_time.resize(num_virtual_nets);
	router.net_router.resize(num_virtual_nets);
	router.net_stats.resize(num_virtual_nets);

	router.iter = 0;

	vector<pthread_t> tids(opts->num_threads);
	for (int i = 1; i < opts->num_threads; ++i) {
		pthread_attr_t attr;
		assert(pthread_attr_init(&attr) == 0);

#ifdef __linux__
		cpu_set_t cpuset;
		CPU_ZERO(&cpuset);
		CPU_SET(4+i, &cpuset);
		assert(pthread_attr_setaffinity_np(&attr, sizeof(cpu_set_t), &cpuset) == 0);
#endif

		//assert(pthread_create(&tids[i], &attr, virtual_worker_thread, (void *)&router) == 0);
		assert(pthread_create(&tids[i], &attr, virtual_topo_worker_thread, &router) == 0);
	}

#ifdef PTHREAD_BARRIER
	pthread_barrier_wait(&router.sync.barrier);
#else
	router.sync.barrier->wait();
#endif

	timer::duration total_route_time = timer::duration::zero();

	/* just to make sure that we are allocated on the correct NUMA node */
	DeltaSteppingRouter *net_router = new DeltaSteppingRouter(router.g, router.state.congestion, router.pres_fac);

	sprintf(buffer, "%d", net_router->get_instance());
	zlog_put_mdc("tid", buffer);

	bool routed = false;
	bool has_resched = false;
	for (; router.iter < opts->max_router_iterations && !routed; ++router.iter) {
		sprintf(buffer, "%d", router.iter);
		zlog_put_mdc("iter", buffer);

		auto resched_start = timer::now();
		if (router.iter == 1 && has_resched) {
			//for (auto &p : router.nets) {
				//p.clear();
			//}
			//schedule_virtual_nets_3(overlap, all_virtual_nets_ptr, order, opts->num_threads, &router.net_route_time, router.nets);
			schedule_virtual_nets_6(overlap, all_virtual_nets_ptr, [&router] (int v) -> float { return router.net_stats[v].num_heap_pushes; }, opts->num_threads, router.directed, router.num_incoming_edges);

			//int num_directed_edges = 0;
			//for (int i = 0; i < router.directed.size(); ++i) {
				//num_directed_edges += router.directed[i].size();
			//}

			//assert(((num_edges-num_orig_directed_edges) % 2) == 0);
			//assert(num_directed_edges-num_orig_directed_edges == (num_edges-num_orig_directed_edges) / 2);
		}
		auto resched_time = timer::now()-resched_start;

		printf("Routing iteration: %d\n", router.iter);

		//assert(router.sync.num_incoming_edges.empty());

		//router.sync.handles.resize(router.topo_nets.size());

		//for (const auto &vnet : router.topo_nets) {
			//router.sync.handles[vnet->global_index] = router.sync.num_incoming_edges.emplace(make_pair(router.num_incoming_edges[vnet->global_index], vnet));
		//}
		
		for (int i = 0; i < num_virtual_nets; ++i) {
			router.sync.current_num_incoming_edges[i].val = router.num_incoming_edges[i];
		}

		for (const auto &vnet : router.topo_nets) {
			if (router.num_incoming_edges[vnet->global_index] == 0) {
				router.sync.global_pending_nets.push(vnet);
				//router.sync.fast_global_pending_nets.enqueue(vnet);
			}
		}

		assert(!router.sync.global_pending_nets.empty());
		//assert(router.sync.fast_global_pending_nets.size_approx() > 0);

		router.sync.net_index = 0;

		for (int i = 0; i < num_virtual_nets; ++i) {
			router.net_stats[i].num_heap_pops = std::numeric_limits<int>::max();
			router.net_stats[i].num_heap_pushes = std::numeric_limits<int>::max();
			router.net_stats[i].num_neighbor_visits = std::numeric_limits<int>::max();
		}

		for (int i = 0; i < opts->num_threads; ++i) {
			router.stats[i].num_heap_pops = 0;
			router.stats[i].num_heap_pushes = 0;
			router.stats[i].num_neighbor_visits = 0;
			router.route_time[i] = timer::duration::zero();
			router.wait_time[i] = timer::duration::zero();
			router.pure_route_time[i] = timer::duration::zero();
			router.update_time[i] = timer::duration::zero();
			router.last_barrier_wait_time[i] = timer::duration::zero();
		}

		for (auto &net : nets) {
			update_sink_criticalities(net, router.state.net_timing[net.vpr_id], router.params);
		}

		//do_virtual_work(&router, *net_router);
		do_virtual_work_topo_4(&router, *net_router);

		/* checking */
		for (int i = 0; i < num_virtual_nets; ++i) {
			assert(router.sync.current_num_incoming_edges[i].val == 0);
		}

		for (int i = 0; i < num_vertices(router.g); ++i) {
			router.state.congestion[i].recalc_occ = 0; 
		}

		for (const auto &net : nets) {
			check_route_tree(router.state.route_trees[net.local_id], net, router.g);
			recalculate_occ(router.state.route_trees[net.local_id], router.g, router.state.congestion);
		}

		printf("Resched time [has_resched %d]: %g\n", has_resched ? 1 : 0, duration_cast<nanoseconds>(resched_time).count()/1e9);

		printf("Route time: ");
		for (int i = 0; i < opts->num_threads; ++i) {
			printf("%g ", duration_cast<nanoseconds>(router.route_time[i]).count()/1e9);
		}
		printf("\n");

		//printf("Route time debug: ");
		//for (auto iter = begin(router.route_time); iter != end(router.route_time); ++iter) {
			//printf("%g ", duration_cast<nanoseconds>(*iter).count()/1e9);
		//}
		//printf("\n");

		auto max_route_time_i = std::max_element(begin(router.route_time), end(router.route_time));
		assert(max_route_time_i != end(router.route_time));
		timer::duration max_route_time = *max_route_time_i;
		total_route_time += max_route_time;

		printf("Max route time: %lu %g\n", max_route_time.count(), max_route_time.count() / 1e9);

		printf("Wait time: ");
		for (int i = 0; i < opts->num_threads; ++i) {
			printf("%g (%g) ", duration_cast<nanoseconds>(router.wait_time[i]).count()/1e9, 100.0*router.wait_time[i].count()/router.route_time[i].count());
		}
		printf("\n");

		printf("Pure route time: ");
		for (int i = 0; i < opts->num_threads; ++i) {
			printf("%g (%g) ", duration_cast<nanoseconds>(router.pure_route_time[i]).count()/1e9, 100.0*router.pure_route_time[i].count()/router.route_time[i].count());
		}
		printf("\n");

		printf("Update time: ");
		for (int i = 0; i < opts->num_threads; ++i) {
			printf("%g (%g) ", duration_cast<nanoseconds>(router.update_time[i]).count()/1e9, 100.0*router.update_time[i].count()/router.route_time[i].count());
		}
		printf("\n");

		printf("Last barrier wait time: ");
		for (int i = 0; i < opts->num_threads; ++i) {
			printf("%g (%g) ", duration_cast<nanoseconds>(router.last_barrier_wait_time[i]).count()/1e9, 100.0*router.last_barrier_wait_time[i].count()/router.route_time[i].count());
		}
		printf("\n");

		unsigned long total_num_heap_pops = 0;
		printf("Num heap pops: ");
		for (int i = 0; i < opts->num_threads; ++i) {
			printf("%d ", router.stats[i].num_heap_pops);
			total_num_heap_pops += router.stats[i].num_heap_pops;
		}
		printf("\n");
		printf("Total num heap pops: %lu\n", total_num_heap_pops);

		unsigned long total_num_neighbor_visits = 0;
		printf("Num neighbor visits: ");
		for (int i = 0; i < opts->num_threads; ++i) {
			printf("%d ", router.stats[i].num_neighbor_visits);
			total_num_neighbor_visits += router.stats[i].num_neighbor_visits;
		}
		printf("\n");
		printf("Total num neighbor visits: %lu\n", total_num_neighbor_visits);

		unsigned long total_num_heap_pushes = 0;
		printf("Num heap pushes: ");
		for (int i = 0; i < opts->num_threads; ++i) {
			printf("%d ", router.stats[i].num_heap_pushes);
			total_num_heap_pushes += router.stats[i].num_heap_pushes;
		}
		printf("\n");
		printf("Total num heap pushes: %lu\n", total_num_heap_pushes);

		sprintf(buffer, "%s/stats_%d.txt", dirname, router.iter);
		FILE *stats = fopen(buffer, "w");

		for (int i = 0; i < router.nets.size(); ++i) {
			for (int j = 0; j < opts->num_threads; ++j) {
				fprintf(stats, "%g %g %g ", duration_cast<nanoseconds>(router.pure_phase_time[j][i]).count()/1e9, duration_cast<nanoseconds>(router.phase_time[j][i]).count()/1e9,
						100.0*router.pure_phase_time[j][i].count()/router.phase_time[j][i].count());
			}
			fprintf(stats, "\n");
		}

		fclose(stats);

		sprintf(buffer, "%s/net_route_time_%d.txt", dirname, router.iter);
		FILE *nrt = fopen(buffer, "w");

		for (int i = 0; i < router.nets.size(); ++i) {
			for (int j = 0; j < router.nets[i].size(); ++j) {
				const auto *vnet = router.nets[i][j];
				fprintf(nrt, "%d %d %d %g %d %d %d %g\n", i, router.nets[i].size(), router.net_router[vnet->global_index], duration_cast<nanoseconds>(router.net_route_time[vnet->global_index]).count()/1e9, vnet->index, vnet->net->sinks.size(), vnet->sinks.size(), bg::area(vnet->bounding_box));
			}
		}

		fclose(nrt);

		unsigned long num_overused_nodes = 0;
		vector<int> overused_nodes_by_type(NUM_RR_TYPES, 0);

		if (feasible_routing(router.g, router.state.congestion)) {
			//dump_route(*current_traces_ptr, "route.txt");
			routed = true;

		} else {
			for (int i = 0; i < num_vertices(router.g); ++i) {
				if (router.state.congestion[i].occ > get_vertex_props(router.g, i).capacity) {
					++num_overused_nodes;
					const auto &v_p = get_vertex_props(router.g, i);
					++overused_nodes_by_type[v_p.type];
				}
			}

			//auto update_cost_start = timer::now();

			if (router.iter == 0) {
				router.pres_fac = opts->initial_pres_fac;
				update_costs(router.g, router.state.congestion, router.pres_fac, 0);
			} else {
				router.pres_fac *= opts->pres_fac_mult;

				/* Avoid overflow for high iteration counts, even if acc_cost is big */
				router.pres_fac = std::min(router.pres_fac, static_cast<float>(HUGE_POSITIVE_FLOAT / 1e5));

				update_costs(router.g, router.state.congestion, router.pres_fac, opts->acc_fac);
			}

			//update_cost_time = timer::now()-update_cost_start;
		}

		if (router.iter == 0) {
			sprintf(buffer, "%s/routes_iter_%d.txt", dirname, router.iter);
			write_routes(buffer, router.state.route_trees, nets.size());
		} else if (routed) {
			sprintf(buffer, "%s/routes_final.txt", dirname);
			write_routes(buffer, router.state.route_trees, nets.size());
		}

		if (routed) {
			sprintf(buffer, "%s/net_dijkstra_stats_final.txt", dirname);
			write_net_dijkstra_stats(buffer, router.net_stats);
			sprintf(buffer, "%s/congestion_state_final.txt", dirname);
			write_congestion_state(buffer, router.state.congestion, num_vertices(router.g));
		} else {
			sprintf(buffer, "%s/net_dijkstra_stats_%d.txt", dirname, router.iter);
			write_net_dijkstra_stats(buffer, router.net_stats);
			sprintf(buffer, "%s/congestion_state_%d.txt", dirname, router.iter);
			write_congestion_state(buffer, router.state.congestion, num_vertices(router.g));
		}

		//auto analyze_timing_start = timer::now();

		float crit_path_delay = analyze_timing(router.state.net_timing);

		//analyze_timing_time = timer::now()-analyze_timing_start;
		
		printf("Overused: %lu Crit path delay: %g\n", num_overused_nodes, crit_path_delay);
		printf("\n");

		//for (int i = 0; i < nets.size(); ++i) {
			//route_tree_clear(router.state.back_route_trees[i]);
		//}
		//std::swap(router.state.route_trees, router.state.back_route_trees);
		for (int i = 0; i < num_virtual_nets; ++i) {
			router.state.back_added_rr_nodes[i].clear();
		}
		std::swap(router.state.added_rr_nodes, router.state.back_added_rr_nodes);
		for (int i = 0; i < nets.size(); ++i) {
			route_tree_clear(router.state.route_trees[i]);
		}
	}

	if (routed) {
		printf("Routed in %d iterations\n", router.iter);
		printf("Total route time: %g\n", duration_cast<nanoseconds>(total_route_time).count() / 1e9);
	}

	return false;
}

bool partitioning_delta_stepping_deterministic_route(t_router_opts *opts)
{
	init_logging();
    zlog_set_record("custom_output", concurrent_log_impl);

	router_t<net_t *> router;

	assert(((unsigned long)&router.sync.stop_routing & 63) == 0);
	assert(((unsigned long)&router.sync.net_index & 63) == 0);
#ifdef __linux__
	assert(((unsigned long)&router.sync.barrier & 63) == 0);
#endif

	router.num_threads = opts->num_threads;

	init_graph(router.g);

	init_sprintf_rr_node(&router.g);

	vector<net_t> nets;
	vector<net_t> global_nets;
	init_nets(nets, global_nets, opts->bb_factor, opts->large_bb);

	//best_case(nets, router.nets);
	test_rtree(nets, router.nets);
	//test_misr(nets);
	write_nets(nets);

	congestion_t *congestion_aligned;
//#ifdef __linux__
	//assert(posix_memalign((void **)&congestion_aligned, 64, sizeof(congestion_t)*num_vertices(router.g)) == 0);
	//assert(((unsigned long)congestion_aligned & 63) == 0);
//#else
	congestion_aligned = (congestion_t *)malloc(sizeof(congestion_t)*num_vertices(router.g));
//#endif 
	router.state.congestion = new(congestion_aligned) congestion_t[num_vertices(router.g)];
    for (int i = 0; i < num_vertices(router.g); ++i) {
        router.state.congestion[i].acc_cost = 1;
        router.state.congestion[i].pres_cost = 1;
        router.state.congestion[i].occ = 0;
    }

	router.state.route_trees = new route_tree_t[nets.size()];
	for (int i = 0; i < nets.size(); ++i) {
		route_tree_init(router.state.route_trees[i]);
	}

    router.state.net_timing = new t_net_timing[nets.size()+global_nets.size()];
    init_net_timing(nets, global_nets, router.state.net_timing);

    router.params.criticality_exp = opts->criticality_exp;
    router.params.astar_fac = opts->astar_fac;
    router.params.max_criticality = opts->max_criticality;
    router.params.bend_cost = opts->bend_cost;

	//std::sort(begin(nets), end(nets), [] (const net_t &a, const net_t &b) -> bool {
			//return a.sinks.size() > b.sinks.size();
			//});

	router.pres_fac = opts->first_iter_pres_fac;

	for (const auto &n : get_vertices(router.g)) {
		const auto &prop = get_vertex_props(router.g, n);
		if (prop.type == CHANX) {
			extern t_rr_indexed_data *rr_indexed_data;
			router.delta = rr_indexed_data[prop.cost_index].base_cost;
			break;
		}
	}

#ifdef PTHREAD_BARRIER
	assert(pthread_barrier_init(&router.sync.barrier, nullptr, opts->num_threads) == 0);
#else
	router.sync.barrier = new boost::barrier(opts->num_threads);
	//router.sync.barrier = new SpinningBarrier(opts->num_threads);
#endif
	router.sync.stop_routing = 0;

	router.stats.resize(opts->num_threads);
	router.pure_route_time.resize(opts->num_threads);
	router.route_time.resize(opts->num_threads);

	router.iter = 0;

	vector<pthread_t> tids(opts->num_threads);
	for (int i = 1; i < opts->num_threads; ++i) {
		pthread_attr_t attr;
		assert(pthread_attr_init(&attr) == 0);

#ifdef __linux__
		cpu_set_t cpuset;
		CPU_ZERO(&cpuset);
		CPU_SET(i, &cpuset);
		assert(pthread_attr_setaffinity_np(&attr, sizeof(cpu_set_t), &cpuset) == 0);
#endif

		assert(pthread_create(&tids[i], &attr, worker_thread, (void *)&router) == 0);
	}

	timer::duration total_route_time = timer::duration::zero();

	/* just to make sure that we are allocated on the correct NUMA node */
	DeltaSteppingRouter *net_router = new DeltaSteppingRouter(router.g, router.state.congestion, router.pres_fac);

	char buffer[256];
	sprintf(buffer, "%d", net_router->get_instance());
	zlog_put_mdc("tid", buffer);

	bool routed = false;
	for (; router.iter < opts->max_router_iterations && !routed; ++router.iter) {
		printf("Routing iteration: %d\n", router.iter);

		do_work(&router, *net_router);

		auto max_duration = *std::max_element(begin(router.route_time), end(router.route_time));
		total_route_time += max_duration;

		printf("Route time: ");
		for (int i = 0; i < opts->num_threads; ++i) {
			printf("%g ", duration_cast<nanoseconds>(router.route_time[i]).count()/1e9);
		}
		printf("\n");

		printf("Max route time: %lu %g\n", max_duration.count(), max_duration.count() / 1e9);

		printf("Pure route time: ");
		for (int i = 0; i < opts->num_threads; ++i) {
			printf("%g (%g) ", duration_cast<nanoseconds>(router.pure_route_time[i]).count()/1e9, 100.0*router.pure_route_time[i].count()/router.route_time[i].count());
		}
		printf("\n");

		printf("Num heap pops: ");
		for (int i = 0; i < opts->num_threads; ++i) {
			printf("%d ", router.stats[i].num_heap_pops);
		}
		printf("\n");

		/* checking */
		for (int i = 0; i < num_vertices(router.g); ++i) {
			router.state.congestion[i].recalc_occ = 0; 
		}

		for (const auto &net : nets) {
			check_route_tree(router.state.route_trees[net.local_id], net, router.g);
			recalculate_occ(router.state.route_trees[net.local_id], router.g, router.state.congestion);
		}

		unsigned long num_overused_nodes = 0;
		vector<int> overused_nodes_by_type(NUM_RR_TYPES, 0);

		if (feasible_routing(router.g, router.state.congestion)) {
			//dump_route(*current_traces_ptr, "route.txt");
			routed = true;
		} else {
			for (int i = 0; i < num_vertices(router.g); ++i) {
				if (router.state.congestion[i].occ > get_vertex_props(router.g, i).capacity) {
					++num_overused_nodes;
					const auto &v_p = get_vertex_props(router.g, i);
					++overused_nodes_by_type[v_p.type];
				}
			}

			//auto update_cost_start = timer::now();

			if (router.iter == 0) {
				router.pres_fac = opts->initial_pres_fac;
				update_costs(router.g, router.state.congestion, router.pres_fac, 0);
			} else {
				router.pres_fac *= opts->pres_fac_mult;

				/* Avoid overflow for high iteration counts, even if acc_cost is big */
				router.pres_fac = std::min(router.pres_fac, static_cast<float>(HUGE_POSITIVE_FLOAT / 1e5));

				update_costs(router.g, router.state.congestion, router.pres_fac, opts->acc_fac);
			}

			//update_cost_time = timer::now()-update_cost_start;
		}

		//auto analyze_timing_start = timer::now();

		float crit_path_delay = analyze_timing(router.state.net_timing);

		//analyze_timing_time = timer::now()-analyze_timing_start;
		
		printf("Overused: %lu Crit path delay: %g\n", num_overused_nodes, crit_path_delay);
		printf("\n");
	}

	if (routed) {
		printf("Routed in %d iterations\n", router.iter);
		printf("Total route time: %g\n", duration_cast<nanoseconds>(total_route_time).count() / 1e9);
	}

	return false;
}
