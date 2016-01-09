#ifndef SCHEDULER_H
#define SCHEDULER_H

template<typename Interval, typename Object>
void max_independent_intervals(vector<pair<Interval, Object>> &intervals, vector<pair<Interval, Object> *> &chosen)
{
	if (intervals.size() == 0) {
		return;
	}

	using namespace boost::numeric;
	std::sort(intervals.begin(), intervals.end(), [] (const pair<Interval, Object> &a, const pair<Interval, Object> &b) -> bool {
			return a.first.lower() < b.first.lower();
			});

	int current_upper = intervals[0].first.upper();
	chosen.push_back(&intervals[0]);
	for (int i = 1; i < intervals.size(); ++i) {
		if (intervals[i].first.lower() > current_upper) {
			chosen.push_back(&intervals[i]);
			current_upper = intervals[i].first.upper();
		}
	}
}

template<typename Box, typename Object>
void max_independent_rectangles_simple(vector<pair<const Box *, Object> *> &bounding_boxes, vector<pair<const Box *, Object> *> &chosen)
{
	using namespace boost::numeric;

	typedef struct interval_wrapper_t {
		const Box *box;

		int lower() const { return box->ymin; }
		int upper() const { return box->ymax; }
	} interval_wrapper_t;

	using IntervalWrapperPair = pair<interval_wrapper_t, pair<const Box *, Object> *>;

	vector<IntervalWrapperPair> intervals(bounding_boxes.size());
	for (int i = 0; i < bounding_boxes.size(); ++i) {
		intervals[i].first.box = bounding_boxes[i]->first;
		intervals[i].second = bounding_boxes[i];
	}

	vector<IntervalWrapperPair *> chosen_intervals;
	max_independent_intervals(intervals, chosen_intervals);

	for (const auto &c : chosen_intervals) {
		chosen.push_back(c->second);
	}
}

template<typename Box, typename Object>
void max_independent_rectangles(vector<pair<const Box *, Object> *> &bounding_boxes, vector<pair<const Box *, Object> *> &chosen)
{
	extern zlog_category_t *independent_log;
	/* base case */
	if (bounding_boxes.size() <= 2) {
		if (bounding_boxes.size() == 2) {
			if (!box_overlap(*bounding_boxes[0]->first, *bounding_boxes[1]->first)) {
				chosen.push_back(bounding_boxes[0]);
				chosen.push_back(bounding_boxes[1]);
			} else {
				if (get_bounding_box_area(*bounding_boxes[0]->first) < get_bounding_box_area(*bounding_boxes[1]->first)) {
					chosen.push_back(bounding_boxes[0]);
				} else {
					chosen.push_back(bounding_boxes[1]);
				}
			}
		} else {
			for (int i = 0; i < bounding_boxes.size(); ++i) {
				chosen.push_back(bounding_boxes[i]);
			}
		}

		return;
	}

	using namespace boost::numeric;
	typedef struct vertical_edge_t {
		int x;
		int box_id;
		bool operator<(const vertical_edge_t &other) const {
			return x < other.x;
		}
	} vertical_edge_t;

	vector<vertical_edge_t> ver_edges;

	for (int i = 0; i < bounding_boxes.size(); ++i) {
		ver_edges.push_back({ bounding_boxes[i]->first->xmin, i });
		ver_edges.push_back({ bounding_boxes[i]->first->xmax, i });
	}

	std::sort(ver_edges.begin(), ver_edges.end());

	float x_median;
	bool odd = ver_edges.size() % 2 ? true : false;
	if (odd) {
		x_median = ver_edges[ver_edges.size()/2].x;
	} else {
		int left = ver_edges[ver_edges.size()/2 - 1].x;
		int right = ver_edges[ver_edges.size()/2].x;
		x_median = (float)(left + right) / 2;
	}

	vector<pair<const Box *, Object> *> left;
	vector<pair<const Box *, Object> *> right;
	vector<pair<const Box *, Object> *> middle;
	for (int i = 0; i < bounding_boxes.size(); ++i) {
		if (bounding_boxes[i]->first->xmin < x_median && bounding_boxes[i]->first->xmax < x_median) {
			left.push_back(bounding_boxes[i]);
			zlog_debug(independent_log, "Box x %d-%d y %d-%d is left of median %g\n", bounding_boxes[i]->first->xmin, bounding_boxes[i]->first->xmax, bounding_boxes[i]->first->ymin, bounding_boxes[i]->first->ymax, x_median);
		} else if (bounding_boxes[i]->first->xmin > x_median && bounding_boxes[i]->first->xmax > x_median) {
			right.push_back(bounding_boxes[i]);
			zlog_debug(independent_log, "Box x %d-%d y %d-%d is right of median %g\n", bounding_boxes[i]->first->xmin, bounding_boxes[i]->first->xmax, bounding_boxes[i]->first->ymin, bounding_boxes[i]->first->ymax, x_median);
		} else {
			zlog_debug(independent_log, "Box x %d-%d y %d-%d is middle of median %g\n", bounding_boxes[i]->first->xmin, bounding_boxes[i]->first->xmax, bounding_boxes[i]->first->ymin, bounding_boxes[i]->first->ymax, x_median);
			assert(bounding_boxes[i]->first->xmin <= x_median && bounding_boxes[i]->first->xmax >= x_median);
			middle.push_back(bounding_boxes[i]);
		}
	}

	vector<pair<const Box *, Object> *> left_chosen;
	vector<pair<const Box *, Object> *> right_chosen;
	vector<pair<const Box *, Object> *> middle_chosen;
	max_independent_rectangles(left, left_chosen);
	max_independent_rectangles(right, right_chosen);
	max_independent_rectangles_simple(middle, middle_chosen);

	if (left_chosen.size() + right_chosen.size() > middle_chosen.size()) {
		for (const auto &c : left_chosen) {
			chosen.push_back(c);
		}
		for (const auto &c : right_chosen) {
			chosen.push_back(c);
		}
	} else {
		for (const auto &c : middle_chosen) {
			chosen.push_back(c);
		}
	}
}

template<typename Box, typename Object>
void max_independent_rectangles_2(vector<pair<const Box *, Object> *> &bounding_boxes, vector<pair<const Box *, Object> *> &chosen)
{
	extern zlog_category_t *independent_log;
	/* base case */
	if (bounding_boxes.size() <= 2) {
		if (bounding_boxes.size() == 2) {
			if (!box_overlap(*bounding_boxes[0]->first, *bounding_boxes[1]->first)) {
				chosen.push_back(bounding_boxes[0]);
				chosen.push_back(bounding_boxes[1]);
			} else {
				if (get_bounding_box_area(*bounding_boxes[0]->first) < get_bounding_box_area(*bounding_boxes[1]->first)) {
					chosen.push_back(bounding_boxes[0]);
				} else {
					chosen.push_back(bounding_boxes[1]);
				}
			}
		} else {
			for (int i = 0; i < bounding_boxes.size(); ++i) {
				chosen.push_back(bounding_boxes[i]);
			}
		}

		return;
	}

	using namespace boost::numeric;
	typedef struct vertical_edge_t {
		int x;
		int box_id;
		bool operator<(const vertical_edge_t &other) const {
			return x < other.x;
		}
	} vertical_edge_t;

	vector<vertical_edge_t> ver_edges;

	for (int i = 0; i < bounding_boxes.size(); ++i) {
		ver_edges.push_back({ bounding_boxes[i]->first->xmin, i });
		ver_edges.push_back({ bounding_boxes[i]->first->xmax, i });
	}

	std::sort(ver_edges.begin(), ver_edges.end());

	float x_median;
	bool odd = ver_edges.size() % 2 ? true : false;
	if (odd) {
		x_median = ver_edges[ver_edges.size()/2].x;
	} else {
		int left = ver_edges[ver_edges.size()/2 - 1].x;
		int right = ver_edges[ver_edges.size()/2].x;
		x_median = (float)(left + right) / 2;
	}

	vector<pair<const Box *, Object> *> left;
	vector<pair<const Box *, Object> *> right;
	vector<pair<const Box *, Object> *> middle;
	for (int i = 0; i < bounding_boxes.size(); ++i) {
		if (bounding_boxes[i]->first->xmin < x_median && bounding_boxes[i]->first->xmax < x_median) {
			left.push_back(bounding_boxes[i]);
			zlog_debug(independent_log, "Box x %d-%d y %d-%d is left of median %g\n", bounding_boxes[i]->first->xmin, bounding_boxes[i]->first->xmax, bounding_boxes[i]->first->ymin, bounding_boxes[i]->first->ymax, x_median);
		} else if (bounding_boxes[i]->first->xmin > x_median && bounding_boxes[i]->first->xmax > x_median) {
			right.push_back(bounding_boxes[i]);
			zlog_debug(independent_log, "Box x %d-%d y %d-%d is right of median %g\n", bounding_boxes[i]->first->xmin, bounding_boxes[i]->first->xmax, bounding_boxes[i]->first->ymin, bounding_boxes[i]->first->ymax, x_median);
		} else {
			zlog_debug(independent_log, "Box x %d-%d y %d-%d is middle of median %g\n", bounding_boxes[i]->first->xmin, bounding_boxes[i]->first->xmax, bounding_boxes[i]->first->ymin, bounding_boxes[i]->first->ymax, x_median);
			assert(bounding_boxes[i]->first->xmin <= x_median && bounding_boxes[i]->first->xmax >= x_median);
			middle.push_back(bounding_boxes[i]);
		}
	}

	vector<pair<const Box *, Object> *> left_chosen;
	vector<pair<const Box *, Object> *> right_chosen;
	vector<pair<const Box *, Object> *> middle_chosen;
	max_independent_rectangles(left, left_chosen);
	max_independent_rectangles(right, right_chosen);
	max_independent_rectangles_simple(middle, middle_chosen);

	if (left_chosen.size() + right_chosen.size() > middle_chosen.size()) {
		for (const auto &c : left_chosen) {
			chosen.push_back(c);
		}
		for (const auto &c : right_chosen) {
			chosen.push_back(c);
		}
	} else {
		for (const auto &c : middle_chosen) {
			chosen.push_back(c);
		}
	}
}

void schedule_nets_ind(vector<net_t *> &nets, vector<vector<const net_t *>> &net_scheduled_at_time);

void verify_ind(const vector<const net_t *> &nets);

void schedule_nets_bounding_box(vector<net_t *> &nets, vector<pair<sink_t *, net_t *>> &scheduled_sinks);

#endif
