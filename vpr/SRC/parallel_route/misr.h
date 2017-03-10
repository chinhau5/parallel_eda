#ifndef MISR_H
#define MISR_H

//using namespace std;

//#define DEBUG_MISR

#ifdef DEBUG_MISR
#define PRINT(msg, ...) printf((msg), __VA_ARGS__)
#else
#define PRINT(msg, ...) 
#endif

template<typename Item, typename ItemToInterval>
void max_independent_intervals_from_sorted(const std::vector<Item *> &items, const ItemToInterval &to_interval, std::vector<Item *> &chosen)
{
	if (items.size() == 0) {
		return;
	}

	int current_finishing = to_interval(*items[0]).second;
	chosen.push_back(items[0]);
	for (int i = 1; i < items.size(); ++i) {
		const auto &interval = to_interval(*items[i]);
		if (interval.first > current_finishing) {
			chosen.push_back(items[i]);
			current_finishing = interval.second;
		}
	}
}

template<typename Item, typename ItemToBox>
void split_ver_edges(const std::vector<std::pair<int, Item *>> &ver_edges, const ItemToBox &to_box,
		float x_median,
		std::vector<std::pair<int, Item *>> &left, std::vector<std::pair<int, Item *>> &right)
{
	for (const auto &ver_edge : ver_edges) {
		const auto &box = to_box(*ver_edge.second);

		int l = bg::get<bg::min_corner, 0>(box);
		int r = bg::get<bg::max_corner, 0>(box);

		assert((ver_edge.first == l && ver_edge.first != r) || (ver_edge.first != l && ver_edge.first == r));

		if (r < x_median) {
			left.push_back(ver_edge);
		} else if (l > x_median) {
			right.push_back(ver_edge);
		} else {
			assert((l < x_median && r >= x_median) || (l <= x_median && r > x_median));
		}
	}
}

template<typename Item, typename ItemToBox>
void split_hor(const std::vector<Item *> &sorted_items, const ItemToBox &to_box,
		float x_median,
		std::vector<Item *> &left, std::vector<Item *> &right, std::vector<Item *> &middle)
{
	for (const auto &item : sorted_items) {
		const auto &box = to_box(*item);

		int l = bg::get<bg::min_corner, 0>(box);
		int r = bg::get<bg::max_corner, 0>(box);

		if (r < x_median) {
			left.push_back(item);
		} else if (l > x_median) {
			right.push_back(item);
		} else {
			assert((l < x_median && r >= x_median) || (l <= x_median && r > x_median));
			middle.push_back(item);
		}
	}
}

template<typename Item, typename ItemToBox, typename ItemToInterval>
void max_independent_rectangles_internal(const std::vector<std::pair<int, Item *>> &ver_edges, const std::vector<Item *> &sorted_items, const ItemToBox &to_box, const ItemToInterval &to_interval, std::vector<Item *> &chosen)
{
	assert((ver_edges.size() % 2) == 0);
	assert(sorted_items.size() == ver_edges.size() / 2);
	/* base case */
	if (sorted_items.size() <= 2) {
		if (sorted_items.size() == 2) {
			const auto &a_box = to_box(*sorted_items[0]);
			const auto &b_box = to_box(*sorted_items[1]);

			if (bg::disjoint(a_box, b_box)) {
				chosen.push_back(sorted_items[0]);
				chosen.push_back(sorted_items[1]);
			} else {
				if (bg::area(a_box) < bg::area(b_box)) {
					chosen.push_back(sorted_items[0]);
				} else {
					chosen.push_back(sorted_items[1]);
				}
			}
		} else {
			for (int i = 0; i < sorted_items.size(); ++i) {
				chosen.push_back(sorted_items[i]);
			}
		}

		return;
	}

	for (const auto &ver_e : ver_edges) {
		PRINT("ver edge %d box %X\n", ver_e.first, ver_e.second);
	}

	float x_median = (float)(ver_edges[ver_edges.size()/2 - 1].first + ver_edges[ver_edges.size()/2].first) / 2;

	PRINT("x median %g\n", x_median);

//#ifdef DEBUG_MISR
	int prev_ver_edge = -1;
	for (const auto &e : ver_edges) {
		assert(e.first >= prev_ver_edge);
		prev_ver_edge = e.first;
	}

	int previous_finishing = -1;
	for (const auto &item : sorted_items) {
		const auto &interval = to_interval(*item);
		assert(interval.second >= previous_finishing);
		previous_finishing = interval.second;
	}
//#endif

	std::vector<std::pair<int, Item *>> left_ver;
	std::vector<std::pair<int, Item *>> right_ver;
	split_ver_edges(ver_edges, to_box,
			x_median,
			left_ver, right_ver);

	std::vector<Item *> left_sorted;
	std::vector<Item *> right_sorted;
	std::vector<Item *> middle_sorted;
	split_hor(sorted_items, to_box,
			x_median,
			left_sorted, right_sorted, middle_sorted);

	std::vector<Item *> left_right_chosen;
	std::vector<Item *> middle_chosen;
	max_independent_rectangles_internal(left_ver, left_sorted, to_box, to_interval, left_right_chosen);
	max_independent_rectangles_internal(right_ver, right_sorted, to_box, to_interval, left_right_chosen);
	max_independent_intervals_from_sorted(middle_sorted, to_interval, middle_chosen);

	if (left_right_chosen.size() > middle_chosen.size()) {
		//chosen = std::move(left_right_chosen);
		chosen.insert(end(chosen), begin(left_right_chosen), end(left_right_chosen));
	} else {
		//chosen = std::move(middle_chosen);
		chosen.insert(end(chosen), begin(middle_chosen), end(middle_chosen));
	}
}

template<typename Item, typename ItemToBox, typename ItemToInterval>
void max_independent_rectangles(std::vector<Item> &items, const ItemToBox &to_box, const ItemToInterval &to_interval, std::vector<Item *> &chosen)
{
	std::vector<std::pair<int, Item *>> ver_edges;
	std::vector<Item *> sorted_items;

	for (int i = 0; i < items.size(); ++i) {
		const auto &box = to_box(items[i]);

		ver_edges.push_back(std::make_pair(bg::get<bg::min_corner, 0>(box), &items[i]));
		ver_edges.push_back(std::make_pair(bg::get<bg::max_corner, 0>(box), &items[i]));

		sorted_items.push_back(&items[i]);
	}

	std::sort(begin(ver_edges), end(ver_edges));
	std::sort(begin(sorted_items), end(sorted_items), [&] (const Item *a, const Item *b) -> bool { return to_interval(*a).second < to_interval(*b).second; });

	max_independent_rectangles_internal(ver_edges, sorted_items, to_box, to_interval, chosen);
}

template<typename Item, typename ItemToBox>
void verify_ind(const std::vector<Item> &chosen, const ItemToBox &to_box)
{
	for (int i = 0; i < chosen.size(); ++i) {
		for (int j = i + 1; j < chosen.size(); ++j) {
			const auto &box_a = to_box(chosen[i]);
			const auto &box_b = to_box(chosen[j]);
			assert(bg::disjoint(box_a, box_b));
		}	
	}
}

#endif
