#ifndef QUADTREE_H
#define QUADTREE_H

#include "geometry.h"
#include <zlog.h>

template<typename Tag>
class QuadTree {
	private:
		int quadrant;
		bounding_box_t box;
		QuadTree *parent;
		QuadTree *children;
		std::vector<std::pair<bounding_box_t, Tag>> items;
		int max_num_items;

	public:
		QuadTree(const bounding_box_t &box, int max_num_items = 16)
			: quadrant(-1), box(box), parent(nullptr), children(nullptr), max_num_items(max_num_items)
		{
		}

	protected:
		QuadTree()
			: quadrant(-1), parent(nullptr), children(nullptr)
		{
		}

	public:
		int num_levels() const
		{
			int num_levels_;
			if (children) {
				num_levels_ = std::numeric_limits<int>::min();
				for (int i = 0; i < 4; ++i) {
					num_levels_ = std::max(num_levels_, 1+children[i].num_levels());
				}
			} else {
				num_levels_ = 1;
			}
			return num_levels_;
		}

		template<typename Func>
		void print_items(const Func &func) const
		{
			for (const auto &item : items) {
				func(item);
			}
			if (children) {
				for (int i = 0; i < 4; ++i) {
					children[i].print_items(func);
				}
			}
		}
		bool empty() const
		{
			bool e = items.empty();
			if (children) {
				for (int i = 0; i < 4; ++i) {
					e = e && children[i].empty();
				}
			}
			return e; 
		}

		void verify(int &num_items) const
		{
			for (const auto &item : items) {
				assert(contains(item.first));
			}
			num_items += items.size();
			if (children) {
				for (int i = 0; i < 4; ++i) {
					children[i].verify(num_items);
				}
			}
		}

		void print_num_items(int level)
		{
			extern zlog_category_t *delta_log;
			for (int i = 0; i < level; ++i) {
				zlog_info(delta_log, "\t");
			}
			zlog_info(delta_log, "num_items level %d quadrant %d: %lu\n", level, quadrant, items.size());
			if (children) {
				for (int i = 0; i < 4; ++i) {
					children[i].print_num_items(level+1);
				}
			}
		}

		void get_non_overlapping_quadrant(QuadTree *current, std::vector<QuadTree *> &non_overlapping_quadrants) {
			if (current->parent) {
				for (int i = 0; i < 4; ++i) {
					if (i != current->quadrant) {
						non_overlapping_quadrants.push_back(current->parent->children[i]);
					}
				}
			}
		}

		QuadTree *insert(const std::pair<bounding_box_t, Tag> &item)
		{
			QuadTree *res = nullptr;

			int q = get_quadrant(item.first);
			if (q >= 0) {
				res = children[q].insert(item);
			} else {
				assert(contains(item.first));
				items.push_back(item);
				res = this;

				if (items.size() > max_num_items && split()) {
					for (auto iter = begin(items); iter != end(items); ) {
						int q = get_quadrant(iter->first);
						if (q >= 0) {
							res = children[q].insert(*iter);
							iter = items.erase(iter);
						} else {
							assert(contains(iter->first));
							++iter;
						}
					}
				}
			}
			assert(res);
			return res;
		}

		void remove(const std::pair<bounding_box_t, Tag> &item)
		{
			int q = get_quadrant(item.first);
			if (q >= 0) {
				children[q].remove(item);
			} else {
				auto iter = find(begin(items), end(items), item);
				assert(iter != end(items));
				items.erase(iter);
			}
		}

		bool has_overlapping_boxes(const std::pair<bounding_box_t, Tag> &item) const
		{
			int q = get_quadrant(item.first);

			bool res; 
			if (q >= 0) {
				//assert(!contains(item));
				bool has_overlap = any_of(begin(items), end(items),
						[&item] (const std::pair<bounding_box_t, Tag> &other_item) -> bool {
						return box_overlap(other_item.first, item.first) && item.second != other_item.second;
						});

				res = has_overlap || children[q].has_overlapping_boxes(item);
			} else {
				assert(contains(item.first));

				bool has_overlap = any_of(begin(items), end(items),
						[&item] (const std::pair<bounding_box_t, Tag> &other_item) -> bool {
						return box_overlap(other_item.first, item.first) && item.second != other_item.second;
						});

				if (children && !has_overlap) {
					for (int i = 0; i < 4 && !has_overlap; ++i) {
						if (box_overlap(item.first, children[i].box)) {
							bounding_box_t overlap_box;
							overlap_box.xmin = std::max(item.first.xmin, children[i].box.xmin);
							overlap_box.xmax = std::min(item.first.xmax, children[i].box.xmax);
							overlap_box.ymin = std::max(item.first.ymin, children[i].box.ymin);
							overlap_box.ymax = std::min(item.first.ymax, children[i].box.ymax);
							has_overlap = has_overlap || children[i].has_overlapping_boxes({ overlap_box, item.second });
						}
					}
				}
				res = has_overlap;
			}
			return res;
		}

		void get_overlapping_boxes(const bounding_box_t &item, std::vector<const std::pair<bounding_box_t, Tag> *> &overlapping_boxes) const
		{
			int q = get_quadrant(item);

			struct ptr_back_inserter {
				std::vector<const std::pair<bounding_box_t, Tag> *> &res;
				ptr_back_inserter(std::vector<const std::pair<bounding_box_t, Tag> *> &res)
					: res(res) 
				{
				}

				ptr_back_inserter &operator++()
				{
					return *this;
				}

				ptr_back_inserter &operator*()
				{
					return *this;
				}

				ptr_back_inserter &operator=(const std::pair<bounding_box_t, Tag> &item)
				{
					res.push_back(&item);
					return *this;
				}
			};
			
			if (q >= 0) {
				//assert(!contains(item));
				copy_if(begin(items), end(items), ptr_back_inserter(overlapping_boxes),
						[&item] (const std::pair<bounding_box_t, Tag> &other_box) -> bool {
						return box_overlap(other_box.first, item);
						});

				children[q].get_overlapping_boxes(item, overlapping_boxes);
			} else {
				assert(contains(item));

				copy_if(begin(items), end(items), ptr_back_inserter(overlapping_boxes),
						[&item] (const std::pair<bounding_box_t, Tag> &other_box) -> bool {
						return box_overlap(other_box.first, item);
						});

				if (children) {
					for (int i = 0; i < 4; ++i) {
						if (box_overlap(item, children[i].box)) {
							bounding_box_t overlap_box;
							overlap_box.xmin = std::max(item.xmin, children[i].box.xmin);
							overlap_box.xmax = std::min(item.xmax, children[i].box.xmax);
							overlap_box.ymin = std::max(item.ymin, children[i].box.ymin);
							overlap_box.ymax = std::min(item.ymax, children[i].box.ymax);
							children[i].get_overlapping_boxes(overlap_box, overlapping_boxes);
						}
					}
				}
			}
		}

	private:
		bool contains(const bounding_box_t &item) const
		{
			return item.xmin >= box.xmin && item.xmax <= box.xmax
				&& item.ymin >= box.ymin && item.ymax <= box.ymax;
		}

		int get_quadrant(const bounding_box_t &item) const 
		{
			int q = -1;
			if (children) {
				int num_matches = 0;
				for (int i = 0; i < 4; ++i) {
					if (children[i].contains(item)) {
						q = i;
						++num_matches;
					}
				}
				assert(num_matches <= 1);
			} 
			return q;
		}

		bool split()
		{
			bool did_split = false;
			if (!children) {
				children = new QuadTree[4];

				for (int i = 0; i < 4; ++i) {
					children[i].max_num_items = max_num_items;
				}

				int xsplit = (box.xmin + box.xmax) / 2;
				int ysplit = (box.ymin + box.ymax) / 2;

				/* top right */
				children[0].box.xmin = xsplit+1; children[0].box.xmax = box.xmax;
				children[0].box.ymin = ysplit+1; children[0].box.ymax = box.ymax;

				/* bottom right */
				children[1].box.xmin = xsplit+1; children[1].box.xmax = box.xmax;
				children[1].box.ymin = box.ymin; children[1].box.ymax = ysplit;

				/* bottom left */
				children[2].box.xmin = box.xmin; children[2].box.xmax = xsplit;
				children[2].box.ymin = box.ymin; children[2].box.ymax = ysplit;

				/* top left */
				children[3].box.xmin = box.xmin; children[3].box.xmax = xsplit;
				children[3].box.ymin = ysplit+1; children[3].box.ymax = box.ymax;

				for (int i = 0; i < 4; ++i) {
					children[i].quadrant = i;
					children[i].parent = this;
				}

				for (int i = 0; i < 4; ++i) {
					assert(children[i].box.xmin < children[i].box.xmax);
					assert(children[i].box.ymin < children[i].box.ymax);
				}

				did_split = true;
			}
			return did_split;
		}
};

#endif
