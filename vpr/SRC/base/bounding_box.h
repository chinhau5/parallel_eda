#ifndef BOUNGING_BOX_H
#define BOUNGING_BOX_H

int get_bounding_box_area(const struct s_bb *bb);

bool bounding_box_overlap(const struct s_bb &bb_a, const struct s_bb &bb_b);

bool bounding_box_overlap(int net_a, int net_b);

int get_bounding_box_overlap_area(const struct s_bb &bb_a, const struct s_bb &bb_b);

#endif
