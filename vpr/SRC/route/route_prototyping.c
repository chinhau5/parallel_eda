#include <stdio.h>
#include <math.h>
#include <time.h>
#include <assert.h>
#include <pthread.h>
#include "util.h"
#include "vpr_types.h"
#include "globals.h"
#include "route_export.h"
#include "route_common.h"
#include "route_tree_timing.h"
#include "route_timing.h"
#include "heapsort.h"
#include "path_delay.h"
#include "net_delay.h"
#include "stats.h"
#include "ReadOptions.h"

void get_number_of_heap_pushes()
{
	float pres_fac;
	float max_criticality;
	float criticality_exp;
	float astar_fac;
	float bend_cost;
	float *pin_criticality;
	int *sink_order;
	t_rt_node ** rt_node_of_sink;
	float *net_delay;
	t_slack * slacks;
	int num_heap_pushes;

	alloc_timing_driven_route_structs(&pin_criticality, &sink_order,
			&rt_node_of_sink);
	timing_driven_route_net(0, pres_fac, max_criticality,
			criticality_exp,  astar_fac,  bend_cost,
			pin_criticality, sink_order,
			rt_node_of_sink,  net_delay, slacks,
			&num_heap_pushes);
}
