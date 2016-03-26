#include "log.h"
#include "utility.h"
#include "vpr_types.h"
#include "route.h"
#include "geometry.h"
#include "log.h"

using namespace std;

void delete_graph(RRGraph &g)
{
	//for (int i = 0; i < num_vertices(g); ++i) {
		//auto &v = get_vertex(g, i);
	//}
}

int get_track_domain(int track_number, int num_partitions)
{
	return track_number % num_partitions;
}

void init_graph(RRGraph &g, vector<vector<RRNode>> &sink_in_nodes, vector<vector<RRNode>> &ipin_in_nodes)
{
	extern s_rr_node *rr_node;
	extern int num_rr_nodes;
	extern t_rr_indexed_data *rr_indexed_data;
	extern struct s_switch_inf *switch_inf;

	sink_in_nodes.resize(num_rr_nodes);
	ipin_in_nodes.resize(num_rr_nodes);

	add_vertex(g, num_rr_nodes);
	for (int i = 0; i < num_vertices(g); ++i) {
		auto &v = get_vertex_props(g, i);
		v.type = rr_node[i].type;
		v.inc_direction = rr_node[i].direction == INC_DIRECTION;
		v.xlow = rr_node[i].xlow;
		v.ylow = rr_node[i].ylow;
		v.xhigh = rr_node[i].xhigh;
		v.yhigh = rr_node[i].yhigh;

		extern struct s_grid_tile **grid;
		auto type = &grid[v.xlow][v.ylow];

		//v.real_xlow = rr_node[i].xlow;
		//v.real_ylow = rr_node[i].ylow;
		//v.real_xhigh = rr_node[i].xhigh;
		//v.real_yhigh = rr_node[i].ylow + type->offset;
		v.R = rr_node[i].R;
		v.C = rr_node[i].C;
		v.cost_index = rr_node[i].cost_index;
		v.capacity = rr_node[i].capacity;

		char buffer[256];
		sprintf_rr_node(i, buffer);
		//zlog_debug(rr_log, "%s: real_xlow: %d real_xhigh: %d real_ylow: %d real_yhigh: %d\n", buffer, v.real_xlow, v.real_xhigh, v.real_ylow, v.real_yhigh);

		for (int j = 0; j < rr_node[i].num_edges; ++j) {
			int neighbor = rr_node[i].edges[j];
			auto &e = add_edge(g, i, neighbor);

			if (rr_node[neighbor].type == SINK) {
				assert(find(begin(sink_in_nodes[neighbor]), end(sink_in_nodes[neighbor]), i) == end(sink_in_nodes[neighbor]));
				sink_in_nodes[neighbor].push_back(i);
			} else if (rr_node[neighbor].type == IPIN) {
				assert(find(begin(ipin_in_nodes[neighbor]), end(ipin_in_nodes[neighbor]), i) == end(ipin_in_nodes[neighbor]));
				ipin_in_nodes[neighbor].push_back(i);
			}

			int si = rr_node[i].switches[j];

			auto &e_p = get_edge_props(g, e);
			
			e_p.buffered = switch_inf[si].buffered; 
			e_p.switch_delay = switch_inf[si].Tdel; 
			e_p.R = switch_inf[si].R; 
		}
	}
	zlog_info(delta_log, "RR graph num vertices: %d\n", num_vertices(g));
	zlog_info(delta_log, "RR graph num edges: %d\n", num_edges(g));
}

void init_graph(RRGraph &g)
{
	extern s_rr_node *rr_node;
	extern int num_rr_nodes;
	extern t_rr_indexed_data *rr_indexed_data;
	extern struct s_switch_inf *switch_inf;

	add_vertex(g, num_rr_nodes);
	for (int i = 0; i < num_vertices(g); ++i) {
		auto &v = get_vertex_props(g, i);
		v.type = rr_node[i].type;
		v.inc_direction = rr_node[i].direction == INC_DIRECTION;
		v.xlow = rr_node[i].xlow;
		v.ylow = rr_node[i].ylow;
		v.xhigh = rr_node[i].xhigh;
		v.yhigh = rr_node[i].yhigh;

		extern struct s_grid_tile **grid;
		auto type = &grid[v.xlow][v.ylow];

		//v.real_xlow = rr_node[i].xlow;
		//v.real_ylow = rr_node[i].ylow;
		//v.real_xhigh = rr_node[i].xhigh;
		//v.real_yhigh = rr_node[i].ylow + type->offset;
		v.R = rr_node[i].R;
		v.C = rr_node[i].C;
		v.cost_index = rr_node[i].cost_index;
		v.capacity = rr_node[i].capacity;

		char buffer[256];
		sprintf_rr_node(i, buffer);
		//zlog_debug(rr_log, "%s: real_xlow: %d real_xhigh: %d real_ylow: %d real_yhigh: %d\n", buffer, v.real_xlow, v.real_xhigh, v.real_ylow, v.real_yhigh);

		for (int j = 0; j < rr_node[i].num_edges; ++j) {
			auto &e = add_edge(g, i, rr_node[i].edges[j]);

			int si = rr_node[i].switches[j];

			auto &e_p = get_edge_props(g, e);
			
			e_p.buffered = switch_inf[si].buffered; 
			e_p.switch_delay = switch_inf[si].Tdel; 
			e_p.R = switch_inf[si].R; 
		}
	}
	zlog_info(delta_log, "RR graph num vertices: %d\n", num_vertices(g));
	zlog_info(delta_log, "RR graph num edges: %d\n", num_edges(g));
}

void delete_net_timing(const vector<net_t> &nets, const vector<net_t> &global_nets, t_net_timing *net_timing)
{
	for (const auto &net : nets) {
	   	delete [] net_timing[net.vpr_id].delay;
	   	delete [] net_timing[net.vpr_id].timing_criticality;
	   	delete [] net_timing[net.vpr_id].slack;
	}
	for (const auto &net : global_nets) {
	   	delete [] net_timing[net.vpr_id].delay;
	   	delete [] net_timing[net.vpr_id].timing_criticality;
	   	delete [] net_timing[net.vpr_id].slack;
	}
}

void init_net_timing(const vector<net_t> &nets, const vector<net_t> &global_nets, t_net_timing *net_timing)
{
	for (const auto &net : nets) {
	   	net_timing[net.vpr_id].delay = new float[net.sinks.size() + 1];
	   	net_timing[net.vpr_id].timing_criticality = new float[net.sinks.size() + 1];
	   	net_timing[net.vpr_id].slack = new float[net.sinks.size() + 1];
	}
	for (const auto &net : global_nets) {
	   	net_timing[net.vpr_id].delay = new float[net.sinks.size() + 1];
	   	net_timing[net.vpr_id].timing_criticality = new float[net.sinks.size() + 1];
	   	net_timing[net.vpr_id].slack = new float[net.sinks.size() + 1];
	}

	float init_timing_criticality_val = 1;

	for (const auto &net : nets) {
		for (int ipin = 1; ipin <= net.sinks.size(); ipin++) {
			net_timing[net.vpr_id].timing_criticality[ipin] = init_timing_criticality_val;
#ifdef PATH_COUNTING
			net_timing[net.vpr_id].path_criticality[ipin] = init_timing_criticality_val;
#endif		
		}
	} 
	
	for (const auto &net : global_nets) {
		/* Set delay of global signals to zero. Non-global net 
		 * 			delays are set by update_net_delays_from_route_tree() 
		 * 						inside timing_driven_route_net(), which is only called
		 * 									for non-global nets. */
		for (int ipin = 1; ipin <= net.sinks.size(); ipin++) {
			net_timing[net.vpr_id].delay[ipin] = 0.;
		}
	}
}

void init_nets(vector<net_t> &nets, vector<net_t> &global_nets, int bb_factor, const vector<vector<RRNode>> &sink_in_nodes)
{
	extern struct s_net *clb_net;
	extern int num_nets;
	extern int **net_rr_terminals; /* [0..num_nets-1][0..num_pins-1] */
	extern struct s_rr_node *rr_node;
	extern struct s_bb *route_bb;
	extern struct s_block *block;

	point_t<int> bl, tr;

	vector<int> ipin_allocation(sink_in_nodes.size(), 0);

	int local_id = 0;
	for (int i = 0; i < num_nets; ++i) {
		net_t net;

		net.vpr_id = i;
		net.local_id = local_id;

		net.source.rr_node = net_rr_terminals[i][0];
		int b = clb_net[i].node_block[0];
		int p = clb_net[i].node_block_pin[0];
		net.source.x = block[b].x;
		net.source.y = block[b].y + block[b].type->pin_height[p];

		bl.x = net.source.x;
		bl.y = net.source.y;
		tr.x = net.source.x;
		tr.y = net.source.y;

		char buffer[256];
		sprintf_rr_node(net.source.rr_node, buffer);
		zlog_debug(net_log, "Net %d source %s\n", net.vpr_id, buffer);

		net.current_source = net.source;
		/*net.previous_source.rr_node = -1;*/
		/*net.previous_source.x = -1;*/
		/*net.previous_source.y = -1;*/
		/*net.previous_source_valid = false;*/
		 
		for (int j = 1; j <= clb_net[i].num_sinks; j++) {
			sink_t sink;

			int sink_rr_node = net_rr_terminals[i][j];
			int sink_alloc = ipin_allocation[sink_rr_node];
			assert(sink_alloc < sink_in_nodes[sink_rr_node].size());
			sink.rr_node = sink_in_nodes[sink_rr_node][sink_alloc];
			++ipin_allocation[sink_rr_node];
			sink.id = j-1;
			int b = clb_net[i].node_block[j];
			int p = clb_net[i].node_block_pin[j];
			sink.x = block[b].x;
			sink.y = block[b].y + block[b].type->pin_height[p];
			sink.criticality_fac = std::numeric_limits<float>::max();
			/* initially all sink will be reached from the net's source */
			/* this is later updated during scheduling */
			sink.source = net.source;
			sink.bb_factor = bb_factor;
			sink.current_bounding_box = get_bounding_box(sink.source, sink, sink.bb_factor); 

			sprintf_rr_node(sink.rr_node, buffer);
			zlog_debug(net_log, "Net %d sink %d (%s) bounding box %d-%d %d-%d\n", net.vpr_id, sink.id, buffer, sink.current_bounding_box.xmin, sink.current_bounding_box.xmax, sink.current_bounding_box.ymin, sink.current_bounding_box.ymax);
			sink.congested_iterations = 0;

			net.sinks.push_back(sink);

			bl.x = std::min(bl.x, sink.x);
			bl.y = std::min(bl.y, sink.y);
			tr.x = std::max(tr.x, sink.x);
			tr.y = std::max(tr.y, sink.y);

			/*int inode = net_rr_terminals[i][j];*/
			/*extern struct s_block *block;*/
			/*assert(clb_net[i].node_block_pin[j] == rr_node[inode].ptc_num);*/
		}

		net.sink_routed.resize(net.sinks.size(), false);

		assert(net.sinks.size() > 0);
		net.has_sink = true;

		int rank = 0;
		for (auto &sink : net.sinks) {
			sink.distance_to_source_rank = rank++;
		}

		net.current_sink_index = 0;

		/*net.current_bounding_box = get_bounding_box(net.current_source, net.sinks[net.current_sink_index]);*/
		/*net.previous_sink_index = -1;*/
		/*net.previous_sink_valid = false;*/

		/*net.previous_bounding_box_valid = false;*/

		zlog_debug(net_log, "Net %d sorted\n", net.vpr_id);
		sprintf_rr_node(net.source.rr_node, buffer);
		zlog_debug(net_log, "%s x: %d y: %d\n", buffer, net.source.x, net.source.y);
		zlog_debug(net_log, "Sorted sinks:\n");
		for (const auto &s : net.sinks) {
			sprintf_rr_node(s.rr_node, buffer);
			zlog_debug(net_log, "%s x: %d y: %d\n", buffer, s.x, s.y);
		}

		bounding_box_t bb = get_bounding_box(bl, tr,  bb_factor);
		assert(bb.xmin == route_bb[i].xmin);
		assert(bb.ymin == route_bb[i].ymin);
		assert(bb.xmax == route_bb[i].xmax);
		assert(bb.ymax == route_bb[i].ymax);

		net.bounding_box = bb;

		//for (auto &sink : net.sinks) {
			//sink.current_bounding_box = bb;
		//}

		/*net.box.xmin = route_bb[i].xmin;*/
		/*net.box.ymin = route_bb[i].ymin;*/
		/*net.box.xmax = route_bb[i].xmax;*/
		/*net.box.ymax = route_bb[i].ymax;*/

		if (clb_net[i].is_global) {
			global_nets.push_back(net);

			zlog_debug(net_log, "Global net %d\n", i);
		} else {
			nets.push_back(net);
			++local_id;
		}
	}

	/* update pointers */
	for (auto &net : nets) {
		net.source.net = &net;
		for (auto &sink : net.sinks) {
			sink.net = &net;
		}
	}
	/*int num_local_nets = local_id;*/
	/*for (auto &net : nets) {*/
		/*net.num_local_nets = num_local_nets;*/
		/*net.overlapping_nets = new bool[num_local_nets];*/
		/*net.non_overlapping_nets = new bool[num_local_nets];*/
	/*}*/
}

void init_nets(vector<net_t> &nets, vector<net_t> &global_nets, int bb_factor)
{
	extern struct s_net *clb_net;
	extern int num_nets;
	extern int **net_rr_terminals; /* [0..num_nets-1][0..num_pins-1] */
	extern struct s_rr_node *rr_node;
	extern struct s_bb *route_bb;
	extern struct s_block *block;

	point_t<int> bl, tr;

	int local_id = 0;
	for (int i = 0; i < num_nets; ++i) {
		net_t net;

		net.vpr_id = i;
		net.local_id = local_id;

		net.source.rr_node = net_rr_terminals[i][0];
		int b = clb_net[i].node_block[0];
		int p = clb_net[i].node_block_pin[0];
		net.source.x = block[b].x;
		net.source.y = block[b].y + block[b].type->pin_height[p];

		bl.x = net.source.x;
		bl.y = net.source.y;
		tr.x = net.source.x;
		tr.y = net.source.y;

		char buffer[256];
		sprintf_rr_node(net.source.rr_node, buffer);
		zlog_debug(net_log, "Net %d source %s\n", net.vpr_id, buffer);

		net.current_source = net.source;
		/*net.previous_source.rr_node = -1;*/
		/*net.previous_source.x = -1;*/
		/*net.previous_source.y = -1;*/
		/*net.previous_source_valid = false;*/
		 
		for (int j = 1; j <= clb_net[i].num_sinks; j++) {
			sink_t sink;

			sink.rr_node = net_rr_terminals[i][j];
			sink.id = j-1;
			int b = clb_net[i].node_block[j];
			int p = clb_net[i].node_block_pin[j];
			sink.x = block[b].x;
			sink.y = block[b].y + block[b].type->pin_height[p];
			sink.criticality_fac = std::numeric_limits<float>::max();
			/* initially all sink will be reached from the net's source */
			/* this is later updated during scheduling */
			sink.source = net.source;
			sink.bb_factor = bb_factor;
			sink.current_bounding_box = get_bounding_box(sink.source, sink, sink.bb_factor); 

			sprintf_rr_node(sink.rr_node, buffer);
			zlog_debug(net_log, "Net %d sink %d (%s) bounding box %d-%d %d-%d\n", net.vpr_id, sink.id, buffer, sink.current_bounding_box.xmin, sink.current_bounding_box.xmax, sink.current_bounding_box.ymin, sink.current_bounding_box.ymax);
			sink.congested_iterations = 0;

			net.sinks.push_back(sink);

			bl.x = std::min(bl.x, sink.x);
			bl.y = std::min(bl.y, sink.y);
			tr.x = std::max(tr.x, sink.x);
			tr.y = std::max(tr.y, sink.y);

			/*int inode = net_rr_terminals[i][j];*/
			/*extern struct s_block *block;*/
			/*assert(clb_net[i].node_block_pin[j] == rr_node[inode].ptc_num);*/
		}

		net.sink_routed.resize(net.sinks.size(), false);

		assert(net.sinks.size() > 0);
		net.has_sink = true;

		int rank = 0;
		for (auto &sink : net.sinks) {
			sink.distance_to_source_rank = rank++;
		}

		net.current_sink_index = 0;

		/*net.current_bounding_box = get_bounding_box(net.current_source, net.sinks[net.current_sink_index]);*/
		/*net.previous_sink_index = -1;*/
		/*net.previous_sink_valid = false;*/

		/*net.previous_bounding_box_valid = false;*/

		zlog_debug(net_log, "Net %d sorted\n", net.vpr_id);
		sprintf_rr_node(net.source.rr_node, buffer);
		zlog_debug(net_log, "%s x: %d y: %d\n", buffer, net.source.x, net.source.y);
		zlog_debug(net_log, "Sorted sinks:\n");
		for (const auto &s : net.sinks) {
			sprintf_rr_node(s.rr_node, buffer);
			zlog_debug(net_log, "%s x: %d y: %d\n", buffer, s.x, s.y);
		}

		bounding_box_t bb = get_bounding_box(bl, tr,  bb_factor);
		assert(bb.xmin == route_bb[i].xmin);
		assert(bb.ymin == route_bb[i].ymin);
		assert(bb.xmax == route_bb[i].xmax);
		assert(bb.ymax == route_bb[i].ymax);

		net.bounding_box = bb;

		//for (auto &sink : net.sinks) {
			//sink.current_bounding_box = bb;
		//}

		/*net.box.xmin = route_bb[i].xmin;*/
		/*net.box.ymin = route_bb[i].ymin;*/
		/*net.box.xmax = route_bb[i].xmax;*/
		/*net.box.ymax = route_bb[i].ymax;*/

		if (clb_net[i].is_global) {
			global_nets.push_back(net);

			zlog_debug(net_log, "Global net %d\n", i);
		} else {
			nets.push_back(net);
			++local_id;
		}
	}

	/* update pointers */
	for (auto &net : nets) {
		net.source.net = &net;
		for (auto &sink : net.sinks) {
			sink.net = &net;
		}
	}
	/*int num_local_nets = local_id;*/
	/*for (auto &net : nets) {*/
		/*net.num_local_nets = num_local_nets;*/
		/*net.overlapping_nets = new bool[num_local_nets];*/
		/*net.non_overlapping_nets = new bool[num_local_nets];*/
	/*}*/
}

