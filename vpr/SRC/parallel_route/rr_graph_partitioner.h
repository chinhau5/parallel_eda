#ifndef RR_GRAPH_PARTITIONER_H
#define RR_GRAPH_PARTITIONER_H

#include "bfs.h"

int get_rr_node_index(int x, int y, t_rr_type rr_type, int ptc, t_ivec *** L_rr_node_indices);
bool starts_at(const rr_node_property_t &rr_node, int x, int y);
std::pair<int, int> get_node_start(int inode);
void dump_rr_graph(const RRGraph &g, const char *filename);

template<typename Predicate>
struct subgraph_vertex_predicate_t {
	const Predicate &vertex_valid;

	subgraph_vertex_predicate_t(const Predicate &vertex_valid) :
		vertex_valid(vertex_valid)
	{
	}

	bool operator()(unsigned long v) const
	{
		return vertex_valid(v);
	}
};

template<typename Graph, typename Predicate>
struct subgraph_edge_predicate_t {
	const Graph &g;
	const Predicate &vertex_valid;

	subgraph_edge_predicate_t(const Graph &g, const Predicate &vertex_valid) : g(g), vertex_valid(vertex_valid) {}

	bool operator()(const RREdge &e) const
	{
		int from = get_source(g, e);
		int to = get_target(g, e);

		return vertex_valid(from) && vertex_valid(to);
	}
};

template<typename Graph>
struct partition_vertex_predicate_t {
	const Graph &g;
	const vector<int> &pid;
	int p;
	
	partition_vertex_predicate_t(const Graph &g, const vector<int> &pid, int p)
		: g(g), pid(pid), p(p)
	{
	}

	bool operator()(unsigned long v) const {
		const auto &ver = get_vertex_props(g, v);
		return (ver.type == CHANX || ver.type == CHANY) && pid[v] == p;
	}
};

class rr_graph_partitioner {
	private:
		int num_tracks;

	public:
		RRGraph orig_g;
		RRGraph undirected_orig_g;
		RRGraph channel_with_interior_g;
		RRGraph channel_without_interior_g;
		int num_levels;
		vector<int> result_pid;
        vector<vector<int>> result_pid_by_level;
		vector<vector<RRNode>> sink_in_nodes;
		vector<vector<RRNode>> ipin_in_nodes;

	private:
		vector<vector<vector<int>>> num_starting_tracks; // [type][x][y]
		vector<int> real_track_num; // [rr_node]

		int num_partitions;
		//vector<int> initial_partition;
		//int partition_index_offset;
		vector<vector<vector<vector<int>>>> num_tracks_in_partition; //[type][part][x][y]

		int total_num_chans;

	public:
		void init_graphs(t_router_opts *opts, struct s_det_routing_arch *det_routing_arch, t_direct_inf *directs, int num_directs, t_segment_inf *segment_inf, t_timing_inf *timing_inf)
		{
			extern int num_types;
			extern struct s_type_descriptor *type_descriptors;
			extern int nx, ny;
			extern struct s_grid_tile **grid;

			//free_rr_graph();

			//int warnings;

			//build_rr_graph(GRAPH_UNIDIR, num_types, type_descriptors, nx, ny, grid,
					//opts->fixed_channel_width, NULL, det_routing_arch->switch_block_type,
					//det_routing_arch->Fs, det_routing_arch->num_segment,
					//det_routing_arch->num_switch, segment_inf,
					//det_routing_arch->global_route_switch,
					//det_routing_arch->delayless_switch, *timing_inf,
					//det_routing_arch->wire_to_ipin_switch, opts->base_cost_type,
					//directs, num_directs, FALSE,
					//&warnings);

			//init_channel_only_graph(channel_with_interior_g);

			//dump_rr_graph(channel_with_interior_g, "/Volumes/DATA/rr_graph_channel_with_interior.txt");

			//init_graph(orig_g, sink_in_nodes, ipin_in_nodes);
			init_graph(orig_g);

			init_sprintf_rr_node(&orig_g);

			dump_rr_graph(orig_g, "rr_graph.txt");

			//undirected_orig_g = orig_g;

			//for (const auto &e : get_edges(orig_g)) {
				//int from = get_source(orig_g, e);
				//int to = get_target(orig_g, e);
				//add_edge(undirected_orig_g, to, from);
			//}

			//free_rr_graph();
			//for (int i = 0; i < det_routing_arch->num_segment; ++i) {
				//for (int j = 0; j < segment_inf[i].sb_len; ++j) {
					//if (j != 0 && j != segment_inf[i].sb_len-1) {
						//segment_inf[i].sb[j] = FALSE;
					//}
				//}
			//}
			//build_rr_graph(GRAPH_UNIDIR, num_types, type_descriptors, nx, ny, grid,
					//opts->fixed_channel_width, NULL, det_routing_arch->switch_block_type,
					//det_routing_arch->Fs, det_routing_arch->num_segment,
					//det_routing_arch->num_switch, segment_inf,
					//det_routing_arch->global_route_switch,
					//det_routing_arch->delayless_switch, *timing_inf,
					//det_routing_arch->wire_to_ipin_switch, opts->base_cost_type,
					//directs, num_directs, FALSE,
					//&warnings);

			//init_channel_only_graph(channel_without_interior_g);

			//dump_rr_graph(channel_without_interior_g, "/Volumes/DATA/rr_graph_channel_without_interior.txt");
		}

		void load_real_track_numbers()
		{
			extern int nx, ny;
			extern t_ivec ***rr_node_indices;

			real_track_num.resize(num_vertices(orig_g), -1);

			num_starting_tracks.resize(2);
			for (int type = 0; type < 2; ++type) {
				num_starting_tracks[type].resize(nx+2);
				for (int x = 0; x < nx+1; ++x) {
					num_starting_tracks[type][x].resize(ny+2, 0);
					for (int y = 0; y < ny+1; ++y) {
						int num_tracks;
						if (type == 0) {
							num_tracks = rr_node_indices[CHANX][y][x].nelem;
						} else {
							num_tracks = rr_node_indices[CHANY][x][y].nelem;
						}
						for (int i = 0; i < num_tracks; ++i) {
							int rr_node = get_rr_node_index(x, y, type == 0 ? CHANX : CHANY, i, rr_node_indices);
							const auto &source = get_vertex_props(channel_without_interior_g, rr_node);

							if (starts_at(source, x, y)) {
								assert(real_track_num[rr_node] == -1);
								real_track_num[rr_node] = num_starting_tracks[type][x][y]++;
							} 
						}
					}
				}
			}
			int total_num_tracks = 0;
			for (int type = 0; type < 2; ++type) {
				for (int x = 0; x < nx+1; ++x) {
					for (int y = 0; y < ny+1; ++y) {
						total_num_tracks += num_starting_tracks[type][x][y];
					}
				}
			}
			assert(total_num_tracks == total_num_chans);
		}

		void init(t_router_opts *opts, struct s_det_routing_arch *det_routing_arch, t_direct_inf *directs, int num_directs, t_segment_inf *segment_inf, t_timing_inf *timing_inf)
		{
			extern int *chan_width_x, *chan_width_y;
			assert(chan_width_x[0] == chan_width_y[0]);
			num_tracks = chan_width_x[0];
			init_graphs(opts, det_routing_arch, directs, num_directs, segment_inf, timing_inf);
			total_num_chans = 0;
			for (const auto &v : get_vertices(orig_g)) {
				const auto &ver = get_vertex_props(orig_g, v);
				if (ver.type == CHANX || ver.type == CHANY) {
					++total_num_chans;
				}
			}
			//load_real_track_numbers();
			//alloc_num_tracks_in_partition();
			result_pid.resize(num_vertices(orig_g));
		}

		//void load_initial_partitions()
		//{
			//initial_partition.resize(num_vertices(orig_g));
			//for (const auto &v : get_vertices(orig_g)) {
				//const auto &ver = get_vertex_props(orig_g, v);
				//if (ver.type == CHANX || ver.type == CHANY) {
					//assert(real_track_num[v] >= 0 && real_track_num[v] < num_tracks);
					//initial_partition[v] = real_track_num[v] % 2;
					//int type = ver.type == CHANX ? 0 : 1;
					//int x, y;
					//std::tie(x, y) = get_node_start(v);

					//++num_tracks_in_partition[type][initial_partition[v]][x][y];
				//}
			//}
		//}

		void alloc_num_tracks_in_partition()
		{
			extern int nx, ny;

			num_tracks_in_partition.resize(2);
			for (int type = 0; type < 2; ++type) {
				num_tracks_in_partition[type].resize(2);
				for (int part = 0; part < 2; ++part) {
					num_tracks_in_partition[type][part].resize(nx+2);
					for (int x = 0; x < nx+2; ++x) {
						num_tracks_in_partition[type][part][x].resize(ny+2);
					}
				}
			}
		}

		void init_num_tracks_in_partition()
		{
			extern int nx, ny;

			for (int type = 0; type < 2; ++type) {
				for (int part = 0; part < 2; ++part) {
					for (int x = 0; x < nx+2; ++x) {
						for (int y = 0; y < ny+2; ++y) {
							num_tracks_in_partition[type][part][x][y] = 0;
						}
					}
				}
			}
		}

		//void load_initial_partitions(int num_partitions)
		//{
		//initial_partition.resize(num_vertices(orig_g));
		//for (const auto &v : get_vertices(orig_g)) {
		//const auto &ver = get_vertex_props(orig_g, v);
		//if (ver.type == CHANX || ver.type == CHANY) {
		//assert(real_track_num[v] >= 0 && real_track_num[v] < num_tracks);
		//initial_partition[v] = real_track_num[v] % num_partitions;
		//int type = ver.type == CHANX ? 0 : 1;
		//int x, y;
		//std::tie(x, y) = get_node_start(v);

		//++num_tracks_in_partition[type][initial_partition[v]][x][y];
		//}
		//}
		//}

		int type_index(int v) const
		{
			const auto &ver = get_vertex_props(orig_g, v);
			assert(ver.type == CHANX || ver.type == CHANY);
			return ver.type == CHANX ? 0 : 1;
		}

		bool is_imbalanced(int v, int to, int &imba) const
		{
			int from = (to+1) % 2;
			int x, y;
			std::tie(x, y) = get_node_start(v);

			assert(from == 0 || from == 1);
			assert(to == 0 || to == 1);

			int average = num_starting_tracks[type_index(v)][x][y] / num_partitions;
			int after = num_tracks_in_partition[type_index(v)][from][x][y] - 1; 
			imba = abs(after-average);

			assert(after >= 0);

			//return after == 0;
			return imba > 3;
		}

		void move(int v, int to)
		{
			int from = (to+1) % 2;
			int x, y;
			std::tie(x, y) = get_node_start(v);

			assert(from == 0 || from == 1);
			assert(to == 0 || to == 1);

			num_tracks_in_partition[type_index(v)][to][x][y]++;
			num_tracks_in_partition[type_index(v)][from][x][y]--;
		}

		void add(int v, int to)
		{
			int x, y;
			std::tie(x, y) = get_node_start(v);

			assert(to == 0 || to == 1);
			num_tracks_in_partition[type_index(v)][to][x][y]++;
		}

		//void recursive_bipartition(int start, int end)
		//{
		//int current_num_partitions = end-start+1;
		//assert(current_num_partitions % 2 == 0);
		//int temp = current_num_partitions / 2;
		//if (temp > 1) {
		//recursive_bipartition(start, start+temp-1);
		//recursive_bipartition(start+temp, end);
		//} else {
		//partition_index_offset = start;
		//assert(current_num_partitions == 2);
		//auto current_g = make_filtered_graph(orig_g, 
		//[&] (unsigned long v) { 
		//const auto &ver = get_vertex_props(orig_g, v);
		//return (ver.type == CHANX || ver.type == CHANY) && (initial_partition[v] == start || initial_partition[v] == end);
		//},
		//[&] (unsigned long e) { 
		//const auto &from = get_vertex_props(orig_g, get_source(orig_g, e));
		//const auto &to = get_vertex_props(orig_g, get_target(orig_g, e));
		//return (from.type == CHANX || from.type == CHANY) && (initial_partition[id(from)] == start || initial_partition[id(from)] == end) &&
		//(to.type == CHANX || to.type == CHANY) && (initial_partition[id(to)] == start || initial_partition[id(to)] == end);
		//});

		//fm<decltype(current_g), rr_graph_partitioner> fm;
		//vector<int> adjusted_initial_partition(num_vertices(current_g), -1);
		//for (const auto &v : get_vertices(current_g)) {
		//adjusted_initial_partition[v] = initial_partition[v]-partition_index_offset;
		//}

		//fm.init(current_g, adjusted_initial_partition, *this);
		//fm.run();
		//}
		//}

		template<typename Graph, typename Valid>
			filtered_graph_t<Graph, subgraph_vertex_predicate_t<Valid>, subgraph_edge_predicate_t<typename Graph::base, Valid>> make_subgraph(const Graph &g, const Valid &vertex_valid)
			{
				return make_filtered_graph(g, subgraph_vertex_predicate_t<Valid>(vertex_valid), subgraph_edge_predicate_t<typename Graph::base, Valid>(g, vertex_valid));
			}

		template<typename Graph>
			void grow_initial_partition_2(const Graph &g, vector<int> &initial_partition)
			{
				struct visitor_t {
					set<int> visited_edges;
					set<int> visited_nodes;
					set<tuple<int, int, t_rr_type, bool>> visited_channels;
					int num_nodes_discovered;

					visitor_t() :
						num_nodes_discovered(0)
					{
					}

					bool tree_edge(const RREdge &e, const Graph &g)
					{
						bool expand = false;
						assert(false);
						//if (visited_nodes.find(get_source(g, e)) != end(visited_nodes)) {
							//bool recorded_to = record(g, get_target(g, e));
							//if (recorded_to) {
								//assert(visited_edges.find(e) == end(visited_edges));
								//visited_edges.insert(e);
								//expand = true;
							//}
						//}
						return expand;
					}

					void examine_edge(const RREdge &e, const Graph &g)
					{
					}

					bool discover_vertex(int v, const Graph &g)
					{
						char buffer[256];
						sprintf_rr_node(v, buffer);
						zlog_level(delta_log, ROUTER_V3, "Current: %s\n", buffer);
						++num_nodes_discovered;
					}

					void examine_vertex(int v, const Graph &g)
					{
					}

					bool record_impl(const Graph &g, int v)
					{
						const auto &ver = get_vertex_props(g, v);
						const auto &start = get_node_start(v);
						const auto &key = make_tuple(start.first, start.second, ver.type, ver.inc_direction);
						bool recorded;
						if (visited_channels.find(key) == visited_channels.end()) {
							visited_channels.insert(key);
							recorded = true;
						} else {
							recorded = false;
						}
						return recorded;
					}

					bool record(const Graph &g, int v)
					{
						bool recorded = record_impl(g, v);
						if (recorded) {
							assert(visited_nodes.find(v) == end(visited_nodes));
							visited_nodes.insert(v);
						}
						return recorded;
					}

					bool record_source(const Graph &g, int v)
					{	
						record_impl(g, v);

						assert(visited_nodes.find(v) == end(visited_nodes));
						visited_nodes.insert(v);

						return true;
					}
				};

				using clock = std::chrono::high_resolution_clock;

				vector<VertexColor> color(num_vertices(g), VertexColor::WHITE);
				int num_nodes = 0;
				set<int> all_visited_nodes;
				vector<vector<int>> visited_nodes;
				for (const auto &v : get_vertices(g)) {
					if (color[v] == VertexColor::WHITE) {
						visitor_t visitor;
						visitor.record_source(g, v);

						auto start = clock::now();
						bfs(g, { static_cast<unsigned long>(v) }, color, visitor);
						auto time = clock::now()-start;

						start = clock::now();
						assert(visitor.visited_nodes.size() == visitor.visited_edges.size()+1);
						assert(visitor.visited_nodes.find(v) != end(visitor.visited_nodes));
						for (const auto &v : visitor.visited_nodes) {
							//initial_partition[v] = current_pid;
							assert(color[v] == VertexColor::BLACK);
							assert(all_visited_nodes.find(v) == all_visited_nodes.end());
							all_visited_nodes.insert(v);
						}
						visited_nodes.emplace_back(begin(visitor.visited_nodes), end(visitor.visited_nodes));
						//for (const auto &v : get_vertices(g)) {
							//if (all_visited_nodes.find(v) == end(all_visited_nodes)) {
								//assert(color[v] == VertexColor::WHITE);
							//} else {
								//assert(color[v] == VertexColor::BLACK);
							//}
						//}
						auto update_time = clock::now()-start;

						//num_nodes_in_partition[current_pid] += visitor.visited_nodes.size();
						//current_pid = (current_pid+1) % 2;

						//printf("ver: %d time = %lld ms update time = %lld ms visited nodes = %lu\n", v, std::chrono::duration_cast<std::chrono::milliseconds>(time).count(), std::chrono::duration_cast<std::chrono::milliseconds>(update_time).count(), visitor.visited_nodes.size());
						//if (visitor.visited_nodes.size() == 1) {
							//printf("
						//}
					}
					++num_nodes;
				}
				std::sort(begin(visited_nodes), end(visited_nodes), [] (const vector<int> &a, const vector<int> &b) -> bool { return a.size()>b.size(); });
				vector<int> num_nodes_in_partition(2, 0);
				int current_pid = 0;
				for (const auto &nodes : visited_nodes) {
					num_nodes_in_partition[current_pid] += nodes.size();
					for (const auto &v : nodes) {
						initial_partition[v] = current_pid;
					}
					current_pid = (current_pid+1) % 2;
				}
				assert(num_nodes_in_partition[0] + num_nodes_in_partition[1] == num_nodes);
				printf("P0 Size = %d P1 Size = %d Imbalance: %d\n", num_nodes_in_partition[0], num_nodes_in_partition[1], num_nodes_in_partition[0]-num_nodes_in_partition[1]);
			}

		template<typename Graph>
			void grow_initial_partition_3(const Graph &g, vector<int> &initial_partition)
			{
				struct visitor_t {
					set<RREdge> visited_edges;
					set<int> visited_nodes;
					set<tuple<int, int, t_rr_type, bool>> visited_channels;
					int num_nodes_discovered;
					const vector<int> &pid;
					bool reachable_partition[2];

					visitor_t(const vector<int> &pid) :
						num_nodes_discovered(0) , pid(pid)
					{
						std::fill(begin(reachable_partition), end(reachable_partition), false);
						assert(!reachable_partition[0]);
						assert(!reachable_partition[1]);
					}

					bool tree_edge(const RREdge &e, const Graph &g)
					{
						bool expand = false;
						if (visited_nodes.find(get_source(g, e)) != end(visited_nodes)) {
							bool recorded_to = record(g, get_target(g, e));
							if (recorded_to) {
								assert(visited_edges.find(e) == end(visited_edges));
								visited_edges.insert(e);
								expand = true;
							}
						}
						return expand;
					}

					void examine_edge(const RREdge &e, const Graph &g)
					{
						int to = get_target(g, e);
						int to_pid = pid[to];
						if (to_pid != -1) {
							assert(to_pid == 0 || to_pid == 1);
							if (!reachable_partition[to_pid]) {
								reachable_partition[to_pid] = true;
							}
						}
					}

					void discover_vertex(int v, const Graph &g)
					{
						char buffer[256];
						sprintf_rr_node(v, buffer);
						zlog_level(delta_log, ROUTER_V3, "Current: %s\n", buffer);
						++num_nodes_discovered;
					}

					void examine_vertex(int v, const Graph &g)
					{
					}

					bool record_impl(const Graph &g, int v)
					{
						const auto &ver = get_vertex_props(g, v);
						const auto &start = get_node_start(v);
						const auto &key = make_tuple(start.first, start.second, ver.type, ver.inc_direction);
						bool recorded;
						if (visited_channels.find(key) == visited_channels.end()) {
							visited_channels.insert(key);
							recorded = true;
						} else {
							recorded = false;
						}
						return recorded;
					}

					bool record(const Graph &g, int v)
					{
						bool recorded = record_impl(g, v);
						if (recorded) {
							assert(visited_nodes.find(v) == end(visited_nodes));
							visited_nodes.insert(v);
						}
						return recorded;
					}

					bool record_source(const Graph &g, int v)
					{	
						assert(record_impl(g, v));

						assert(visited_nodes.find(v) == end(visited_nodes));
						visited_nodes.insert(v);

						return true;
					}
				};

				using clock = std::chrono::high_resolution_clock;

				vector<VertexColor> color(num_vertices(g), VertexColor::WHITE);
				int num_nodes = 0;
				set<int> all_visited_nodes;
				vector<int> num_nodes_in_partition(2, 0);
				//vector<int> sorted;
				//topological_sort(g, sorted);
				for (const auto &v : get_vertices(g)) {
					if (color[v] == VertexColor::WHITE) {
						visitor_t visitor(initial_partition);
						visitor.record_source(g, v);

						auto start = clock::now();
						bfs(g, { static_cast<unsigned long>(v) }, color, visitor);
						auto time = clock::now()-start;

						assert(visitor.visited_nodes.size() == visitor.visited_edges.size()+1);
						assert(visitor.visited_nodes.find(v) != end(visitor.visited_nodes));
						int target_partition;
						if (num_nodes_in_partition[0] == 0) {
							target_partition = 0;
						} else if (num_nodes_in_partition[1] == 0) {
							target_partition = 1;
						} else {
							if ((visitor.reachable_partition[0] && visitor.reachable_partition[1]) ||
									(!visitor.reachable_partition[0] && !visitor.reachable_partition[1])) {
								if (num_nodes_in_partition[0] < num_nodes_in_partition[1]) {
									target_partition = 0;
								} else {
									target_partition = 1;
								}
							} else {
								if (visitor.reachable_partition[0]) {
									target_partition = 0;
								} else {
									assert(visitor.reachable_partition[1]); 
									target_partition = 1;
								} 
							}
						}
						for (const auto &v : visitor.visited_nodes) {
							assert(color[v] == VertexColor::BLACK);
							assert(all_visited_nodes.find(v) == all_visited_nodes.end());
							initial_partition[v] = target_partition;
							all_visited_nodes.insert(v);
						}
						num_nodes_in_partition[target_partition] += visitor.visited_nodes.size();
					}
					++num_nodes;
				}
				assert(num_nodes_in_partition[0] + num_nodes_in_partition[1] == num_nodes);
				printf("P0 Size = %d P1 Size = %d Imbalance: %d\n", num_nodes_in_partition[0], num_nodes_in_partition[1], num_nodes_in_partition[0]-num_nodes_in_partition[1]);
			}

		template<typename Graph>
			void grow_initial_partition(const Graph &g, vector<int> &initial_partition)
			{
				struct visitor_t {
					set<RREdge> visited_edges;
					set<int> visited_nodes;
					set<tuple<int, int, t_rr_type, bool>> visited_channels;
					int num_nodes_discovered;
					const vector<int> &pid;
					bool reachable_partition[2];

					visitor_t(const vector<int> &pid) :
						num_nodes_discovered(0) , pid(pid)
					{
						std::fill(begin(reachable_partition), end(reachable_partition), false);
						assert(!reachable_partition[0]);
						assert(!reachable_partition[1]);
					}

					bool tree_edge(const RREdge &e, const Graph &g)
					{
						bool expand = false;
						if (visited_nodes.find(get_source(g, e)) != end(visited_nodes)) {
							bool recorded_to = record(g, get_target(g, e));
							if (recorded_to) {
								assert(visited_edges.find(e) == end(visited_edges));
								visited_edges.insert(e);
								expand = true;
							}
						}
						return expand;
					}

					void examine_edge(const RREdge &e, const Graph &g)
					{
						int to = get_target(g, e);
						int to_pid = pid[to];
						if (to_pid != -1) {
							assert(to_pid == 0 || to_pid == 1);
							if (!reachable_partition[to_pid]) {
								reachable_partition[to_pid] = true;
							}
						}
					}

					void discover_vertex(int v, const Graph &g)
					{
						char buffer[256];
						sprintf_rr_node(v, buffer);
						zlog_level(delta_log, ROUTER_V3, "Current: %s\n", buffer);
						++num_nodes_discovered;
					}

					void examine_vertex(int v, const Graph &g)
					{
					}

					bool record_impl(const Graph &g, int v)
					{
						const auto &ver = get_vertex_props(g, v);
						const auto &start = get_node_start(v);
						const auto &key = make_tuple(start.first, start.second, ver.type, ver.inc_direction);
						bool recorded;
						if (visited_channels.find(key) == visited_channels.end()) {
							visited_channels.insert(key);
							recorded = true;
						} else {
							recorded = false;
						}
						return recorded;
					}

					bool record(const Graph &g, int v)
					{
						bool recorded = record_impl(g, v);
						if (recorded) {
							assert(visited_nodes.find(v) == end(visited_nodes));
							visited_nodes.insert(v);
						}
						return recorded;
					}

					bool record_source(const Graph &g, int v)
					{	
						assert(record_impl(g, v));

						assert(visited_nodes.find(v) == end(visited_nodes));
						visited_nodes.insert(v);

						return true;
					}
				};

				using clock = std::chrono::high_resolution_clock;

				vector<VertexColor> color(num_vertices(g), VertexColor::WHITE);
				int num_nodes = 0;
				set<int> all_visited_nodes;
				vector<int> num_nodes_in_partition(2, 0);
				//vector<int> sorted;
				//topological_sort(g, sorted);
				for (const auto &v : get_vertices(g)) {
					if (color[v] == VertexColor::WHITE) {
						visitor_t visitor(initial_partition);
						visitor.record_source(g, v);

						auto start = clock::now();
						bfs(g, { static_cast<unsigned long>(v) }, color, visitor);
						auto time = clock::now()-start;

						assert(visitor.visited_nodes.size() == visitor.visited_edges.size()+1);
						assert(visitor.visited_nodes.find(v) != end(visitor.visited_nodes));
						int target_partition;
						if (num_nodes_in_partition[0] == 0) {
							target_partition = 0;
						} else if (num_nodes_in_partition[1] == 0) {
							target_partition = 1;
						} else {
							if ((visitor.reachable_partition[0] && visitor.reachable_partition[1]) ||
									(!visitor.reachable_partition[0] && !visitor.reachable_partition[1])) {
								if (num_nodes_in_partition[0] < num_nodes_in_partition[1]) {
									target_partition = 0;
								} else {
									target_partition = 1;
								}
							} else {
								if (visitor.reachable_partition[0]) {
									target_partition = 0;
								} else {
									assert(visitor.reachable_partition[1]); 
									target_partition = 1;
								} 
							}
						}
						for (const auto &v : visitor.visited_nodes) {
							assert(color[v] == VertexColor::BLACK);
							assert(all_visited_nodes.find(v) == all_visited_nodes.end());
							initial_partition[v] = target_partition;
							all_visited_nodes.insert(v);
						}
						num_nodes_in_partition[target_partition] += visitor.visited_nodes.size();
					}
					++num_nodes;
				}
				assert(num_nodes_in_partition[0] + num_nodes_in_partition[1] == num_nodes);
				printf("P0 Size = %d P1 Size = %d Imbalance: %d\n", num_nodes_in_partition[0], num_nodes_in_partition[1], num_nodes_in_partition[0]-num_nodes_in_partition[1]);
			}

		template<typename Graph>
			void recursive_bipartition(const Graph &g, const Graph &undirected_g, int level, int node)
			{
				init_num_tracks_in_partition();

				vector<int> initial_partition(num_vertices(g), -1);

				/* get initial partition by BFS */
				grow_initial_partition(g, initial_partition);

				/* refine partition */
				//fm<Graph, rr_graph_partitioner> fm;
				//fm.init(undirected_g, initial_partition, *this);
				//fm.run();

				//const vector<int> &pid = fm.get_pid();
				const vector<int> &pid = initial_partition;

                for (const auto &v : get_vertices(g)) {
                    assert(result_pid_by_level[level][v] == -1);
                    result_pid_by_level[level][v] = 2*node+pid[v];
                }

				if (level > 0) {
					//auto valid_vertex_0 = [] (unsigned long v) -> bool { return true; };
					//auto valid_vertex_1 = [] (unsigned long v) -> bool { return true; };

					for (int i = 0; i < 2; ++i) {
						recursive_bipartition(
								make_subgraph(g.base_graph(), partition_vertex_predicate_t<typename Graph::base>(g.base_graph(), pid, i)),
								make_subgraph(undirected_g.base_graph(), partition_vertex_predicate_t<typename Graph::base>(undirected_g.base_graph(), pid, i)),
								level-1, 2*node+i);
					}
				} else {
					for (const auto &v : get_vertices(g)) {
						assert(result_pid[v] == -1);
						assert(pid[v] == 0 || pid[v] == 1);
						result_pid[v] = 2*node+pid[v];
					}
					//printf("last level index: %d %d\n", 2*node, 2*node+1);
				}
			}

		void partition_without_ipin(int _num_partitions, vector<RRGraph *> &graphs)
		{
			assert(std::pow(2, std::log2(_num_partitions)) == _num_partitions);

			num_partitions = _num_partitions;
            num_levels = std::log2(num_partitions)+1;

            result_pid_by_level.resize(num_levels, vector<int>(num_vertices(orig_g), -1));

			for (const auto &v : get_vertices(orig_g)) {
				const auto &props = get_vertex_props(orig_g, v);
				if (props.type == CHANX || props.type == CHANY) {
					result_pid_by_level[num_levels-1][v] = 0;
				} else {
					assert(result_pid_by_level[num_levels-1][v] == -1);
				}
			}

			std::fill(begin(result_pid), end(result_pid), -1);

			if (num_partitions == 1) {
				assert(num_levels == 1);
				for (const auto &v : get_vertices(orig_g)) {
					const auto &props = get_vertex_props(orig_g, v);
					if (props.type == CHANX || props.type == CHANY) {
						result_pid[v] = 0;
					} else {
						assert(result_pid[v] == -1);
					}
				}
			} else {
				vector<int> pid(num_vertices(orig_g), 0);

				recursive_bipartition(
						make_subgraph(orig_g.base_graph(), partition_vertex_predicate_t<typename decltype(orig_g)::base>(orig_g.base_graph(), pid, 0)),
						make_subgraph(undirected_orig_g.base_graph(), partition_vertex_predicate_t<typename decltype(undirected_orig_g)::base>(undirected_orig_g.base_graph(), pid, 0)), num_levels-2, 0);
			}

			for (int level = 0; level < result_pid_by_level.size(); ++level) {
				int max_pid = 0;
				for (int i = 0; i < result_pid_by_level[level].size(); ++i) {
					max_pid = std::max(max_pid, result_pid_by_level[level][i]);
				}
				assert(max_pid == pow(2, result_pid_by_level.size()-level-1)-1);
				vector<int> num_nodes(max_pid+1, 0);
				for (const auto &v : get_vertices(orig_g)) {
					const auto &props = get_vertex_props(orig_g, v);
					if (props.type == CHANX || props.type == CHANY) {
						++num_nodes[result_pid_by_level[level][v]];
					} else {
						assert(result_pid_by_level[level][v] == -1);
					}
				}
				for (int i = 0; i < num_nodes.size(); ++i) {
					printf("level %d p%d num nodes: %d\n", level, i, num_nodes[i]);
				}
			}

			//int num_nodes_partitioned = 0;
			//vector<int> num_nodes_in_partition(num_partitions, 0);
			//for (const auto &p : result_pid) {
				//if (p != -1) {
					//assert(p >= 0 && p < num_partitions);
					//++num_nodes_in_partition[p];
					//++num_nodes_partitioned;
				//}
			//}
			//int num_ipins = 0;
			//for (const auto &v : get_vertices(orig_g)) {
				//const auto &props = get_vertex_props(orig_g, v);
				//if (props.type == IPIN) {
					//++num_ipins;
				//}
			//}
			//assert(num_nodes_partitioned == total_num_chans + num_ipins_partitioned);

			//for (int i = 0; i < num_partitions; ++i) {
				//printf("Num nodes in partition %d = %d\n", i, num_nodes_in_partition[i]);
				//char filename[256];
				//sprintf(filename, "/Volumes/DATA/graph_part_%d.txt", i);
				//FILE *file = fopen(filename, "w");
				//for (const auto &v : get_vertices(orig_g)) {
					//if (result_pid[v] != i) {
						//continue;
					//}
					//const auto &ver = get_vertex_props(orig_g, v);
					//if ((ver.type == CHANX || ver.type == CHANY)) {
						//const auto &start = get_node_start(v);
						//fprintf(file, "%d %d\n", start.first, start.second);
					//}
				//}
				//fclose(file);

				//RRGraph &new_g = *new RRGraph;

				//add_vertex(new_g, num_vertices(orig_g));
				//for (const auto &v : get_vertices(orig_g)) {
					//get_vertex_props(new_g, v) = get_vertex_props(orig_g, v);
				//}

				//for (const auto &e : get_edges(orig_g)) {
					//int from = get_source(orig_g, e);
					//int to = get_target(orig_g, e);

					//const auto &from_ver = get_vertex_props(orig_g, from);
					//const auto &to_ver = get_vertex_props(orig_g, to);

					//if ((from_ver.type == CHANX || from_ver.type == CHANY)
							//&& (to_ver.type == CHANX || to_ver.type == CHANY)) {
						//assert(result_pid[from] != -1 && result_pid[to] != -1);
						//if (result_pid[from] == i && result_pid[from] == result_pid[to]) {
							//const RREdge &new_e = add_edge(new_g, from, to);
							//get_edge_props(new_g, new_e) = get_edge_props(orig_g, e);
						//}
					//} else {
						//assert(result_pid[from] == -1 || result_pid[to] == -1);
						//const RREdge &new_e = add_edge(new_g, from, to);
						//get_edge_props(new_g, new_e) = get_edge_props(orig_g, e);
					//}
				//}

				//graphs.push_back(&new_g);
			//}
		}

		void partition(int _num_partitions, vector<RRGraph *> &graphs)
		{
			assert(std::pow(2, std::log2(_num_partitions)) == _num_partitions);

			num_partitions = _num_partitions;
            num_levels = std::log2(num_partitions)+1;

            result_pid_by_level.resize(num_levels, vector<int>(num_vertices(orig_g), -1));

			for (const auto &v : get_vertices(orig_g)) {
				const auto &props = get_vertex_props(orig_g, v);
				if (props.type == CHANX || props.type == CHANY) {
					result_pid_by_level[num_levels-1][v] = 0;
				} else {
					assert(result_pid_by_level[num_levels-1][v] == -1);
				}
			}

			std::fill(begin(result_pid), end(result_pid), -1);

			if (num_partitions == 1) {
				assert(num_levels == 1);
				for (const auto &v : get_vertices(orig_g)) {
					const auto &props = get_vertex_props(orig_g, v);
					if (props.type == CHANX || props.type == CHANY || props.type == IPIN) {
						result_pid[v] = 0;
					} else {
						assert(result_pid[v] == -1);
					}
				}
			} else {
				vector<int> pid(num_vertices(orig_g), 0);

				recursive_bipartition(
						make_subgraph(orig_g.base_graph(), partition_vertex_predicate_t<typename decltype(orig_g)::base>(orig_g.base_graph(), pid, 0)),
						make_subgraph(undirected_orig_g.base_graph(), partition_vertex_predicate_t<typename decltype(undirected_orig_g)::base>(undirected_orig_g.base_graph(), pid, 0)), num_levels-2, 0);
			}

			int num_ipins_partitioned = 0;
			for (int level = 0; level < result_pid_by_level.size(); ++level) {
				int cur_num_partitions = pow(2, result_pid_by_level.size()-level-1);

				vector<vector<int>> sink_partition_connect_count(num_vertices(orig_g), vector<int>(cur_num_partitions, 0));

				for (const auto &v : get_vertices(orig_g)) {
					const auto &rr_node_p = get_vertex_props(orig_g, v);

					if (rr_node_p.type == IPIN) {
						assert(num_out_edges(orig_g, v) == 1);

						const auto &e = *begin(get_out_edges(orig_g, v));

						const auto &sink_rr_node = get_target(orig_g, e);

						assert(get_vertex_props(orig_g, sink_rr_node).type == SINK);

						if (!ipin_in_nodes[v].empty()) {
							int min_partition = -1;
							int min_count = std::numeric_limits<int>::max();
							vector<int> num_in_nodes_per_partition(cur_num_partitions, 0);
							for (const auto &in_node : ipin_in_nodes[v]) {
								const auto &in_node_p = get_vertex_props(orig_g, in_node);

								assert(in_node_p.type == CHANX || in_node_p.type == CHANY);

								int part = result_pid_by_level[level][in_node];

								assert(part >= 0 && part < cur_num_partitions);

								++num_in_nodes_per_partition[part];

								if (sink_partition_connect_count[sink_rr_node][part] < min_count) {
									min_partition = part;
									min_count = sink_partition_connect_count[sink_rr_node][part];
								}
							}
							assert(min_partition != -1);

							sink_partition_connect_count[sink_rr_node][min_partition] += num_in_nodes_per_partition[min_partition];
							assert(result_pid_by_level[level][v] == -1);
							result_pid_by_level[level][v] = min_partition;

							++num_ipins_partitioned;
						}
					}
				}
			}

			for (int level = 0; level < result_pid_by_level.size(); ++level) {
				int max_pid = 0;
				for (int i = 0; i < result_pid_by_level[level].size(); ++i) {
					max_pid = std::max(max_pid, result_pid_by_level[level][i]);
				}
				assert(max_pid == pow(2, result_pid_by_level.size()-level-1)-1);
				vector<int> num_nodes(max_pid+1, 0);
				for (const auto &v : get_vertices(orig_g)) {
					const auto &props = get_vertex_props(orig_g, v);
					if (props.type == CHANX || props.type == CHANY || props.type == IPIN) {
						int part = result_pid_by_level[level][v];
						if (part == -1) {
							assert(props.type == IPIN);
						} else {
							++num_nodes[part];
						}
					} else {
						assert(result_pid_by_level[level][v] == -1);
					}
				}
				for (int i = 0; i < num_nodes.size(); ++i) {
					printf("level %d p%d num nodes: %d\n", level, i, num_nodes[i]);
				}
			}

			//int num_nodes_partitioned = 0;
			//vector<int> num_nodes_in_partition(num_partitions, 0);
			//for (const auto &p : result_pid) {
				//if (p != -1) {
					//assert(p >= 0 && p < num_partitions);
					//++num_nodes_in_partition[p];
					//++num_nodes_partitioned;
				//}
			//}
			//int num_ipins = 0;
			//for (const auto &v : get_vertices(orig_g)) {
				//const auto &props = get_vertex_props(orig_g, v);
				//if (props.type == IPIN) {
					//++num_ipins;
				//}
			//}
			//assert(num_nodes_partitioned == total_num_chans + num_ipins_partitioned);

			//for (int i = 0; i < num_partitions; ++i) {
				//printf("Num nodes in partition %d = %d\n", i, num_nodes_in_partition[i]);
				//char filename[256];
				//sprintf(filename, "/Volumes/DATA/graph_part_%d.txt", i);
				//FILE *file = fopen(filename, "w");
				//for (const auto &v : get_vertices(orig_g)) {
					//if (result_pid[v] != i) {
						//continue;
					//}
					//const auto &ver = get_vertex_props(orig_g, v);
					//if ((ver.type == CHANX || ver.type == CHANY)) {
						//const auto &start = get_node_start(v);
						//fprintf(file, "%d %d\n", start.first, start.second);
					//}
				//}
				//fclose(file);

				//RRGraph &new_g = *new RRGraph;

				//add_vertex(new_g, num_vertices(orig_g));
				//for (const auto &v : get_vertices(orig_g)) {
					//get_vertex_props(new_g, v) = get_vertex_props(orig_g, v);
				//}

				//for (const auto &e : get_edges(orig_g)) {
					//int from = get_source(orig_g, e);
					//int to = get_target(orig_g, e);

					//const auto &from_ver = get_vertex_props(orig_g, from);
					//const auto &to_ver = get_vertex_props(orig_g, to);

					//if ((from_ver.type == CHANX || from_ver.type == CHANY)
							//&& (to_ver.type == CHANX || to_ver.type == CHANY)) {
						//assert(result_pid[from] != -1 && result_pid[to] != -1);
						//if (result_pid[from] == i && result_pid[from] == result_pid[to]) {
							//const RREdge &new_e = add_edge(new_g, from, to);
							//get_edge_props(new_g, new_e) = get_edge_props(orig_g, e);
						//}
					//} else {
						//assert(result_pid[from] == -1 || result_pid[to] == -1);
						//const RREdge &new_e = add_edge(new_g, from, to);
						//get_edge_props(new_g, new_e) = get_edge_props(orig_g, e);
					//}
				//}

				//graphs.push_back(&new_g);
			//}
		}
};

#endif
