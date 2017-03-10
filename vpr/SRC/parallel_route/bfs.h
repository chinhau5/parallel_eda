#ifndef BFS_H
#define BFS_H

enum class VertexColor {
	WHITE,
	GRAY,
	BLACK
};

template<typename Graph, typename Visitor>
//void bfs(const RRGraph &g, const RRNode &rr_node, bool horizontal, bool inc_direction, vector<bool> &visited, vector<int> &visited_edges) 
void bfs(const Graph &g, const std::vector<unsigned long> &rr_nodes, std::vector<VertexColor> &color, Visitor &visitor)
{
	//struct queue_item {
		//const RRNode *node;
		//const RREdge *edge;
	//};

	std::queue<int> s;
	for (auto rr_node : rr_nodes) {
		assert(color[rr_node] == VertexColor::WHITE);
		visitor.discover_vertex(rr_node, g);
		s.push(rr_node);
		color[rr_node] = VertexColor::GRAY;
	}

	//extern int nx, ny;
	//vector<vector<map<pair<int, int>, int>>> visit_state(nx+2, vector<map<pair<int, int>, int>>(ny+2));
	//for (int x = 0; x < nx+2; ++x) {
		//for (int y = 0; y < ny+2; ++y) {
			//visit_state[x][y]
		//}
	//}
	//
	//char buffer[256];
	//extern zlog_category_t *delta_log;

	while (!s.empty()) {
		int current = s.front();
		s.pop();
		visitor.examine_vertex(current, g);
		//int rr_node_id = id(*current);
		//sprintf_rr_node(rr_node_id, buffer);
		//zlog_level(delta_log, ROUTER_V3, "Current: %s Add: %d\n", buffer, item.add ? 1 : 0);

		//if (get_vertex_props(g, rr_node_id).type != CHANX && get_vertex_props(g, rr_node_id).type != CHANY) {assert(false);}
		//
		//bool visit_current = visit(*current);
		//if (item.add && visit_current) {
		//visited_nodes.push_back(rr_node_id);

		//if (item.edge) {
		//visited_edges.push_back(id(*item.edge));
		//}
		//}

		for (auto e : get_out_edges(g, current)) {
			visitor.examine_edge(e, g);
			int neighbor = get_target(g, e);
			//if (!((current.type == CHANX || current.type == CHANY) &&
			//(neighbor.type == CHANX || neighbor.type == CHANY))
			//|| current.inc_direction == neighbor.inc_direction) {
			//s.push({ &neighbor, &e });
			//}

			VertexColor vc = color[neighbor];
			if (vc == VertexColor::WHITE) {
				if (visitor.tree_edge(e, g)) {
					visitor.discover_vertex(neighbor, g);
					s.push(neighbor);
					color[neighbor] = VertexColor::GRAY;
				}
				//sprintf_rr_node(id(neighbor), buffer);
				//zlog_level(delta_log, ROUTER_V3, "\tNeighbor: %s Parent add: %d Parent visit: %d\n", buffer, item.add ? 1 : 0, visit_current ? 1 : 0);
			} else if (vc == VertexColor::GRAY) {
			} else {
				assert(vc == VertexColor::BLACK);
				//visitor.black_target(
			}

			//if (current.type == CHANX || current.type == CHANY) {
			//if ((horizontal && neighbor.type == CHANX) || (!horizontal && neighbor.type == CHANY)) {
			//if (neighbor.inc_direction == inc_direction) {
			//s.push({ &neighbor, &e });
			//} else {
			//printf("Not pushing");
			//}
			//} else {
			//s.push({ &neighbor, &e });
			//}
			//} else {
			//s.push({ &neighbor, &e });
			//}
		} 
		color[current] = VertexColor::BLACK;
	}
}

#endif
