void postorder_dfs(int inode, NodeColor *node_color)
{
	struct stack_item {
		int inode;
	};

	/*struct node_state {*/
		/*node_state() {*/
			/*visited = false;*/
			/*neighbor_visited = 0;*/
			/*parent_inode = -1;*/
		/*}*/
		/*bool visited;*/
		/*int neighbor_visited;*/
		/*int parent_inode;*/
	/*};*/
	/*node_state *ns = new node_state[num_rr_nodes];*/

	std::deque<stack_item> s;
	s.push_back({ inode });
	node_color[inode] = NodeColor::GRAY;

	while (!s.empty()) {
		const stack_item item = s.back();
		
		t_rr_node *current = &rr_node[item.inode];

		if (node_color[item.inode] == NodeColor::GRAY) {
			/*dzlog_debug("Already visited %d from %d\n", item.inode, ns[item.inode].parent_inode);*/
			dzlog_debug("Current: %d\n", item.inode);
			for (int iedge = 0; iedge < current->num_edges; ++iedge) {
				int neighbor_inode = current->edges[iedge];
				
				current->reachable_nodes.insert(neighbor_inode);

				dzlog_debug("Neighbor %d is color %d\n", neighbor_inode, node_color[neighbor_inode]);

				if (node_color[neighbor_inode] == NodeColor::WHITE) {
					dzlog_debug("Pushing Neighbor: %d\n", neighbor_inode);

					/*if (ns[neighbor_inode].parent_inode == -1) {*/
					/*ns[neighbor_inode].parent_inode = item.inode;*/
					/*}*/
					node_color[neighbor_inode] = NodeColor::GRAY;
					s.push_back({ neighbor_inode });
				}
			}
		} else {
			assert(node_color[item.inode] == NodeColor::BLACK);
			dzlog_debug("Already visited %d\n", item.inode);
		}

		/*int parent_inode = ns[item.inode].parent_inode;*/
		/*if (parent_inode != -1) {*/
			/*dzlog_debug("Parent: %d\n", parent_inode);*/

			/*t_rr_node *parent = &rr_node[parent_inode];*/
			/*if (++ns[parent_inode].neighbor_visited == parent->num_edges) {*/
				/*dzlog_debug("Done visiting children of %d\n", parent_inode);*/
				/*for (int iedge = 0; iedge < parent->num_edges; ++iedge) {*/
					/*int sibling_inode = parent->edges[iedge];*/
					/*dzlog_debug("Sibling: %d\n", sibling_inode);*/
					/*t_rr_node *sibling = &rr_node[sibling_inode];*/
					/*parent->reachable_nodes.insert(*/
							/*sibling->reachable_nodes.begin(), sibling->reachable_nodes.end());*/
				/*}*/
			/*}*/
		/*}*/
		assert(node_color[item.inode] == NodeColor::GRAY || node_color[item.inode] == NodeColor::BLACK);
		if (node_color[item.inode] == NodeColor::GRAY) {
			bool all_children_visited = true;
			for (int iedge = 0; iedge < current->num_edges && all_children_visited; ++iedge) {
				int child_inode = rr_node[item.inode].edges[iedge];
				bool child_visited = node_color[child_inode] == NodeColor::BLACK;
				if (!child_visited) {
					s.push_back({ child_inode });
				}
				all_children_visited = child_visited;
			}
			if (all_children_visited) {
				for (int iedge = 0; iedge < current->num_edges; ++iedge) {
					int child_inode = rr_node[item.inode].edges[iedge];
					current->reachable_nodes.insert(
							rr_node[child_inode].reachable_nodes.begin(),
							rr_node[child_inode].reachable_nodes.end());
				}
				dzlog_debug("Visited all children of %d\n", item.inode);
				node_color[item.inode] = NodeColor::BLACK;
				s.pop_back();
			}
		}
	}

	/*delete [] ns;*/
}

void postorder_dfs(int inode, NodeColor *node_color)
{
	struct stack_item {
		int inode;
	};

	/*struct node_state {*/
		/*node_state() {*/
			/*visited = false;*/
			/*neighbor_visited = 0;*/
			/*parent_inode = -1;*/
		/*}*/
		/*bool visited;*/
		/*int neighbor_visited;*/
		/*int parent_inode;*/
	/*};*/
	/*node_state *ns = new node_state[num_rr_nodes];*/

	std::deque<stack_item> s;
	s.push_back({ inode });

	while (!s.empty()) {
		const stack_item item = s.back();
		
		t_rr_node *current = &rr_node[item.inode];

		if (node_color[item.inode] != NodeColor::WHITE) {
			/*dzlog_debug("Already visited %d from %d\n", item.inode, ns[item.inode].parent_inode);*/
			dzlog_debug("Already visited %d\n", item.inode);
		} else {
			node_color[item.inode] = NodeColor::GRAY;

			dzlog_debug("Current: %d\n", item.inode);
			for (int iedge = 0; iedge < current->num_edges; ++iedge) {
				int neighbor_inode = current->edges[iedge];
				
				current->reachable_nodes.insert(neighbor_inode);

				dzlog_debug("Neighbor %d is color %d\n", neighbor_inode, node_color[neighbor_inode]);

				if (node_color[neighbor_inode] == NodeColor::WHITE) {
					dzlog_debug("Pushing Neighbor: %d\n", neighbor_inode);

					/*if (ns[neighbor_inode].parent_inode == -1) {*/
					/*ns[neighbor_inode].parent_inode = item.inode;*/
					/*}*/
					s.push_back({ neighbor_inode });
				}
			}
		}

		/*int parent_inode = ns[item.inode].parent_inode;*/
		/*if (parent_inode != -1) {*/
			/*dzlog_debug("Parent: %d\n", parent_inode);*/

			/*t_rr_node *parent = &rr_node[parent_inode];*/
			/*if (++ns[parent_inode].neighbor_visited == parent->num_edges) {*/
				/*dzlog_debug("Done visiting children of %d\n", parent_inode);*/
				/*for (int iedge = 0; iedge < parent->num_edges; ++iedge) {*/
					/*int sibling_inode = parent->edges[iedge];*/
					/*dzlog_debug("Sibling: %d\n", sibling_inode);*/
					/*t_rr_node *sibling = &rr_node[sibling_inode];*/
					/*parent->reachable_nodes.insert(*/
							/*sibling->reachable_nodes.begin(), sibling->reachable_nodes.end());*/
				/*}*/
			/*}*/
		/*}*/
		assert(node_color[item.inode] == NodeColor::GRAY || node_color[item.inode] == NodeColor::BLACK);
		if (node_color[item.inode] == NodeColor::GRAY) {
			bool all_children_visited = true;
			for (int iedge = 0; iedge < current->num_edges && all_children_visited; ++iedge) {
				int child_inode = rr_node[item.inode].edges[iedge];
				all_children_visited = node_color[child_inode] == NodeColor::BLACK;
			}
			if (all_children_visited) {
				for (int iedge = 0; iedge < current->num_edges; ++iedge) {
					int child_inode = rr_node[item.inode].edges[iedge];
					current->reachable_nodes.insert(
							rr_node[child_inode].reachable_nodes.begin(),
							rr_node[child_inode].reachable_nodes.end());
				}
				dzlog_debug("Visited all children of %d\n", item.inode);
				node_color[item.inode] = NodeColor::BLACK;
				s.pop_back();
			}
		}
	}

	/*delete [] ns;*/
}
