#ifndef FM_H
#define FM_H

//class bucket {
	//private:
		//vector<int> nodes;
	//public:
		//void add(int v)
		//{
			//nodes.push_back(v);
		//}
		
		//int remove();

		//bool empty() const;
//};
//

//#define fm_debug(...) printf(__VA_ARGS__)
#define fm_debug(...)

template<typename Graph, typename Imbalance>
class fm_partition_builder;

template<typename Graph, typename Imbalance>
class fm_partition {
	private:
	int this_pid;

	const Graph &g;

	int max_possible_gain;

	vector<int> &pid;
	int num_nodes;

	vector<bool> &locked;
	int num_locked;

	vector<pair<int, int>> node_bucket_lookup;
	vector<vector<int>> buckets;
	int highest_gain_bucket;
	int highest_gain_bucket_item_id;

	const Imbalance &imbalance;

	friend class fm_partition_builder<Graph, Imbalance>;

	fm_partition(int this_pid, const Graph &g, int max_possible_gain, vector<int> &pid, vector<bool> &locked, const Imbalance &imbalance)
		: this_pid(this_pid), g(g), max_possible_gain(max_possible_gain), pid(pid), num_nodes(0), locked(locked), num_locked(0), node_bucket_lookup(num_vertices(g), make_pair(-1, -1)), buckets(max_possible_gain*2+1), highest_gain_bucket(-1), highest_gain_bucket_item_id(-1), imbalance(imbalance)
	{
	}

	public:

	bool has_unlocked() const
	{
		return num_locked != num_nodes;
	}

	void lock(int v)
	{
		assert(!locked[v]);
		locked[v] = true;
		++num_locked;
	}

	void unlock_all()
	{
		for (const auto &bucket : buckets) {
			for (const auto &item : bucket) {
				assert(pid[item] == this_pid);
				locked[item] = false;
			}
		}
		num_locked = 0;
	}

	pair<int, int> get_highest_gain_node()
	{
		int b = get_highest_gain_bucket();

		int node;
		int gain;

		if (b >= 0) {
			node = -1;
			int min_imba = std::numeric_limits<int>::max();
			for (int i = 0; i < buckets[b].size(); ++i) {
				int v = buckets[b][i];
				int imba;
				bool is_imba = imbalance.is_imbalanced(v, (this_pid+1) % 2, imba);
				if (locked[v] || is_imba) {
					continue;
				}
				if (imba < min_imba) {
					node = v;
					min_imba = imba;
					highest_gain_bucket_item_id = i;
				}
			}
			if (node != -1) {
				gain = bucket_id_to_gain(b);
			} else {
				gain = 0;
			}
		} else {
			node = -1;
			gain = 0;
		}

		return make_pair(node, gain);
	}

	void remove_from_bucket(int v)
	{
		int bid, iid;
		std::tie(bid, iid) = node_bucket_lookup[v];
		assert(node_bucket_lookup[v] != make_pair(-1, -1));
		assert(find(begin(buckets[bid]), end(buckets[bid]), v) != end(buckets[bid]));
		buckets[bid].erase(buckets[bid].begin()+iid);
		for (const auto &item : buckets[bid]) {
			if (node_bucket_lookup[item].second > iid) {
				--node_bucket_lookup[item].second;
			}
		}
		node_bucket_lookup[v] = make_pair(-1, -1);
		/* allows the highest_gain_bucket variable to be updated */
		get_highest_gain_bucket();
	}

	void remove(int v)
	{
		assert(pid[v] == this_pid);
		pid[v] = -1;
		remove_from_bucket(v);
		--num_nodes;
		assert(num_nodes >= 0);
	}

	void add_to_bucket(int v)
	{
		int gain = get_gain(v);
		int bid = gain_to_bucket_id(gain);
		int iid = buckets[bid].size();
		assert(find(begin(buckets[bid]), end(buckets[bid]), v) == end(buckets[bid]));
		buckets[bid].push_back(v);
		assert(node_bucket_lookup[v] == make_pair(-1, -1));
		node_bucket_lookup[v] = make_pair(bid, iid);
		if (bid > highest_gain_bucket) {
			highest_gain_bucket = bid;
		}
		locked[v] = false;
	}

	void raw_add(int v)
	{
		assert(pid[v] == -1);
		pid[v] = this_pid;
		++num_nodes;
	}

	void add(int v)
	{
		assert(pid[v] == -1);
		/* set part first so that v is added to correct bucket */
		pid[v] = this_pid;
		add_to_bucket(v);
		++num_nodes;
	}

	void update_gain(int v)
	{
		assert(pid[v] == this_pid);
		int new_b = gain_to_bucket_id(get_gain(v));
		int old_b = node_bucket_lookup[v].first;
		if (new_b != old_b) {
			remove_from_bucket(v);
			add_to_bucket(v);
		}
	}

	int get_highest_gain_bucket()
	{
		while (highest_gain_bucket >= 0 && buckets[highest_gain_bucket].empty()) {
			--highest_gain_bucket;
		}
		return highest_gain_bucket;
	}

	int bucket_id_to_gain(int b) const 
	{
		return b-max_possible_gain;
	}

	int gain_to_bucket_id(int gain) const
	{
		return gain+max_possible_gain;
	}

	int get_gain(int from) const
	{
		int gain = 0;
		for (const auto &e : get_out_edges(g, from)) {
			int to = get_target(g, e);
			if (pid[from] != pid[to]) {
				++gain;
			} else {
				--gain;
			}
		}
		assert(gain >= -max_possible_gain && gain <= max_possible_gain);
		return gain;
	}

	int get_imbalance(int v)
	{
		return imbalance.get_imbalance(v);
	}
};

template<typename Graph, typename Imbalance>
class fm_partition_builder {
	private:
		const Graph *_g;
		const Imbalance *_imbalance;
		int _this_pid;
		int _max_possible_gain;
		vector<int> *_pid;
		vector<bool> *_locked;

	public:
		fm_partition_builder()
			: _g(nullptr), _imbalance(nullptr), _this_pid(-1), _max_possible_gain(-1), _pid(nullptr), _locked(nullptr)
		{
		}
		void set_this_pid(int this_pid)
		{
			_this_pid = this_pid;
		}
		void set_graph(const Graph *g)
		{
			_g = g;
		}
		void set_max_possible_gain(int max_possible_gain)
		{
			_max_possible_gain = max_possible_gain;
		}
		void set_locked(vector<bool> *locked)
		{
			_locked = locked;
		}
		void set_pid(vector<int> *pid)
		{
			_pid = pid;
		}
		void set_imbalance(const Imbalance *imbalance)
		{
			_imbalance = imbalance;
		}
		fm_partition<Graph, Imbalance> *build()
		{
			return new fm_partition<Graph, Imbalance>(_this_pid, *_g, _max_possible_gain, *_pid, *_locked, *_imbalance);
		}
};


template<typename Graph, typename Imbalance>
class fm {
	private:
		const Graph *g;
		vector<int> initial_pid;
		vector<int> pid;
		vector<bool> locked;
		vector<pair<int, int>> moves;
		fm_partition<Graph, Imbalance> *part[2];
		int max_possible_gain;
		Imbalance *imbalance;

	public:	
		~fm()
		{
			delete part[0];
			delete part[1];
		}

		bool valid_gain(int gain) const
		{
			return gain >= -max_possible_gain && gain <= max_possible_gain;
		}

		const vector<int> &get_pid() const
		{
			return pid;
		}

		void move(int v)
		{
			int from = pid[v];
			int to = (from+1) % 2;

			part[from]->remove(v);
			part[to]->add(v);
			part[to]->lock(v);

			for (const auto &e : get_out_edges(*g, v)) {
				int n = get_target(*g, e);
				part[pid[n]]->update_gain(n);
			}
		}

		void init(const Graph &_g, const vector<int> &_initial_pid, Imbalance &_imbalance)
		{
			g = &_g;
			imbalance = &_imbalance;
			initial_pid = _initial_pid;
			pid.resize(num_vertices(_g), -1);
			locked.resize(num_vertices(_g), false);
			moves.clear();

			max_possible_gain = std::numeric_limits<int>::min();
			for (const auto &v : get_vertices(_g)) {
				int num_out_edges = 0;
				for (const auto &e : get_out_edges(_g, v)) {
					++num_out_edges;
				}
				max_possible_gain = std::max(max_possible_gain, num_out_edges);
			}

			fm_partition_builder<Graph, Imbalance> builder;

			builder.set_graph(&_g);
			builder.set_locked(&locked);
			builder.set_pid(&pid);
			builder.set_max_possible_gain(max_possible_gain);
			builder.set_imbalance(&_imbalance);

			for (int i = 0; i < 2; ++i) {
				builder.set_this_pid(i);

				part[i] = builder.build();
			}

			for (const auto &v : get_vertices(_g)) {
				int p = initial_pid[v];
				assert(p == 0 || p == 1);
				part[p]->raw_add(v);	
				imbalance->add(v, p);
			}

			for (const auto &v : get_vertices(_g)) {
				part[initial_pid[v]]->add_to_bucket(v);	
			}
		}

		int get_best_prefix()
		{
			int max_gain_sum = 0;
			int gain_sum = 0;
			int best_prefix = -1;
			for (int i = 0; i < moves.size(); ++i) {
				int node, gain;
				std::tie(node, gain) = moves[i];
				gain_sum += gain;
				if (gain_sum > max_gain_sum) {
					max_gain_sum = gain_sum;
					best_prefix = i;
				}
			}
			return best_prefix;
		}

		int get_node_to_move(int &gain)
		{
			int p0_node, p0_gain, p1_node, p1_gain;

			std::tie(p0_node, p0_gain) = part[0]->get_highest_gain_node();
			std::tie(p1_node, p1_gain) = part[1]->get_highest_gain_node(); 

			int node;

			if (p0_node == -1 && p1_node == -1) {
				/* no more valid nodes */
				fm_debug("No more valid nodes to move\n");
				node = -1;
			} else if (p0_node != -1 && p1_node != -1) {
				if (p0_gain > p1_gain) {
					node = p0_node;
					gain = p0_gain;
				} else {
					node = p1_node;
					gain = p1_gain;
				}
			} else if (p0_node != -1) {
				assert(p1_node == -1);
				node = p0_node;
				gain = p0_gain;
			} else {
				assert(p0_node == -1);
				assert(p1_node != -1);
				node = p1_node;
				gain = p1_gain;
			}

			return node;
		}

		bool has_node_to_move()
		{
			return part[0]->has_unlocked() || part[1]->has_unlocked();
		}

		void move_all()
		{
			int node, gain;
			while (has_node_to_move() && (node = get_node_to_move(gain)) != -1) {
				move(node);
				imbalance->move(node, pid[node]);
				moves.emplace_back(node, gain);
				fm_debug("Moved node %d from %d to %d with gain %d\n", node, (pid[node]+1)%2, pid[node], gain);
			}
		}

		void commit(int best_prefix)
		{
			for (int i = 0; i <= best_prefix; ++i) {
				int node, gain;
				std::tie(node, gain) = moves[i];
				initial_pid[node] = (initial_pid[node]+1) % 2;
			}
			for (const auto &v : get_vertices(*g)) {
				if (pid[v] != initial_pid[v]) {
					move(v);
					imbalance->move(v, pid[v]);
				}
			}
			moves.clear();
			for (int i = 0; i < 2; ++i) {
				part[i]->unlock_all();
			}
		}

		int get_cut_size()
		{
			int cut_size = 0;
			for (const auto &e : get_edges(*g)) {
				int from = get_source(*g, e);
				int to = get_target(*g, e);

				if (pid[from] != pid[to]) {
					++cut_size;
				}
			}
			assert(cut_size % 2 == 0);
			/* dividing cut size by 2 because of duplicate edges created during conversion
			 * from directed to undirected graph */
			return cut_size/2;
		}

		int get_num_nodes(int p)
		{
			int num_nodes = 0;
			for (const auto &v : get_vertices(*g)) {
				if (pid[v] == p) {
					++num_nodes;
				}
			}
			return num_nodes;
		}

		void print_stats()
		{
			int cut_size = get_cut_size();
			int p0_size = get_num_nodes(0);
			int p1_size = get_num_nodes(1);
			int num_vertices = 0;
			for (const auto &v : get_vertices(*g)) {
				++num_vertices;
			}

			assert(p0_size + p1_size == num_vertices);

			printf("Cut size = %d Partition sizes = (%d,%d) Imbalance = %d\n", cut_size, p0_size, p1_size, abs(p0_size-p1_size));
		}

		void run()
		{
			int best_prefix;
			int num_rounds = 0;
			printf("Initial stats\n");
			print_stats();
			do {
				move_all();
				best_prefix = get_best_prefix();
				commit(best_prefix);
				printf("Round %d num moves = %d\n", num_rounds, best_prefix);
				print_stats();
				++num_rounds;
			} while (best_prefix != -1);
		}
};

#endif
