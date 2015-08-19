template<typename T>
class Interval {
	private:
		T low;
		T high;
	public:
		Interval(T _low, T _high) {
			if (_low < _high) {
				low = _low;
				high = _high;
			} else {
				low = _high;
				high = _low;
			}
		}
		bool contains(const T &val) const {
			return (val >= low && val <= high);
		}

		T intersect_size(const Interval &other) const {
			return std::max(0, std::min(high, other.high) - std::max(low, other.low));
		}

		bool intersects(const Interval &other) const {
			return low <= other.high && other.low <= high; 
		}
};

int get_num_threads();
boolean try_parallel_timing_driven_route(struct s_router_opts router_opts,
		t_net_timing *net_timing, t_ivec ** clb_opins_used_locally, boolean timing_analysis_enabled);

