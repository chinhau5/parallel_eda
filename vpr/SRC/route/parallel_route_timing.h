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
};

int get_num_threads();

boolean try_parallel_timing_driven_route(struct s_router_opts router_opts,
		float **net_delay, t_slack * slacks, t_ivec ** clb_opins_used_locally, boolean timing_analysis_enabled);
