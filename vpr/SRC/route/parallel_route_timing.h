template<typename T>
class Interval {
	private:
		T low;
		T high;
	public:
		Interval(T _low, T _high) {
			set(_low, _high);
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

		void set(const T &_low, const T &_high) {
			if (_low < _high) {
				low = _low;
				high = _high;
			} else {
				low = _high;
				high = _low;
			}
		}

		void set(const Interval &other) {
			low = other.low;
			high = other.high;
		}

		T size() const {
			return abs(high-low);
		}

		T getLow() const { return low; }
		T getHigh() const { return high; }
};

class BoundingBox {
	private:
		Interval<int> horizontal;
		Interval<int> vertical;
	public:
		BoundingBox(const Interval<int> &horizontal, const Interval<int> &vertical)
			: horizontal(horizontal), vertical(vertical)
		{
		}

		void set_horizontal(int low, int high) {
			horizontal.set(low, high);
		}
		void set_horizontal(const Interval<int> &hor) {
			horizontal.set(hor);
		}
		void set_vertical(int low, int high) {
			vertical.set(low, high);
		}
		void set_vertical(const Interval<int> &vert) {
			vertical.set(vert);
		}

		const Interval<int> &get_horizontal() const {
			return horizontal;
		}

		const Interval<int> &get_vertical() const {
			return vertical;
		}

		int area() const {
			return (horizontal.size()+1) * (vertical.size()+1);
		}
};

int get_num_threads();
boolean try_parallel_timing_driven_route_top(struct s_router_opts router_opts,
		t_net_timing *net_timing, t_ivec ** clb_opins_used_locally, boolean timing_analysis_enabled);

boolean try_new_cost_parallel_timing_driven_route_top(struct s_router_opts router_opts,
		t_net_timing *net_timing, t_ivec ** clb_opins_used_locally, boolean timing_analysis_enabled);

boolean try_fine_grained_parallel_timing_driven_route_top(struct s_router_opts router_opts,
		t_net_timing *net_timing, t_ivec ** clb_opins_used_locally, boolean timing_analysis_enabled);

boolean try_fine_grained_parallel_timing_driven_route_top_2(struct s_router_opts router_opts,
		t_net_timing *net_timing, t_ivec ** clb_opins_used_locally, boolean timing_analysis_enabled);
