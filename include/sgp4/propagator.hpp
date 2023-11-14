#pragma once

#include "tle.hpp"
#include "vec3.hpp"

#include <chrono>

namespace sgp4 {
	struct state_vecs {
		vec3 position;
		vec3 velocity;
	};

	class propagator {
	public:
		// loading tle data
		propagator(const tle_set& set);
		void load(const tle_set& set);
	
		state_vecs run(std::chrono::utc_clock::time_point time);

	private: // internal state

	};
}