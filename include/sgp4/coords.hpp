#pragma once

#include <chrono>
#include <ostream>

#include "vec3.hpp"

namespace sgp4 {
	struct earth_coords {
		double latitude, longitude;
	};

	earth_coords from_eci_to_coords_sphere(const vec3& pos, std::chrono::utc_clock::time_point time);
	earth_coords from_eci_to_coords_ellipsoid(const vec3& pos, std::chrono::utc_clock::time_point time);
	
	std::ostream& operator<<(std::ostream& os, const earth_coords& coords);
}