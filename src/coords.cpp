#include "sgp4/coords.hpp"
#include "sgp4/time_utils.hpp"

#include <cmath>

static const double PI = 3.14159265359;
static const double EARTH_FLAT = 1.0 / 298.257223563;
static const double EARTH_SEMI_MAJOR = 6378.137; // km

// position in km!
sgp4::earth_coords sgp4::from_eci_to_coords_sphere(const vec3& pos, std::chrono::utc_clock::time_point time) {
	earth_coords coords;

	double pos_xy_mag = sqrt(pos.x * pos.x + pos.y * pos.y);
	double earth_rot_angle = sgp4::time_utils::to_sidereal_secs(time) / 86400 * 2 * PI;

	coords.latitude = atan(pos.z / pos_xy_mag);
	coords.longitude = (atan2(pos.y, pos.x) - earth_rot_angle);
	if (coords.longitude < -PI) {
		coords.longitude += 2 * PI;
	}

	return coords;
}

// position in km!
sgp4::earth_coords sgp4::from_eci_to_coords_ellipsoid(const vec3& pos, std::chrono::utc_clock::time_point time) {
	earth_coords coords = from_eci_to_coords_sphere(pos, time);
	double new_lat = coords.latitude;
	double e_sqrd = 2 * EARTH_FLAT - EARTH_FLAT * EARTH_FLAT;
	double pos_xy_mag = sqrt(pos.x * pos.x + pos.y * pos.y);

	// TODO determine the desired number of iterations
	for (int i = 0; i < 10; i++) {
		double sin_lat = sin(new_lat);
		double c_coef = 1.0 / sqrt(1 - e_sqrd * sin_lat * sin_lat);
		new_lat = atan((pos.z + EARTH_SEMI_MAJOR * c_coef * e_sqrd * sin_lat) / pos_xy_mag);
	}

	coords.latitude = new_lat;

	return coords;
}
