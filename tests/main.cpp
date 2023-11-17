#include <iostream>
#include <cstring>
#include <vector>

#include <cassert>
#include <thread>

#include "SGP4/SGP4.h"
#include "SGP4/SGP4_mod.h"
#include "SGP4/SGP4_mod_near.h"

#include "sgp4/sgp4.hpp"

static const double PI = 3.14159265359;
static const double RAD_TO_DEG = 180.0 / PI;

static const std::string TLE = R""""(
ISS (ZARYA)             
1 25544U 98067A   23321.41134024  .00012713  00000+0  23531-3 0  9990
2 25544  51.6431 289.5970 0000839 305.0931 191.9975 15.49391010425558
)"""";

void split_tle(const std::string& tle, char* line_one, char* line_two) {
	std::istringstream tle_stream(tle);
	std::string line;

	do std::getline(tle_stream, line);
	while (!std::isdigit(line[0]));
	
	// line 1
	strcpy(line_one, line.c_str());

	// line 3
	std::getline(tle_stream, line);
	strcpy(line_two, line.c_str());
}

void print_coords_from_pos(double pos[3], std::chrono::utc_clock::time_point time) {
	sgp4::earth_coords coords =
		sgp4::from_eci_to_coords_ellipsoid({ pos[0], pos[1], pos[2] }, time);

	std::cout << "Latitude longitude: " << std::endl;
	std::cout << (coords.latitude * RAD_TO_DEG) << " " 
			  << (coords.longitude * RAD_TO_DEG) << std::endl;
}

bool assert_d3(double d1[3], double d2[3], const char* msg) {
	static const double epsilon = 1e-4;

	bool a0 = abs(d1[0] - d2[0]) < epsilon;
	bool a1 = abs(d1[1] - d2[1]) < epsilon;
	bool a2 = abs(d1[2] - d2[2]) < epsilon;

	if (!a0 || !a1 || !a2) {
		std::cerr << " --------- ASSERTION FAILED! --------- " << std::endl;
		std::cerr << msg << std::endl;
		return false;
	}
	return true;
}

bool assert_d3v3(double d[3], sgp4::vec3 vec, const char* msg) {
	static const double epsilon = 1e-4;

	bool a0 = abs(d[0] - vec.x) < epsilon;
	bool a1 = abs(d[1] - vec.y) < epsilon;
	bool a2 = abs(d[2] - vec.z) < epsilon;

	if (!a0 || !a1 || !a2) {
		std::cerr << " --------- ASSERTION FAILED! --------- " << std::endl;
		std::cerr << msg << std::endl;
		return false;
	}
	return true;
}


int main() {
	char tle_line_one[130], tle_line_two[130];
	split_tle(TLE, tle_line_one, tle_line_two);


	auto entry = sgp4::parse_tle_entry(TLE);

	SGP4Funcs::elsetrec sgp4_rec;
	double startmfe, stopmfe, deltamin;
	SGP4Funcs::twoline2rv(tle_line_one, tle_line_two,
		'c', 0, 'i', SGP4Funcs::wgs84, startmfe, stopmfe, deltamin, sgp4_rec);
	
	sgp4::propagator iss_propagator(TLE);

	while (true) {
		auto now = std::chrono::utc_clock::now();

		double seconds_offset
			= std::chrono::duration<double>(now - entry.epoch).count();

		auto [pos, vel] = iss_propagator.run(seconds_offset / 60.0);

		sgp4::earth_coords coords =
			sgp4::from_eci_to_coords_ellipsoid(pos, now);

		std::cout << "Latitude longitude: " << std::endl;
		std::cout << (coords.latitude * RAD_TO_DEG) << " "
			<< (coords.longitude * RAD_TO_DEG) << std::endl;

		

		double pos_sgp4[3];
		double vel_sgp4[3];

		SGP4Funcs::sgp4(sgp4_rec, seconds_offset / 60.0, pos_sgp4, vel_sgp4);
		if (sgp4_rec.error != 0) {
			std::cerr << "SGP4 error " << sgp4_rec.error << std::endl;
		}

		if (!assert_d3v3(pos_sgp4, pos, "Positions don't match!")) {
			std::cout << "SGP4: " << std::endl;
			print_coords_from_pos(pos_sgp4, now);
		}

		assert_d3v3(vel_sgp4, vel, "Velocities don't match!");

		std::this_thread::sleep_for(std::chrono::seconds(1));
		std::cout << std::endl << std::endl;
	}
}