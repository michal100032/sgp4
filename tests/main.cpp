#include <iostream>
#include <cstring>
#include <vector>

#include <cassert>
#include <thread>

#include "SGP4.h"
#include "SGP4_mod.h"
#include "sgp4/sgp4.hpp"

static const double PI = 3.14159265359;
static const double RAD_TO_DEG = 180.0 / PI;

static const std::string TLE = R""""(
ISS (ZARYA)             
1 25544U 98067A   23310.23018765  .00031907  00000+0  56109-3 0  9994
2 25544  51.6419 344.9652 0001009  76.3839 283.7262 15.50187395423820
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

void assert_d3(double d1[3], double d2[3], const char* msg) {
	static const double epsilon = 1e-4;

	bool a0 = abs(d1[0] - d2[0]) < epsilon;
	bool a1 = abs(d1[1] - d2[1]) < epsilon;
	bool a2 = abs(d1[2] - d2[2]) < epsilon;

	if (!a0 || !a1 || !a2) {
		std::cerr << "ASSERTION FAILED!" << std::endl;
		std::cerr << msg << std::endl;
	}
}


int main() {
	char tle_line_one[130], tle_line_two[130];
	split_tle(TLE, tle_line_one, tle_line_two);


	auto entry = sgp4::parse_tle_entry(TLE);

	SGP4Funcs::elsetrec sgp4_rec;
	double startmfe, stopmfe, deltamin;
	SGP4Funcs::twoline2rv(tle_line_one, tle_line_two,
		'c', 0, 'i', SGP4Funcs::wgs84, startmfe, stopmfe, deltamin, sgp4_rec);
	
	SGP4Funcs_mod::elsetrec sgp4_mod_rec = 
		SGP4Funcs_mod::twoline2rv(entry, 'i');


	while (true) {
		double pos_sgp4[3], pos_sgp4_mod[3];
		double vel_sgp4[3], vel_sgp4_mod[3];

		auto now = std::chrono::utc_clock::now();

		double seconds_offset
			= std::chrono::duration<double>(now - entry.epoch).count();

		SGP4Funcs::sgp4(sgp4_rec, seconds_offset / 60.0, pos_sgp4, vel_sgp4);
		if (sgp4_rec.error != 0) {
			std::cerr << "SGP4 error " << sgp4_rec.error << std::endl;
		}
		std::cout << "SGP4: " << std::endl;
		print_coords_from_pos(pos_sgp4, now);

		SGP4Funcs_mod::sgp4(sgp4_mod_rec, seconds_offset / 60.0, pos_sgp4_mod, vel_sgp4_mod);
		if (sgp4_mod_rec.error != 0) {
			std::cerr << "SGP4_mod error " << sgp4_mod_rec.error << std::endl;
		}
		std::cout << "SGP4_mod: " << std::endl;
		print_coords_from_pos(pos_sgp4_mod, now);
		
		assert_d3(pos_sgp4, pos_sgp4_mod,
			"SGP4 and SGP4_mod yield different positions");
	
		std::this_thread::sleep_for(std::chrono::seconds(1));
		std::cout << std::endl << std::endl;
	}
}