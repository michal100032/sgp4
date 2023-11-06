#include <iostream>
#include <cstring>
#include <vector>

#include <thread>

#include "SGP4_mod.h"
#include "sgp4/sgp4.hpp"

const double PI = 3.14159265359;
const double RAD_TO_DEG = 180 / PI;

const char* TLE = R""""(
ISS (ZARYA)             
1 25544U 98067A   23310.23018765  .00031907  00000+0  56109-3 0  9994
2 25544  51.6419 344.9652 0001009  76.3839 283.7262 15.50187395423820
)"""";

int main() {
	double startmfe, stopmfe, deltamin;
	std::cout << "NOAA 15 position: " << std::endl;
	elsetrec record;
	auto entry = sgp4::parse_tle_entry(TLE);
	SGP4Funcs::twoline2rv(entry,
		/* opsmode: improved */   'i',
		/* whichconst:       */   wgs84,
		/* returned elserec  */   record
	);
	double pos[3];
	double vel[3];

	while (true) {
		// system("cls");

		auto now = std::chrono::utc_clock::now();

		double seconds_offset
			= std::chrono::duration<double>(now - entry.epoch).count();

		SGP4Funcs::sgp4(record, seconds_offset / 60.0, pos, vel);
		if (record.error != 0) {
			std::cerr << "ERROR " << record.error << std::endl;
		}

		sgp4::earth_coords coords = sgp4::from_eci_to_coords_ellipsoid({ pos[0], pos[1], pos[2] }, now);

		std::cout << "Latitude longitude: " << std::endl;
		std::cout << (coords.latitude * RAD_TO_DEG) << " " << (coords.longitude * RAD_TO_DEG) << std::endl;

		using namespace std::chrono_literals;
		std::this_thread::sleep_for(1s);
	}
}