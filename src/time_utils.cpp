#include "sgp4/time_utils.hpp"

#include <iostream>

std::tm sgp4::time_utils::to_tm(std::chrono::utc_clock::time_point time) {
	using namespace std::chrono;

	system_clock::time_point sys_time = utc_clock::to_sys(time);
	time_t sys_time_t = system_clock::to_time_t(sys_time);

	tm utc_tm = *gmtime(&sys_time_t);

	return utc_tm;
}

double sgp4::time_utils::to_julian(std::chrono::utc_clock::time_point time) {
	tm utc_tm = to_tm(time);
	return to_julian_from_tm(utc_tm);
}

inline double sgp4::time_utils::to_julian_from_tm(const tm& utc_tm) {
	int day_number = 
		(1461 * (utc_tm.tm_year + 6700 + (utc_tm.tm_mon - 13) / 12)) / 4 
		+ (367 * (utc_tm.tm_mon - 1 - 12 * ((utc_tm.tm_mon - 13) / 12))) / 12 
		- (3 * ((utc_tm.tm_year + 6700 + (utc_tm.tm_mon - 13) / 12) / 100)) / 4 
		+ utc_tm.tm_mday - 32075;
	double julian = day_number 
		+ ((double)utc_tm.tm_hour - 12.0) / 24 
		+ (double)utc_tm.tm_min / 1440 
		+ (double)utc_tm.tm_sec / 86400;

	return julian;
}

double sgp4::time_utils::to_sidereal_secs(std::chrono::utc_clock::time_point time) {
	double julian = to_julian(time);
	
	// julian fraction since midnight (<day>.5)
	double julian_since_midnight = julian + 0.5f - truncf(julian + 0.5f);
	double julian_until_midnight = julian - julian_since_midnight;
	
	double julian_since_2000 = julian_until_midnight - 2451545.0;

	// from the formula
	double year_fraction = julian_since_2000 / 36525;
	double midnight_sidereal_secs = 
		  24110.54841 
		+ 8640184.812866 * year_fraction
		+ 0.093104       * year_fraction * year_fraction
		- 6.2e-6         * year_fraction * year_fraction * year_fraction;
	
	// Why?????
	double sidereal_secs =
		fmod(midnight_sidereal_secs + 86636.555367 * julian_since_midnight, 86400);

	return sidereal_secs;
}
