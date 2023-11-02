#include "sgp4/time_utils.hpp"

#include <iostream>

// rounds seconds down!
std::tm sgp4::time_utils::to_tm(std::chrono::utc_clock::time_point time) {
	using namespace std::chrono;

	system_clock::time_point sys_time = utc_clock::to_sys(time);
	time_t sys_time_t = system_clock::to_time_t(sys_time);

	tm utc_tm = *std::gmtime(&sys_time_t);

	return utc_tm;
}

std::chrono::utc_clock::time_point sgp4::time_utils::to_utc(const tm& time) {
	using namespace std::chrono;

	std::time_t local_time = std::mktime(const_cast<std::tm*>(&time));
	system_clock::time_point sys_time = system_clock::from_time_t(local_time);

	sys_time += get_utc_offset();
	
	utc_clock::time_point utc = utc_clock::from_sys(sys_time);

	return utc;
}

std::chrono::utc_clock::duration sgp4::time_utils::get_utc_offset() {
	// cache utc offset to avoid doing the same calculations multiple times
	static bool offset_calculated = false;
	static std::chrono::utc_clock::duration utc_offset;

	if (!offset_calculated) {
		std::time_t now = std::time(0);
		std::tm t = *gmtime(&now);

		std::time_t local_time = std::mktime(const_cast<std::tm*>(&t));
		std::time_t utc_time = std::mktime(std::gmtime(&local_time));

		utc_offset = std::chrono::seconds(local_time - utc_time);
		offset_calculated = true;
	}

	return utc_offset;
}

double sgp4::time_utils::to_julian(std::chrono::utc_clock::time_point time) {
	std::tm utc_tm = to_tm(time);
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
// sidereal time from 0 to 84600!
// 1 sidereal second = 86400 / 86164.0905 seconds = 1.00273790971 seconds
double sgp4::time_utils::to_sidereal_secs(std::chrono::utc_clock::time_point time) {
	double julian = to_julian(time);
	
	// julian fraction since midnight (<day>.5)
	double julian_since_midnight = julian + 0.5f - truncf(julian + 0.5f);
	double julian_until_midnight = julian - julian_since_midnight;
	
	double julian_since_2000 = julian_until_midnight - 2451545.0;

	// from the formula
	double julian_centuries_since_2000 = julian_since_2000 / 36525;

	// GMST at UT1 = 0
	double midnight_sidereal_secs = 
		  24110.54841 
		+ 8640184.812866 * julian_centuries_since_2000
		+ 0.093104       * julian_centuries_since_2000 * julian_centuries_since_2000
		- 6.2e-6         * julian_centuries_since_2000 * julian_centuries_since_2000 * julian_centuries_since_2000;
	
	double sidereal_secs_since_midnight = 84600 * 1.00273790971 * julian_since_midnight;

	double sidereal_secs =
		fmod(midnight_sidereal_secs + sidereal_secs_since_midnight, 86400);
		
	return sidereal_secs;
}
