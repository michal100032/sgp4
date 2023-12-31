#pragma once

#include <chrono>

namespace sgp4 {
	namespace time_utils {
		inline std::tm to_tm(std::chrono::utc_clock::time_point time);
		std::chrono::utc_clock::time_point to_utc(const tm& time);
		std::chrono::utc_clock::duration get_utc_offset();

		double to_julian(std::chrono::utc_clock::time_point time);
		double to_julian_from_tm(const tm& utc_tm);

		double to_sidereal_secs(std::chrono::utc_clock::time_point time);
	}
}