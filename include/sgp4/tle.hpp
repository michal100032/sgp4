#pragma once

#include <vector>
#include <string>
#include <string_view>
#include <chrono>

namespace sgp4 {

	struct tle_set {
		// Line 0
		std::string name;

		// Line 1
		int catalog_number;

		// <should i keep it?>
		enum class classification_type {
			unclassified, classified, secret
		};
		classification_type classification;

		struct cospar_id {
			int year;
			int launch_number;
			char piece[4];
		};
		cospar_id int_designator;
		// </should i keep it?>


		std::chrono::utc_clock::time_point epoch;

		// first derivative of mean motion [rad/min^2]
		double d_mean_motion;

		// second derivative of mean motion [rad/min^3]
		double dd_mean_motion;

		double rad_press_coef;

		int ephemeris_type;

		int set_num;
		// Line 2

		// all angles in radians
		double inclination;
		double right_ascension; // of the ascending node
		double eccentricity;

		double arg_of_perigee;
		double mean_anomaly;

		// mean motion [rad/min]
		double mean_motion;

		int rev_num;
	};

	tle_set parse_tle_entry(std::string_view str);
	std::vector<tle_set> parse_tle_entries(std::string_view str);
}