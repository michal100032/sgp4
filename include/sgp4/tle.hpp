#pragma once

#include <vector>
#include <string>
#include <string_view>

namespace sgp4 {

	struct tle_entry {
		// Line 0
		std::string name;

		// Line 1
		int catalog_number;

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

		int epoch_year;
		double epoch_day_frac;

		double d_mean_motion;
		double dd_mean_motion;

		double rad_press_coef;

		int ephemeris_type;

		int set_num;
		// Line 2

		double inclination;
		double right_ascension; // of the ascending node
		double eccentricity;

		double arg_of_perigee;
		double mean_anomaly;
		double mean_motion;

		int rev_num;
	};

	std::vector<tle_entry> parse_tle_entries(std::string_view str);
}