#include "sgp4/tle.hpp"

#include <sstream>
#include <charconv>

static const double PI = 3.14159265359;
static const double RAD_TO_DEG = 180.0 / PI;
static const double DEG_TO_RAD = PI / 180.0;

static std::string_view trim_end(std::string_view str) {
	size_t last = str.find_last_not_of(" \n\r");
	return str.substr(0, last + 1);
}

template<typename T>
static inline T parse_num(std::string_view line, size_t start, size_t end) {
	while (line[start] == ' ')
		start++;
	T data;
	std::from_chars_result res = 
		std::from_chars(line.data() + start, line.data() + end, data);
	if (res.ec != std::errc())
		return T();
	return data;
}

static double parse_double_pow(std::string_view line, size_t start, size_t end) {
	size_t plus_min_pos = line.find_first_of("+-", start);
	if (plus_min_pos == std::string_view::npos || plus_min_pos >= end - 1)
		return parse_num<double>(line, start, end);
	double without_pow = parse_num<double>(line, start, plus_min_pos);

	double power = parse_num<double>(line, plus_min_pos, end)
		+ start - plus_min_pos + 1;

	return without_pow * pow(10.0, power);
}

static std::chrono::utc_clock::time_point tle_epoch_to_utc(int ep_year, double ep_day_frac) {
	int actual_year = ep_year < 57 ? ep_year + 2000 : ep_year + 1900;
	int day_of_year = static_cast<int>(ep_day_frac);

	double day_frac = ep_day_frac - day_of_year;
	int milisecs = static_cast<int>(86400000.0 * day_frac);

	using namespace std::chrono;

	system_clock::time_point system_time = sys_days{ year(actual_year) / month(1) / day(0) }
	+ hours{ 24 * day_of_year } + milliseconds{ milisecs };
	
	utc_clock::time_point utc = utc_clock::from_sys(system_time);
	
	return utc;
}

static void parse_line1(std::string_view line, sgp4::tle_set& out_entry) {
	out_entry.catalog_number = parse_num<int>(line, 2, 7);
	out_entry.classification = line[7] == 'U'
		? sgp4::tle_set::classification_type::unclassified
		: (line[8] == 'C'
			? sgp4::tle_set::classification_type::classified
			: sgp4::tle_set::classification_type::secret
		);
	
	out_entry.int_designator.year = parse_num<int>(line, 9, 11);
	out_entry.int_designator.launch_number = parse_num<int>(line, 11, 14);
	
	int i;
	for (i = 0; i < 3 && line[i + 14] != ' '; i++)
		out_entry.int_designator.piece[i] = line[i + 14];
	out_entry.int_designator.piece[i] = '\0';

	int epoch_year = parse_num<int>(line, 18, 20);
	double epoch_day_frac = parse_num<double>(line, 20, 32);
	out_entry.epoch = tle_epoch_to_utc(epoch_year, epoch_day_frac);

	out_entry.d_mean_motion = parse_num<double>(line, 33, 43) * 2 * PI / 1440.0 / 1444.0;
	out_entry.dd_mean_motion = parse_double_pow(line, 44, 52) * 2 * PI / 1440.0 / 1440.0 / 1440.0;
	out_entry.rad_press_coef = parse_double_pow(line, 53, 61);

	out_entry.ephemeris_type = line[62] - '0';
	out_entry.set_num = parse_num<int>(line, 64, 68);
}

static void parse_line2(std::string_view line, sgp4::tle_set& out_entry) {
	out_entry.inclination = parse_num<double>(line, 8, 16) * DEG_TO_RAD;
	out_entry.right_ascension = parse_num<double>(line, 17, 25) * DEG_TO_RAD;
	out_entry.eccentricity = parse_num<double>(line, 26, 33) / 10000000;
	out_entry.arg_of_perigee = parse_num<double>(line, 34, 42) * DEG_TO_RAD;
	out_entry.mean_anomaly = parse_num<double>(line, 43, 51) * DEG_TO_RAD;
	out_entry.mean_motion = parse_num<double>(line, 52, 63) * 2 * PI / 1440.0;

	out_entry.rev_num = parse_num<int>(line, 63, 68);
}

sgp4::tle_set sgp4::parse_tle_entry(std::string_view str) {
	tle_set entry;
	std::istringstream stream(str.data());

	std::string line;
	while (std::getline(stream, line)) {
		if (std::isdigit(line[0])) {
			int entry_line = line[0] - '0';

			if (entry_line == 1) {
				parse_line1(line, entry);
			}
			else {
				parse_line2(line, entry);
				return entry;
			}
		}
		else {
			entry.name = trim_end(line);
		}
	}

	return entry;
}

std::vector<sgp4::tle_set> sgp4::parse_tle_entries(std::string_view str) {
	std::vector<tle_set> entries = std::vector<tle_set>();
	std::istringstream stream(str.data());

	std::string line;
	tle_set current_entry;

	while (std::getline(stream, line)) {
		if (std::isdigit(line[0])) {
			int entry_line = line[0] - '0';

			if (entry_line == 1) {
				parse_line1(line, current_entry);
			} else {
				parse_line2(line, current_entry);
				entries.push_back(current_entry);
				current_entry = tle_set();
			}
		} else {
			current_entry.name = trim_end(line);
		}
	}

	return entries;
}