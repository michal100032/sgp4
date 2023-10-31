#include "sgp4/tle.hpp"

#include <sstream>
#include <charconv>
static std::string_view trim_end(std::string_view str) {
	size_t last = str.find_last_not_of(" \n\r");
	return str.substr(0, last + 1);
}

template<typename T>
static inline T parse_num(std::string_view line, size_t start, size_t end) {
	while (line[start] == ' ')
		start++;
	T data;
	std::from_chars(line.data() + start, line.data() + end, data);
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

static void parse_line1(std::string_view line, sgp4::tle_entry& out_entry) {
	out_entry.catalog_number = parse_num<int>(line, 2, 7);
	out_entry.classification = line[7] == 'U'
		? sgp4::tle_entry::classification_type::unclassified
		: (line[8] == 'C'
			? sgp4::tle_entry::classification_type::classified
			: sgp4::tle_entry::classification_type::secret
		);
	
	out_entry.int_designator.year = parse_num<int>(line, 9, 11);
	out_entry.int_designator.launch_number = parse_num<int>(line, 11, 14);
	
	int i;
	for (i = 0; i < 3 && line[i + 14] != ' '; i++)
		out_entry.int_designator.piece[i] = line[i + 14];
	out_entry.int_designator.piece[i] = '\0';

	out_entry.epoch_year = parse_num<int>(line, 18, 20);
	out_entry.epoch_day_frac = parse_num<double>(line, 20, 32);

	out_entry.d_mean_motion = parse_num<double>(line, 33, 43);
	out_entry.dd_mean_motion = parse_double_pow(line, 44, 52);
	out_entry.rad_press_coef = parse_double_pow(line, 53, 61);

	out_entry.ephemeris_type = line[62] - '0';
	out_entry.set_num = parse_num<int>(line, 64, 68);
}

static void parse_line2(std::string_view line, sgp4::tle_entry& out_entry) {
	out_entry.inclination = parse_num<double>(line, 8, 16);
	out_entry.right_ascension = parse_num<double>(line, 17, 25);
	out_entry.eccentricity = parse_num<double>(line, 26, 33) / 10000000;
	out_entry.arg_of_perigee = parse_num<double>(line, 34, 42);
	out_entry.mean_anomaly = parse_num<double>(line, 43, 51);
	out_entry.mean_motion = parse_num<double>(line, 52, 63);

	out_entry.rev_num = parse_num<int>(line, 63, 68);
}

std::vector<sgp4::tle_entry> sgp4::parse_tle_entries(std::string_view str) {
	std::vector<tle_entry> entries = std::vector<tle_entry>();
	std::istringstream stream(str.data());

	std::string line;
	tle_entry current_entry;

	while (std::getline(stream, line)) {
		if (std::isdigit(line[0])) {
			int entry_line = line[0] - '0';

			if (entry_line == 1) {
				parse_line1(line, current_entry);
			} else {
				parse_line2(line, current_entry);
				entries.push_back(current_entry);
				current_entry = tle_entry();
			}
		} else {
			current_entry.name = trim_end(line);
		}
	}

	return entries;
}