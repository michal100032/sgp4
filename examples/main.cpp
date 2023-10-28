#include <iostream>
#include <chrono>

#include "sgp4/sgp4.hpp"

int main() {
	std::cout << "Hello there!" << std::endl;
	auto utc_now = std::chrono::utc_clock::now();
	std::cout << std::setprecision(10) << sgp4::time_utils::to_julian(utc_now) << std::endl;
	double sidereal = sgp4::time_utils::to_sidereal_secs(utc_now);
	std::cout << ((int)sidereal / 3600) << std::endl;
	std::cout << ((int)sidereal % 3600 / 60) << std::endl;
	std::cout << ((int)sidereal % 3600 % 60) << std::endl;
	std::cout << sidereal;
}