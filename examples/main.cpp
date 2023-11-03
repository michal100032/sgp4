#include <iostream>
#include <chrono>
#include <thread>

#include "sgp4/sgp4.hpp"

const char* TLE = R""""(
ATLAS CENTAUR 2
1 00694U 63047A   23306.67866552  .00001938  00000+0  23728-3 0  9993
2 00694  30.3548 289.9788 0571435 111.2764 254.9806 14.06138642 10495
)"""";

/*
ATLAS CENTAUR 2
1 00694U 63047A   23306.67866552  .00001938  00000+0  23728-3 0  9993
2 00694  30.3548 289.9788 0571435 111.2764 254.9806 14.06138642 10495
THOR AGENA D R/B
1 00733U 64002A   23306.86105797  .00000505  00000+0  20134-3 0  9994
2 00733  99.0420 260.9537 0032554 235.3063 124.5051 14.32942425115395
*/

int main() {
	using namespace std::chrono_literals;

	while (true) {
		system("cls");
		auto now = std::chrono::utc_clock::now();
		std::cout << "UTC" << std::endl;
		std::cout << now << std::endl;;
		int sidereal_secs = sgp4::time_utils::to_sidereal_secs(now);
		std::cout << "Sidereal" << std::endl;
		std::cout << (sidereal_secs / 3600) << ":"
			<< (sidereal_secs % 3600 / 60) << ":"
			<< (sidereal_secs % 3600 % 60);

		std::this_thread::sleep_for(1s);
	}
}