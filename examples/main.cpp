#include <iostream>
#include <chrono>

#include "sgp4/sgp4.hpp"

const char* tle_test = R""""(
NOAA 1 [-]              
ATLAS CENTAUR 2         
1 00694U 63047A   23303.69695459  .00002227  00000+0  27508-3 0  9990
2 00694  30.3548 306.4872 0571517  85.2415 281.3177 14.06127117 10072
THOR AGENA D R/B        
1 00733U 64002A   23303.85853536  .00000536  00000+0  21239-3 0  9992
2 00733  99.0413 257.8209 0032468 244.3071 115.4758 14.32938990114968
)"""";

int main() {
	std::cout << tle_test << std::endl;

	auto now = std::chrono::utc_clock::now();
	tm now_tm = sgp4::time_utils::to_tm(now);
	auto now2 = sgp4::time_utils::from_tm(now_tm);

	std::cout << now << std::endl;
	std::cout << now2 << std::endl;

}