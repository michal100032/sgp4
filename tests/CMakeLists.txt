cmake_minimum_required(VERSION 3.12)

project(sgp4_test)
add_executable(
	sgp4_test 
	main.cpp
	"SGP4/SGP4_mod.cpp"
	"SGP4/SGP4.cpp"
 "SGP4/SGP4_mod_ds.cpp" "SGP4/SGP4_mod_ds.h" "SGP4/SGP4_mod_near.cpp" "SGP4/SGP4_mod_near.h")

target_include_directories(sgp4_test PRIVATE ${CMAKE_SOURCE_DIR}/include)
target_link_libraries(sgp4_test sgp4)
