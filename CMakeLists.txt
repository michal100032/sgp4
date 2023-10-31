cmake_minimum_required (VERSION 3.12)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++20")

project (sgp4)
add_library(sgp4
    src/time_utils.cpp
 "src/tle.cpp" "src/propagator.cpp" "src/vec3.cpp")

target_include_directories(sgp4
    PUBLIC
    ${CMAKE_SOURCE_DIR}/include
)

add_subdirectory(examples)