cmake_minimum_required (VERSION 3.12)

project(sgp4)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED On)
set(CMAKE_CXX_EXTENSIONS Off)
# set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++20")

add_library(
    sgp4
    src/time_utils.cpp
    src/tle.cpp
    src/vec3.cpp
    src/coords.cpp
    src/propagator.cpp
)

target_include_directories(sgp4
    PUBLIC
    ${CMAKE_SOURCE_DIR}/include
)

add_subdirectory(tests)
