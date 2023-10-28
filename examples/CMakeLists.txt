cmake_minimum_required(VERSION 3.12)

project(sgp4_examples)
add_executable(sgp4_examples main.cpp)

target_include_directories(sgp4_examples PRIVATE ${CMAKE_SOURCE_DIR}/include)
target_link_libraries(sgp4_examples sgp4)
