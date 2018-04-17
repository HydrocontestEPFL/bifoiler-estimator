cmake_minimum_required(VERSION 2.8.3)
project(bifoiler_state_estimator)


# Compiler flags
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -Wall")
SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wall")

# Libraries
find_package(yaml-cpp 0.6 REQUIRED)

add_executable(
    config_test
    test/config_test.cpp
    )

target_link_libraries(
    config_test
    yaml-cpp
    )
