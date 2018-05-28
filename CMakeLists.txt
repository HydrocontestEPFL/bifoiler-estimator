cmake_minimum_required(VERSION 2.8.3)
project(bifoiler_state_estimator)

# add local cmake/ directory to module path
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH};${PROJECT_SOURCE_DIR}/cmake)

# Compiler flags
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -Wall")
SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wall")

# Libraries
find_package(casadi REQUIRED)
find_package(yaml-cpp 0.6 REQUIRED)

# Eigen library
find_package(Eigen3 3.3 REQUIRED NO_MODULE)


link_directories(
    $ENV{CASADI_PREFIX}
    )

include_directories(
    ./
    src/
    ${CASADI_INCLUDE_DIR}
    )

add_executable(
    casadi_test
    test/casadi_test.cpp
    )

target_link_libraries(
    casadi_test
    casadi
    )

add_executable(
    config_test
    test/config_test.cpp
    )

target_link_libraries(
    config_test
    yaml-cpp
    )

add_executable(
    boat_model_test
    test/boat_model_test.cpp
    src/boat_model.cpp
    src/boat_properties.cpp
    )

target_link_libraries(
    boat_model_test
    casadi
    yaml-cpp
    )

add_executable(
    codegen_test
    test/codegen_test.cpp
    test/codegen_dynamics.c
    )

target_link_libraries(
    codegen_test
    Eigen3::Eigen
    )
