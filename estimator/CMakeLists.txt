cmake_minimum_required(VERSION 3.0 FATAL_ERROR)
project(estimator)

SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -Wall")
SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wall")

# Eigen library
find_package(Eigen3 3.3 REQUIRED NO_MODULE)

include_directories(
    ./
    src/
    test/
    )

add_executable(
    estimator_test
    test/estimator_test.cpp
    src/model.cpp
    src/MEKF.cpp
    codegen/integrator.c
    codegen/estimator_jacobian.c
    codegen/propagation_matrix.c
    codegen/output_map.c
    codegen/output_map_jacobian.c
    )

target_link_libraries(
    estimator_test
    Eigen3::Eigen
    )
