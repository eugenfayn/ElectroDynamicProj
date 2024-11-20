cmake_minimum_required(VERSION 3.20)
project(name cxx)

set(CMAKE_CXX_STANDART 11)
set(CMAKE_CXX_FLAGS "-O3")

find_package(LAPACK REQUIRED)
find_package(BLAS REQUIRED)
find_package(OPENMP REQUIRED)

add_executable(main file1.cpp file2.cpp main.cpp)
target_link_libraries(main PUBLIC LAPACK::LAPACK)
target_link_libraries(main PUBLIC BLAS::BLAS)
target_link_libraries(main PUBLIC OPENMP::OPENMP)