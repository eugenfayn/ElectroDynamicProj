cmake_minimum_required(VERSION 3.10)
project(ElectroDynamicProj VERSION 1.0)

# Add compiler warnings
add_compile_options(-Wall -Wextra -Wpedantic -Werror)

# Specify the C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Optimization options
option(OPTIMIZE_O0 "No optimizations (-O0)" OFF)
option(OPTIMIZE_O1 "Basic optimizations (-O1)" OFF)
option(OPTIMIZE_O2 "Medium optimizations (-O2)" OFF)
option(OPTIMIZE_O3 "Heavy optimizations (-O3)" OFF)

# Set optimization flags based on options
if(OPTIMIZE_O0)
    add_compile_options(-O0)
    set(OPTIMIZATION_LEVEL "O0")
elseif(OPTIMIZE_O1)
    add_compile_options(-O1)
    set(OPTIMIZATION_LEVEL "O1")
elseif(OPTIMIZE_O2)
    add_compile_options(-O2)
    set(OPTIMIZATION_LEVEL "O2")
elseif(OPTIMIZE_O3)
    add_compile_options(-O3)
    set(OPTIMIZATION_LEVEL "O3")
else()
    # Default to O0 if no optimization selected
    add_compile_options(-O0)
    set(OPTIMIZATION_LEVEL "O0")
endif()

# Define project directories
set(PROJECT_ROOT_DIR "${CMAKE_CURRENT_SOURCE_DIR}")
set(PROJECT_DATA_DIR "${PROJECT_ROOT_DIR}/data")
set(PROJECT_INPUT_DIR "${PROJECT_DATA_DIR}/input")
set(PROJECT_OUTPUT_DIR "${PROJECT_DATA_DIR}/output")

# Create directories if they don't exist
file(MAKE_DIRECTORY ${PROJECT_DATA_DIR})
file(MAKE_DIRECTORY ${PROJECT_INPUT_DIR})
file(MAKE_DIRECTORY ${PROJECT_OUTPUT_DIR})

# Generate configuration header
configure_file(
    "${PROJECT_SOURCE_DIR}/include/config.h.in"
    "${PROJECT_BINARY_DIR}/include/config.h"
)

# Set output directories
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)

# Add include directories
include_directories(
    ${PROJECT_SOURCE_DIR}/include
    ${PROJECT_BINARY_DIR}/include  # For generated config.h
)

# Create library from geometry sources
add_library(geometry
    src/geometry.cpp
)

# Create executable
add_executable(integral_solver
    src/main.cpp
)

# Link the executable with the geometry library
target_link_libraries(integral_solver
    PRIVATE geometry
)

# Find Python (optional, if you need to call Python from C++)
find_package(Python3 COMPONENTS Interpreter Development REQUIRED)

# Create a custom target for installing Python dependencies
add_custom_target(python_deps
    COMMAND ${Python3_EXECUTABLE} -m pip install -r ${CMAKE_SOURCE_DIR}/scripts/requirements.txt
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
)

# Print selected optimization information
message(STATUS "")
message(STATUS "Build configuration:")
message(STATUS "  Optimization level: -${OPTIMIZATION_LEVEL}")
message(STATUS "")