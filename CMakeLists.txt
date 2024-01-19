# Minimum version of CMake required to build this project
cmake_minimum_required(VERSION 3.4.1)

# Name of the project
project(EIT)

# Enable C++11 features
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_FLAGS "-Wall -Wextra -m64 -O3 -g -std=c++11 -fprofile-arcs -ftest-coverage")

set(RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/build)
set(EXECUTABLE_OUTPUT_PATH ${CMAKE_SOURCE_DIR}/build)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/libs)

set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)
# Find MPI package
find_package(MPI REQUIRED)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake_modules/morse_cmake/modules" )

# include(MorseInit)

# Add include directories for MPI
include_directories(${MPI_INCLUDE_PATH})

include_directories (${CMAKE_SOURCE_DIR}/src)

#find_package(Eigen3 REQUIRED)
#include(${EIGEN3_USE_FILE})
#include_directories(${EIGEN3_INCLUDE_DIR})
enable_testing()

include_directories (${CMAKE_SOURCE_DIR}/src)
#include_directories (${CMAKE_SOURCE_DIR}/tst)

# Add subdirectory for source code
add_subdirectory (${CMAKE_SOURCE_DIR}/src)
add_subdirectory (${CMAKE_SOURCE_DIR}/tst)


set(CMAKE_BUILD_TYPE "Debug")
set(CMAKE_EXE_LINKER_FLAGS "-fprofile-arcs -ftest-coverage")

add_custom_target(coverage
    COMMAND ${CMAKE_COMMAND} -E chdir ${CMAKE_BINARY_DIR} lcov --capture --directory . --output-file coverage.info
    COMMAND ${CMAKE_COMMAND} -E chdir ${CMAKE_BINARY_DIR} lcov --remove coverage.info '/usr/*' '*googletest*' '*test/*' --output-file coverage.info.cleaned
    COMMAND ${CMAKE_COMMAND} -E chdir ${CMAKE_BINARY_DIR} genhtml coverage.info.cleaned --output-directory coverage-report
    DEPENDS tst
)
