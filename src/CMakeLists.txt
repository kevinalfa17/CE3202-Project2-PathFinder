cmake_minimum_required(VERSION 2.6)
set (CMAKE_CXX_STANDARD 11)

INCLUDE (CheckIncludeFiles)

option(IS_SIMD_ACTIVE "Activate SIMD" ON)

CHECK_INCLUDE_FILES (xmmintrin.h IS_SIMD_AVAILABLE)
CONFIGURE_FILE(${CMAKE_CURRENT_SOURCE_DIR}/config.h.in ${CMAKE_CURRENT_BINARY_DIR}/config.h)
add_executable(pathFinder CE3202-Project2-PathFinder.cpp)
