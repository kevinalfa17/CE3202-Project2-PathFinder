cmake_minimum_required(VERSION 3.5)

project(PathFinder VERSION 1.0.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 11)
set(COMMON_SOURCES)
ENABLE_LANGUAGE(C)


#FIND AND INCLUDE
find_package(OpenCV REQUIRED )
find_package(GTest REQUIRED)

INCLUDE (CheckIncludeFiles)
include_directories(${GTEST_INCLUDE_DIRS})


# SIMD OPTION
option(IS_SIMD_ACTIVE "Activate SIMD" ON)

#VERIFY IF SIMD IS AVAILABLE
CHECK_INCLUDE_FILES (xmmintrin.h IS_SIMD_AVAILABLE)
CONFIGURE_FILE(${CMAKE_CURRENT_SOURCE_DIR}/config.h.in ${CMAKE_CURRENT_BINARY_DIR}/config.h)



#EXECUTABLES
add_executable(pathFinder main.cpp)
add_executable(runTests test.cpp)

#LINK LIBRARIES
target_link_libraries(runTests ${GTEST_LIBRARIES} pthread)
target_link_libraries( pathFinder ${OpenCV_LIBS} python2.7)


