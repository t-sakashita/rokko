cmake_minimum_required(VERSION 2.8.12 FATAL_ERROR)
project(myproject NONE)

find_package(Rokko REQUIRED HINTS ${ROKKO_ROOT_DIR} $ENV{ROKKO_ROOT})
message(STATUS "Found Rokko: ${ROKKO_ROOT_DIR}")
include(${ROKKO_USE_FILE})

enable_testing()
include(CTest)

add_executable(syev syev.cpp)
target_link_libraries(syev ${ROKKO_LIBRARIES})
add_test(syev syev)

add_executable(dsyev dsyev.c)
target_link_libraries(dsyev ${ROKKO_LIBRARIES})
add_test(dsyev dsyev)

