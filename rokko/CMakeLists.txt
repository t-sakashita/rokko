include_directories(${ELEMENTAL_INCLUDE_DIRS})

set(ROKKO_SOURCES solver_factory.cpp)
add_library(rokko ${ROKKO_SOURCES})
