include_directories(${ELEMENTAL_INCLUDE_DIRS})

set(ROKKO_SOURCES solver_factory.cpp)
set(ROKKO_SOURCES ${ROKKO_SOURCES} elemental/elemental.cpp)

add_library(rokko ${ROKKO_SOURCES})
target_link_libraries(rokko ${ELEMENTAL_LIBRARIES})
