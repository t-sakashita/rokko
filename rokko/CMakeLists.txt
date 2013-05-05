include_directories(${ELEMENTAL_INCLUDE_DIRS})

set(ROKKO_SOURCES solver_factory.cpp)
#set(ROKKO_SOURCES ${ROKKO_SOURCES} elemental/elemental.cpp)
#set(ROKKO_SOURCES ${ROKKO_SOURCES} eigen_s/eigen_s.cpp)
set(ROKKO_SOURCES ${ROKKO_SOURCES} scalapack/scalapack.cpp)

add_library(rokko ${ROKKO_SOURCES})
target_link_libraries(rokko ${ELEMENTAL_LIBRARIES})
