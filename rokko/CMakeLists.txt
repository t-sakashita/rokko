include_directories(${ELEMENTAL_INCLUDE_DIRS})

set(ROKKO_SOURCES solver_factory.cpp)

set(DIRECTORY "${BOOST_INCLUDE_DIR}/libs/system")
set(SOURCES error_code.cpp)
foreach (D ${DIRECTORY}/src)
  foreach(S ${SOURCES})
    if(EXISTS ${D}/${S})
      set(ROKKO_SOURCES ${ROKKO_SOURCES} ${D}/${S})
    endif(EXISTS ${D}/${S})
  endforeach(S)
endforeach(D)

if(BUILD_EIGEN_S)
  set(ROKKO_SOURCES ${ROKKO_SOURCES} eigen_s/eigen_s.cpp)
endif(BUILD_EIGEN_S)
  
if(BUILD_EIGEN_SX)
  set(ROKKO_SOURCES ${ROKKO_SOURCES} eigen_sx/eigen_sx.cpp)
endif(BUILD_EIGEN_SX)

if( DEFINED ELEMENTAL_LIBRARIES AND ELEMENTAL_INCLUDE_DIRS )
  set(ROKKO_SOURCES ${ROKKO_SOURCES} elemental/elemental.cpp)
endif( DEFINED ELEMENTAL_LIBRARIES AND ELEMENTAL_INCLUDE_DIRS )

if ( DEFINED SCALAPACK_LIB )
  set(ROKKO_SOURCES ${ROKKO_SOURCES} scalapack/scalapack.cpp)
endif ( DEFINED SCALAPACK_LIB )

add_library(rokko ${ROKKO_SOURCES})
if(BUILD_EIGEN_S)
  target_link_libraries(rokko eigen_s)
endif(BUILD_EIGEN_S)
if(BUILD_EIGEN_SX)
  target_link_libraries(rokko eigen_sx)
endif(BUILD_EIGEN_SX)
