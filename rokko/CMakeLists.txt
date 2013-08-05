include_directories(${ELEMENTAL_INCLUDE_DIRS})

set(ROKKO_SOURCES solver_factory.cpp grid.cpp)

set(ROKKO_SOURCES ${ROKKO_SOURCES} blacs/blacs_barrier.c blacs/blacs_exit.c blacs/blacs_get.c blacs/blacs_gridexit.c blacs/blacs_gridinfo.c blacs/blacs_gridinit.c blacs/blacs_pinfo.c blacs/descinit.c)
set(ROKKO_SOURCES ${ROKKO_SOURCES} pblas/pdgemm.c)

set(DIRECTORY "${BOOST_INCLUDE_DIR}/libs/system")
if (EXISTS ${DIRECTORY})
  set(SOURCES error_code.cpp)
  foreach (D ${DIRECTORY}/src)
    foreach(S ${SOURCES})
      if(EXISTS ${D}/${S})
        set(ROKKO_SOURCES ${ROKKO_SOURCES} ${D}/${S})
      endif(EXISTS ${D}/${S})
    endforeach(S)
  endforeach(D)
else (EXISTS ${DIRECTORY})
  find_package(Boost 1.48.0 COMPONENTS system)
endif (EXISTS ${DIRECTORY})

if(BUILD_EIGEN_S)
  set(ROKKO_SOURCES ${ROKKO_SOURCES} eigen_s/eigen_s.cpp)
endif(BUILD_EIGEN_S)

if(BUILD_EIGEN_SX)
  set(ROKKO_SOURCES ${ROKKO_SOURCES} eigen_sx/eigen_sx.cpp)
endif(BUILD_EIGEN_SX)

if(BUILD_EIGEN_EXA)
  set(ROKKO_SOURCES ${ROKKO_SOURCES} eigen_exa/eigen_exa.cpp)
endif(BUILD_EIGEN_EXA)

if(BUILD_ELPA)
  set(ROKKO_SOURCES ${ROKKO_SOURCES} elpa/elpa.cpp)
endif(BUILD_ELPA)

if( DEFINED ELEMENTAL_LIBRARIES AND ELEMENTAL_INCLUDE_DIRS )
  set(ROKKO_SOURCES ${ROKKO_SOURCES} elemental/elemental.cpp)
endif( DEFINED ELEMENTAL_LIBRARIES AND ELEMENTAL_INCLUDE_DIRS )

if ( DEFINED SCALAPACK_LIB )
  set(ROKKO_SOURCES ${ROKKO_SOURCES} scalapack/pdelget.c scalapack/pdelset.c scalapack/pdlamch.c scalapack/pdlaprnt.c scalapack/pdsyev.c scalapack/pdsyevd.c scalapack/pdsyevx.c)
  if ( USE_PDSYEVR )
    set(ROKKO_SOURCES ${ROKKO_SOURCES} scalapack/pdsyevr.c)
  endif ( USE_PDSYEVR )
  set(ROKKO_SOURCES ${ROKKO_SOURCES} scalapack/scalapack.cpp)
  set(ROKKO_SOURCES ${ROKKO_SOURCES} scalapack/scalapack_pdsyevd.cpp)
#  set(ROKKO_SOURCES ${ROKKO_SOURCES} scalapack/scalapack_pdsyevx.cpp)
endif ( DEFINED SCALAPACK_LIB )

add_library(rokko ${ROKKO_SOURCES})
target_link_libraries(rokko ${ELEMENTAL_LIBRARIES})
if(BUILD_LAPACKE)
  target_link_libraries(rokko lapacke)
endif(BUILD_LAPACKE)
if(BUILD_EIGEN_S)
  target_link_libraries(rokko eigen_s)
endif(BUILD_EIGEN_S)
if(BUILD_EIGEN_SX)
  target_link_libraries(rokko eigen_sx)
endif(BUILD_EIGEN_SX)
if(BUILD_EIGEN_EXA)
  target_link_libraries(rokko ${EIGEN_EXA_LIB})
endif(BUILD_EIGEN_EXA)
if(BUILD_ELPA)
  target_link_libraries(rokko ${ELPA_LIB})
endif(BUILD_ELPA)

if (Boost_FOUND)
  target_link_libraries(rokko ${Boost_LIBRARIES})
endif (Boost_FOUND)

install(TARGETS rokko COMPONENT libraries
        ARCHIVE DESTINATION lib
        LIBRARY DESTINATION lib
        RUNTIME DESTINATION bin)
install(DIRECTORY . COMPONENT headers DESTINATION include/rokko
        FILES_MATCHING PATTERN "*.h" PATTERN "*.hpp" PATTERN "*.ipp"
        PATTERN ".svn" EXCLUDE PATTERN "CVS" EXCLUDE)
