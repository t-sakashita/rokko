set(ROKKO_SOURCES serial_solver_factory.cpp solver_factory.cpp grid.cpp)

# Eigen3
set(ROKKO_SOURCES ${ROKKO_SOURCES} eigen3/eigen3.cpp)

# LAPACK
set(ROKKO_SOURCES ${ROKKO_SOURCES} lapack/lapack.cpp)

# BLACS
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

if(BUILD_SCALAPACK)
  set(ROKKO_SOURCES ${ROKKO_SOURCES} scalapack/pdelget.c scalapack/pdelset.c scalapack/pdlamch.c scalapack/pdlaprnt.c scalapack/pdsyev.c scalapack/pdsyevd.c scalapack/pdsyevx.c)
  if(USE_PDSYEVR)
    set(ROKKO_SOURCES ${ROKKO_SOURCES} scalapack/pdsyevr.c)
  endif(USE_PDSYEVR)
  set(ROKKO_SOURCES ${ROKKO_SOURCES} scalapack/scalapack.cpp)
  set(ROKKO_SOURCES ${ROKKO_SOURCES} scalapack/scalapack_pdsyevd.cpp)
  # set(ROKKO_SOURCES ${ROKKO_SOURCES} scalapack/scalapack_pdsyevx.cpp)
endif(BUILD_SCALAPACK)

if(BUILD_EIGEN_EXA)
  include_directories(${EIGENEXA_INCLUDE_DIR})
  set(ROKKO_SOURCES ${ROKKO_SOURCES} eigen_exa/eigen_exa.cpp)
endif(BUILD_EIGEN_EXA)

if(BUILD_ELEMENTAL)
  include_directories(${ELEMENTAL_INCLUDE_DIR})
  set(ROKKO_SOURCES ${ROKKO_SOURCES} elemental/elemental.cpp)
endif(BUILD_ELEMENTAL)

if(BUILD_ELPA)
  set(ROKKO_SOURCES ${ROKKO_SOURCES} elpa/elpa.cpp)
endif(BUILD_ELPA)

add_library(rokko ${ROKKO_SOURCES})
target_link_libraries(rokko lapacke)
if(BUILD_SCALAPACK)
  target_link_libraries(rokko ${SCALAPACK_LIBRARIES})
endif(BUILD_SCALAPACK)
if(BUILD_EIGENEXA)
  target_link_libraries(rokko ${EIGENEXA_LIBRARIES})
endif(BUILD_EIGENEXA)
if(BUILD_ELEMENTAL)
  target_link_libraries(rokko ${ELEMENTAL_LIBRARIES})
endif(BUILD_ELEMENTAL)
if(BUILD_ELPA)
  target_link_libraries(rokko ${ELPA_LIBRARIES})
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
