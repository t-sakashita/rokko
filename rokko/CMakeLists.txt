set(ROKKO_SOURCES serial_solver_factory.cpp solver_factory.cpp grid.cpp)

# Eigen3
set(ROKKO_SOURCES ${ROKKO_SOURCES} eigen3/eigen3.cpp)

# LAPACK
set(ROKKO_SOURCES ${ROKKO_SOURCES} lapack/lapack.cpp)

# BLACS
set(ROKKO_SOURCES ${ROKKO_SOURCES} blacs/blacs_barrier.c blacs/blacs_exit.c blacs/blacs_get.c blacs/blacs_gridexit.c blacs/blacs_gridinfo.c blacs/blacs_gridinit.c blacs/blacs_pinfo.c blacs/descinit.c)
set(ROKKO_SOURCES ${ROKKO_SOURCES} pblas/pdgemm.c)

if(Boost_SORUCE_DIR)
  set(SOURCES error_code.cpp)
  foreach (D ${Boost_SOURCE_DIR}/libs/system/src)
    foreach(S ${SOURCES})
      if(EXISTS ${D}/${S})
        set(ROKKO_SOURCES ${ROKKO_SOURCES} ${D}/${S})
      endif(EXISTS ${D}/${S})
    endforeach(S)
  endforeach(D)
endif(Boost_SORUCE_DIR)

if(BUILD_SCALAPACK)
  set(ROKKO_SOURCES ${ROKKO_SOURCES} scalapack/pdelget.c scalapack/pdelset.c scalapack/pdlamch.c scalapack/pdlaprnt.c scalapack/pdsyev.c scalapack/pdsyevd.c scalapack/pdsyevx.c)
  if(USE_PDSYEVR)
    set(ROKKO_SOURCES ${ROKKO_SOURCES} scalapack/pdsyevr.c)
  endif(USE_PDSYEVR)
  set(ROKKO_SOURCES ${ROKKO_SOURCES} scalapack/scalapack.cpp)
  set(ROKKO_SOURCES ${ROKKO_SOURCES} scalapack/scalapack_pdsyevd.cpp)
  # set(ROKKO_SOURCES ${ROKKO_SOURCES} scalapack/scalapack_pdsyevx.cpp)
endif(BUILD_SCALAPACK)

if(BUILD_EIGENEXA)
  include_directories(${EIGENEXA_INCLUDE_DIR})
  set(ROKKO_SOURCES ${ROKKO_SOURCES} eigen_exa/eigen_exa.cpp eigen_exa/eigen_exa_wrap.f90 eigen_exa/eigen_init.c eigen_exa/eigen_free.c eigen_exa/cstab_get_optdim.c eigen_exa/eigen_sx.c)
endif(BUILD_EIGENEXA)

if(BUILD_ELEMENTAL)
  include_directories(${ELEMENTAL_INCLUDE_DIR})
  set(ROKKO_SOURCES ${ROKKO_SOURCES} elemental/elemental.cpp)
endif(BUILD_ELEMENTAL)

if(BUILD_ELPA)
  set(ROKKO_SOURCES ${ROKKO_SOURCES} elpa/elpa.cpp)
endif(BUILD_ELPA)

add_library(rokko ${ROKKO_SOURCES})
target_link_libraries(rokko lapacke ${MPI_LIBRARY})
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

if(Boost_LIBRARIES)
  target_link_libraries(rokko ${Boost_LIBRARIES})
endif(Boost_LIBRARIES)

install(TARGETS rokko COMPONENT libraries
        ARCHIVE DESTINATION lib
        LIBRARY DESTINATION lib
        RUNTIME DESTINATION bin)
install(DIRECTORY . COMPONENT headers DESTINATION include/rokko
        FILES_MATCHING PATTERN "*.h" PATTERN "*.hpp" PATTERN "*.ipp"
        PATTERN ".svn" EXCLUDE PATTERN "CVS" EXCLUDE)
