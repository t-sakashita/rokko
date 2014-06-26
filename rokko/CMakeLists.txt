#
# Rokko: Integrated Interface for libraries of eigenvalue decomposition
#
# Copyright (C) 2012-2013 by Tatsuya Sakashita <t-sakashita@issp.u-tokyo.ac.jp>,
#                            Synge Todo <wistaria@comp-phys.org>
#
# Distributed under the Boost Software License, Version 1.0. (See accompanying
# file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
#

set(ROKKO_SOURCES grid.cpp serial_dense_solver.cpp parallel_dense_solver.cpp)
if(BUILD_C_BINDING)
  set(ROKKO_SOURCES ${ROKKO_SOURCES} grid_c.cpp solver_c.cpp localized_vector_c.cpp distributed_matrix_c.cpp)
  set(ROKKO_SOURCES ${ROKKO_SOURCES} utility/frank_matrix_c.cpp)
endif(BUILD_C_BINDING)
if(BUILD_Fortran_BINDING)
  set(ROKKO_SOURCES ${ROKKO_SOURCES} rokko.f90 solver_f.f90)
  set(ROKKO_SOURCES ${ROKKO_SOURCES} utility/frank_matrix_f.f90)
endif(BUILD_Fortran_BINDING)

# Eigen3
set(ROKKO_SOURCES ${ROKKO_SOURCES} eigen3/eigen3.cpp)

# LAPACK
set(ROKKO_SOURCES ${ROKKO_SOURCES} lapack/lapack.cpp)

# BLACS
set(ROKKO_SOURCES ${ROKKO_SOURCES} blacs/blacs_barrier.c blacs/blacs_exit.c blacs/blacs_get.c blacs/blacs_gridexit.c blacs/blacs_gridinfo.c blacs/blacs_gridinit.c blacs/blacs_pinfo.c blacs/descinit.c)
set(ROKKO_SOURCES ${ROKKO_SOURCES} pblas/pdgemm.c)

if(Boost_SOURCE_DIR)
  set(SOURCES error_code.cpp)
  foreach (D ${Boost_SOURCE_DIR}/libs/system/src)
    foreach(S ${SOURCES})
      if(EXISTS ${D}/${S})
        set(ROKKO_SOURCES ${ROKKO_SOURCES} ${D}/${S})
      endif(EXISTS ${D}/${S})
    endforeach(S)
  endforeach(D)
endif(Boost_SOURCE_DIR)

if(BUILD_SCALAPACK)
  set(ROKKO_SOURCES ${ROKKO_SOURCES} scalapack/pdelget.c scalapack/pdelset.c scalapack/pdlamch.c scalapack/pdlaprnt.c scalapack/pdsyev.c scalapack/pdsyevd.c scalapack/pdsyevx.c)
  if(USE_PDSYEVR)
    set(ROKKO_SOURCES ${ROKKO_SOURCES} scalapack/pdsyevr.c)
  endif(USE_PDSYEVR)
  if(BUILD_SHARED_LIBS)
    set(ROKKO_SOURCES ${ROKKO_SOURCES} scalapack/scalapack.cpp)
    set(ROKKO_SOURCES ${ROKKO_SOURCES} scalapack/scalapack_pdsyevd.cpp)
    # set(ROKKO_SOURCES ${ROKKO_SOURCES} scalapack/scalapack_pdsyevx.cpp)
  endif(BUILD_SHARED_LIBS)
endif(BUILD_SCALAPACK)

if(BUILD_EIGENEXA)
  include_directories(${EIGENEXA_INCLUDE_DIR})
  set(ROKKO_SOURCES ${ROKKO_SOURCES} eigen_exa/eigen_exa_wrap.f90 eigen_exa/eigen_init.c eigen_exa/eigen_free.c eigen_exa/cstab_get_optdim.c eigen_exa/eigen_sx.c)
  if(BUILD_SHARED_LIBS)
    set(ROKKO_SOURCES ${ROKKO_SOURCES} eigen_exa/eigen_exa.cpp)
  endif(BUILD_SHARED_LIBS)
endif(BUILD_EIGENEXA)

if(BUILD_ELEMENTAL)
  include_directories(${ELEMENTAL_INCLUDE_DIR})
  if(BUILD_SHARED_LIBS)
    set(ROKKO_SOURCES ${ROKKO_SOURCES} elemental/elemental.cpp)
  endif(BUILD_SHARED_LIBS)
endif(BUILD_ELEMENTAL)

if(BUILD_ELPA)
  if(BUILD_SHARED_LIBS)
    set(ROKKO_SOURCES ${ROKKO_SOURCES} elpa/elpa.cpp)
  endif(BUILD_SHARED_LIBS)
endif(BUILD_ELPA)

add_library(rokko ${ROKKO_SOURCES})
if (BUILD_SHARED_LIBS)
  set_target_properties(rokko PROPERTIES COMPILE_DEFINITIONS ROKKO_BUILD_SHARED_LIBS=1)
endif (BUILD_SHARED_LIBS)
target_link_libraries(rokko lapacke)

if (MPI_CXX_FOUND)
  target_link_libraries(rokko ${MPI_CXX_LIBRARIES})
endif(MPI_CXX_FOUND)
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
if(BUILD_PETSC)
  target_link_libraries(rokko ${PETSC_LIBRARIES})
endif(BUILD_PETSC)
if(BUILD_SLEPC)
  target_link_libraries(rokko ${SLEPC_LIBRARIES})
endif(BUILD_SLEPC)
if(BUILD_ANASAZI)
  target_link_libraries(rokko ${ANASAZI_LIBRARIES})
endif(BUILD_ANASAZI)

if(Boost_FOUND)
  target_link_libraries(rokko ${Boost_LIBRARIES})
endif(Boost_FOUND)

install(TARGETS rokko COMPONENT libraries
        ARCHIVE DESTINATION lib
        LIBRARY DESTINATION lib
        RUNTIME DESTINATION bin)
install(DIRECTORY . COMPONENT headers DESTINATION include/rokko
        FILES_MATCHING PATTERN "*.h" PATTERN "*.hpp" PATTERN "*.ipp"
        PATTERN ".svn" EXCLUDE PATTERN "CVS" EXCLUDE)
