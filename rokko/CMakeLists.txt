#
# Rokko: Integrated Interface for libraries of eigenvalue decomposition
#
# Copyright (C) 2012-2015 Rokko Developers https://github.com/t-sakashita/rokko
#
# Distributed under the Boost Software License, Version 1.0. (See accompanying
# file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
#

# For serial solvers
set(ROKKO_SOURCES serial_dense_ev.cpp utility/timer.cpp)
if(BUILD_C_BINDING)
  set(ROKKO_SOURCES ${ROKKO_SOURCES} localized_vector_c.cpp localized_matrix_c.cpp
    serial_dense_ev_c.cpp utility/frank_matrix_c.cpp parameters_c.cpp solver_name_c.cpp)
endif(BUILD_C_BINDING)
if(BUILD_Fortran_BINDING)
  set(ROKKO_SOURCES ${ROKKO_SOURCES} string.F90 parameters.F90 solver_name.F90 rokko_serial_dense.F90)
endif(BUILD_Fortran_BINDING)

if(BUILD_PARALLEL_DENSE_SOLVER)
  set(ROKKO_SOURCES ${ROKKO_SOURCES} grid.cpp matrix_major.cpp)
  set(ROKKO_SOURCES ${ROKKO_SOURCES} parallel_dense_ev.cpp)
  if(BUILD_C_BINDING)
    set(ROKKO_SOURCES ${ROKKO_SOURCES} grid_c.cpp mapping_bc_c.cpp distributed_matrix_c.cpp)
    set(ROKKO_SOURCES ${ROKKO_SOURCES} parallel_dense_ev_c.cpp)
  endif(BUILD_C_BINDING)
  if(BUILD_Fortran_BINDING)
    set(ROKKO_SOURCES ${ROKKO_SOURCES} distributed_matrix.f90 collective.f90)
    set(ROKKO_SOURCES ${ROKKO_SOURCES} rokko_parallel_dense_classes.F90 rokko_parallel_dense.F90 parallel_dense_ev.f90)
  endif(BUILD_Fortran_BINDING)
endif(BUILD_PARALLEL_DENSE_SOLVER)

if(BUILD_PARALLEL_SPARSE_SOLVER)
  set(ROKKO_SOURCES ${ROKKO_SOURCES} parallel_sparse_ev.cpp)
  if(BUILD_C_BINDING)
    set(ROKKO_SOURCES ${ROKKO_SOURCES} parallel_sparse_ev_c.cpp distributed_crs_matrix_c.cpp distributed_mfree_c.cpp)
  endif(BUILD_C_BINDING)
  if(BUILD_Fortran_BINDING)
    set(ROKKO_SOURCES ${ROKKO_SOURCES} rokko_sparse.F90 parallel_sparse_ev.f90 distributed_mfree_f.cpp)
  endif(BUILD_Fortran_BINDING)
endif(BUILD_PARALLEL_SPARSE_SOLVER)
# For all solvers
if(BUILD_Fortran_BINDING)
  set(ROKKO_SOURCES ${ROKKO_SOURCES} rokko_dense.F90 rokko.F90)
endif(BUILD_Fortran_BINDING)

# Boost
if(Boost_SOURCE_DIR)
  set(SOURCES error_code.cpp)
  foreach (D ${Boost_SOURCE_DIR}/libs/system/src)
    foreach(S ${SOURCES})
      if(EXISTS ${D}/${S})
        set(ROKKO_SOURCES ${ROKKO_SOURCES} ${D}/${S})
      endif(EXISTS ${D}/${S})
    endforeach(S)
  endforeach(D)
  if(BUILD_Python_BINDING)
    set(SOURCES dict.cpp errors.cpp exec.cpp import.cpp list.cpp long.cpp
        module.cpp numeric.cpp object_operators.cpp object_protocol.cpp slice.cpp
        str.cpp tuple.cpp wrapper.cpp converter/arg_to_python_base.cpp
        converter/builtin_converters.cpp converter/from_python.cpp
        converter/registry.cpp converter/type_id.cpp object/class.cpp
        object/enum.cpp object/function.cpp object/function_doc_signature.cpp
        object/inheritance.cpp object/iterator.cpp object/life_support.cpp
        object/pickle_support.cpp object/stl_iterator.cpp)
    foreach(D ${Boost_SOURCE_DIR}/libs/python/src)
      foreach(S ${SOURCES})
        if(EXISTS ${D}/${S})
          set(ROKKO_SOURCES ${ROKKO_SOURCES} ${D}/${S})
        endif(EXISTS ${D}/${S})
      endforeach(S)
    endforeach(D)
    add_definitions(-DBOOST_PYTHON_SOURCE)
  endif(BUILD_Python_BINDING)
endif(Boost_SOURCE_DIR)

# Eigen3
set(ROKKO_SOURCES ${ROKKO_SOURCES} eigen3/eigen3.cpp)

# LAPACK
set(ROKKO_SOURCES ${ROKKO_SOURCES} lapack/lapack.cpp lapack/geqrf.cpp lapack/gesvd.cpp lapack/heev.cpp lapack/orgqr.cpp lapack/syev.cpp lapack/sygv.cpp lapack/ungqr.cpp lapack/diagonalize_get_parameters.cpp)

# BLACS, PBLAS, ScaLAPACK
if(BUILD_SCALAPACK)
  set(ROKKO_SOURCES ${ROKKO_SOURCES} blacs/blacs_barrier.c blacs/blacs_exit.c blacs/blacs_get.c
      blacs/blacs_gridexit.c blacs/blacs_gridinfo.c blacs/blacs_gridinit.c blacs/blacs_pinfo.c
      blacs/descinit.c)
  set(ROKKO_SOURCES ${ROKKO_SOURCES} pblas/pcopy.cpp pblas/pdot.cpp pblas/pdotc.cpp pblas/pgemv.cpp pblas/pgemm.cpp)
  set(ROKKO_SOURCES ${ROKKO_SOURCES} scalapack/pdelget.c scalapack/pdelset.c scalapack/pdlamch.c
      scalapack/pdlaprnt.c scalapack/pdsyev_work.c scalapack/pdsyev.c scalapack/pdsyevd_work.c scalapack/pdsyevd.c scalapack/pdsyevx_work.c scalapack/pdsyevx.c)
  if(ROKKO_HAVE_PDSYEVR)
    set(ROKKO_SOURCES ${ROKKO_SOURCES} scalapack/pdsyevr.c scalapack/pdsyevr_work.c)
  endif(ROKKO_HAVE_PDSYEVR)
  if(BUILD_SHARED_LIBS)
    set(ROKKO_SOURCES ${ROKKO_SOURCES} scalapack/scalapack.cpp)
  endif(BUILD_SHARED_LIBS)
endif(BUILD_SCALAPACK)

# EigenExa
if(BUILD_EIGENEXA)
  set(ROKKO_SOURCES ${ROKKO_SOURCES} eigen_exa/eigen_exa_wrap.f90 eigen_exa/eigen_exa_init.c
      eigen_exa/eigen_exa_free.c eigen_exa/eigen_exa_get_optdim.c eigen_exa/eigen_exa_get_matdims.c
      eigen_exa/eigen_exa_s.c eigen_exa/eigen_exa_sx.c)
  if(BUILD_SHARED_LIBS)
    set(ROKKO_SOURCES ${ROKKO_SOURCES} eigen_exa/eigen_exa.cpp)
  endif(BUILD_SHARED_LIBS)
endif(BUILD_EIGENEXA)

# Elemental
if(BUILD_ELEMENTAL)
  if(BUILD_SHARED_LIBS)
    set(ROKKO_SOURCES ${ROKKO_SOURCES} elemental/elemental.cpp)
  endif(BUILD_SHARED_LIBS)
endif(BUILD_ELEMENTAL)

# ELPA
if(BUILD_ELPA)
  if(BUILD_SHARED_LIBS)
    set(ROKKO_SOURCES ${ROKKO_SOURCES} elpa/elpa.cpp)
  endif(BUILD_SHARED_LIBS)
endif(BUILD_ELPA)

# SLEPc
if(BUILD_SLEPC)
  if(BUILD_SHARED_LIBS)
    set(ROKKO_SOURCES ${ROKKO_SOURCES} slepc/slepc.cpp)
  endif(BUILD_SHARED_LIBS)
endif(BUILD_SLEPC)

# Anasazi
if(BUILD_ANASAZI)
  if(BUILD_SHARED_LIBS)
    set(ROKKO_SOURCES ${ROKKO_SOURCES} anasazi/anasazi.cpp)
  endif(BUILD_SHARED_LIBS)
endif(BUILD_ANASAZI)

add_library(rokko ${ROKKO_SOURCES})
if (BUILD_SHARED_LIBS)
  set_target_properties(rokko PROPERTIES COMPILE_DEFINITIONS ROKKO_BUILD_SHARED_LIBS=1)
endif (BUILD_SHARED_LIBS)
target_link_libraries(rokko lapacke)

if(MPI_CXX_FOUND)
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

if(BUILD_Python_BINDING)
  target_link_libraries(rokko ${PYTHON_LIBRARIES})
  python_add_module(rokko_ext rokko_ext.cpp)
  target_link_libraries(rokko_ext rokko)
endif(BUILD_Python_BINDING)

# install

install(TARGETS rokko COMPONENT libraries
        ARCHIVE DESTINATION lib
        LIBRARY DESTINATION lib
        RUNTIME DESTINATION bin)
install(DIRECTORY . COMPONENT headers DESTINATION include/rokko
        FILES_MATCHING PATTERN "*.h" PATTERN "*.hpp" PATTERN "*.ipp"
        PATTERN ".svn" EXCLUDE PATTERN "CVS" EXCLUDE)
if(BUILD_Fortran_BINDING)
  install(FILES ${CMAKE_CURRENT_BINARY_DIR}/rokko.mod DESTINATION include)
endif(BUILD_Fortran_BINDING)
if(BUILD_Python_BINDING)
  install(TARGETS rokko_ext
          ARCHIVE DESTINATION lib
          LIBRARY DESTINATION lib
          RUNTIME DESTINATION bin)
  install(DIRECTORY . DESTINATION lib
          FILES_MATCHING PATTERN "*.py"
	  PATTERN ".svn" EXCLUDE PATTERN "CVS" EXCLUDE)
endif(BUILD_Python_BINDING)
