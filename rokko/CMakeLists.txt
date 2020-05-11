#
# Rokko: Integrated Interface for libraries of eigenvalue decomposition
#
# Copyright (C) 2012-2017 Rokko Developers https://github.com/t-sakashita/rokko
#
# Distributed under the Boost Software License, Version 1.0. (See accompanying
# file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
#

set(ROKKO_SOURCES)

# For serial solvers
if(BUILD_SERIAL_DENSE_SOLVER)
  set(ROKKO_SOURCES ${ROKKO_SOURCES} serial_dense_ev.cpp utility/timer.cpp)
  # LAPACK
  set(ROKKO_SOURCES ${ROKKO_SOURCES} lapack/lapack.cpp lapack/diagonalize_get_parameters.cpp)
  # Eigen3
  set(ROKKO_SOURCES ${ROKKO_SOURCES} eigen3/eigen3.cpp)


  if(BUILD_C_BINDING)
    set(ROKKO_SOURCES ${ROKKO_SOURCES} eigen_vector_c.cpp eigen_matrix_c.cpp
      serial_dense_ev_c.cpp utility/minij_matrix_c.cpp utility/frank_matrix_c.cpp utility/helmert_matrix_c.cpp utility/matrix012_c.cpp copy_string.cpp parameters_c.cpp solver_name_c.cpp)
  endif(BUILD_C_BINDING)
  if(BUILD_Fortran_BINDING)
    set(ROKKO_SOURCES ${ROKKO_SOURCES} rokko_dense.F90 rokko.F90)
    set(ROKKO_SOURCES ${ROKKO_SOURCES} string.F90 command_line.F90 parameters.F90 solver_name.F90 eigen_vector.F90 eigen_matrix.F90 serial_dense_ev.F90)
    set(ROKKO_SOURCES ${ROKKO_SOURCES} serial_frank_matrix.F90 serial_laplacian_matrix.F90 serial_minij_matrix.F90 serial_helmert_matrix.F90 serial_matrix012.F90)
  endif(BUILD_Fortran_BINDING)
endif(BUILD_SERIAL_DENSE_SOLVER)

if(BUILD_PARALLEL_DENSE_SOLVER)
  set(ROKKO_SOURCES ${ROKKO_SOURCES} grid.cpp matrix_major.cpp)
  set(ROKKO_SOURCES ${ROKKO_SOURCES} mapping_1d.cpp distributed_crs_matrix.cpp parallel_dense_ev.cpp)
  set(ROKKO_SOURCES ${ROKKO_SOURCES} parallel_frank_matrix.F90 parallel_laplacian_matrix.F90 parallel_minij_matrix.F90 parallel_helmert_matrix.F90 parallel_matrix012.F90 various_mpi_comm.F90)

  # PBLAS, ScaLAPACK
  if(ROKKO_HAVE_SCALAPACK)
    set(ROKKO_SOURCES ${ROKKO_SOURCES} pblas/pcopy.c pblas/pdot.c pblas/pdotu.c pblas/pdotc.c pblas/pgemv.c pblas/pgemm.c)
    set(ROKKO_SOURCES ${ROKKO_SOURCES} scalapack/descinit.c scalapack/indxg2p.c scalapack/numroc.c scalapack/pdelget.c scalapack/pdelset.c scalapack/pdlamch.c scalapack/pdlange.c scalapack/pdlange_work.c scalapack/pdlaprnt.c scalapack/psyev_work.c scalapack/psyev.c scalapack/psyevd_work.c scalapack/psyevd.c scalapack/psyevx_work.c scalapack/psyevx.c)
    set(ROKKO_SOURCES ${ROKKO_SOURCES} scalapack/pselget.c scalapack/pselset.c scalapack/pslamch.c scalapack/pslange.c scalapack/pslange_work.c scalapack/pslaprnt.c)
    set(ROKKO_SOURCES ${ROKKO_SOURCES} scalapack/pstebz.c scalapack/pstebz_work.c scalapack/pstein.c scalapack/pstein_work.c scalapack/pstedc.c scalapack/pstedc_work.c)
    if(ROKKO_HAVE_SCALAPACK_PDSYEVR)
      set(ROKKO_SOURCES ${ROKKO_SOURCES} scalapack/psyevr.c scalapack/psyevr_work.c)
    endif(ROKKO_HAVE_SCALAPACK_PDSYEVR)
    if(BUILD_SHARED_LIBS)
      set(ROKKO_SOURCES ${ROKKO_SOURCES} scalapack/scalapack.cpp)
    endif(BUILD_SHARED_LIBS)
  endif(ROKKO_HAVE_SCALAPACK)

  # EigenExa
  if(ROKKO_HAVE_EIGENEXA)
    set(ROKKO_SOURCES ${ROKKO_SOURCES} eigenexa/eigenexa_wrap.f90 eigenexa/eigen_init.c eigenexa/eigen_free.c eigenexa/eigen_get_matdims.c eigenexa/eigen_get_procs.c eigenexa/eigen_get_id.c eigenexa/eigen_loop_start.c eigenexa/eigen_loop_end.c eigenexa/eigen_translate_l2g.c eigenexa/eigen_translate_g2l.c eigenexa/eigen_s.c eigenexa/eigen_sx.c)
    if(BUILD_SHARED_LIBS)
      set(ROKKO_SOURCES ${ROKKO_SOURCES} eigenexa/eigenexa.cpp)
    endif(BUILD_SHARED_LIBS)
  endif(ROKKO_HAVE_EIGENEXA)

  # Elemental
  if(ROKKO_HAVE_ELEMENTAL)
    if(BUILD_SHARED_LIBS)
      set(ROKKO_SOURCES ${ROKKO_SOURCES} elemental/elemental.cpp)
    endif(BUILD_SHARED_LIBS)
  endif(ROKKO_HAVE_ELEMENTAL)

  # ELPA
  if(ROKKO_HAVE_ELPA)
    include_directories(${ELPA_INCLUDE_DIR})
    if(BUILD_SHARED_LIBS)
      set(ROKKO_SOURCES ${ROKKO_SOURCES} elpa/elpa.cpp)
    endif(BUILD_SHARED_LIBS)
  endif(ROKKO_HAVE_ELPA)

  if(BUILD_C_BINDING)
    set(ROKKO_SOURCES ${ROKKO_SOURCES} grid_c.cpp mapping_bc_c.cpp distributed_matrix_c.cpp collective_c.cpp)
    set(ROKKO_SOURCES ${ROKKO_SOURCES} parallel_dense_ev_c.cpp)
  endif(BUILD_C_BINDING)
  if(BUILD_Fortran_BINDING)
    set(ROKKO_SOURCES ${ROKKO_SOURCES} grid.F90 mapping_bc.F90 distributed_matrix.F90 collective.F90)
    set(ROKKO_SOURCES ${ROKKO_SOURCES} parallel_dense_ev.F90 rokko_parallel_dense.F90)
  endif(BUILD_Fortran_BINDING)
endif(BUILD_PARALLEL_DENSE_SOLVER)

if(BUILD_PARALLEL_SPARSE_SOLVER)
  set(ROKKO_SOURCES ${ROKKO_SOURCES} parallel_sparse_ev.cpp)
  if(BUILD_C_BINDING)
    set(ROKKO_SOURCES ${ROKKO_SOURCES} mapping_1d_c.cpp parallel_sparse_ev_c.cpp distributed_crs_matrix_c.cpp distributed_mfree_c.cpp)
  endif(BUILD_C_BINDING)
  if(BUILD_Fortran_BINDING)
    set(ROKKO_SOURCES ${ROKKO_SOURCES} mapping_1d.F90 rokko_parallel_sparse.F90 parallel_sparse_ev.F90 distributed_crs_matrix.F90 distributed_mfree.F90 distributed_mfree_f.cpp)
  endif(BUILD_Fortran_BINDING)

  # original Lanczos
  if(BUILD_SHARED_LIBS)
    set(ROKKO_SOURCES ${ROKKO_SOURCES} original/original.cpp)
  endif(BUILD_SHARED_LIBS)

  # SLEPc
  if(ROKKO_HAVE_SLEPC)
    if(BUILD_SHARED_LIBS)
      set(ROKKO_SOURCES ${ROKKO_SOURCES} slepc/slepc.cpp)
    endif(BUILD_SHARED_LIBS)
  endif(ROKKO_HAVE_SLEPC)

  # Anasazi
  if(ROKKO_HAVE_ANASAZI)
    if(BUILD_SHARED_LIBS)
      set(ROKKO_SOURCES ${ROKKO_SOURCES} anasazi/anasazi.cpp)
    endif(BUILD_SHARED_LIBS)
  endif(ROKKO_HAVE_ANASAZI)
endif(BUILD_PARALLEL_SPARSE_SOLVER)

if(ROKKO_SOURCES)
  add_library(rokko ${ROKKO_SOURCES})
  if (BUILD_SHARED_LIBS)
    set_target_properties(rokko PROPERTIES COMPILE_DEFINITIONS ROKKO_BUILD_SHARED_LIBS=1)
  endif (BUILD_SHARED_LIBS)

  if(BUILD_CBLAS)
    target_link_libraries(rokko cblas)
  else(BUILD_CBLAS)
    target_link_libraries(rokko ${BLAS_LIBRARIES})
  endif(BUILD_CBLAS)

  if(BUILD_LAPACKE)
    target_link_libraries(rokko lapacke)
  else(BUILD_LAPACKE)
    target_link_libraries(rokko ${LAPACK_LIBRARIES})
  endif(BUILD_LAPACKE)

  if(ROKKO_HAVE_MPI)
    if(BUILD_Fortran_BINDING)
      target_link_libraries(rokko ${MPI_Fortran_LIBRARIES})
    endif(BUILD_Fortran_BINDING)
    target_link_libraries(rokko ${MPI_CXX_LIBRARIES})
  endif(ROKKO_HAVE_MPI)

  if(ROKKO_HAVE_SCALAPACK)
    target_link_libraries(rokko ${SCALAPACK_LIBRARIES})
  endif(ROKKO_HAVE_SCALAPACK)

  if(ROKKO_HAVE_EIGENEXA)
    target_link_libraries(rokko ${EIGENEXA_LIBRARIES})
  endif(ROKKO_HAVE_EIGENEXA)

  if(ROKKO_HAVE_ELEMENTAL)
    target_link_libraries(rokko ${ELEMENTAL_LIBRARIES})
  endif(ROKKO_HAVE_ELEMENTAL)

  if(ROKKO_HAVE_ELPA)
    target_link_libraries(rokko ${ELPA_LIBRARIES})
  endif(ROKKO_HAVE_ELPA)

  if(ROKKO_HAVE_PETSC)
    target_link_libraries(rokko ${PETSC_LIBRARIES})
  endif(ROKKO_HAVE_PETSC)

  if(ROKKO_HAVE_SLEPC)
    target_link_libraries(rokko ${SLEPC_LIBRARIES})
  endif(ROKKO_HAVE_SLEPC)

  if(ROKKO_HAVE_ANASAZI)
    target_link_libraries(rokko ${ANASAZI_LIBRARIES})
  endif(ROKKO_HAVE_ANASAZI)

  if(BUILD_Python_BINDING)
    pybind11_add_module(pyrokko pyrokko.cpp)
    target_include_directories(pyrokko PUBLIC ${MPI4PY_DIR}/include)
    target_link_libraries(pyrokko PRIVATE rokko)
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
    install(TARGETS pyrokko
            ARCHIVE DESTINATION lib
            LIBRARY DESTINATION lib
            RUNTIME DESTINATION bin)
    install(DIRECTORY . DESTINATION lib
            FILES_MATCHING PATTERN "*.py"
            PATTERN ".svn" EXCLUDE PATTERN "CVS" EXCLUDE)
  endif(BUILD_Python_BINDING)

endif(ROKKO_SOURCES)
