if(BUILD_LAPACKE)
  include(FortranCInterface)
  FortranCInterface_HEADER(${CMAKE_CURRENT_BINARY_DIR}/include/lapacke_mangling.h 
                           MACRO_NAMESPACE "LAPACK_"
                           SYMBOL_NAMESPACE "LAPACK_")
  set(SOURCES src/lapacke_dsyev.c src/lapacke_dsyev_work.c src/lapacke_dsyevd.c src/lapacke_dsyevd_work.c src/lapacke_dsyevr.c src/lapacke_dsyevr_work.c src/lapacke_dsyevx.c src/lapacke_dsyevx_work.c)
  add_library(lapacke ${SOURCES})
  install(TARGETS lapacke ARCHIVE DESTINATION lib LIBRARY DESTINATION lib RUNTIME DESTINATION bin)
endif(BUILD_LAPACKE)
