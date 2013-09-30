include(FortranCInterface)
FortranCInterface_HEADER(${CMAKE_CURRENT_BINARY_DIR}/include/lapacke_mangling.h 
                         MACRO_NAMESPACE "LAPACK_"
                         SYMBOL_NAMESPACE "LAPACK_")
set(SOURCES src/lapacke_dsyev.c src/lapacke_dsyev_work.c src/lapacke_dsyevd.c src/lapacke_dsyevd_work.c src/lapacke_dsyevr.c src/lapacke_dsyevr_work.c src/lapacke_dsyevx.c src/lapacke_dsyevx_work.c
  utils/lapacke_d_nancheck.c utils/lapacke_dge_nancheck.c utils/lapacke_dge_trans.c utils/lapacke_dsy_nancheck.c utils/lapacke_dsy_trans.c utils/lapacke_dtr_nancheck.c utils/lapacke_dtr_trans.c utils/lapacke_lsame.c utils/lapacke_xerbla.c)
add_library(lapacke ${SOURCES})
if (NOT STATIC_SCALAPACK)
  target_link_libraries(lapacke ${SCALAPACK_LIB})
endif (NOT STATIC_SCALAPACK)
install(TARGETS lapacke ARCHIVE DESTINATION lib LIBRARY DESTINATION lib RUNTIME DESTINATION bin)
