#ifndef ROKKO_EIGEN_SX_DIAGONALIZE_HPP
#define ROKKO_EIGEN_SX_DIAGONALIZE_HPP

#include <mpi.h>
#include <rokko/distributed_matrix.hpp>
#include <rokko/localized_vector.hpp>
#include <rokko/eigen_sx/eigen_sx.hpp>

namespace rokko {
namespace eigen_sx {

template <typename MATRIX_MAJOR>
void diagonalize(rokko::distributed_matrix<MATRIX_MAJOR>& mat, double* eigvals, rokko::distributed_matrix<MATRIX_MAJOR>& eigvecs, timer& timer_in) {
  if(mat.get_grid().is_row_major()) throw "eigen_sx doesn't support grid_row_major.  Use eigen_sx with grid_col_major.";
  if(mat.is_row_major()) throw "eigen_sx doesn't support matrix_row_major.  Use eigen_sx with matrix_col_major.";

  int m = 32;  // block_size

  int iflag;
  //if (flag_eigvecs)
  iflag = 0;
  //else
  //iflag = 1;

  int nme = ((mat.get_m_global() - 1) / 2 + 1) * 2;
  double* d = new double[nme];
  if (d == 0) {
    std::cerr << "failed to allocate d." << std::endl;
    MPI_Abort(MPI_COMM_WORLD, 3);
  }
  double* e = new double[2 * nme];
  if (e == 0) {
    std::cerr << "failed to allocate e." << std::endl;
    MPI_Abort(MPI_COMM_WORLD, 3);
  }

  int dim = mat.get_m_global();
  int lld = mat.get_lld();
  double* mat_array = mat.get_array_pointer();
  double* eigvecs_array= eigvecs.get_array_pointer();

  timer_in.start(1);
  eigen_sx_(dim, mat_array, lld, eigvals, eigvecs_array, d, e, nme, m, iflag);
  timer_in.stop(1);

  delete[] d;
  delete[] e;
}

template<typename MATRIX_MAJOR>
void diagonalize(rokko::distributed_matrix<MATRIX_MAJOR>& mat, localized_vector& eigvals, rokko::distributed_matrix<MATRIX_MAJOR>& eigvecs, timer& timer_in) {
  if(mat.get_grid().is_row_major()) throw "eigen_sx doesn't support grid_row_major.  Use eigen_sx with grid_col_major.";
  if(mat.is_row_major()) throw "eigen_sx doesn't support matrix_row_major.  Use eigen_sx with matrix_col_major.";

  diagonalize(mat, &eigvals[0], eigvecs, timer_in);
}

} // namespace eigen_sx
} // namespace rokko

#endif // ROKKO_EIGEN_SX_DIAGONALIZE_HPP
