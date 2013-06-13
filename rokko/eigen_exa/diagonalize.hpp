#ifndef ROKKO_EIGEN_EXA_DIAGONALIZE_HPP
#define ROKKO_EIGEN_EXA_DIAGONALIZE_HPP

#include <mpi.h>
#include <rokko/distributed_matrix.hpp>
#include <rokko/localized_vector.hpp>
#include <rokko/eigen_exa/eigen_exa.hpp>

namespace rokko {
namespace eigen_exa {

template <typename MATRIX_MAJOR>
void diagonalize(rokko::distributed_matrix<MATRIX_MAJOR>& mat, double* eigvals, rokko::distributed_matrix<MATRIX_MAJOR>& eigvecs, timer& timer_in) {
  //if(mat.g.is_row_major()) throw "eigen_sx doesn't support grid_row_major.  Use eigen_sx with grid_col_major.";
  if(mat.is_row_major()) throw "eigen_sx doesn't support matrix_row_major.  Use eigen_sx with matrix_col_major.";

  MPI_Fint comm = MPI_Comm_c2f(MPI_COMM_WORLD);
  char order = 'C'; 
  eigen_init(&comm, &order);

  int m = 32;  // block_size

  int dim = mat.get_m_global();
  int lld = mat.get_lld();
  double* mat_array = mat.get_array_pointer();
  double* eigvecs_array= eigvecs.get_array_pointer();

  std::cout << "dim=" << dim << std::endl;
  std::cout << "lld=" << lld << std::endl;

  int m_forward = 8;
  int m_backward = 128;

  timer_in.start(1);
  eigen_sx(dim, dim, mat_array, lld, eigvals, eigvecs_array, lld, m_forward, m_backward);
  timer_in.stop(1);

  eigen_free();
}

template<typename MATRIX_MAJOR>
void diagonalize(rokko::distributed_matrix<MATRIX_MAJOR>& mat, localized_vector& eigvals, rokko::distributed_matrix<MATRIX_MAJOR>& eigvecs, timer& timer_in) {
  //if(mat.g.is_row_major()) throw "eigen_sx doesn't support grid_row_major.  Use eigen_sx with grid_col_major.";
  if(mat.is_row_major()) throw "eigen_sx doesn't support matrix_row_major.  Use eigen_sx with matrix_col_major.";

  diagonalize(mat, &eigvals[0], eigvecs, timer_in);
}

} // namespace eigen_exa
} // namespace rokko

#endif // ROKKO_EIGEN_EXA_DIAGONALIZE_HPP
