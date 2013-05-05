#ifndef ROKKO_EIGEN_S_DIAGONALIZE_HPP
#define ROKKO_EIGEN_S_DIAGONALIZE_HPP

#include <mpi.h>

#include <rokko/eigen_s/eigen_s.hpp>


namespace rokko {

namespace eigen_s {

template<typename MATRIX_MAJOR>
void diagonalize(rokko::distributed_matrix<MATRIX_MAJOR>& mat, Eigen::VectorXd& eigvals, rokko::distributed_matrix<MATRIX_MAJOR>& eigvecs)
{
}

template<typename MATRIX_MAJOR>
void diagonalize(rokko::distributed_matrix<MATRIX_MAJOR>& mat, double* eigvals, rokko::distributed_matrix<MATRIX_MAJOR>& eigvecs)
{
}

template<>
void diagonalize(rokko::distributed_matrix<rokko::matrix_row_major>& mat, double* eigvals, rokko::distributed_matrix<rokko::matrix_row_major>& eigvecs)
{
  int m = 32;  // block_size

  int iflag;
  //if (flag_eigvecs)
  iflag = 0;
  //else
  //  iflag = 1;

  double start, end;
  MPI_Barrier(MPI_COMM_WORLD);
  start = MPI_Wtime();

  //eigen_s_(mat.m_global, mat.array, mat.lld, &eigvals[0], double_null_ptr, mat.lld, m, iflag);
  eigen_s_(mat.m_global, mat.array, mat.lld, eigvals, eigvecs.array, mat.lld, m, iflag);

  MPI_Barrier(MPI_COMM_WORLD);
  end = MPI_Wtime();
}

template<>
void diagonalize(rokko::distributed_matrix<rokko::matrix_row_major>& mat, Eigen::VectorXd& eigvals, rokko::distributed_matrix<rokko::matrix_row_major>& eigvecs)
{
  return diagonalize(mat, &eigvals[0], eigvecs);
}

} // namespace eigen_s

} // namespace rokko

#endif // ROKKO_EIGEN_S_DIAGONALIZE_HPP

