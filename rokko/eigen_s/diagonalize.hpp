#ifndef ROKKO_EIGEN_S_DIAGONALIZE_HPP
#define ROKKO_EIGEN_S_DIAGONALIZE_HPP

#include <mpi.h>

#include <rokko/eigen_s/eigen_s.hpp>


namespace rokko {

template<typename T>
void diagonalize(rokko::distributed_matrix<T>& mat, Eigen::VectorXd& eigvals, rokko::distributed_matrix<T>& eigvecs, bool flag_eigvecs)
{
}

template<typename T>
void diagonalize(rokko::distributed_matrix<T>& mat, double* eigvals, rokko::distributed_matrix<T>& eigvecs, bool flag_eigvecs)
{
}

template<>
void diagonalize(rokko::distributed_matrix<rokko::eigen_s>& mat, double* eigvals, rokko::distributed_matrix<rokko::eigen_s>& eigvecs, bool flag_eigvecs)
{
  int m = 32;  // block_size

  int iflag;
  if (flag_eigvecs)
    iflag = 0;
  else
    iflag = 1;

  double start, end;
  MPI_Barrier(MPI_COMM_WORLD);
  start = MPI_Wtime();

  //eigen_s_(mat.m_global, mat.array, mat.lld, &eigvals[0], double_null_ptr, mat.lld, m, iflag);
  eigen_s_(mat.m_global, mat.array, mat.lld, eigvals, eigvecs.array, mat.lld, m, iflag);

  MPI_Barrier(MPI_COMM_WORLD);
  end = MPI_Wtime();
}

template<>
void diagonalize(rokko::distributed_matrix<rokko::eigen_s>& mat, Eigen::VectorXd& eigvals, rokko::distributed_matrix<rokko::eigen_s>& eigvecs, bool flag_eigvecs)
{
  return diagonalize(mat, &eigvals[0], eigvecs, flag_eigvecs);
}


} // namespace rokko

#endif // ROKKO_EIGEN_S_DIAGONALIZE_HPP

