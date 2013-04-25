#ifndef ROKKO_EIGEN_S_DIAGONALIZE_H
#define ROKKO_EIGEN_S_DIAGONALIZE_H

#include <mpi.h>

#include <rokko/eigen_s/eigen_s.hpp>


namespace rokko {

template<typename T>
void diagonalize(rokko::distributed_matrix<T>& mat, Eigen::VectorXd& eigvals, rokko::distributed_matrix<T>& eigvecs, bool flag_eigvecs)
{
}

template<>
void diagonalize(rokko::distributed_matrix<rokko::eigen_s>& mat, Eigen::VectorXd& eigvals, rokko::distributed_matrix<rokko::eigen_s>& eigvecs, bool flag_eigvecs)
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

  eigen_s_(mat.m_global, mat.array, mat.lld, &eigvals[0], eigvecs.array, mat.lld, m, iflag);

  MPI_Barrier(MPI_COMM_WORLD);
  end = MPI_Wtime();
}

} // namespace rokko

#endif // ROKKO_EIGEN_S_DIAGONALIZE_H

