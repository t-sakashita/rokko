#ifndef ROKKO_EIGEN_SX_DIAGONALIZE_HPP
#define ROKKO_EIGEN_SX_DIAGONALIZE_HPP

#include <mpi.h>

#include <rokko/eigen_sx/eigen_sx.hpp>

namespace rokko {

template<typename T>
void diagonalize(rokko::distributed_matrix<T>& mat, Eigen::VectorXd& eigvals, rokko::distributed_matrix<T>& eigvecs, bool flag_eigvecs)
{
}


template <>
void diagonalize(rokko::distributed_matrix<rokko::eigen_sx>& mat, Eigen::VectorXd& eigvals, rokko::distributed_matrix<rokko::eigen_sx>& eigvecs, bool flag_eigvecs)
{
  int m = 32;  // block_size

  int iflag;
  if (flag_eigvecs)
    iflag = 0;
  else
    iflag = 1;

  double* d = new double[mat.nme];
  if (d == NULL) {
    cerr << "failed to allocate d." << endl;
    MPI_Abort(MPI_COMM_WORLD, 3);
  }
  double* e = new double[2 * mat.nme];
  if (e == NULL) {
    cerr << "failed to allocate e." << endl;
    MPI_Abort(MPI_COMM_WORLD, 3);
  }

  double start, end;
  MPI_Barrier(MPI_COMM_WORLD);
  start = MPI_Wtime();

  eigen_sx_(mat.m_global, mat.array, mat.lld, &eigvals[0], eigvecs.array, d, e, mat.nme, m, iflag);

  MPI_Barrier(MPI_COMM_WORLD);
  end = MPI_Wtime();

  delete[] d;
  delete[] e;
}

} // namespace rokko

#endif // ROKKO_EIGEN_SX_DIAGONALIZE_HPP
