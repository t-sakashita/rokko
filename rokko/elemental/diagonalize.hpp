#ifndef ROKKO_ELEMENTAL_DIAGONALIZE_HPP
#define ROKKO_ELEMENTAL_DIAGONALIZE_HPP

#include <Eigen/Dense>

#include "elemental.hpp"

namespace rokko {

template<typename T>
void diagonalize(rokko::distributed_matrix<T>& mat, Eigen::VectorXd& eigvals, rokko::distributed_matrix<rokko::elemental>& eigvecs)
{
}

template<typename T>
void diagonalize(rokko::distributed_matrix<T>& mat, Eigen::VectorXd& eigvals)
{
}

template<>
void diagonalize<rokko::elemental>(rokko::distributed_matrix<rokko::elemental>& mat, Eigen::VectorXd& eigvals, rokko::distributed_matrix<rokko::elemental>& eigvecs)
{
  //int m = 32;  // block_size

  elem::DistMatrix<double,elem::VR,elem::STAR> w(*(mat.g.get_elem_grid()) );

  double start, end;
  MPI_Barrier(MPI_COMM_WORLD);
  start = MPI_Wtime();

  elem::HermitianEig(elem::LOWER, mat.mat, w, eigvecs.mat); // only access lower half of H

  MPI_Barrier(MPI_COMM_WORLD);
  end = MPI_Wtime();

  for (int i=0; i<eigvals.size(); ++i)
    eigvals(i) = w.Get(i, 0);
}

template<>
void diagonalize<rokko::elemental>(rokko::distributed_matrix<rokko::elemental>& mat, Eigen::VectorXd& eigvals)
{
  //int m = 32;  // block_size

  elem::DistMatrix<double,elem::VR,elem::STAR> w(*(mat.g.get_elem_grid()) );

  double start, end;
  MPI_Barrier(MPI_COMM_WORLD);
  start = MPI_Wtime();

  elem::HermitianEig(elem::LOWER, mat.mat, w); // only access lower half of H

  MPI_Barrier(MPI_COMM_WORLD);
  end = MPI_Wtime();

  for (int i=0; i<eigvals.size(); ++i)
    eigvals(i) = w.Get(i, 0);
}


} // namespace rokko

#endif // ROKKO_ELEMENTAL_DIAGONALIZE_HPP

