#ifndef ROKKO_ELEMENTAL_DIAGONALIZE_HPP
#define ROKKO_ELEMENTAL_DIAGONALIZE_HPP

#include <Eigen/Dense>

#include "elemental.hpp"

namespace rokko {
  //namespace elemental {

template<typename T>
void diagonalize(rokko::distributed_matrix<T>& mat, Eigen::VectorXd& eigvals, rokko::distributed_matrix<rokko::elemental>& eigvecs)
//void diagonalize<T>(MATRIX& mat, VECTOR& eigvals, MATRIX& eigvecs)
{
}

  template<>
  //template <class MATRIX, class VECTOR>
  void diagonalize<rokko::elemental>(rokko::distributed_matrix<rokko::elemental>& mat, Eigen::VectorXd& eigvals, rokko::distributed_matrix<rokko::elemental>& eigvecs)

 //template <class MATRIX, class VECTOR>
 //void diagonalize<rokko::elemental>(MATRIX& mat, VECTOR& eigvals, MATRIX& eigvecs)
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

  //} // namespace elemental
} // namespace rokko

#endif // ROKKO_ELEMENTAL_DIAGONALIZE_HPP

