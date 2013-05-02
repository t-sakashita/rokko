#ifndef ROKKO_ELEMENTAL_DIAGONALIZE_HPP
#define ROKKO_ELEMENTAL_DIAGONALIZE_HPP

#include <Eigen/Dense>

#include "elemental.hpp"
#include <rokko/elemental/core.hpp>

#include <rokko/grid.hpp>
#include <rokko/distributed_matrix.hpp>

namespace rokko {

void solver_elemental::diagonalize(rokko::distributed_matrix& mat, Eigen::VectorXd& eigvals, rokko::distributed_matrix& eigvecs)
{
  //int m = 32;  // block_size

  MPI_Comm comm = MPI_COMM_WORLD;
  cout << "mat.npcol = " << mat.npcol << endl;
  cout << "mat.lld = " << mat.lld << endl;

  elem::Grid   elem_grid(comm, mat.nprow, mat.npcol);

  elem::DistMatrix<double> elem_mat(mat.m_global, mat.n_global, 0, 0, mat.get_array(), mat.lld, elem_grid);

  //elem::DistMatrix<double> elem_eigvecs0(eigvecs.m_global, eigvecs.n_local, 0, 0, eigvecs.get_array(), eigvecs.lld, elem_grid);
  //elem::DistMatrix<double> elem_eigvecs;  //(eigvecs.m_global, eigvecs.n_local, elem_grid);
  //elem::View(elem_eigvecs, elem_eigvecs0);
  elem::DistMatrix<double> elem_eigvecs(mat.m_global, mat.n_global, elem_grid);

  elem::DistMatrix<double,elem::VR,elem::STAR> elem_w(elem_grid);

  double start, end;
  MPI_Barrier(MPI_COMM_WORLD);
  start = MPI_Wtime();

  elem::HermitianEig(elem::LOWER, elem_mat, elem_w, elem_eigvecs); // only access lower half of H

  eigvecs.array = elem_eigvecs.Buffer();


  MPI_Barrier(MPI_COMM_WORLD);
  end = MPI_Wtime();

  for (int i=0; i<eigvals.size(); ++i)
    eigvals(i) = elem_w.Get(i, 0);
}

/*
template<>
void diagonalize<rokko::elemental>(rokko::distributed_matrix<rokko::elemental>& mat, Eigen::VectorXd& eigvals)
{
  //int m = 32;  // block_size

  elem::DistMatrix<double,elem::VR,elem::STAR> w(*(mat.get_elem_grid()) );

  double start, end;
  MPI_Barrier(MPI_COMM_WORLD);
  start = MPI_Wtime();

  elem::HermitianEig(elem::LOWER, mat.mat, w); // only access lower half of H

  MPI_Barrier(MPI_COMM_WORLD);
  end = MPI_Wtime();

  for (int i=0; i<eigvals.size(); ++i)
    eigvals(i) = w.Get(i, 0);
}
*/


} // namespace rokko

#endif // ROKKO_ELEMENTAL_DIAGONALIZE_HPP

