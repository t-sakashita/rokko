#ifndef ROKKO_ELEMENTAL_DIAGONALIZE_HPP
#define ROKKO_ELEMENTAL_DIAGONALIZE_HPP

#include <mpi.h>
#include <elemental.hpp>
#include <rokko/grid.hpp>
#include <rokko/distributed_matrix.hpp>

namespace rokko {
namespace elemental {
    
template<typename MATRIX_MAJOR>
void diagonalize(rokko::distributed_matrix<MATRIX_MAJOR>& mat, Eigen::VectorXd& eigvals,
  rokko::distributed_matrix<MATRIX_MAJOR>& eigvecs) {

  MPI_Comm comm = MPI_COMM_WORLD;
  elem::Grid elem_grid(comm, mat.nprow, mat.npcol);
  elem::DistMatrix<double> elem_mat(mat.m_global, mat.n_global, 0, 0, mat.get_array(),
    mat.lld, elem_grid);
  elem::DistMatrix<double> elem_eigvecs(mat.m_global, mat.n_global, elem_grid);
  elem::DistMatrix<double, elem::VR, elem::STAR> elem_w(elem_grid);

  double start, end;
  MPI_Barrier(MPI_COMM_WORLD);
  start = MPI_Wtime();

  elem::HermitianEig(elem::LOWER, elem_mat, elem_w, elem_eigvecs); // only access lower half of H
  eigvecs.array = elem_eigvecs.Buffer();

  MPI_Barrier(MPI_COMM_WORLD);
  end = MPI_Wtime();

  for (int i = 0; i < eigvals.size(); ++i)
    eigvals(i) = elem_w.Get(i, 0);
}

} // namespace elemental
} // namespace rokko

#endif // ROKKO_ELEMENTAL_DIAGONALIZE_HPP
