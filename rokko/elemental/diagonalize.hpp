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

  MPI_Comm comm = mat.g.get_comm();
  if(mat.g.is_row_major()) throw "Elemental doesn't support grid_row_major.  Use Elemental with grid_col_major.";
  if(mat.is_row_major()) throw "Elemental doesn't support matrix_row_major.  Use Elemental with matrix_col_major.";

  elem::Grid elem_grid(comm, mat.get_nprow(), mat.get_npcol());
  elem::DistMatrix<double> elem_mat(mat.get_m_global(), mat.get_n_global(), 0, 0, mat.get_array_pointer(),
    mat.lld, elem_grid);
  //elem::DistMatrix<double> elem_eigvecs(mat.m_global, mat.n_global, elem_grid);
  elem::DistMatrix<double> elem_eigvecs(0, 0, elem_grid);

  elem::DistMatrix<double, elem::VR, elem::STAR> elem_w(elem_grid);

  double start, end;
  MPI_Barrier(MPI_COMM_WORLD);
  start = MPI_Wtime();

  elem::HermitianEig(elem::LOWER, elem_mat, elem_w, elem_eigvecs); // only access lower half of H

  MPI_Barrier(MPI_COMM_WORLD);
  end = MPI_Wtime();

  for (int i = 0; i < eigvals.size(); ++i)
    eigvals(i) = elem_w.Get(i, 0);

  double* result_mat = elem_eigvecs.Buffer();
  for(int local_i=0; local_i<mat.m_local; ++local_i) {
    for(int local_j=0; local_j<mat.n_local; ++local_j) {
      eigvecs.set_local(local_i, local_j, result_mat[local_j * mat.get_lld() + local_i]);
    }
  }

  /*
  elem_eigvecs.Print("elem_eigvecs=");
  for (int proc=0; proc<mat.nprocs; ++proc) {
    if (proc == mat.myrank) {
      printf("Rank = %d  myrow=%d mycol=%d\n", mat.myrank, mat.myrow, mat.mycol);
      elem::Matrix<double> Y = elem_eigvecs.Matrix();
      std::cout << "Matrix():" << std::endl;
      Y.Print();
      printf("Local Matrix after diagonalize:\n");
      for (int local_i=0; local_i<mat.m_local; ++local_i) {
        for (int local_j=0; local_j<mat.n_local; ++local_j) {
          printf("%e ", (eigvecs.array)[local_i + local_j * mat.lld]);
        }
        printf("\n");
      }
    }
    MPI_Barrier(MPI_COMM_WORLD);
  }
  */

}

/*
template<typename MATRIX_MAJOR>
void diagonalize(rokko::distributed_matrix<MATRIX_MAJOR>& mat, Eigen::VectorXd& eigvals)
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

} // namespace elemental
} // namespace rokko

#endif // ROKKO_ELEMENTAL_DIAGONALIZE_HPP
