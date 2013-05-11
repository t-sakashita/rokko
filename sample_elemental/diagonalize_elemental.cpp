#include <mpi.h>
#include <iostream>

#include <rokko/solver.hpp>
#include <rokko/grid.hpp>
#include <rokko/distributed_matrix.hpp>
#include <rokko/collective.hpp>

#include <rokko/utility/frank_matrix.hpp>
#include <rokko/utility/sort_eigenpairs.hpp>

int main(int argc, char *argv[]) {
  MPI_Init(&argc, &argv);
  typedef rokko::grid_col_major grid_major;
  //typedef rokko::matrix_row_major matrix_major;
  typedef rokko::matrix_col_major matrix_major;

  rokko::solver solver("elemental");
  solver.initialize(argc, argv);

  MPI_Comm comm = MPI_COMM_WORLD;
  rokko::grid<grid_major> g(comm);
  int myrank = g.myrank;

  const int root = 0;
  const int dim = 10;

  rokko::distributed_matrix<matrix_major> mat(dim, dim, g, solver);
  rokko::generate_frank_matrix(mat);
  mat.print();

  Eigen::MatrixXd global_mat;
  rokko::gather(mat, global_mat, root);
  if (myrank == root)
    std::cout << "global_mat:" << std::endl << global_mat << std::endl;


  Eigen::VectorXd w(dim);
  rokko::distributed_matrix<matrix_major> Z(dim, dim, g, solver);

  try {
    solver.diagonalize(mat, w, Z);
  }
  catch (const char *e) {
    std::cout << "Exception : " << e << std::endl;
    MPI_Abort(MPI_COMM_WORLD, 22);
  }

  // gather of eigenvectors
  Eigen::MatrixXd eigvec_global;
  Eigen::MatrixXd eigvec_sorted(dim, dim);
  Eigen::VectorXd eigval_sorted(dim);
  rokko::gather(Z, eigvec_global, root);
  Z.print();
  if (myrank == root) {
    std::cout << "eigvec:" << std::endl << eigvec_global << std::endl;
  }

  std::cout.precision(20);
  /*
  std::cout << "w=" << std::endl;
  for (int i=0; i<dim; ++i) {
    std::cout << w[i] << " ";
  }
  std::cout << std::endl;
  */

  if (myrank == root) {
    rokko::sort_eigenpairs(w, eigvec_global, eigval_sorted, eigvec_sorted);
    std::cout << "Computed Eigenvalues= " << eigval_sorted.transpose() << std::endl;

    std::cout.precision(3);
    std::cout << "Check the orthogonality of Eigenvectors:" << std::endl
	 << eigvec_sorted * eigvec_sorted.transpose() << std::endl;   // Is it equal to indentity matrix?
    //<< eigvec_global.transpose() * eigvec_global << std::endl;   // Is it equal to indentity matrix?

    std::cout << "residual := A x - lambda x = " << std::endl
         << global_mat * eigvec_sorted.col(1)  -  eigval_sorted(1) * eigvec_sorted.col(1) << std::endl;
    std::cout << "Are all the following values equal to some eigenvalue = " << std::endl
	 << (global_mat * eigvec_sorted.col(0)).array() / eigvec_sorted.col(0).array() << std::endl;
    //cout << "global_matrix=" << std::endl << global_matrix << std::endl;
  }


  /*
  double time;
  if (rank == 0) {
    time = end - start;
    std::cout << "time=" << time << std::endl;
    ofs << "time=" << time << std::endl;
    //cout << "iter=" << iter << std::endl;
    //ofs << "iter=" << iter << std::endl;
  }
  */

  solver.finalize();
  MPI_Finalize();
  return 0;
}
