#include <mpi.h>

// Eigen3に関するヘッダファイル
#include <Eigen/Dense>

#include <iostream>

using namespace std;

#include <rokko/eigen_s/eigen_s.hpp>
#include <rokko/eigen_s/core.hpp>
#include <rokko/eigen_s/grid.hpp>
#include <rokko/eigen_s/distributed_matrix.hpp>
#include <rokko/eigen_s/diagonalize.hpp>

//#include <rokko/pblas/pblas.hpp>

#include <rokko/collective.hpp>
#include <rokko/utility/frank_matrix.hpp>
#include <rokko/utility/sort_eigenpairs.hpp>


int main (int argc, char *argv[])
{
  typedef rokko::eigen_s solver;
  MPI_Init(&argc, &argv);
  rokko::Initialize<solver>(argc, argv);
  MPI_Comm comm = MPI_COMM_WORLD;
  rokko::grid<solver> g(comm);
  int myrank = g.myrank, nprocs = g.nprocs;

  const int root = 0;
  const int dim = 10;

  rokko::distributed_matrix<solver> mat(dim, dim, g);
  rokko::distributed_matrix<solver> Z(dim, dim, g);

  rokko::generate_frank_matrix(mat);
  //rokko::generate_frank_matrix_global(mat);
  Eigen::MatrixXd global_mat;
  //Eigen::MatrixXd global_mat(dim, dim);

  //rokko::scatter(frank_mat, mat_global, root);
  rokko::gather(mat, global_mat, root);
  mat.print();
  //rokko::print_matrix(mat_frank);
  if (myrank == root)
    cout << "global_mat:" << endl << global_mat << endl;

  Eigen::VectorXd w(dim);
  rokko::diagonalize<solver>(mat, w, Z, true);

  Z.print();
  // gather of eigenvectors
  Eigen::MatrixXd global_eigvec;
  Eigen::MatrixXd eigvec_sorted(dim, dim);
  Eigen::VectorXd eigval_sorted(dim);
  rokko::gather(Z, global_eigvec, root);
  rokko::print_matrix(mat);

  if (myrank == root) {
    cout << "eigvec:" << endl << global_eigvec << endl;
    rokko::sort_eigenpairs(w, global_eigvec, eigval_sorted, eigvec_sorted);
    //cout.precision(20);
    //cout << "w=" << endl;
    //for (int i=0; i<dim; ++i) {
    //  cout << w[i] << " ";
    //}
    cout << endl;

    cout << "Computed Eigenvalues= " << eigval_sorted.transpose() << endl;

    cout.precision(3);
    cout << "Check the orthogonality of Eigenvectors:" << endl
	 << eigvec_sorted * eigvec_sorted.transpose() << endl;   // Is it equal to indentity matrix?

    cout << "residual := A x - lambda x = " << endl
         << global_mat * eigvec_sorted.col(1)  -  eigval_sorted(1) * eigvec_sorted.col(1) << endl;
    cout << "Are all the following values equal to some eigenvalue = " << endl
      << (global_mat * eigvec_sorted.col(0)).array() / eigvec_sorted.col(0).array() << endl;
  }

  MPI_Finalize();
  return 0;
}

