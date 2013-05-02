#include <mpi.h>

#include <Eigen/Dense>

#include <iostream>

using namespace std;

#include <mpi.h>
#include <rokko/elemental/core.hpp>
#include <rokko/grid.hpp>
#include <rokko/distributed_matrix.hpp>
#include <rokko/elemental/diagonalize.hpp>
#include <rokko/collective.hpp>

#include <rokko/utility/frank_matrix.hpp>
#include <rokko/utility/sort_eigenpairs.hpp>



int main(int argc, char *argv[])
{
  typedef rokko::grid_col_major grid_major;
  typedef rokko::matrix_col_major matrix_major;

  MPI_Init(&argc, &argv);
  //rokko::solver solver("elemental");
  rokko::solver_elemental solver;
  solver.initialize(argc, argv);


  MPI_Comm comm = MPI_COMM_WORLD;
  //rokko::grid  g(comm); //, solver);
  rokko::grid<grid_major>  g(comm); //, solver);
  int myrank = g.myrank, nprocs = g.nprocs;

  const int root = 0;
  const int dim = 10;

  rokko::distributed_matrix<matrix_major> mat(dim, dim, g); //, solver);

  rokko::generate_frank_matrix(mat);

  //mat.mat.Print("elemental_matrix");
  Eigen::MatrixXd global_mat;

  rokko::gather(mat, global_mat, root);
  mat.print();
  //rokko::print_matrix(mat);
  if (myrank == root)
    cout << "global_mat:" << endl << global_mat << endl;


  Eigen::VectorXd w(dim);
  rokko::distributed_matrix<matrix_major> Z(dim, dim, g); //, true);

  solver.diagonalize(mat, w, Z);

  // gather of eigenvectors
  Eigen::MatrixXd eigvec_global;  //(dim, dim);
  Eigen::MatrixXd eigvec_sorted(dim, dim);
  Eigen::VectorXd eigval_sorted(dim);
  rokko::gather(Z, eigvec_global, root);
  Z.print();
  if (myrank == root) {
    cout << "eigvec:" << endl << eigvec_global << endl;
  }

  cout.precision(20);
  cout << "w=" << endl;
  for (int i=0; i<dim; ++i) {
    cout << w[i] << " ";
  }
  cout << endl;


  if (myrank == root) {
    rokko::sort_eigenpairs(w, eigvec_global, eigval_sorted, eigvec_sorted);
    cout << "Computed Eigenvalues= " << eigval_sorted.transpose() << endl;

    cout.precision(3);
    cout << "Check the orthogonality of Eigenvectors:" << endl
	 << eigvec_sorted * eigvec_sorted.transpose() << endl;   // Is it equal to indentity matrix?
    //<< eigvec_global.transpose() * eigvec_global << endl;   // Is it equal to indentity matrix?

    cout << "residual := A x - lambda x = " << endl
         << global_mat * eigvec_sorted.col(1)  -  eigval_sorted(1) * eigvec_sorted.col(1) << endl;
    cout << "Are all the following values equal to some eigenvalue = " << endl
	 << (global_mat * eigvec_sorted.col(0)).array() / eigvec_sorted.col(0).array() << endl;
    //cout << "global_matrix=" << endl << global_matrix << endl;
  }


  /*
  double time;
  if (rank == 0) {
    time = end - start;
    cout << "time=" << time << endl;
    ofs << "time=" << time << endl;
    //cout << "iter=" << iter << endl;
    //ofs << "iter=" << iter << endl;
  }
  */

  solver.finalize();
  MPI_Finalize();
  return 0;
}

