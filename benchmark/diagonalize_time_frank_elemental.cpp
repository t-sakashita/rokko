#include <mpi.h>

// Eigen3に関するヘッダファイル
#include <Eigen/Dense>

#include <iostream>
#include <fstream>

using namespace std;

#include <rokko/elemental/elemental.hpp>
#include <rokko/elemental/core.hpp>
#include <rokko/elemental/grid.hpp>
#include <rokko/elemental/distributed_matrix.hpp>
#include <rokko/elemental/diagonalize.hpp>

#include <rokko/pblas/pblas.hpp>

#include <rokko/collective.hpp>
#include <rokko/utility/frank_matrix.hpp>
#include <rokko/utility/sort_eigenpairs.hpp>


int main (int argc, char *argv[])
{
  typedef rokko::elemental  solver;
  MPI_Init(&argc, &argv);
  rokko::initialize<solver>(argc, argv);
  MPI_Comm comm = MPI_COMM_WORLD;
  //rokko::grid<solver, rokko::grid_row_major<solver> > g(comm);
  rokko::grid<solver> g(comm);
  int myrank = g.myrank, nprocs = g.nprocs;

  ofstream ofs;
  if (myrank == 0) {
   ofs.open("elemental_time.txt");
   if (!ofs) {
     MPI_Abort(MPI_COMM_WORLD, 22) ;
   }
  }

  const int root = 0;
  const int dim = 100;

  rokko::distributed_matrix<solver> mat(dim, dim, g);
  rokko::distributed_matrix<solver> Z(dim, dim, g);

  rokko::generate_frank_matrix(mat);
  //Eigen::MatrixXd global_mat;
  //Eigen::MatrixXd global_mat(dim, dim);

  //rokko::scatter(frank_mat, mat_global, root);
  //rokko::gather(mat, global_mat, root);
  //mat.print();
  //if (myrank == root)
  //cout << "global_mat:" << endl << global_mat << endl;

  Eigen::VectorXd w(dim);

  //rokko::diagonalize<solver, rokko::grid_row_major>(mat, w, Z);
  //rokko::diagonalize<rokko::grid_row_major>(mat, w, Z);

  // Solve the problem
  double start, end;

  MPI_Barrier(MPI_COMM_WORLD);
  start = MPI_Wtime();
  rokko::diagonalize(mat, w, Z);
  MPI_Barrier(MPI_COMM_WORLD);
  end = MPI_Wtime();

  /*
  Z.print();
  // gather eigenvectors
  Eigen::MatrixXd global_eigvec;
  Eigen::MatrixXd eigvec_sorted(dim, dim);
  Eigen::VectorXd eigval_sorted(dim);
  rokko::gather(Z, global_eigvec, root);

  if (myrank == root) {
    cout << "eigvec:" << endl << global_eigvec << endl;
    rokko::sort_eigenpairs(w, global_eigvec, eigval_sorted, eigvec_sorted);
    cout.precision(20);
    //cout << "w=" << endl;
    //for (int i=0; i<dim; ++i) {
    //  cout << w[i] << " ";
    //}
    //cout << endl;
    cout << "Computed Eigenvalues= " << eigval_sorted.transpose() << endl;

    cout.precision(3);
    cout << "Check the orthogonality of Eigenvectors:" << endl
	 << eigvec_sorted * eigvec_sorted.transpose() << endl;   // Is it equal to indentity matrix?

    cout << "residual := A x - lambda x = " << endl
         << global_mat * eigvec_sorted.col(1)  -  eigval_sorted(1) * eigvec_sorted.col(1) << endl;
    cout << "Are all the following values equal to some eigenvalue = " << endl
      << (global_mat * eigvec_sorted.col(0)).array() / eigvec_sorted.col(0).array() << endl;
  }
  */

  double time;
  if (myrank == 0) {
    time = end - start;
    cout << "time=" << time << endl;
    ofs << "time=" << time << endl;
    //cout << "iter=" << iter << endl;
    //ofs << "iter=" << iter << endl;
  }

  rokko::finalize<solver>();
  MPI_Finalize();
  return 0;
}

