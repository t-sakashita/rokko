#include <mpi.h>

// Eigen3に関するヘッダファイル
#include <Eigen/Dense>

#include <iostream>

using namespace std;

#include <rokko/grid.hpp>
#include <rokko/distributed_matrix_eigenK.hpp>
#include <rokko/eigen_s.hpp>
#include <rokko/pblas.hpp>

#include <rokko/collective_eigenK.hpp>
#include <rokko/frank_matrix.hpp>

//typedef Eigen::MatrixXd matrix_type;


#undef __FUNCT__
#define __FUNCT__ "main"
int main (int argc, char *argv[])
{
  MPI_Init(&argc, &argv);
  MPI_Comm comm = MPI_COMM_WORLD;
  rokko::grid g(comm);
  int myrank = g.myrank, nprocs = g.nprocs;

  const int root = 0;
  const int dim = 10;

  rokko::distributed_matrix frank_mat(dim, dim, g);
  rokko::distributed_matrix Z(dim, dim, g);

  //rokko::generate_frank_matrix_local(frank_mat);
  rokko::generate_frank_matrix_global(frank_mat);
  Eigen::MatrixXd frank_mat_global;
  //Eigen::MatrixXd frank_mat_global(dim, dim);

  //rokko::scatter(frank_mat, frank_mat_global, root);
  rokko::gather(frank_mat, frank_mat_global, root);
  frank_mat.print();
  //rokko::print_matrix(mat_frank);
  if (myrank == root)
    cout << "global_mat_123:" << endl << frank_mat_global << endl;

  /*
  double* w = new double[n];
  if (w==NULL) {
    cerr << "error: w" << endl;
    return 1;
  }
  */

  Eigen::VectorXd w(dim);

  rokko::eigen_s::diagonalize(frank_mat, w, Z, true);

  Z.print();
  // gather of eigenvectors
  Eigen::MatrixXd eigvec_global;  //(dim, dim);
  Eigen::MatrixXd eigvec_sorted(dim, dim);
  Eigen::VectorXd eigval_sorted(dim);
  rokko::gather(Z, eigvec_global, root);
  //rokko::print_matrix(mat_frank);
  if (myrank == root) {
    cout << "eigvec:" << endl << eigvec_global << endl;
  }

  int* q = new int[dim];
  if (myrank == root) {
    // 固有値を（絶対値のではなく）昇順に並べる
    if (q==NULL) {
      cerr << "error: q" << endl;
      return 1;
    }

    double emax;
    for (int i=0; i<dim; ++i) q[i] = i;
    for (int k=0; k<dim; ++k) {
      emax = w[q[k]];
      for (int i=k+1; i<dim; ++i) {
	if (emax < w[q[i]]) {       // 昇順になっていないとき、交換
	  emax = w[q[i]];
	  int qq = q[k];
	  q[k] = q[i];
	  q[i] = qq;
	}
      }
      eigval_sorted(k) = w[q[k]];
      eigvec_sorted.col(k) = eigvec_global.col(q[k]);
      //eigvec_sorted.row(k) = eigvec_global.col(q[k]);
    }

    cout.precision(20);
    cout << "w=" << endl;
    for (int i=0; i<dim; ++i) {
      cout << w[i] << " ";
    }
    cout << endl;

    cout << "Computed Eigenvalues= " << eigval_sorted.transpose() << endl;

    cout.precision(3);
    cout << "Check the orthogonality of Eigenvectors:" << endl
	 << eigvec_sorted * eigvec_sorted.transpose() << endl;   // Is it equal to indentity matrix?
      //<< eigvec_global.transpose() * eigvec_global << endl;   // Is it equal to indentity matrix?

    Eigen::MatrixXd A_global_matrix = frank_mat_global;
    cout << "residual := A x - lambda x = " << endl
         << A_global_matrix * eigvec_sorted.col(1)  -  eigval_sorted(1) * eigvec_sorted.col(1) << endl;
    cout << "Are all the following values equal to some eigenvalue = " << endl
      << (A_global_matrix * eigvec_sorted.col(0)).array() / eigvec_sorted.col(0).array() << endl;
    cout << "A_global_matrix=" << endl << A_global_matrix << endl;
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

  delete [] w;
  delete [] q;
  */

  MPI_Finalize();
  return 0;
}

