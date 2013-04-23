#include <mpi.h>

// Eigen3に関するヘッダファイル
#include <Eigen/Dense>

#include <iostream>
using namespace std;


#include <rokko/grid.hpp>
#include <rokko/distributed_matrix.hpp>
#include <rokko/scalapack.hpp>
#include <rokko/pblas.hpp>

#include <rokko/collective.hpp>
#include <rokko/frank_matrix.hpp>


int main(int argc, char* argv[])
{
  MPI_Init(&argc, &argv);
  MPI_Comm comm = MPI_COMM_WORLD;
  rokko::grid g(comm);

  int myrank, nprocs;
  MPI_Comm_rank(comm, &myrank);
  MPI_Comm_size(comm, &nprocs);

  int dim = 10;
  int root = 0;

  rokko::distributed_matrix frank_mat(dim, dim, g);
  //rokko::generate_frank_matrix_local(frank_mat);
  rokko::generate_frank_matrix_global(frank_mat);
  Eigen::MatrixXd frank_mat_global;
  rokko::gather(frank_mat, frank_mat_global, root);
  rokko::scatter(frank_mat, frank_mat_global, root);

  frank_mat.print();
  //rokko::print_matrix(mat_frank);
  if (myrank == root)
    cout << "global_mat_123:" << endl << frank_mat_global << endl;


  Eigen::VectorXd eigvals(dim);
  rokko::distributed_matrix eigvecs(dim, dim, g);
  rokko::scalapack::diagonalize(frank_mat, eigvals, eigvecs);


  Eigen::MatrixXd eigvecs_global;
  rokko::gather(eigvecs, eigvecs_global, root);
  eigvecs.print();
  //rokko::print_matrix(eigvecs);
  if (myrank == root) {
    std::cout << eigvecs_global << std::endl;
  }

  // 固有値の絶対値の降順に固有値(と対応する固有ベクトル)をソート
  // ソート後の固有値の添字をベクトルqに求める．
  double absmax;
  int qq;

  Eigen::VectorXi q(dim);

  Eigen::VectorXd eigval_sorted(dim);
  Eigen::MatrixXd eigvec_sorted(dim,dim);

  if(myrank == root) {
    // 固有値・固有ベクトルを絶対値の降順にソート
    for (int i=0; i<eigvals.size(); ++i) q[i] = i;
    for (int m=0; m<eigvals.size(); ++m) {
      absmax = abs(eigvals[q[m]]);
      for (int i=m+1; i<eigvecs_global.rows(); ++i) {
	if (absmax < abs(eigvals[q[i]])) {
	  absmax = eigvals[q[i]];
	  qq = q[m];
	  q[m] = q[i];
	  q[i] = qq;
	}
      }
      eigval_sorted(m) = eigvals(q[m]);
      eigvec_sorted.col(m) = eigvecs_global.col(q[m]);
    }

    cout << "Computed Eigenvalues= " << eigval_sorted.transpose() << endl;

    cout << "Eigenvector:" << endl << eigvec_sorted << endl << endl;
    cout << "Check the orthogonality of Eigenvectors:" << endl
      //<< eigvec_sorted * eigvec_sorted.transpose() << endl;   // Is it equal to indentity matrix?
         << eigvecs_global.transpose() * eigvecs_global << endl;   // Is it equal to indentity matrix?

    //cout << "residual := A x - lambda x = " << endl
    //	 << A_global_matrix * eigvec_sorted.col(0)  -  eigval_sorted(0) * eigvec_sorted.col(0) << endl;
    //cout << "Are all the following values equal to some eigenvalue = " << endl
    //	 << (A_global_matrix * eigvec_sorted.col(0)).array() / eigvec_sorted.col(0).array() << endl;
    //cout << "A_global_matrix=" << endl << A_global_matrix << endl;

  }

  //g.~grid();
  //rokko::scalapack::Finzalize();
  MPI_Finalize();
}
