#include <mpi.h>

// Eigen3に関するヘッダファイル
#include <Eigen/Dense>

#include <iostream>
using namespace std;

#define TIMER
#define DETAIL

#define eigen_init_wrapper eigen_init_wrapper_
#define eigen_free_wrapper eigen_free_wrapper_
#define CSTAB_get_optdim cstab_get_optdim_  // オブジェクトファイルでは小文字に変換する
#define matrix_set matrix_set_
#define matrix_adjust_s matrix_adjust_s_
#define eigen_s eigen_s_
#define ev_test_2D ev_test_2d_

//extern "C" void eigen_init(int&);
//extern "C" void eigen_free(int&);

extern "C" void eigen_init_wrapper(int&, int&, int&);
extern "C" void eigen_free_wrapper(int&);

extern "C" void CSTAB_get_optdim(int&, int&, int&, int&, int&);
extern "C" void matrix_set(int&, double*);
extern "C" void matrix_adjust_s(int&, double*, double*, int&);
extern "C" void eigen_s(int&, double*, int&, double*, double*, int&, int&, int&);
extern "C" void ev_test_2D(int&, double*, int&, double*, double*, int&);

extern "C" struct
{
  int   my_col, size_of_col, mpi_comm_col,
    my_row, size_of_row, mpi_comm_row,
    p0_      ,q0_      , n_common,
    diag_0, diag_1;
} cycl2d_;



#include <rokko/grid.hpp>
#include <rokko/distributed_matrix_eigenK.hpp>
#include <rokko/scalapack.hpp>
#include <rokko/pblas.hpp>

#include <rokko/collective_eigenK.hpp>
#include <rokko/frank_matrix.hpp>


typedef Eigen::MatrixXd matrix_type;


#undef __FUNCT__
#define __FUNCT__ "main"
int main (int argc, char *argv[])
{
  int para_int;

  MPI_Init(&argc, &argv);

  MPI_Comm comm = MPI_COMM_WORLD;
  rokko::grid g(comm);

  int myrank = g.myrank, nprocs = g.nprocs;

  int n, dim;
  n = dim = 10;
  int root = 0;

  int m = 32;

  int size_of_col_local, size_of_row_local;
  para_int = 2;   eigen_init_wrapper(para_int, size_of_col_local, size_of_row_local);

  int NPROW = cycl2d_.size_of_row;
  int NPCOL = cycl2d_.size_of_col;
  int MYROW = cycl2d_.my_row;
  int MYCOL = cycl2d_.my_col;

  cout << "NPROW=" << NPROW << "  NPCOL=" << NPCOL << endl;
  int nx = ((n-1)/NPROW+1);

  // main.fでは整数型変数の宣言はなかった。
  int i1 = 6, i2 = 16*4, i3 = 16*4*2, nm;
  CSTAB_get_optdim(nx, i1, i2, i3, nm);
  para_int = 0;   eigen_free_wrapper(para_int);

  // nm = 8; // 上書き
  //int NB = 3;
  int NB  = 64+32;
  int nmz = ((n-1)/NPROW+1);
  nmz = ((nmz-1)/NB+1)*NB+1;
  int nmw = ((n-1)/NPCOL+1);
  nmw = ((nmw-1)/NB+1)*NB+1;

  cout << "nm=" << nm << endl;
  cout << "nmz=" << nmz << endl;
  cout << "nmw=" << nmw << endl;

  int larray = std::max(nmz,nm)*nmw;
  cout << "larray=" << larray << endl;


  rokko::distributed_matrix frank_mat(dim, dim, g, nm, larray);
  rokko::distributed_matrix Z(dim, dim, g, nm, larray);

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

  double* w = new double[n];
  if (w==NULL) {
    cerr << "error: w" << endl;
    return 1;
  }

  if (myrank == 0) {
    cout << "n=" << n << endl;
    cout << "nm=" << nm << endl;
  }

  int ZERO = 0, ONE = 1;
  // index 1 -> 0
  double start, end;
  MPI_Barrier(MPI_COMM_WORLD);
  start = MPI_Wtime();

  eigen_s(n, frank_mat.array, nm, w, Z.array, nm, m, ZERO);

  MPI_Barrier(MPI_COMM_WORLD);
  end = MPI_Wtime();

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

  int* q = new int[n];
  if (myrank == root) {
    // 固有値を（絶対値のではなく）昇順に並べる
    if (q==NULL) {
      cerr << "error: q" << endl;
      return 1;
    }

    double emax;
    for (int i=0; i<n; ++i) q[i] = i;
    for (int k=0; k<n; ++k) {
      emax = w[q[k]];
      for (int i=k+1; i<n; ++i) {
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
    for (int i=0; i<n; ++i) {
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

