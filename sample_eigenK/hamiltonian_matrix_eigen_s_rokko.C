#include <iostream>
#include <Eigen/Dense>
#include <mpi.h>

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

/*
extern "C" struct
{
  int   nprocs, myrank;
} usempi_;
*/

extern "C" struct
{
  int   my_col, size_of_col, mpi_comm_col,
    my_row, size_of_row, mpi_comm_row,
    p0_      ,q0_      , n_common,
    diag_0, diag_1;
} cycl2d_;

using namespace std;

typedef Eigen::MatrixXd matrix_type;


#undef __FUNCT__
#define __FUNCT__ "main"
int main (int argc, char *argv[])
{
  int rank, size;
  int para_int;

  MPI_Init(&argc, &argv);/* starts MPI */
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);/* get current process id */
  MPI_Comm_size(MPI_COMM_WORLD, &size);/* get number of processes */

  printf( "Hello world from process %d of %d\n", rank, size );

  std::ifstream  ifs(argv[1]);
  alps::Parameters params(ifs);
  //params["L"] = 7;
  cout << "L=" << params["L"] << endl;
  barista::Hamiltonian<> hamiltonian(params);
  matrix_type matrix(hamiltonian.dimension(), hamiltonian.dimension());
  hamiltonian.fill<double>(matrix);
  //std::cout << matrix << std::endl;
  int n = hamiltonian.dimension();

  std::ofstream ofs;
  if (rank == 0) {
    ofs.open("eigen_s_time.txt");
    if (!ofs) {
      MPI_Finalize() ;
      return -1;
    }
  }

  double* a = new double[n*n];
  if (a==NULL) {
    cerr << "error: a" << endl;
    return 1;
  }
  for (int i=0; i<n; ++i) {
    for (int j=0; j<n; ++j)
      a[i*n+j] = matrix(i,j);
  }


  int m = 32;

  int size_of_col_local, size_of_row_local;
  para_int = 2;   eigen_init_wrapper(para_int, size_of_col_local, size_of_row_local);

  int NPROW = cycl2d_.size_of_col;
  int NPCOL = cycl2d_.size_of_row;
  int nx = ((n-1)/NPROW+1);

  // main.fでは整数型変数の宣言はなかった。
  int i1 = 6, i2 = 16*4, i3 = 16*4*2, nm;
  CSTAB_get_optdim(nx, i1, i2, i3, nm);
  para_int = 0;   eigen_free_wrapper(para_int);

  int NB  = 64+32;
  int nmz = ((n-1)/NPROW+1);
  nmz = ((nmz-1)/NB+1)*NB+1;
  int nmw = ((n-1)/NPCOL+1);
  nmw = ((nmw-1)/NB+1)*NB+1;

  cout << "nmz=" << nmz << endl;
  cout << "nm=" << nm << endl;
  cout << "nmw=" << nmw << endl;

  int larray = std::max(nmz,nm)*nmw;
  cout << "larray=" << larray << endl;

  double* b = new double[larray];
  if (b==NULL) {
    cerr << "error: b" << endl;
    return 1;
  }

  double* z = new double[larray];
  if (z==NULL) {
    cerr << "error: z" << endl;
    return 1;
  }
  double* w = new double[n];
  if (w==NULL) {
    cerr << "error: w" << endl;
    return 1;
  }

  //matrix_set(n, a);
  matrix_adjust_s(n, a, b, nm);

  if (rank == 0) {
    cout << "n=" << n << endl;
    //    for (int i=0; i<n; ++i) {
    //      for (int j=0; j<n; ++j)
    //	std::cout << a[i*n+j] << " ";
    //      cout << endl;
    //}
    cout << "nm=" << nm << endl;
  }

  int zero = 0;
  // index 1 -> 0
  double start, end;
  MPI_Barrier(MPI_COMM_WORLD);
  start = MPI_Wtime();
  eigen_s(n, b, nm, w, z, nm, m, zero);
  MPI_Barrier(MPI_COMM_WORLD);
  end = MPI_Wtime();

  int* q = new int[n];
  if (rank == 0) {
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
	if (emax > w[q[i]]) {       // 昇順になっていないとき、交換
	  emax = w[q[i]];
	  int qq = q[k];
	  q[k] = q[i];
	  q[i] = qq;
	}
      }
    }
  }

  cout.precision(20);
  if (rank == 0) {
    cout << "w=" << endl;
    for (int i=0; i<n; ++i) {
      cout << w[i] << " ";
    }
    cout << endl;
  }

  double time;
  if (rank == 0) {
    time = end - start;
    cout << "time=" << time << endl;
    ofs << "time=" << time << endl;
    //cout << "iter=" << iter << endl;
    //ofs << "iter=" << iter << endl;
  }

  delete [] a;
  delete [] b;
  delete [] z;
  delete [] w;
  delete [] q;

  MPI_Finalize();
  return 0;
}

