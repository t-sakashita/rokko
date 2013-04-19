#ifndef ROKKO_EIGEN_S_H
#define ROKKO_EIGEN_S_H

#include <mpi.h>

#define eigen_init_wrapper eigen_init_wrapper_
#define eigen_free_wrapper eigen_free_wrapper_
#define CSTAB_get_optdim cstab_get_optdim_  // オブジェクトファイルでは小文字に変換される
#define matrix_set matrix_set_
#define matrix_adjust_s matrix_adjust_s_
//#define eigen_s eigen_s_
#define ev_test_2D ev_test_2d_

extern "C" {
  void eigen_init_wrapper(int&, int&, int&);
  void eigen_free_wrapper(int&);
  void CSTAB_get_optdim(int&, int&, int&, int&, int&);
  void eigen_sx_(int&, double*, int&, double*, double*, double*, double*, int&, int&, int&);
}

/*
extern "C" struct
{
  int   my_col, size_of_col, mpi_comm_col,
    my_row, size_of_row, mpi_comm_row,
    p0_      ,q0_      , n_common,
    diag_0, diag_1;
} cycl2d_;
*/

namespace rokko {
namespace eigen_sx {

/*
int eigen_s_init()
{
  int size_of_col_local, size_of_row_local;
  int ndims = 2;
  eigen_init_wrapper(ndims, size_of_col_local, size_of_row_local);

  int NPROW = cycl2d_.size_of_row;
  int NPCOL = cycl2d_.size_of_col;
  int MYROW = cycl2d_.my_row;
  int MYCOL = cycl2d_.my_col;

  //cout << "NPROW=" << NPROW << "  NPCOL=" << NPCOL << endl;
  int nx = ((n-1)/NPROW+1);
  int i1 = 6, i2 = 16*4, i3 = 16*4*2, nm;
  CSTAB_get_optdim(nx, i1, i2, i3, nm);  // return an optimized (possiblly) leading dimension of local block-cyclic matrix to nm.
  int para_int = 0;   eigen_free_wrapper(para_int);
}
*/

template <class MATRIX, class VECTOR>
void diagonalize(MATRIX& mat, VECTOR& eigvals, MATRIX& eigvecs, bool flag_eigvecs)
{
  int m = 32;  // block_size

  int iflag;
  if (flag_eigvecs)
    iflag = 0;
  else
    iflag = 1;

  double* d = new double[mat.nme];
  if (d == NULL) {
    cerr << "failed to allocate d." << endl;
    MPI_Abort(MPI_COMM_WORLD, 3);
  }
  double* e = new double[2 * mat.nme];
  if (e == NULL) {
    cerr << "failed to allocate e." << endl;
    MPI_Abort(MPI_COMM_WORLD, 3);
  }

  double start, end;
  MPI_Barrier(MPI_COMM_WORLD);
  start = MPI_Wtime();

  eigen_sx_(mat.m_global, mat.array, mat.lld, &eigvals[0], eigvecs.array, d, e, mat.nme, m, iflag);

  MPI_Barrier(MPI_COMM_WORLD);
  end = MPI_Wtime();

  delete[] d;
  delete[] e;
}

} // namespace eigen_sx
} // namespace rokko

#endif // ROKKO_EIGEN_S_H
