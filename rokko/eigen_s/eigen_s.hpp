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
  void matrix_set(int&, double*);
  void matrix_adjust_s(int&, double*, double*, int&);
  void eigen_s_(int&, double*, int&, double*, double*, int&, int&, int&);
  void ev_test_2D(int&, double*, int&, double*, double*, int&);
}


namespace rokko {
namespace eigen_s {

template <class MATRIX, class VECTOR>
void diagonalize(MATRIX& mat, VECTOR& eigvals, MATRIX& eigvecs, bool flag_eigvecs)
{
  int m = 32;  // block_size

  int iflag;
  if (flag_eigvecs)
    iflag = 0;
  else
    iflag = 1;

  double start, end;
  MPI_Barrier(MPI_COMM_WORLD);
  start = MPI_Wtime();

  eigen_s_(mat.m_global, mat.array, mat.lld, &eigvals[0], eigvecs.array, mat.lld, m, iflag);

  MPI_Barrier(MPI_COMM_WORLD);
  end = MPI_Wtime();
}

} // namespace eigen_s
} // namespace rokko

#endif // ROKKO_EIGEN_S_H
