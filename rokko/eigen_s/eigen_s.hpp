#ifndef ROKKO_EIGEN_S_H
#define ROKKO_EIGEN_S_H


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


extern "C" struct
{
  int   my_col, size_of_col, mpi_comm_col,
    my_row, size_of_row, mpi_comm_row,
    p0_      ,q0_      , n_common,
    diag_0, diag_1;
} cycl2d_;


#endif // ROKKO_EIGEN_S_H
