#ifndef ROKKO_EIGEN_EXA_HPP
#define ROKKO_EIGEN_EXA_HPP

#define eigen_init eigen_init_
#define eigen_free eigen_free_
#define CSTAB_get_optdim cstab_get_optdim_  // オブジェクトファイルでは小文字に変換される
#define eigen_sx eigen_sx_

extern "C" {
  //void eigen_init();
  void eigen_init(const MPI_Fint&, char*);
  void eigen_free();
  void CSTAB_get_optdim(int&, int&, int&, int&, int&);
  void eigen_sx(int&, int&, double*, int&, double*, double*, int&, int&, int&);
}

#endif // ROKKO_EIGEN_EXA_HPP

