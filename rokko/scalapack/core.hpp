#ifndef ROKKO_SCALAPACK_CORE_HPP
#define ROKKO_SCALAPACK_CORE_HPP

#include <rokko/scalapack/scalapack.hpp>
#include <rokko/scalapack/diagonalize.hpp>
#include <iostream>

namespace rokko {
namespace scalapack {

class solver {
public:
  void initialize(int& argc, char**& argv) { }

  void finalize() { }

  void optimized_grid_size() {}

  //void optimized_matrix_size(int dim, int nprow, int npcol, int& mb, int& nb, int& lld, int& len_array) {
  template <typename MATRIX_MAJOR>
  void optimized_matrix_size(distributed_matrix<MATRIX_MAJOR>& mat) {
    // ローカル行列の形状を指定
    //m_global = n_global = dim;
    mat.mb = mat.m_global / mat.nprow;
    if (mat.mb == 0) mat.mb = 1;
    mat.nb = mat.n_global / mat.npcol;
    if (mat.nb == 0) mat.nb = 1;
    // mbとnbを最小値にそろえる．（注意：pdsyevではmb=nbでなければならない．）
    mat.mb = std::min(mat.mb, mat.nb);
    mat.nb = mat.mb;

    // determine m_local, n_local from m_global, n_global, mb, nb
    mat.m_local = mat.get_row_size();
    mat.n_local = mat.get_col_size();
    mat.lld = mat.get_lld();
    mat.length_array = mat.m_local * mat.n_local;

    /*
    for (int proc=0; proc<nprocs; ++proc) {
      if (proc == myrank) {
	std::cout << "proc=" << proc << std::endl;
	std::cout << "  mb=" << mb << "  nb=" << nb << std::endl;
	std::cout << "  mA=" << m_local << "  nprow=" << nprow << std::endl;
	std::cout << "  nA=" << n_local << "  npcol=" << npcol << std::endl;
	std::cout << " m_local=" << m_local << " n_local=" << n_local << std::endl;
        std::cout << " myrow=" << myrow << " mycol=" << mycol << std::endl;
      }
      MPI_Barrier(MPI_COMM_WORLD);
    }
    */


  }

  template<typename MATRIX_MAJOR>
  void diagonalize(rokko::distributed_matrix<MATRIX_MAJOR>& mat, Eigen::VectorXd& eigvals,
                   rokko::distributed_matrix<MATRIX_MAJOR>& eigvecs) {
    rokko::scalapack::diagonalize(mat, eigvals, eigvecs);
  }
};

} // namespace sclapack
} // namespace rokko


#endif // ROKKO_SCALAPACK_CORE_HPP

