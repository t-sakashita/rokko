#ifndef ROKKO_SCALAPACK_CORE_HPP
#define ROKKO_SCALAPACK_CORE_HPP

#include <rokko/scalapack/scalapack.hpp>
#include <rokko/scalapack/diagonalize.hpp>

namespace rokko {
namespace scalapack {

class solver {
public:
  void initialize(int& argc, char**& argv) { }

  void finalize() { }

  void optimized_grid_size() {}

  void optimized_matrix_size(int dim, int nprow, int npcol, int& mb, int& nb, int& lld, int& len_array) {
    // ローカル行列の形状を指定
    m_global = n_global = dim;
    mb = m_global / nprow;
    if (mb == 0) mb = 1;
    nb = n_global / npcol;
    if (nb == 0) nb = 1;
    // mbとnbを最小値にそろえる．（注意：pdsyevではmb=nbでなければならない．）
    mb = min(mb, nb);
    nb = mb;

    lld = get_lld();

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

    len_array = m_local * n_local;
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

