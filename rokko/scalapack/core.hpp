#ifndef ROKKO_SCALAPACK_CORE_HPP
#define ROKKO_SCALAPACK_CORE_HPP

#include <rokko/scalapack/scalapack.hpp>
#include <rokko/scalapack/diagonalize.hpp>
#include <iostream>

namespace rokko {
namespace scalapack {

class solver {
public:
  void initialize(int& argc, char**& argv) {}

  void finalize() {}

  void optimized_grid_size() {}

  template <typename MATRIX_MAJOR>
  void optimized_matrix_size(distributed_matrix<MATRIX_MAJOR>& mat) {
    // Determine mb, nb, lld, larray
    int mb = mat.get_m_global() / mat.get_nprow();
    if (mb == 0)  mb = 1;
    int nb = mat.get_n_global() / mat.get_npcol();
    if (nb == 0)  nb = 1;
    // Note: it should be that mb = nb in pdsyev.
    int tmp = std::min(mb, nb);
    mat.set_block_size(tmp, tmp);

    // Determine m_local, n_local from m_global, n_global, mb, nb
    mat.set_default_local_size();
    mat.set_default_lld();
    mat.set_default_length_array();
  }

  template<typename MATRIX_MAJOR>
  void diagonalize(distributed_matrix<MATRIX_MAJOR>& mat, localized_vector& eigvals,
                   distributed_matrix<MATRIX_MAJOR>& eigvecs) {
    rokko::scalapack::diagonalize(mat, eigvals, eigvecs);
  }
};

} // namespace sclapack
} // namespace rokko


#endif // ROKKO_SCALAPACK_CORE_HPP

