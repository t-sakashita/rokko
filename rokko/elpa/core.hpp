#ifndef ROKKO_ELPA_CORE_HPP
#define ROKKO_ELPA_CORE_HPP

#include <rokko/elpa/elpa.hpp>
#include <rokko/elpa/diagonalize.hpp>
#include <iostream>

namespace rokko {
namespace elpa {

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
    //mat.set_block_size(tmp, tmp);
    mat.set_block_size(1, 1);

    // Determine m_local, n_local from m_global, n_global, mb, nb
    mat.set_default_local_size();
    mat.set_default_lld();
    mat.set_default_length_array();
  }

  void diagonalize(distributed_matrix<rokko::matrix_row_major>& mat, localized_vector& eigvals,
                   distributed_matrix<rokko::matrix_row_major>& eigvecs, timer& timer_in);

  void diagonalize(distributed_matrix<rokko::matrix_col_major>& mat, localized_vector& eigvals,
                   distributed_matrix<rokko::matrix_col_major>& eigvecs, timer& timer_in);
};

/*
template<typename MATRIX_MAJOR>
inline void solver<rokko::elpa::pdsyev>::diagonalize(distributed_matrix<MATRIX_MAJOR>& mat, localized_vector& eigvals,
                                                   distributed_matrix<MATRIX_MAJOR>& eigvecs, timer& timer_in) {
  rokko::elpa::diagonalize(mat, eigvals, eigvecs, timer_in);
}
*/

inline void solver::diagonalize(distributed_matrix<rokko::matrix_row_major>& mat, localized_vector& eigvals,
                                                          distributed_matrix<rokko::matrix_row_major>& eigvecs, timer& timer_in) {
  rokko::elpa::diagonalize(mat, eigvals, eigvecs, timer_in);
}


inline void solver::diagonalize(distributed_matrix<rokko::matrix_col_major>& mat, localized_vector& eigvals,
                                                          distributed_matrix<rokko::matrix_col_major>& eigvecs, timer& timer_in) {
  rokko::elpa::diagonalize(mat, eigvals, eigvecs, timer_in);
}

} // namespace elpa
} // namespace rokko


#endif // ROKKO_ELPA_CORE_HPP

