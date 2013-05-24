#ifndef ROKKO_ELEMENTAL_CORE_HPP
#define ROKKO_ELEMENTAL_CORE_HPP

#include <elemental.hpp>
#include <rokko/elemental/diagonalize.hpp>

namespace rokko {
namespace elemental {

class solver {
public:
  void initialize(int& argc, char**& argv) { elem::Initialize(argc, argv); }

  void finalize() { elem::Finalize(); }

  void optimized_grid_size() {}

  template <typename MATRIX_MAJOR>
  void optimized_matrix_size(distributed_matrix<MATRIX_MAJOR>& mat) {
    mat.set_block_size(1, 1); // set mb = nb = 1

    // Determine m_local, n_local from m_global, n_global, mb, nb
    int m_local = mat.calculate_row_size();
    int n_local = mat.calculate_col_size();
    mat.set_local_size(m_local, n_local);
    mat.set_default_lld();
    mat.set_default_length_array();
  }

  template<typename MATRIX_MAJOR>
  void diagonalize(distributed_matrix<MATRIX_MAJOR>& mat, localized_vector& eigvals,
                   distributed_matrix<MATRIX_MAJOR>& eigvecs, timer& timer_in) {
    rokko::elemental::diagonalize(mat, eigvals, eigvecs, timer_in);
  }
};

} // namespace elemental
} // namespace rokko

#endif // ROKKO_ELEMENTAL_CORE_HPP
