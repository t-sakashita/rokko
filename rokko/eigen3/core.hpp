#ifndef ROKKO_EIGEN3_CORE_HPP
#define ROKKO_EIGEN3_CORE_HPP

#include <rokko/eigen3/diagonalize.hpp>
#include <iostream>

namespace rokko {
namespace eigen3 {

class solver {
public:
  void initialize(int& argc, char**& argv) {}

  void finalize() {}

  void diagonalize(localized_matrix<rokko::matrix_row_major>& mat, localized_vector& eigvals,
                   localized_matrix<rokko::matrix_row_major>& eigvecs, timer& timer_in);

  void diagonalize(localized_matrix<rokko::matrix_col_major>& mat, localized_vector& eigvals,
                   localized_matrix<rokko::matrix_col_major>& eigvecs, timer& timer_in);
};

inline void solver::diagonalize(localized_matrix<rokko::matrix_row_major>& mat, localized_vector& eigvals,
                                                          localized_matrix<rokko::matrix_row_major>& eigvecs, timer& timer_in) {
  rokko::eigen3::diagonalize(mat, eigvals, eigvecs, timer_in);
}

inline void solver::diagonalize(localized_matrix<rokko::matrix_col_major>& mat, localized_vector& eigvals,
                                                          localized_matrix<rokko::matrix_col_major>& eigvecs, timer& timer_in) {
  rokko::eigen3::diagonalize(mat, eigvals, eigvecs, timer_in);
}

} // namespace eigen3
} // namespace rokko


#endif // ROKKO_EIGEN3_CORE_HPP

