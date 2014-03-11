#ifndef ROKKO_LAPACK_CORE_HPP
#define ROKKO_LAPACK_CORE_HPP

#include <rokko/lapack/diagonalize.hpp>
//#include <rokko/lapack/diagonalize_dsyevd.hpp>
//#include <rokko/lapack/diagonalize_dsyevx.hpp>
#include <iostream>

namespace rokko {
namespace lapack {

struct dsyev {};
struct dsyevd {};
struct dsyevx {};

template<typename ROUTINE>
class solver {
public:
  void initialize(int& argc, char**& argv) {}

  void finalize() {}

  void diagonalize(localized_matrix<rokko::matrix_row_major>& mat, localized_vector& eigvals,
                   localized_matrix<rokko::matrix_row_major>& eigvecs, timer& timer_in);

  void diagonalize(localized_matrix<rokko::matrix_col_major>& mat, localized_vector& eigvals,
                   localized_matrix<rokko::matrix_col_major>& eigvecs, timer& timer_in);
};

template<>
inline void solver<rokko::lapack::dsyev>::diagonalize(localized_matrix<rokko::matrix_row_major>& mat, localized_vector& eigvals,
                                                          localized_matrix<rokko::matrix_row_major>& eigvecs, timer& timer_in) {
  rokko::lapack::diagonalize(mat, eigvals, eigvecs, timer_in);
}

template<>
inline void solver<rokko::lapack::dsyev>::diagonalize(localized_matrix<rokko::matrix_col_major>& mat, localized_vector& eigvals,
                                                          localized_matrix<rokko::matrix_col_major>& eigvecs, timer& timer_in) {
  rokko::lapack::diagonalize(mat, eigvals, eigvecs, timer_in);
}

/*
template<>
inline void solver<rokko::lapack::dsyevd>::diagonalize(localized_matrix<rokko::matrix_row_major>& mat, localized_vector& eigvals,
                                                           localized_matrix<rokko::matrix_row_major>& eigvecs, timer& timer_in) {
  rokko::lapack::diagonalize_d(mat, eigvals, eigvecs, timer_in);
}

template<>
inline void solver<rokko::lapack::dsyevd>::diagonalize(localized_matrix<rokko::matrix_col_major>& mat, localized_vector& eigvals,
                                                           localized_matrix<rokko::matrix_col_major>& eigvecs, timer& timer_in) {
  rokko::lapack::diagonalize_d(mat, eigvals, eigvecs, timer_in);
}

template<>
inline void solver<rokko::lapack::dsyevx>::diagonalize(localized_matrix<rokko::matrix_row_major>& mat, localized_vector& eigvals,
                                                           localized_matrix<rokko::matrix_row_major>& eigvecs, timer& timer_in) {
  rokko::lapack::diagonalize_x(mat, eigvals, eigvecs, timer_in);
}

template<>
inline void solver<rokko::lapack::dsyevx>::diagonalize(localized_matrix<rokko::matrix_col_major>& mat, localized_vector& eigvals,
                                                           localized_matrix<rokko::matrix_col_major>& eigvecs, timer& timer_in) {
  rokko::lapack::diagonalize_x(mat, eigvals, eigvecs, timer_in);
}
*/

} // namespace lapack
} // namespace rokko


#endif // ROKKO_LAPACK_CORE_HPP

