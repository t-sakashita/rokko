#ifndef ROKKO_LOCALIZED_MATRIX_HPP
#define ROKKO_LOCALIZED_MATRIX_HPP

#include <Eigen/Dense>

namespace rokko {

template<typename MATRIX_MAJOR = rokko::matrix_row_major>
struct eigen3_matrix_major;

template<>
struct eigen3_matrix_major<rokko::matrix_row_major> {
  static const int value= Eigen::RowMajor;
};

template<>
struct eigen3_matrix_major<rokko::matrix_col_major> {
  static const int value = Eigen::ColMajor;
};

template<typename MATRIX_MAJOR = rokko::matrix_row_major>
struct localized_matrix : public Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, eigen3_matrix_major<MATRIX_MAJOR>::value > {
public:
  localized_matrix() : Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, eigen3_matrix_major<MATRIX_MAJOR>::value >() {};
  localized_matrix(int rows, int cols) : Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, eigen3_matrix_major<MATRIX_MAJOR>::value >(rows, cols) {};
};

} // namespace rokko

#endif // ROKKO_LOCALIZED_MATRIX_HPP
