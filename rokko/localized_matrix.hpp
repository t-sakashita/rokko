#ifndef ROKKO_LOCALIZED_H
#define ROKKO_LOCALIZED_H

#include <iostream>
#include <cstdlib>
#include <boost/type_traits/is_same.hpp>

#include <Eigen/Core>
#include <Eigen/Dense>


namespace rokko {

class solver;

/*
struct matrix_row_major {};

struct matrix_col_major {};
*/

template<typename MATRIX_MAJOR = rokko::matrix_row_major>
struct eigen3_matrix_major;


struct eigen3_matrix_major<rokko::matrix_row_major>
{
  static const int value= Eigen::RowMajor;
};

template<>
struct eigen3_matrix_major<rokko::matrix_col_major>
{
  static const int value = Eigen::ColMajor;
};

template<typename MATRIX_MAJOR = rokko::matrix_row_major>
struct localized_matrix : public Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, eigen3_matrix_major<MATRIX_MAJOR>::value > {
public:
  localized_matrix(int rows, int cols) : Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, eigen3_matrix_major<MATRIX_MAJOR>::value >(rows, cols) {};
};


} // namespace rokko

#endif // ROKKO_LOCALIZED_H
