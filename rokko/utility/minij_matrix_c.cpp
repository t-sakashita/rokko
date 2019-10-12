/*****************************************************************************
*
* Rokko: Integrated Interface for libraries of eigenvalue decomposition
*
* Copyright (C) 2012-2019 Rokko Developers https://github.com/t-sakashita/rokko
*
* Distributed under the Boost Software License, Version 1.0. (See accompanying
* file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
*
*****************************************************************************/

#include <rokko/utility/minij_matrix.hpp>
#include <rokko/utility/minij_matrix.h>

#if defined(ROKKO_HAVE_PARALLEL_DENSE_SOLVER)
void rokko_minij_matrix_generate_distributed_matrix(rokko_distributed_matrix matrix) {
  if (matrix.major == rokko_matrix_col_major)
    rokko::minij_matrix::generate(*static_cast<rokko::distributed_matrix<double, rokko::matrix_col_major>*>(matrix.ptr));
  else
    rokko::minij_matrix::generate(*static_cast<rokko::distributed_matrix<double, rokko::matrix_row_major>*>(matrix.ptr));
}
#endif

void rokko_minij_matrix_generate_eigen_matrix(rokko_eigen_matrix matrix) {
  if (matrix.major == rokko_matrix_col_major)
    rokko::minij_matrix::generate(*static_cast<Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::ColMajor>*>(matrix.ptr));
  else
    rokko::minij_matrix::generate(*static_cast<Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>*>(matrix.ptr));
}
