/*****************************************************************************
*
* Rokko: Integrated Interface for libraries of eigenvalue decomposition
*
* Copyright (C) 2012-2015 Rokko Developers https://github.com/t-sakashita/rokko
*
* Distributed under the Boost Software License, Version 1.0. (See accompanying
* file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
*
*****************************************************************************/

#pragma once

#include <rokko/distributed_matrix.hpp>

namespace rokko {

template<typename T, typename MATRIX_MAJOR>
int check_orthogonality(rokko::distributed_matrix<T, MATRIX_MAJOR>& mat) {

  rokko::distributed_matrix<T, MATRIX_MAJOR> mat2(mat.get_mapping());
  product(1, mat, false, mat, true, 0, mat2);
  std::cout << "check_orthogonality:" << std::endl;
  mat2.print(std::cout);
  return 0;
}

} // namespace rokko
