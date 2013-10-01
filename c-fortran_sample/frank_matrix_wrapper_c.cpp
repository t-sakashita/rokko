/*****************************************************************************
*
* Rokko: Integrated Interface for libraries of eigenvalue decomposition
*
* Copyright (C) 2012-2013 by Tatsuya Sakashita <t-sakashita@issp.u-tokyo.ac.jp>,
*                            Synge Todo <wistaria@comp-phys.org>,
*                            Tsuyoshi Okubo <t-okubo@issp.u-tokyo.ac.jp>
*    
* Distributed under the Boost Software License, Version 1.0. (See accompanying
* file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
*
*****************************************************************************/
#include <rokko/utility/frank_matrix.hpp>
#include "frank_matrix_wrapper_c.h"

namespace rokko {
  void frank_generate_distributed_matrix_col_major(void* mat){ 
      distributed_matrix<matrix_col_major>* mat_ = static_cast<distributed_matrix<matrix_col_major>*>(mat);
      frank_matrix::generate(*mat_);
  }
  void frank_generate_distributed_matrix_row_major(void* mat){ 
      distributed_matrix<matrix_row_major>* mat_ = static_cast<distributed_matrix<matrix_row_major>*>(mat);
      frank_matrix::generate(*mat_);
  }
  void frank_generate_localized_matrix_col_major(void* mat){ 
      localized_matrix<matrix_col_major>* mat_ = static_cast<localized_matrix<matrix_col_major>*>(mat);
      frank_matrix::generate(*mat_);
  }
  void frank_generate_localized_matrix_row_major(void* mat){ 
      localized_matrix<matrix_row_major>* mat_ = static_cast<localized_matrix<matrix_row_major>*>(mat);
      frank_matrix::generate(*mat_);
  }
}
