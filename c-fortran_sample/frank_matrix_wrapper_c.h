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
#ifdef __cplusplus
namespace rokko {
  extern "C"{
#endif
    
    void frank_generate_distributed_matrix_col_major(void*);
    void frank_generate_distributed_matrix_row_major(void*);
    void frank_generate_localized_matrix_col_major(void*);
    void frank_generate_localized_matrix_row_major(void*);
    void generate_distributed_matrix_col_major(void* mat, double (*func)(int i, int j));
#ifdef __cplusplus
  }
}
#endif
