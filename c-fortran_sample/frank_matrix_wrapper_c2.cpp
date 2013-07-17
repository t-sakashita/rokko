#include <rokko/utility/frank_matrix.hpp>
#include "frank_matrix_wrapper_c2.h"

namespace rokko {
  void generate_distributed_matrix_col_major(void* mat){ 
      distributed_matrix<matrix_col_major>* mat_ = static_cast<distributed_matrix<matrix_col_major>*>(mat);
      frank_matrix::generate(*mat_);
  }
  void generate_distributed_matrix_row_major(void* mat){ 
      distributed_matrix<matrix_row_major>* mat_ = static_cast<distributed_matrix<matrix_row_major>*>(mat);
      frank_matrix::generate(*mat_);
  }
  void generate_localized_matrix_col_major(void* mat){ 
      localized_matrix<matrix_col_major>* mat_ = static_cast<localized_matrix<matrix_col_major>*>(mat);
      frank_matrix::generate(*mat_);
  }
  void generate_localized_matrix_row_major(void* mat){ 
      localized_matrix<matrix_row_major>* mat_ = static_cast<localized_matrix<matrix_row_major>*>(mat);
      frank_matrix::generate(*mat_);
  }
  

}
