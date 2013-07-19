#include "frank_matrix_wrapper_c.h"
  
void generate_distributed_matrix_col_major_(void** mat){ 
  generate_distributed_matrix_col_major(*mat);
}
void generate_distributed_matrix_row_major_(void** mat){ 
  generate_distributed_matrix_row_major(*mat);
}
void generate_localized_matrix_col_major_(void** mat){ 
  generate_localized_matrix_col_major(*mat);
}
void generate_localized_matrix_row_major_(void** mat){ 
  generate_localized_matrix_row_major(*mat);
}
