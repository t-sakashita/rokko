#include "frank_matrix_wrapper_c.h"
  
void frank_generate_distributed_matrix_col_major_(void** mat){ 
  frank_generate_distributed_matrix_col_major(*mat);
}
void frank_generate_distributed_matrix_row_major_(void** mat){ 
  frank_generate_distributed_matrix_row_major(*mat);
}
void frank_generate_localized_matrix_col_major_(void** lmat){ 
  frank_generate_localized_matrix_col_major(*lmat);
}
void frank_generate_localized_matrix_row_major_(void** lmat){ 
  frank_generate_localized_matrix_row_major(*lmat);
}
