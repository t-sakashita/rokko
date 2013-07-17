#ifdef __cplusplus
namespace rokko {
  extern "C"{
#endif
    
    void generate_distributed_matrix_col_major(void*);
    void generate_distributed_matrix_row_major(void*);
    void generate_localized_matrix_col_major(void*);
    void generate_localized_matrix_row_major(void*);


#ifdef __cplusplus
  }
}
#endif
