//#include <rokko/solver.hpp>
//#include <rokko/grid.hpp>
//#include <rokko/distributed_matrix.hpp>
//#include <rokko/localized_vector.hpp>

#ifdef __cplusplus
namespace rokko {
  
  extern "C"{
#endif


    void* initialize_distributed_matrix_row_major(int , int, void*, void* );
    void* initialize_distributed_matrix_col_major(int, int, void*, void*);
    void delete_distributed_matrix_col_major(void* );
    void delete_distributed_matrix_row_major(void* );

    void* initialize_localized_vector( int);
    void delete_localized_vector(void*);
    double localized_vector_get_element(void*, int);



    void* initialize_grid_col_major(MPI_Comm);
    void* initialize_grid_row_major(MPI_Comm); 
    int grid_get_myrank(void*);
    int grid_get_nprocs(void*);
    void delete_grid(void*);

    void* initialize_solver(char* ,int, char*[]);
    void delete_solver(void* );

    void solver_diagonalize_matrix_col_major(void*, void*, void*, void*, void*);
    void solver_diagonalize_matrix_row_major(void*, void*, void*, void*, void*);



#ifdef __cplusplus
  }
}
#endif
