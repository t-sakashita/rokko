#include "wrapper_c.h";
#include <mpi.h>

void initialize_distributed_matrix_row_major_(void** mat, int* dim1, int* dim2, void** g, void** solver_ ){ 
  
  *mat = initialize_distributed_matrix_row_major(*dim1, *dim2, *g, *solver_); 
}

void initialize_distributed_matrix_col_major_(void** mat, int* dim1, int* dim2, void** g, void** solver_ ){ 
  *mat =  initialize_distributed_matrix_col_major(*dim1, *dim2, *g, *solver_);
}

void delete_distributed_matrix_col_major_(void**  mat){      
  delete_distributed_matrix_col_major(*mat);
}

void delete_distributed_matrix_row_major_(void**  mat){
  delete_distributed_matrix_row_major(*mat);
}

void initialize_localized_vector_(void** vec, int* dim){
  *vec = initialize_localized_vector(*dim);
}

void delete_localized_vector_(void** w ){
  delete_localized_vector(*w);
}

double localized_vector_get_element_(void** w, int* i){
  return localized_vector_get_element(*w, *i);
}

void initialize_grid_col_major_(void** grid_, MPI_Fint* comm){ 
  MPI_Comm comm_c = MPI_Comm_f2c(*comm);
  *grid_ = initialize_grid_col_major(comm_c);
}

void initialize_grid_row_major_(void** grid_, MPI_Fint* comm){ 
  MPI_Comm comm_c = MPI_Comm_f2c(*comm);
  *grid_ = initialize_grid_row_major(comm_c);
}

int grid_get_myrank_(void** g){
  return grid_get_myrank(*g);
}

int grid_get_nprocs_(void** g){
    return grid_get_nprocs(*g);
}

void delete_grid_(void** g){
  delete_grid(*g);
}

void initialize_solver_(void** solver_, char* solver_name, int* argc, char* argv[], long length_name, long length_arg ){ 
  *solver_ = initialize_solver(solver_name, *argc, argv);
}

void delete_solver_(void** solver_){
  delete_solver(*solver_);
}

void solver_diagonalize_matrix_col_major_(void** solver_ ,void** mat, void** w, void** Z, void** timer){
  solver_diagonalize_matrix_col_major(*solver_, *mat, *w, *Z, *timer);
}

void solver_diagonalize_matrix_row_major_(void** solver_ ,void** mat, void** w, void** Z, void** timer){
  solver_diagonalize_matrix_row_major(*solver_, *mat, *w, *Z, *timer);
}

