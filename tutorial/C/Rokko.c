#include <rokko/rokko.h>
#include <stdio.h>
#include <stdlib.h>

int main(int argc, char* argv[]) {
  struct rokko_serial_dense_solver solver;
  struct rokko_localized_matrix frank, eigvecs;
  struct rokko_localized_vector eigvals;
  const char* solver_name = argc>1 ? argv[1] : "lapack";
  int dim = argc>2 ? atoi(argv[2]) : 4;

  int i,j;

  rokko_serial_dense_solver_construct(&solver, solver_name, argc, argv);

  rokko_localized_matrix_construct(&frank, dim, dim, rokko_matrix_col_major);
  rokko_localized_matrix_construct(&eigvecs, dim, dim, rokko_matrix_col_major);
  rokko_localized_vector_construct(&eigvals, dim);

  for(i=0; i<dim; ++i){
    for(j=0; j<dim; ++j){
      rokko_localized_matrix_set_local(&frank, i, j, dim - (i>j ? i:j));
    }
  }

  printf("Frank matrix: \n");
  rokko_localized_matrix_print(frank);
  printf("\n");

  rokko_serial_dense_solver_diagonalize_localized_matrix(&solver, &frank, &eigvals, &eigvecs);

  printf("Eigenvalues: \n");
  for (i = 0; i < dim; ++i){
    printf("%30.20f ", rokko_localized_vector_get(eigvals, i));
  }
  printf("\n");

  printf("Eigenvectors: \n");
  rokko_localized_matrix_print(eigvecs);
  printf("\n");

  rokko_localized_vector_destruct(&eigvals);
  rokko_localized_matrix_destruct(&eigvecs);
  rokko_localized_matrix_destruct(&frank);
  rokko_serial_dense_solver_destruct(&solver);

  return 0;
}
