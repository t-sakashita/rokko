#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "mpi.h"

#define MAX(a,b) (a>b? a:b)

extern void Cblacs_pinfo(int *my_pnum, int *nprocs);
extern void Cblacs_get(int context, int what, int* val);
extern void Cblacs_gridinit(int* context, char* order, int np_row, int np_col);
extern void Cblacs_gridinfo(int context, int* np_row, int* np_col, int* my_row, int* my_col);
extern void Cblacs_gridexit(int context);
extern void Cblacs_exit(int notdone);

extern int numroc_(int*, int*, int*, int*, int*);
extern void descinit_(int* desc, int* M, int* N, int* mB, int* nB, int* irsrc, int* icsrc, int* ictxt, int* lld, int* info);
extern void pdgeadd_(char* trans, int* M, int* N, double* alpha, double* A, int* iA, int* jA, int* descA,
                     double* beta, double* C, int* iC, int* jC, int* descC);
extern void pdsyev_(char* jobZ, char* uplo, int* N, double* A, int* iA, int* jA, int* descA, double* W,
                    double* Z, int* iZ, int* jZ, int* descZ, double* work, int* lwork, int* info);

int main(int argc, char** argv)
{
  double zero = 0.0;
  double one = 1.0;
  int i_zero = 0;
  int i_one = 1;

  int info=0, context;

  int nprocs, np_row, np_col;
  int my_id, my_row, my_col;

  int N = 4;
  int block_size = 2;

  int local_nrow, local_ncol;

  double *local_frank = 0, *distributed_frank = 0;
  double *eigvals = 0;
  double *local_eigvecs = 0, *distributed_eigvecs = 0;
  int desc_local[11], desc_distributed[11];

  double *work = 0;
  int lwork = -1;

  int i,j;

  /*
   * initialize BLACS
   */
  Cblacs_pinfo( &my_id, &nprocs);
  np_row = sqrt((double)nprocs);
  np_col = nprocs / np_row;

  Cblacs_get(-1, 0, &context);
  Cblacs_gridinit(&context, "R", np_row, np_col);
  Cblacs_gridinfo(context, &np_row, &np_col, &my_row, &my_col);

  /*
   * initialize problem
   */
  if(argc>1)
    N = atoi(argv[1]);

  local_nrow = numroc_(&N, &block_size, &my_row, &i_zero, &np_row);
  local_ncol = numroc_(&N, &block_size, &my_col, &i_zero, &np_col);

  local_frank = (double*)malloc(N*N*sizeof(double));
  distributed_frank = (double*)malloc(local_nrow*local_ncol*sizeof(double));
  distributed_eigvecs = (double*)malloc(local_nrow*local_ncol*sizeof(double));
  local_eigvecs = (double*)malloc(N*N*sizeof(double));
  eigvals = (double*)malloc(N*sizeof(double));

  descinit_(desc_local, &N, &N, &N, &N, &i_zero, &i_zero, &context, &N, &info);
  descinit_(desc_distributed, &N, &N, &block_size, &block_size, &i_zero, &i_zero, &context, &local_nrow, &info);

  /*
   * initialize matrix
   */
  if(my_id == 0){
    for(i=0; i<N; ++i){
      for(j=0; j<N; ++j){
        local_frank[i*N+j] = N - MAX(i,j);
      }
    }
  }

  /*
   * scatter
   */
  pdgeadd_("N", &N, &N, &one, local_frank, &i_one, &i_one, desc_local,
           &zero, distributed_frank, &i_one, &i_one, desc_distributed);

  /*
   * query size for work
   */
  pdsyev_("V", "U", &N, distributed_frank, &i_one, &i_one, desc_distributed, 
          eigvals, distributed_eigvecs, &i_one, &i_one, desc_distributed, eigvals, &lwork, &info);

  lwork = eigvals[0];
  work = (double*)malloc(lwork*sizeof(double));

  /*
   * diagonalize
   */
  pdsyev_("V", "U", &N, distributed_frank, &i_one, &i_one, desc_distributed, 
          eigvals, distributed_eigvecs, &i_one, &i_one, desc_distributed, work, &lwork, &info);

  /*
   * gather
   */
  pdgeadd_("N", &N, &N, &one, distributed_eigvecs, &i_one, &i_one, desc_distributed,
           &zero, local_eigvecs, &i_one, &i_one, desc_local);

  /*
   * print result
   */
  if(my_id == 0){
    printf("Frank Matrix:\n");
    for(i=0; i<N; ++i){
      for(j=0; j<N; ++j){
        printf("%lf ", local_frank[i+N*j]);
      }
      printf("\n");
    }
    printf("\nEigenvalues:");
    for(i=0; i<N; ++i){
      printf(" %lf", eigvals[i]);
    }
    printf("\n\nEigenvectors\n");

    for(i=0; i<N; ++i){
      for(j=0; j<N; ++j){
        printf("%lf ", local_eigvecs[i+N*j]);
      }
      printf("\n");
    }
  }

  /*
   * finalize
   */
  free(work);
  free(eigvals);
  free(local_eigvecs);
  free(distributed_eigvecs);
  free(distributed_frank);
  free(local_frank);

  Cblacs_gridexit(context);
  Cblacs_exit(0);

  return 0;
}
