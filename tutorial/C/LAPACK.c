/*****************************************************************************
 *
 * Rokko: Integrated Interface for libraries of eigenvalue decomposition
 *
 * Copyright (C) 2012-2014 by Tatsuya Sakashita <t-sakashita@issp.u-tokyo.ac.jp>,
 *                            Synge Todo <wistaria@comp-phys.org>,
 *                            Tsuyoshi Okubo <t-okubo@issp.u-tokyo.ac.jp>
 *                            Yuichi Motoyama <y-motoyama@issp.u-tokyo.ac.jp>
 *    
 * Distributed under the Boost Software License, Version 1.0. (See accompanying
 * file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
 *
 *****************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define N 4
#define LDA N

extern void dsyev_(char* jobz,
                  char* uplo, int* n, double* a, int* lda, 
                  double* w, double* work, int* lwork, int* info );

void print_matrix(double* a, int nrow, int ncol, int lda)
{
  int i,j;
  for(i=0; i<nrow; ++i){
    for(j=0; j<ncol; ++j){
      int index = i+j*lda;
      printf(" %f", a[index]);
    }
    printf("\n");
  }
}

int main()
{
  int i,j;
  int n = N, lda = LDA, info;
  int lwork = -1;
  double wkopt;
  double* work;
  double eigvals[N];
  double frank[LDA*N];
  for(i=0; i<N; ++i){
    for(j=0; j<N; ++j){
      int index = i*N+j;
      frank[index] = N - (i > j ? i : j);
    }
  }

  printf("Frank matrix: \n");
  print_matrix(frank, N, N, N);
  printf("\n");

  dsyev_("V", "U", &n, frank, &lda, eigvals, &wkopt, &lwork, &info);
  lwork = (int)wkopt;
  work = (double*)malloc( lwork*sizeof(double));
  dsyev_("V", "U", &n, frank, &lda, eigvals, work, &lwork, &info);

  printf("Eigenvalues: \n");
  print_matrix(eigvals, 1, N, 1);
  printf("\n");

  printf("Eigenstates: \n");
  print_matrix(frank, N, N, N);
  printf("\n");

  free(work);

  return 0;
}

