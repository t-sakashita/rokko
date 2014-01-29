#include <stdlib.h>
#include <stdio.h>
#include <lapacke.h>

#define N 5

int main() {
  int i, info;
  double w[N];
  double a[N * N] = {
    5, 4, 3, 2, 1,
    4, 4, 3, 2, 1,
    3, 3, 3, 2, 1,
    2, 2, 2, 2, 1,
    1, 1, 1, 1, 1 
  };
  info = LAPACKE_dsyev(LAPACK_COL_MAJOR, 'V', 'U', N, a, N, w);
  printf("Eigenvalues of Frank matrix: ");
  for (i = 0; i < N; ++i) {
    printf(" %6.2f", w[i]);
  }
  printf("\n");
  return 0;
}
