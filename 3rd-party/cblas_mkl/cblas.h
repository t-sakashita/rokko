#ifndef __cplusplus
# include <complex.h>
# undef I
# define MKL_Complex8 complex float
# define MKL_Complex16 complex double
#endif

#include <mkl.h>
