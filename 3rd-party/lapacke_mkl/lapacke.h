#ifdef __cplusplus
# include <complex>
# define lapack_complex_float std::complex<float>
# define lapack_complex_double std::complex<double>
#else
# include <complex.h>
# undef I
# define MKL_Complex8 complex float
# define MKL_Complex16 complex double
#endif

#include <mkl.h>
#include <mkl_lapacke.h>
