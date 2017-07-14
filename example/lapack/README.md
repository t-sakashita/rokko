# Examples of LAPACK, LAPACKE, and rokko::lapack

## LAPACK (Fortran API) examples

  * dlange_f.f90 : matrix norm (double precision)
  * slange_f.f90 : matrix norm (single precision)

## LAPACKE (C API) examples

  * column-major matrix examples

    * dgesvd.c : singular value decomposition (double precision)
    * dgetri.c : matrix inversion (double precision)
    * dgetrs.c : LU decomposition (double precision)
    * dlange.c : matrix norm (double precision)
    * dsyev.c : eigenvalue decomposition of symmetric matrix (double precision)
    * slange.c : matrix norm (single precision)
    * zheev.c : eigenvalue decomposition of Hermitian matrix (double precision)
  
  * row-major matrix examples
  
    * dsyev_r.c : eigenvalue decomposition of symmetric matrix (double precision)

## rokko::lapack (C++ API) examples

  * column-major matrix examples

    * gesvd.cpp : singular value decomposition of real matrix (double precision)
    * gesvd_z.cpp : singular value decomposition of complex matrix (double precision)
    * getrs.cpp : LU decomposition of real matrix (double precision)
    * getrs_z.cpp : LU decomposition of complex matrix (double precision)
    * heev.cpp : eigenvalue decomposition of Hermitian matrix (double precision)
    * heev_c.cpp : eigenvalue decomposition of Hermitian matrix (single precision)
    * lange.cpp : matrix norms
    * orgqr.cpp : QR factorization of real matrix (double precision)
    * syev.cpp : eigenvalue decomposition of symmetric matrix (double precision)
    * syev_s.cpp : eigenvalue decomposition of symmetric matrix (single precision)
    * ungqr.cpp : QR factorization of complex matrix (double precision)
