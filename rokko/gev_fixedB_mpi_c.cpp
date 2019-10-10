/*****************************************************************************
*
* Rokko: Integrated Interface for libraries of eigenvalue decomposition
*
* Copyright (C) 2012-2019 Rokko Developers https://github.com/t-sakashita/rokko
*
* Distributed under the Boost Software License, Version 1.0. (See accompanying
* file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
*
*****************************************************************************/

#include <mpi.h>
#include <rokko/rokko.h>
#include <rokko/utility/frank_matrix.h>
#include <stdio.h>
#include <stdlib.h>

using matrix_major = rokko::matrix_col_major;

void diagonalize_fixedB(struct rokko_parallel_dense_ev* solver_in, struct rokko_distributed_matrix* A_in, struct rokko_distributed_matrix* B, struct rokko_eigen_vector* eigval_in, struct rokko_distributed_matrix* eigvec_in, double tol) {
  rokko::mapping_bc<matrix_major> map = A->ptr->get_mapping();
  rokko::distributed_matrix<double, matrix_major> tmp(map), Binvroot(map), mat(map);
  rokko::parameters params = *(params->ptr);
  int myrank = A.get_myrank();
  params.set(params, "routine", "");
  solver.diagonalize(B, eigval, eigvec, params);
  // computation of B^{-1/2}
  for(int i=0; i<eigval.size(); ++i)
    eigval(i) = (eigval(i) > tol) ? sqrt(1/eigval(i)) : 0;
  function_matrix(eigval, eigvec, Binvroot, tmp);
  
  // computation of B^{-1/2} A B^{-1/2}
  product(1, Binvroot, false, A, false, 0, tmp);
  product(1, tmp, false, Binvroot, false, 0, mat);
  // diagonalization of B^{-1/2} A B^{-1/2}
  solver.diagonalize(mat, eigval, tmp, params);

  // computation of {eigvec of Ax=lambda Bx} = B^{-1/2} {eigvec of B^{-1/2} A B^{-1/2}}
  product(1, Binvroot, false, tmp, false, 0, eigvec);
}

