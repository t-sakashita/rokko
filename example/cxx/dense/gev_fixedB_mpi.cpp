/*****************************************************************************
*
* Rokko: Integrated Interface for libraries of eigenvalue decomposition
*
* Copyright (C) 2012-2015 Rokko Developers https://github.com/t-sakashita/rokko
*
* Distributed under the Boost Software License, Version 1.0. (See accompanying
* file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
*
*****************************************************************************/

#include <rokko/rokko.hpp>
#include <rokko/mapping_bc.hpp>
#include <rokko/collective.hpp>
#include <rokko/utility/solver_name.hpp>
#include <rokko/utility/frank_matrix.hpp>
#include <rokko/utility/check_orthogonality.hpp>
#include <boost/lexical_cast.hpp>
#include <iostream>


typedef rokko::matrix_col_major matrix_major;

template<typename T, typename MATRIX_MAJOR>
void function_matrix(rokko::localized_vector<double> const& eigval_tmp, rokko::distributed_matrix<T, MATRIX_MAJOR> const& eigvec, rokko::distributed_matrix<T, MATRIX_MAJOR>& result, rokko::distributed_matrix<T, MATRIX_MAJOR>& tmp) {
  for (int local_j=0; local_j<eigvec.get_n_local(); ++local_j) {
    int global_j = eigvec.translate_l2g_col(local_j);
    double coeff = eigval_tmp(global_j);
    for (int local_i=0; local_i<eigvec.get_m_local(); ++local_i) {
      double value = eigvec.get_local(local_i, local_j);
      tmp.set_local(local_i, local_j, coeff * value); 
    }
  }
  product(1, tmp, false, eigvec, true, 0, result);
}

template<typename T, typename MATRIX_MAJOR>
void diagonalize_fixedB(rokko::parallel_dense_solver& solver, rokko::distributed_matrix<T, MATRIX_MAJOR>& A, rokko::distributed_matrix<T, MATRIX_MAJOR>& B,
			rokko::localized_vector<double>& eigval, rokko::distributed_matrix<T, MATRIX_MAJOR>& eigvec, T tol = 0) {
  rokko::distributed_matrix<double, matrix_major> tmp(A.get_mapping()), Binvroot(A.get_mapping()), mat(A.get_mapping());
  rokko::parameters params;
  int myrank = A.get_myrank();
  params.set("routine", "");
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

template<typename T, typename MATRIX_MAJOR>
void set_A_B(rokko::localized_matrix<T, MATRIX_MAJOR>& locA, rokko::localized_matrix<T, MATRIX_MAJOR>& locB) {
  if ((locA.rows() != 4) || (locA.cols() != 4) || (locB.rows() != 4) || (locB.cols() != 4)) {
    std::cerr << "error: size must be 4!" << std::endl;
    throw;
  }
  locA << 0.24, 0.39, 0.42, -0.16,
          0.39, -0.11, 0.79, 0.63,
          0.42, 0.79, -0.25, 0.48,
         -0.16, 0.63, 0.48, -0.03;

  locB << 4.16, -3.12, 0.56, -0.10,
         -3.12, 5.03, -0.83, 1.09,
          0.56, -0.83, 0.76, 0.34,
         -0.10, 1.09, 0.34, 1.18;
}

  
int main(int argc, char *argv[]) {
  int provided;
  MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &provided);
  MPI_Comm comm = MPI_COMM_WORLD;
  std::string library_routine(rokko::parallel_dense_solver::default_solver());
  std::string library, routine;
  int dim = 4;
  if (argc >= 2) library_routine = argv[1];
  rokko::split_solver_name(library_routine, library, routine);

  rokko::grid g(comm);
  int myrank = g.get_myrank();

  std::cout.precision(5);

  rokko::parallel_dense_solver solver(library);
  solver.initialize(argc, argv);
  if (myrank == 0)
    std::cout << "Eigenvalue decomposition of Frank matrix" << std::endl
	      << "library:routine = " << library_routine << std::endl
              << "num_procs = " << g.get_nprocs() << std::endl
              #ifdef _OPENMP
              << "num_threads per process = " << omp_get_max_threads() << std::endl
              #endif
	      << "routine = " << routine << std::endl
              << "dimension = " << dim << std::endl;

  rokko::localized_matrix<double, matrix_major> locA(dim, dim), locB(dim, dim);
  set_A_B(locA, locB);
  if (myrank == 0) std::cout << "locA:" << std::endl << locA << std::endl;  

  rokko::mapping_bc<matrix_major> map(dim, g, solver);
  rokko::distributed_matrix<double, matrix_major> A(map), B(map), eigvec(map);
  rokko::localized_vector<double> eigval(dim);
  rokko::scatter(locA, A, 0);
  rokko::scatter(locB, B, 0);
  MPI_Barrier(comm);
  if (myrank == 0) std::cout << "A:" << std::endl;
  std::cout << A << std::endl;

  diagonalize_fixedB(solver, A, B, eigval, eigvec);

  rokko::localized_matrix<double, matrix_major> eigvec_loc(dim, dim);
  rokko::gather(eigvec, eigvec_loc, 0);

  set_A_B(locA, locB);

  if (myrank == 0) {
    bool sorted = true;
    for (unsigned int i = 1; i < dim; ++i) sorted &= (eigval(i-1) <= eigval(i));
    if (!sorted) std::cout << "Warning: eigenvalues are not sorted in ascending order!\n";
    std::cout << "largest eigenvalues:";
    for (int i = 0; i < std::min(dim, 10); ++i) std::cout << ' ' << eigval(dim - 1 - i);
    std::cout << std::endl;

    std::cout << "eigenvalues:\n" << eigval.transpose() << std::endl
	      << "eigvectors:\n" << eigvec_loc << std::endl;
    std::cout << "orthogonality of eigenvectors:" << std::endl
	      << eigvec_loc.transpose() * locB * eigvec_loc << std::endl;
    std::cout << "residual of the smallest eigenvalue/vector (A x - lambda B x):" << std::endl
	      << (locA * eigvec_loc.col(0) - eigval(0) * locB * eigvec_loc.col(0)).transpose() << std::endl;
  }

  solver.finalize();
  MPI_Finalize();
}
