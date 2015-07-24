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
void diagonalize_fixedB(rokko::parallel_dense_solver& solver, rokko::distributed_matrix<T, MATRIX_MAJOR>& A, rokko::distributed_matrix<T, MATRIX_MAJOR>& B, rokko::localized_vector<double>& eigval, rokko::distributed_matrix<T, MATRIX_MAJOR>& eigvec) {
  rokko::localized_vector<double> eigval_inv(eigval.size());
  rokko::distributed_matrix<double, matrix_major> tmp(A.get_mapping()), Binv(B.get_mapping());

  int myrank = A.get_myrank();
  // diagonalization of B
  std::string routine = "";
  try {
    solver.diagonalize(A, eigval, eigvec);
  }
  catch (const char *e) {
    if (myrank == 0) std::cout << "Exception : " << e << std::endl;
    MPI_Abort(MPI_COMM_WORLD, 22);
  }

  // computation of B^{-1}
  for(int i=0; i<eigval.size(); ++i) {
    eigval_inv(i) = 1;// / eigval(i);
  }
  for (int local_j=0; local_j<eigvec.get_n_local(); ++local_j) {
    int global_j = eigvec.translate_l2g_col(local_j);
    double coeff = eigval_inv(global_j);  //1 / eigval(global_j);
    for (int local_i=0; local_i<eigvec.get_m_local(); ++local_i) {
      double value = eigvec.get_local(local_i, local_j);
      tmp.set_local(local_i, local_j, coeff * value); 
    }
  }
  product(1, tmp, false, eigvec, true, 0, Binv);
  std::cout << "inverseB:" << std::endl << Binv << std::endl;
  //Binv.print();

  // computation of inverse B^{-1} A
  product(1, Binv, false, A, false, 0, tmp);
  // diagonalization of B^{-1} A
  try {
    solver.diagonalize(tmp, eigval, eigvec);
  }
  catch (const char *e) {
    if (myrank == 0) std::cout << "Exception : " << e << std::endl;
    MPI_Abort(MPI_COMM_WORLD, 22);
  }
}

int main(int argc, char *argv[]) {
  int provided;
  MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &provided);
  MPI_Comm comm = MPI_COMM_WORLD;
  std::string library_routine(rokko::parallel_dense_solver::default_solver());
  std::string library, routine;
  int dim = 10;
  if (argc >= 2) library_routine = argv[1];
  if (argc >= 3) dim = boost::lexical_cast<int>(argv[2]);
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

  rokko::mapping_bc<matrix_major> map(dim, g, solver);

  rokko::distributed_matrix<double, matrix_major> A(map);
  rokko::distributed_matrix<double, matrix_major> B(map);
  rokko::frank_matrix::generate(A);
  rokko::localized_vector<double> eigval(dim);
  rokko::distributed_matrix<double, matrix_major> eigvec(map);
  diagonalize_fixedB(solver, A, B, eigval, eigvec);
  
  /*
  rokko::localized_matrix<double, matrix_major> eigvec_loc(dim, dim);
  rokko::gather(eigvec, eigvec_loc, 0);
  if (myrank == 0) {
    bool sorted = true;
    for (unsigned int i = 1; i < dim; ++i) sorted &= (eigval(i-1) <= eigval(i));
    if (!sorted) std::cout << "Warning: eigenvalues are not sorted in ascending order!\n";

    std::cout << "largest eigenvalues:";
    for (int i = 0; i < std::min(dim, 10); ++i) std::cout << ' ' << eigval(dim - 1 - i);
    std::cout << std::endl;
    std::cout << "residual of the largest eigenvalue/vector: |x A x - lambda| = "
              << std::abs(eigvec_loc.col(dim - 1).transpose() * mat_loc * eigvec_loc.col(dim - 1)
                          - eigval(dim - 1))
              << std::endl;
  }
  */

  solver.finalize();
  MPI_Finalize();
}
