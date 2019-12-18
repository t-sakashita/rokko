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
#include <rokko/collective.hpp>
#include <rokko/utility/solver_name.hpp>
#include <rokko/utility/helmert_matrix.hpp>
#include <boost/lexical_cast.hpp>
#include <iostream>

using matrix_major = rokko::matrix_col_major;


double func(int i, int j) {
  Eigen::MatrixXd mat_loc(9, 9);
  //  mat_loc << 2.828968253968254132,  0.8289682539682539097,  0.3289682539682539097,  -0.004365079365079349571,  -0.2543650793650793496,  -0.4543650793650793052,  -0.6210317460317460458,  -0.7638888888888888395,  -0.8888888888888888395;
  
  mat_loc << 2.828968253968254132,  0.8289682539682539097,  0.3289682539682539097,  -0.004365079365079349571,  -0.2543650793650793496,  -0.4543650793650793052,  -0.6210317460317460458,  -0.7638888888888888395,  -0.8888888888888888395,
    0.8289682539682539097,  2.828968253968254132,  0.3289682539682539097,  -0.004365079365079349571,  -0.2543650793650793496,  -0.4543650793650793052,  -0.6210317460317460458,  -0.7638888888888888395,  -0.8888888888888888395,
    0.3289682539682539097,  0.3289682539682539097,  3.328968253968254132,  -0.004365079365079349571,  -0.2543650793650793496,  -0.4543650793650793052,  -0.6210317460317460458,  -0.7638888888888888395,  -0.8888888888888888395,
    -0.004365079365079349571,  -0.004365079365079349571,  -0.004365079365079349571,  3.99563492063492065,  -0.2543650793650793496,  -0.4543650793650793052,  -0.6210317460317460458,  -0.7638888888888888395,  -0.8888888888888888395,
    -0.2543650793650793496,  -0.2543650793650793496,  -0.2543650793650793496,  -0.2543650793650793496,  4.745634920634920206,  -0.4543650793650793052,  -0.6210317460317460458,  -0.7638888888888888395,  -0.8888888888888888395,
    -0.4543650793650793052,  -0.4543650793650793052,  -0.4543650793650793052,  -0.4543650793650793052,  -0.4543650793650793052,  5.545634920634920917,  -0.6210317460317460458,  -0.7638888888888888395,  -0.8888888888888888395,
    -0.6210317460317460458,  -0.6210317460317460458,  -0.6210317460317460458,  -0.6210317460317460458,  -0.6210317460317460458,  -0.6210317460317460458,  6.378968253968253954,  -0.7638888888888888395,  -0.8888888888888888395,
    -0.7638888888888888395,  -0.7638888888888888395,  -0.7638888888888888395,  -0.7638888888888888395,  -0.7638888888888888395,  -0.7638888888888888395,  -0.7638888888888888395,  7.236111111111110716,  -0.8888888888888888395,
    -0.8888888888888888395,  -0.8888888888888888395,  -0.8888888888888888395,  -0.8888888888888888395,  -0.8888888888888888395,  -0.8888888888888888395,  -0.8888888888888888395,  -0.8888888888888888395,  8.111111111111110716;

  return mat_loc(i,j);
}
  
int main(int argc, char *argv[]) {

  int provided;
  MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &provided);
  MPI_Comm comm = MPI_COMM_WORLD;
  std::string library_routine(rokko::parallel_dense_ev::default_solver());
  std::string library, routine;
  int dim = 9;
  if (argc >= 2) library_routine = argv[1];
  if (argc >= 3) dim = boost::lexical_cast<int>(argv[2]);
  rokko::split_solver_name(library_routine, library, routine);

  rokko::grid g(comm);
  int myrank = g.get_myrank();

  std::cout.precision(std::numeric_limits<unsigned long long>::digits10);

  rokko::parallel_dense_ev solver(library);
  solver.initialize(argc, argv);
  if (myrank == 0)
    std::cout << "Eigenvalue decomposition of Helmert matrix" << std::endl
              << "num_procs = " << g.get_nprocs() << std::endl
              #ifdef _OPENMP
              << "num_threads per process = " << omp_get_max_threads() << std::endl
              #endif
      	      << "library:routine = " << library_routine << std::endl
              << "dimension = " << dim << std::endl;

  rokko::mapping_bc<matrix_major> map = solver.default_mapping(dim, g);
  rokko::distributed_matrix<double, matrix_major> mat(map);
  mat.generate(&func);
  mat.print();

  Eigen::VectorXd eigval(dim);
  rokko::distributed_matrix<double, matrix_major> eigvec(map);
  rokko::parameters params;
  params.set("routine", routine);
  solver.diagonalize(mat, eigval, eigvec, params);

  if (myrank == 0) {
    std::cout << "largest eigenvalues:";
    //for (int i = 0; i < std::min(dim, 10); ++i) std::cout << ' ' << eigval(dim - 1 - i);
    for (int i = 0; i < dim; ++i) std::cout << ' ' << eigval(dim - 1 - i);
    std::cout << std::endl;
  }

  solver.finalize();
  MPI_Finalize();
}
