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
#include <rokko/utility/laplacian_matrix.hpp>
#include <rokko/utility/solver_name.hpp>

int main(int argc, char *argv[]) {
  int provided;
  MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &provided);
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  std::string library_routine(rokko::parallel_sparse_ev::default_solver());
  std::string library, routine;
  if (argc >= 2) library_routine = argv[1];
  rokko::split_solver_name(library_routine, library, routine);
  
  int dim = (argc >= 3) ? boost::lexical_cast<int>(argv[2]) : 100;

  rokko::parameters params;
  if (!routine.empty()) params.set("routine", routine);
  params.set("Block Size", 5);
  params.set("Maximum Iterations", 500);
  params.set("Convergence Tolerance", 1.0e-8);
  params.set("num_eigenvalues", 10);
  rokko::parallel_sparse_ev solver(library);
  rokko::distributed_crs_matrix mat(dim, dim, solver);
  std::vector<double> values(3);
  std::vector<int> cols(3);
  int start_row = mat.start_row();
  if (start_row == 0) {
    values.push_back(1.);  values.push_back(-1.);
    cols.push_back(0);  cols.push_back(1);
    ++start_row;
    mat.insert(0, cols, values);
  }
  int end_loop_row = mat.end_row();
  bool has_end_row = (mat.end_row() == dim);
  if (has_end_row) {
    --end_loop_row;
  }
  values.clear();
  values.push_back(-1.);  values.push_back(2.);  values.push_back(-1.);
  for (int row = start_row; row < end_loop_row; ++row) {
    cols.clear();
    cols.push_back(row-1);  cols.push_back(row);  cols.push_back(row+1);
    mat.insert(row, cols, values);
  }
  cols.clear();
  values.clear();
  if (has_end_row) {
    values.push_back(-1.);  values.push_back(2.);
    cols.push_back(dim-2);  cols.push_back(dim-1);
    mat.insert(dim-1, cols, values);
  }
  mat.complete();
  if (rank == 0)
    std::cout << "Eigenvalue decomposition of Laplacian" << std::endl
	      << "solver = " << library << std::endl
	      << "dimension = " << dim << std::endl;
  
  rokko::parameters info = solver.diagonalize(mat, params);
  
  int num_conv = info.get<int>("num_conv");
  if (num_conv == 0) MPI_Abort(MPI_COMM_WORLD, -1);
  std::vector<double> eigvec;
  solver.eigenvector(0, eigvec);
  if (rank == 0) {
    std::cout << "number of converged eigenpairs = " << num_conv << std::endl;
    std::cout << "largest eigenvalues:     ";
    for (int i = 0; i < num_conv; ++i) std::cout << solver.eigenvalue(i) << ' ';
    std::cout << std::endl;
    std::cout << "theoretical eigenvalues: ";
    for (int k = dim - 1; k >= dim - num_conv; --k)
        std::cout << rokko::laplacian_matrix::eigenvalue(dim, k) << ' ';
    std::cout << std::endl;
    std::cout << "largest eigenvector: ";
    for (int j = 0; j < eigvec.size(); ++j) std::cout << eigvec[j] << ' ';
    std::cout << std::endl;
  }
  
  solver.finalize();
  MPI_Finalize();
}
