/*****************************************************************************
*
* Rokko: Integrated Interface for libraries of eigenvalue decomposition
*
* Copyright (C) 2012-2020 Rokko Developers https://github.com/t-sakashita/rokko
*
* Distributed under the Boost Software License, Version 1.0. (See accompanying
* file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
*
*****************************************************************************/

#include <rokko/rokko.hpp>
#include <rokko/utility/solver_name.hpp>

#include <stdexcept>

int main(int argc, char *argv[]) {
  int provided;
  MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &provided);
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  
  std::string library_routine(rokko::parallel_sparse_ev::default_solver());
  std::string library, routine;
  if (argc >= 2) library_routine = argv[1];
  rokko::split_solver_name(library_routine, library, routine);

  int L = (argc >= 3) ? std::stoi(argv[2]) : 10;
  int dim = 1 << L;
  std::vector<std::pair<int, int>> lattice;
  for (int i = 0; i < L; ++i) lattice.emplace_back(std::make_pair(i, (i+1) % L));

  rokko::parameters params;
  params.set("routine", routine);
  params.set("Block Size", 5);
  params.set("Maximum Iterations", 500);
  params.set("Convergence Tolerance", 1.0e-8);
  params.set("num_eigenvalues", 10);

  rokko::parallel_sparse_ev solver(library);
  auto map = solver.default_mapping(dim, rokko::mpi_comm{MPI_COMM_WORLD});
  const int num_entries_per_row = lattice.size() + 1;
  rokko::distributed_crs_matrix mat(map, num_entries_per_row);
  std::vector<int> cols;
  std::vector<double> values;
  cols.reserve(num_entries_per_row);
  values.reserve(num_entries_per_row);

  for (int row = map.start_row(); row < map.end_row(); ++row) {
    cols.clear();
    values.clear();
    double diag = 0;
    for (int l = 0; l < lattice.size(); ++l) {
      int i = lattice[l].first;
      int j = lattice[l].second;
      int m1 = 1 << i;
      int m2 = 1 << j;
      int m3 = m1 + m2;
      if (((row & m3) == m1) || ((row & m3) == m2)) {
        cols.emplace_back(row^m3);
        values.emplace_back(0.5);
        diag += -0.25;
      } else {
        diag += 0.25;
      }
    }
    if (diag != 0.) {
      cols.emplace_back(row);
      values.emplace_back(diag);
    }
    std::cout << "row=" << row << " size=" << cols.size() << std::endl;
    mat.insert(row, cols, values);
  }
  mat.complete();
  mat.print();
  if (rank == 0)
    std::cout << "Eigenvalue decomposition of antiferromagnetic Heisenberg chain" << std::endl
	      << "solver = " << library << std::endl
	      << "L = " << L << std::endl
	      << "dimension = " << dim << std::endl;
  
  rokko::parameters info = solver.diagonalize(mat, params);
  
  int num_conv = info.get<int>("num_conv");
  if (num_conv == 0)
    throw std::runtime_error("num_conv=0: solver did not converge");
  std::vector<double> eigvec;
  solver.eigenvector(0, eigvec);
  if (rank == 0) {
    std::cout << "number of converged eigenpairs = " << num_conv << std::endl;
    std::cout << "smallest eigenvalues: ";
    for (int i = 0; i < num_conv; ++i) std::cout << ' ' << solver.eigenvalue(i);
    std::cout << std::endl;
    std::cout << "smallest eigenvector: ";
    for (int j = 0; j < eigvec.size(); ++j) std::cout << eigvec[j] << ' ';
    std::cout << std::endl;
  }
  
  solver.finalize();
  MPI_Finalize();
}
