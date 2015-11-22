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

int main(int argc, char *argv[]) {
  int provided;
  MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &provided);
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  std::vector<std::string> solvers;
  if (argc >= 2) {
    solvers.push_back(argv[1]);
  } else {
    solvers = rokko::parallel_sparse_ev::solvers();
  }

  int L = (argc >= 3) ? boost::lexical_cast<int>(argv[2]) : 10;
  int dim = 1 << L;
  std::vector<std::pair<int, int> > lattice;
  for (int i = 0; i < L; ++i) lattice.push_back(std::make_pair(i, (i+1) % L));

  rokko::parameters params;
  params.set("Block Size", 5);
  params.set("Maximum Iterations", 500);
  params.set("Convergence Tolerance", 1.0e-8);
  params.set("num_eigenvalues", 10);
  BOOST_FOREACH(std::string const& name, solvers) {
    rokko::parallel_sparse_ev solver(name);
    rokko::distributed_crs_matrix mat(dim, dim, solver);
    std::vector<double> values;
    std::vector<int> cols;
    for (int row = mat.start_row(); row <= mat.end_row(); ++row) {
      cols.clear();
      values.clear();
      double diag = 0;
      for (int l = 0;  l < lattice.size(); ++l) {
        int i = lattice[l].first;
        int j = lattice[l].second;
        int m1 = 1 << i;
        int m2 = 1 << j;
        int m3 = m1 + m2;
        if (((row & m3) == m1) || ((row & m3) == m2)) {
          cols.push_back(row^m3);
          values.push_back(0.5);
          diag += -0.25;
        } else {
          diag += 0.25;
        }
      }
      cols.push_back(row);
      values.push_back(diag);
      mat.insert(row, cols, values);
    }
    mat.complete();
    if (rank == 0)
      std::cout << "Eigenvalue decomposition of antiferromagnetic Heisenberg chain" << std::endl
                << "solver = " << name << std::endl
                << "L = " << L << std::endl
                << "dimension = " << dim << std::endl;

    rokko::parameters info = solver.diagonalize(mat, params);

    int num_conv = info.get<int>("num_conv");
    if (num_conv == 0) MPI_Abort(MPI_COMM_WORLD, -1);
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
  }

  MPI_Finalize();
}
