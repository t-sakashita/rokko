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

#include <rokko/rokko.hpp>
#include <rokko/utility/heisenberg_hamiltonian_mfree.hpp>
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

  int L = (argc >= 3) ? boost::lexical_cast<int>(argv[2]) : 10;
  std::vector<std::pair<int, int> > lattice;
  for (int i = 0; i < L; ++i) lattice.push_back(std::make_pair(i, (i+1) % L));

  rokko::parameters params;
  if (!routine.empty()) params.set("routine", routine);
  params.set("Block Size", 5);
  params.set("Maximum Iterations", 500);
  params.set("Convergence Tolerance", 1.0e-8);
  params.set("num_eigenvalues", 10);
  
  rokko::parallel_sparse_ev solver(library);
  rokko::heisenberg_mfree mat(L, lattice);
  int dim = mat.get_dim();
  if (rank == 0)
    std::cout << "Eigenvalue decomposition of antiferromagnetic Heisenberg chain" << std::endl
	      << "solver = " << library << std::endl
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
  
  solver.finalize();
  MPI_Finalize();
}
