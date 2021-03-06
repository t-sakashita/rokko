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
#include <rokko/utility/heisenberg_hamiltonian_mfree.hpp>
#include <rokko/utility/math.hpp>
#include <rokko/utility/output_matrix_market.hpp>

int main(int argc, char *argv[]) {
  int provided;
  MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &provided);
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  std::string library(rokko::parallel_sparse_ev::default_solver());
  if (argc >= 2) library = argv[1];
  
  int L = (argc >= 3) ? std::stoi(argv[2]) : 10;
  int dim = 1 << L;
  std::vector<std::pair<int, int>> lattice;
  for (int i = 0; i < L; ++i) lattice.emplace_back(std::make_pair(i, (i+1) % L));

  rokko::parallel_sparse_ev solver(library);
  rokko::heisenberg_mfree op(L, lattice);
  rokko::output_matrix_market(op);

  solver.finalize();
  MPI_Finalize();
}
