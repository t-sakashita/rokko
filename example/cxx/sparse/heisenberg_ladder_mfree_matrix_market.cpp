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
#include <rokko/utility/lattice.hpp>
#include <rokko/utility/math.hpp>
#include <rokko/distributed_mfree_to_crs.hpp>

int main(int argc, char *argv[]) {
  int provided;
  MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &provided);
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  const std::string library = (argc >= 2) ? argv[1] : rokko::parallel_sparse_ev::default_solver();

  const int len_ladder = (argc >= 3) ? std::stoi(argv[2]) : 5;
  const auto L = 2 * len_ladder;
  const auto dim = 1 << L;
  const auto lattice = rokko::create_ladder_lattice_1dim(len_ladder);

  rokko::parallel_sparse_ev solver(library);
  rokko::heisenberg_mfree op(L, lattice);
  const auto map = solver.default_mapping(dim, rokko::mpi_comm{MPI_COMM_WORLD});
  const int num_entries_per_row = lattice.size() + 1;
  rokko::distributed_crs_matrix mat(map, num_entries_per_row);
  rokko::distributed_mfree_to_crs(op, mat);
  mat.output_matrix_market();
  //mat.print();

  solver.finalize();
  MPI_Finalize();
}
