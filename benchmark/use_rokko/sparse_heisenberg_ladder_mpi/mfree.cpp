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
#include <rokko/utility/solver_name.hpp>
#include <rokko/utility/heisenberg_hamiltonian_mfree.hpp>
#include <rokko/utility/lattice.hpp>
#include <rokko/utility/machine_info.hpp>

int main(int argc, char *argv[]) {
  int provided;
  MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &provided);
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  const std::string library_routine = (argc >= 2) ? argv[1] : rokko::parallel_sparse_ev::default_solver();
  const auto [library, routine] = rokko::split_solver_name(library_routine);
  const int len_ladder = (argc >= 3) ? std::stoi(argv[2]) : 5;

  const auto L = 2 * len_ladder;
  const auto lattice = rokko::create_ladder_lattice_1dim(len_ladder);
  if (rank == 0)
    rokko::print_lattice(lattice);
  const int dim = 1 << L;
  if (rank == 0)
    std::cout << "Eigenvalue decomposition of antiferromagnetic Heisenberg 1D ladder lattice" << std::endl
      	      << "solver = " << library << std::endl
              << "routine = " << routine << std::endl
              << "L = " << L << std::endl
              << "dimension = " << dim << std::endl;

  MPI_Barrier(MPI_COMM_WORLD);
  const auto init_tick = MPI_Wtime();
  rokko::parallel_sparse_ev solver(library);
  MPI_Barrier(MPI_COMM_WORLD);
  const auto initend_tick = MPI_Wtime();
  
  MPI_Barrier(MPI_COMM_WORLD);
  const auto gen_tick = MPI_Wtime();
  const rokko::heisenberg_mfree mat(L, lattice);
  
  MPI_Barrier(MPI_COMM_WORLD);
  const auto diag_tick = MPI_Wtime();
  rokko::parameters params;
  params.set("routine", routine);
  //params.set("max_block_size", 5);
  params.set("max_iters", 200);
  params.set("conv_tol", 1.0e-8);
  //params.set("num_eigvals", 1)
  const auto params_out = solver.diagonalize(mat, params);
  MPI_Barrier(MPI_COMM_WORLD);
  const auto end_tick = MPI_Wtime();

  const auto num_conv = params_out.get<int>("num_conv");
  if (num_conv == 0) {
    throw std::runtime_error("diagonalize : solver does not converge.");
  }
  std::vector<double> eigvec;
  solver.eigenvector(0, eigvec);
  if (rank == 0) {
    std::cout << "number of converged eigenpairs = " << num_conv << std::endl;
    std::cout << "smallest eigenvalues: ";
    for (int i = 0; i < num_conv; ++i) std::cout << ' ' << solver.eigenvalue(i);
    std::cout << std::endl;
    std::cout << "init_time = " << initend_tick - init_tick << std::endl
              << "gen_time = " << diag_tick - gen_tick << std::endl
              << "diag_time = " << end_tick - diag_tick << std::endl;
    rokko::machine_info();
  }

  solver.finalize();
  MPI_Finalize();
}
