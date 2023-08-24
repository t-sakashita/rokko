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

int main(int argc, char *argv[]) {
  int provided;
  MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &provided);
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  const std::string library = (argc >= 2) ? argv[1] : "anasazi";
  const std::string lattice_file = (argc >= 3) ? argv[2] : "xyz.dat";
  const auto [L, lattice] = rokko::read_lattice_file(lattice_file);
  const auto dim = 1 << L;
  if (rank == 0)
    std::cout << "Eigenvalue decomposition of antiferromagnetic Heisenberg chain" << std::endl
              << "library = " << library << std::endl
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
  //params.set("max_block_size", 5);
  //params.set("max_iters", 500);
  //params.set("conv_tol", 1.0e-8);
  //params.set("num_eigvals", 1);
  const auto info = solver.diagonalize(mat, params);
  MPI_Barrier(MPI_COMM_WORLD);
  const auto end_tick = MPI_Wtime();

  const auto num_conv = info.get<int>("num_conv");
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
  }

  solver.finalize();
  MPI_Finalize();
}
