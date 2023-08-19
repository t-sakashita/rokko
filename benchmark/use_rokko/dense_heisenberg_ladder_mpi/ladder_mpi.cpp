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
#include <rokko/mapping_bc.hpp>
#include <rokko/collective.hpp>
#include <rokko/utility/solver_name.hpp>
#include <rokko/utility/heisenberg_hamiltonian_mpi.hpp>
#include <rokko/utility/lattice.hpp>
#include <rokko/utility/machine_info.hpp>
#include <iostream>


using matrix_major = rokko::matrix_col_major;

int main(int argc, char *argv[]) {
  int provided;
  MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &provided);
  MPI_Comm comm = MPI_COMM_WORLD;
  const std::string library_routine = (argc >= 2) ? argv[1] : rokko::parallel_dense_ev::default_solver();
  const auto [library, routine] = rokko::split_solver_name(library_routine);

  const int len_ladder = (argc >= 3) ? std::stoi(argv[2]) : 5;
  const auto L = 2 * len_ladder;
  const auto lattice = rokko::create_ladder_lattice_1dim(len_ladder);
  const auto dim = 1 << L;

  rokko::grid g(comm);
  const auto myrank = g.get_myrank();
  if (myrank == 0) {
    std::cout << "Eigenvalue decomposition of Frank matrix" << std::endl
              << "library = " << library << std::endl
              << "routine = " << routine << std::endl
              << "dimension = " << dim << std::endl;
    rokko::print_lattice(lattice);
  }

  MPI_Barrier(MPI_COMM_WORLD);
  const auto init_tick = MPI_Wtime();
  rokko::parallel_dense_ev solver(library);
  solver.initialize(argc, argv);
  MPI_Barrier(MPI_COMM_WORLD);
  const auto initend_tick = MPI_Wtime();

  MPI_Barrier(MPI_COMM_WORLD);
  const auto gen_tick = MPI_Wtime();
  const rokko::mapping_bc<matrix_major> map = solver.default_mapping(dim, g);
  rokko::distributed_matrix<double, matrix_major> mat(map);
  rokko::heisenberg_hamiltonian::generate(L, lattice, mat);
  //  Eigen::MatrixXd mat_loc(dim, dim);
  //  rokko::gather(mat, mat_loc, 0);

  MPI_Barrier(MPI_COMM_WORLD);
  const auto diag_tick = MPI_Wtime();
  Eigen::VectorXd eigval(dim);
  rokko::distributed_matrix<double, matrix_major> eigvec(map);
  rokko::parameters params;
  params.set("routine", "tri");
  solver.diagonalize(mat, eigval, eigvec); //, params);
  MPI_Barrier(MPI_COMM_WORLD);
  const auto end_tick = MPI_Wtime();

  if (myrank == 0) {
    std::cout << "init_time = " << initend_tick - init_tick << std::endl
              << "gen_time = " << diag_tick - gen_tick << std::endl
              << "diag_time = " << end_tick - diag_tick << std::endl;
    rokko::machine_info();
    bool sorted = true;
    for (unsigned int i = 1; i < dim; ++i) sorted &= (eigval(i-1) <= eigval(i));
    if (!sorted) std::cout << "Warning: eigenvalues are not sorted in ascending order!\n";
    std::cout << "largest eigenvalues:";
    for (int i = 0; i < dim; ++i) std::cout << ' ' << eigval(dim - 1 - i);
    //for (int i = 0; i < std::min(dim, 10); ++i) std::cout << ' ' << eigval(dim - 1 - i);
    std::cout << std::endl;
  }

  solver.finalize();
  MPI_Finalize();
}
