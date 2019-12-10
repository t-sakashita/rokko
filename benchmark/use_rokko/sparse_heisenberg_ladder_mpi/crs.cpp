/*****************************************************************************
*
* Rokko: Integrated Interface for libraries of eigenvalue decomposition
*
* Copyright (C) 2012-2016 Rokko Developers https://github.com/t-sakashita/rokko
*
* Distributed under the Boost Software License, Version 1.0. (See accompanying
* file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
*
*****************************************************************************/

#include <rokko/rokko.hpp>
#include <rokko/utility/solver_name.hpp>
#include <rokko/utility/lattice.hpp>
#include <rokko/utility/machine_info.hpp>

int main(int argc, char *argv[]) {
  int provided;
  MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &provided);
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  double init_tick, initend_tick, gen_tick, diag_tick, end_tick;

  std::string library_routine(rokko::parallel_sparse_ev::default_solver());
  if (argc >= 2) library_routine = argv[1];
  std::string library, routine;
  rokko::split_solver_name(library_routine, library, routine);
  int len_ladder = 5;
  if (argc >= 3) len_ladder = boost::lexical_cast<int>(argv[2]);

  int L = 2 * len_ladder;
  std::vector<std::pair<int, int>> lattice;
  rokko::create_ladder_lattice_1dim(len_ladder, lattice);
  if (rank == 0)
    rokko::print_lattice(lattice);
  int dim = 1 << L;
  if (rank == 0)
    std::cout << "Eigenvalue decomposition of antiferromagnetic Heisenberg 1D ladder lattice" << std::endl
	      << "solver = " << library << std::endl
	      << "routine = " << routine << std::endl
	      << "L = " << L << std::endl
	      << "dimension = " << dim << std::endl;

  MPI_Barrier(MPI_COMM_WORLD);
  init_tick = MPI_Wtime();
  rokko::parallel_sparse_ev solver(library);
  MPI_Barrier(MPI_COMM_WORLD);
  initend_tick = MPI_Wtime();

  MPI_Barrier(MPI_COMM_WORLD);
  gen_tick = MPI_Wtime();
  int num_entries_per_row = 2 * L;
  rokko::distributed_crs_matrix mat(dim, dim, num_entries_per_row, solver);
  std::vector<double> values;
  std::vector<int> cols;
  for (int row = mat.start_row(); row < mat.end_row(); ++row) {
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
    if (diag != 0.) {
      cols.push_back(row);
      values.push_back(diag);
    }
    //std::cout << "num_cols=" << cols.size() << std::endl;
    mat.insert(row, cols, values);
  }
  mat.complete();
  //mat.print();
  
  MPI_Barrier(MPI_COMM_WORLD);
  diag_tick = MPI_Wtime();
  rokko::parameters params;
  params.set("routine", routine);
  //params.set("max_block_size", 5);
  params.set("max_iters", 200);
  params.set("conv_tol", 1.0e-8);
  //params.set("num_eigvals", 1);
  rokko::parameters params_out = solver.diagonalize(mat, params);
  MPI_Barrier(MPI_COMM_WORLD);
  end_tick = MPI_Wtime();
  
  int num_conv = params_out.get<int>("num_conv");
  if (num_conv == 0) MPI_Abort(MPI_COMM_WORLD, -1);
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
