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
#include <rokko/utility/lattice.hpp>

int main(int argc, char *argv[]) {
  int provided;
  MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &provided);
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  std::string name("anasazi");
  if (argc >= 2) name = argv[1];
  int len_ladder = (argc >= 3) ? boost::lexical_cast<int>(argv[2]) : 5;
  int L = 2 * len_ladder;
  int dim = 1 << L;
  std::vector<std::pair<int, int>> lattice;
  rokko::create_ladder_lattice_1dim(len_ladder, lattice);
  //if (rank == 0)
  //  rokko::print_lattice(lattice);

  rokko::parallel_sparse_ev solver(name);
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
    double diag = 0.;
    for (int l = 0;  l < lattice.size(); ++l) {
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
    mat.insert(row, cols, values);
  }
  mat.complete();
  //mat.print();
  mat.output_matrix_market();

  solver.finalize();
  MPI_Finalize();
}
