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


template <typename T>
auto get_num_entries_per_row(std::vector<std::vector<T>> const& vecs) {
  return std::max_element(vecs.cbegin(), vecs.cend(),
                          [] (auto const& a, auto const& b) {
                            return a.size() < b.size();
                          })->size();
}


int main(int argc, char *argv[]) {
  int provided;
  MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &provided);
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  const std::string library = (argc >= 2) ? argv[1] : rokko::parallel_sparse_ev::default_solver();

  constexpr int dim = 8;
  const std::vector<std::vector<int>> cols = {{0, 4}, {3}, {5}, {1, 7}, {0, 5, 6}, {2, 4}, {4, 7}, {3, 6}};
  const std::vector<std::vector<double>> values = {{7.1, 2.8}, {6.4}, {0.5}, {6.4, 3.5}, {2.8, 0.2, 1.4}, {0.5, 0.2}, {1.4, 4.3}, {3.5, 4.3}};

  const auto num_entries_per_row = get_num_entries_per_row(cols);

  if (rank == 0) std::cout << "[solver = " << library << "]" << std::endl;
  rokko::parallel_sparse_ev solver(library);
  const auto map = solver.default_mapping(dim, rokko::mpi_comm{MPI_COMM_WORLD});
  rokko::distributed_crs_matrix mat(map, num_entries_per_row);
  for (int row = map.start_row(); row < map.end_row(); ++row) {
    mat.insert(row, cols[row], values[row]);
  }
  mat.complete();
  mat.print();

  solver.finalize();
  MPI_Finalize();
}
