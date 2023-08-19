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
  const int L = (argc >= 3) ? std::stoi(argv[2]) : 10;
  const auto dim = 1 << L;
  std::vector<std::pair<int, int>> lattice;
  for (int i = 0; i < L; ++i) lattice.emplace_back(std::make_pair(i, (i+1) % L));

  rokko::parallel_sparse_ev solver("anasazi");
  rokko::distributed_crs_matrix mat(dim, dim, solver);
  std::vector<double> values;
  std::vector<int> cols;
  for (int row = 0; row < dim; ++row) {
    cols.clear();
    values.clear();
    double diag = 0;
    for (size_t l = 0;  l < lattice.size(); ++l) {
      const auto i = lattice[l].first;
      const auto j = lattice[l].second;
      const auto m1 = 1 << i;
      const auto m2 = 1 << j;
      const auto m3 = m1 + m2;
      if (((row & m3) == m1) || ((row & m3) == m2)) {
        cols.emplace_back(row^m3);
        values.emplace_back(0.5);
        diag += -0.25;
      } else {
        diag += 0.25;
      }
    }
    cols.emplace_back(row);
    values.emplace_back(diag);
    mat.insert(row, cols, values);
  }
  mat.complete();
  //mat.print();
  mat.output_matrix_market();

  solver.finalize();
  MPI_Finalize();
}
