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

#ifndef ROKKO_UTILITY_HEISENBERG_HAMILTONIAN_MFREE_HPP
#define ROKKO_UTILITY_HEISENBERG_HAMILTONIAN_MFREE_HPP

#include <rokko/utility/heisenberg_hamiltonian_mpi.hpp>

namespace rokko {

class heisenberg_mfree : public rokko::distributed_mfree {
public:
  heisenberg_mfree(int L, const std::vector<std::pair<int, int>>& lattice)
    : L_(L), lattice_(lattice) {
    comm_ = MPI_COMM_WORLD;
    int size, rank;
    MPI_Comm_size(comm_, &size);
    MPI_Comm_rank(comm_, &rank);
    int p = find_power_of_two(size);
    dim_ = 1 << L;
    num_local_rows_ = 1 << (L-p);
    local_offset_ = num_local_rows_ * rank;
    buffer_.assign(num_local_rows_, 0);
  }
  ~heisenberg_mfree() {}

  int find_power_of_two(int n) {
    int p = -1;
    do {
      n /= 2;
      ++p;
    } while (n > 0);
    return p;
  }

  void multiply(const double *const x, double *const y) const {
    rokko::heisenberg_hamiltonian::multiply(comm_, L_, lattice_, x, y, buffer_.data());
  }
  int get_dim() const { return dim_; }
  int get_local_offset() const { return local_offset_; }
  int get_num_local_rows() const { return num_local_rows_; }

private:
  MPI_Comm comm_;
  int L_;
  std::vector<std::pair<int, int>> lattice_;
  int dim_, local_offset_, num_local_rows_;
  mutable std::vector<double> buffer_;
};

} // namespace rokko

#endif // ROKKO_UTILITY_HEISENBERG_HAMILTONIAN_MFREE_HPP
