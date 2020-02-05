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
#include <rokko/distributed_mfree.hpp>
#include <rokko/utility/math.hpp>

namespace rokko {

class heisenberg_mfree : public rokko::distributed_mfree_default {
public:
  heisenberg_mfree(int L, const std::vector<std::pair<int, int>>& lattice, MPI_Comm comm = MPI_COMM_WORLD)
    : L_(L), lattice_(lattice), distributed_mfree_default{1 << L, rokko::mpi_comm{comm}} {
    buffer_.resize(get_num_local_rows());
  }
  ~heisenberg_mfree() = default;

  void multiply(const double *const x, double *const y) const override {
    rokko::heisenberg_hamiltonian::multiply(get_comm(), L_, lattice_, x, y, buffer_.data());
  }

private:
  int L_;
  std::vector<std::pair<int, int>> lattice_;
  mutable std::vector<double> buffer_;
};

} // namespace rokko

#endif // ROKKO_UTILITY_HEISENBERG_HAMILTONIAN_MFREE_HPP
