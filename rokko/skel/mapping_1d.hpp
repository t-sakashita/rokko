/*****************************************************************************
*
* Rokko: Integrated Interface for libraries of eigenvalue decomposition
*
* Copyright (C) 2012-2020 by Rokko Developers https://github.com/t-sakashita/rokko
*
* Distributed under the Boost Software License, Version 1.0. (See accompanying
* file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
*
*****************************************************************************/

#ifndef ROKKO_SKEL_MAPPING_1D_HPP
#define ROKKO_SKEL_MAPPING_1D_HPP

#include <rokko/mpi_communicator.hpp>
#include <rokko/mapping_1d.hpp>

namespace rokko {

namespace skel {

class mapping_1d : public detail::ps_mapping_1d_base {
public:
  explicit mapping_1d() : mapping_1d(0) {}
  explicit mapping_1d(int dim) : mapping_1d(dim, mpi_comm{MPI_COMM_WORLD}) {}
  explicit mapping_1d(int dim, mpi_comm const& mpi_comm_in)
    : detail::ps_mapping_1d_base(dim, mpi_comm_in) {}

  int get_num_local_rows() const override { return calculate_num_local_rows(); }
  int calculate_num_local_rows() const {
    return calculate_num_local_rows(get_mpi_comm().get_myrank());
  }
  int calculate_num_local_rows(int proc) const {
    return (get_dim() + get_mpi_comm().get_nprocs() - proc - 1) / get_mpi_comm().get_nprocs();
  }
  int start_row() const override {
    const int nprocs = get_mpi_comm().get_nprocs();
    const int myrank = get_mpi_comm().get_myrank();
    return get_dim() / nprocs * myrank + std::min(get_dim() % nprocs, myrank);
  }
  int end_row() const override {
    return start_row() + get_num_local_rows();
  }
};

} // namespace skel

} // namespace rokko

#endif // ROKKO_SKEL_MAPPING_1D_HPP
