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
#include <rokko/eigen3.hpp>

namespace rokko {

namespace skel {

class mapping_1d : public detail::ps_mapping_1d_base {
public:
  explicit mapping_1d() : mapping_1d(0) {}
  explicit mapping_1d(int dim) : mapping_1d(dim, mpi_comm{MPI_COMM_WORLD}) {}
  explicit mapping_1d(int dim, mpi_comm const& mpi_comm_in)
    : detail::ps_mapping_1d_base(dim, mpi_comm_in), num_local_rows_(calculate_num_local_rows()), start_row_(calculate_start_row()), end_row_(start_row_ + num_local_rows_) {}

  explicit mapping_1d(int dim, int num_local_rows, mpi_comm const& mpi_comm_in)
    : detail::ps_mapping_1d_base(dim, mpi_comm_in),
    num_local_rows_(num_local_rows), start_row_(calculate_start_row_by_gather(num_local_rows)), end_row_(start_row_ + num_local_rows_) {}

  int get_num_local_rows() const override { return num_local_rows_; }
  int start_row() const override { return start_row_; }
  int end_row() const override { return end_row_; }

  int calculate_start_row_by_gather(int num_local_rows) const {
    const int myrank = get_mpi_comm().get_myrank();
    Eigen::VectorXi target(get_mpi_comm().get_nprocs());

    MPI_Allgather(&num_local_rows, 1, MPI_INT,
                  target.data(), 1, MPI_INT,
                  get_mpi_comm().get_comm());

    int sum = 0;
    for (int rank = 0; rank < myrank; ++rank)
      sum += target(rank);
    return sum;
  }

  int calculate_num_local_rows() const {
    return calculate_num_local_rows(get_mpi_comm().get_myrank());
  }
  int calculate_num_local_rows(int proc) const {
    return (get_dim() + get_mpi_comm().get_nprocs() - proc - 1) / get_mpi_comm().get_nprocs();
  }

  int calculate_start_row() {
    const int nprocs = get_mpi_comm().get_nprocs();
    const int myrank = get_mpi_comm().get_myrank();
    return get_dim() / nprocs * myrank + std::min(get_dim() % nprocs, myrank);
  }

private:
  int num_local_rows_;
  int start_row_, end_row_;
};

} // namespace skel

} // namespace rokko

#endif // ROKKO_SKEL_MAPPING_1D_HPP
