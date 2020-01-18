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

#ifndef ROKKO_DISTRIBUTED_MFREE_HPP
#define ROKKO_DISTRIBUTED_MFREE_HPP

#include <rokko/mapping_1d.hpp>

namespace rokko {

class distributed_mfree {
public:
  distributed_mfree() {}
  virtual ~distributed_mfree() {}

  virtual void multiply(const double *const x, double *const y) const = 0;
  virtual int get_dim() const = 0;
  virtual int get_local_offset() const = 0;
  virtual int get_num_local_rows() const = 0;
  virtual MPI_Comm get_comm() const = 0;
};


class distributed_mfree_default : public distributed_mfree {
public:
  distributed_mfree_default() {}
  distributed_mfree_default(int dim) : distributed_mfree_default(dim, rokko::mpi_comm{MPI_COMM_WORLD}) {}
  distributed_mfree_default(int dim, rokko::mpi_comm const& mpi_comm) : dim_(dim), mpi_comm_(mpi_comm) {}

  virtual ~distributed_mfree_default() {}

  virtual void multiply(const double *const x, double *const y) const = 0;

  int get_dim() const { return dim_; }
  int get_num_local_rows() const { return calculate_num_local_rows(); }
  int calculate_num_local_rows() const {
    return calculate_num_local_rows(get_mpi_comm().get_myrank());
  }
  int calculate_num_local_rows(int proc) const {
    return (get_dim() + get_mpi_comm().get_nprocs() - proc - 1) / get_mpi_comm().get_nprocs();
  }
  int get_local_offset() const { return start_row(); }
  int start_row() const {
    const int nprocs = get_mpi_comm().get_nprocs();
    const int myrank = get_mpi_comm().get_myrank();
    return get_dim() / nprocs * myrank + std::min(get_dim() % nprocs, myrank);
  }
  int end_row() const {
    return start_row() + get_num_local_rows();
  }

  const rokko::mpi_comm& get_mpi_comm() const { return mpi_comm_; }

  MPI_Comm get_comm() const { return get_mpi_comm().get_comm(); }

protected:
  int dim_;
  rokko::mpi_comm mpi_comm_;
};


class distributed_mfree_slepc : public rokko::distributed_mfree {
public:
  virtual void diagonal(double* x) const = 0;
};

} // end namespace rokko

#endif // ROKKO_DISTRIBUTED_MFREE_HPP
