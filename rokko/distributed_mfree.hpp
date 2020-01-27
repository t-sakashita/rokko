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
#include <rokko/mapping_1d_common.hpp>
#include <rokko/skel/mapping_1d.hpp>
#include <rokko/eigen3.hpp>

namespace rokko {

class distributed_mfree : virtual public detail::mapping_1d_common {
public:
  distributed_mfree() = default;
  virtual ~distributed_mfree() = default;

  virtual void multiply(const double *const x, double *const y) const = 0;
  virtual int get_local_offset() const = 0;
  virtual MPI_Comm get_comm() const = 0;
};

class distributed_mfree_default : virtual public distributed_mfree, public skel::mapping_1d {
public:
  distributed_mfree_default() = default;
  distributed_mfree_default(int dim) : distributed_mfree_default(dim, rokko::mpi_comm{MPI_COMM_WORLD}) {}
  distributed_mfree_default(int dim, rokko::mpi_comm const& mpi_comm) : skel::mapping_1d(dim, mpi_comm) {}
  distributed_mfree_default(skel::mapping_1d const& map) : skel::mapping_1d(map) {}

  virtual ~distributed_mfree_default() = default;

  int get_local_offset() const override { return start_row(); }

  MPI_Comm get_comm() const override { return get_mpi_comm().get_comm(); }
};


class distributed_mfree_slepc : public rokko::distributed_mfree {
public:
  virtual void diagonal(double* x) const = 0;
};

} // end namespace rokko

#endif // ROKKO_DISTRIBUTED_MFREE_HPP
