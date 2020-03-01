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

using MapVec = Eigen::Map<Eigen::Vector<double>>;
using ConstMapVec = const Eigen::Map<const Eigen::Vector<double>>;

class distributed_mfree_holder : public rokko::distributed_mfree_default {
public:
  distributed_mfree_holder(std::function<void(const double *const, double *const)> const& multiply_in, rokko::skel::mapping_1d const& map)
    : multiply_(multiply_in), rokko::distributed_mfree_default(map) {}

  // For function pointer
  distributed_mfree_holder(void (*multiply_in)(const double *const, double *const), rokko::skel::mapping_1d const& map)
    : distributed_mfree_holder(std::function<void(const double *const, double *const)>(multiply_in), map) {}

  // For C binding
  distributed_mfree_holder(std::function<void(const double *const, double *const, void*)> const& multiply_in, void* vars, rokko::skel::mapping_1d const& map)
    : distributed_mfree_holder([multiply_in, vars](const double *const x, double *const y) {
        multiply_in(x, y, vars);
      },
      map) {}

  // For Fortran binding
  distributed_mfree_holder(std::function<void(const int*, const double *const, double *const)> const& multiply_in, rokko::skel::mapping_1d const& map)
    : distributed_mfree_holder([this, multiply_in](const double *const x, double *const y) {
        int num_local_rows = get_num_local_rows();
        multiply_in(&num_local_rows, x, y);
      },
      map) {}

  // Use Eigen3 for Python binding
  distributed_mfree_holder(std::function<void(ConstMapVec,MapVec)> const& multiply_in, rokko::skel::mapping_1d const& map)
    : multiply_([this, multiply_in](const double *const x, double *const y) {
        const int num_local_rows = get_num_local_rows();
        ConstMapVec  X(x, num_local_rows);
        MapVec  Y(y, num_local_rows);
        multiply_in(X, Y);
      }), rokko::distributed_mfree_default(map) {}

  ~distributed_mfree_holder() = default;

  void multiply(const double *const x, double *const y) const override {
    multiply_(x, y);
  }

private:
  std::function<void(const double *const, double *const)> multiply_;
};

} // end namespace rokko

#endif // ROKKO_DISTRIBUTED_MFREE_HPP
