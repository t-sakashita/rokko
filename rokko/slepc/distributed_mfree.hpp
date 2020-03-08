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

#ifndef ROKKO_SLEPC_DISTRIBUTED_MFREE_HPP
#define ROKKO_SLEPC_DISTRIBUTED_MFREE_HPP

#include <rokko/distributed_mfree.hpp>

#include <rokko/slepc.hpp>

namespace rokko {

namespace slepc {

class distributed_mfree_holder : public rokko::distributed_mfree_holder {
  using rokko::distributed_mfree_holder::distributed_mfree_holder;
public:
  distributed_mfree_holder(std::function<void(const double *const, double *const)> const& multiply_in,
                           std::function<void(double *const)> const& fill_diagonal_in,
                           rokko::skel::mapping_1d const& map)
    : rokko::distributed_mfree_holder(multiply_in, map), fill_diagonal_(fill_diagonal_in) {}

  // For function pointer
  distributed_mfree_holder(void (*multiply_in)(const double *const, double *const),
                           void (*fill_diagonal_in)(double *const),
                           rokko::skel::mapping_1d const& map)
    : rokko::distributed_mfree_holder(multiply_in, map), fill_diagonal_(fill_diagonal_in) {}

  // For C binding
  distributed_mfree_holder(std::function<void(const double *const, double *const, void*)> const& multiply_in,
                           std::function<void(double *const, void*)> const& fill_diagonal_in,
                           void* vars, rokko::skel::mapping_1d const& map)
    : rokko::distributed_mfree_holder(multiply_in, vars, map), fill_diagonal_([fill_diagonal_in, vars](double *const x) {
        fill_diagonal_in(x, vars);
      }) {}

  // For Fortran binding
  distributed_mfree_holder(std::function<void(const int*, const double *const, double *const)> const& multiply_in,
                           std::function<void(const int*, double *const)> const& fill_diagonal_in,
                           rokko::skel::mapping_1d const& map)
    : rokko::distributed_mfree_holder(multiply_in, map), fill_diagonal_([this, fill_diagonal_in](double *const x) {
        int num_local_rows = get_num_local_rows();
        fill_diagonal_in(&num_local_rows, x);
      }) {}

  // Use Eigen3 for Python binding
  distributed_mfree_holder(std::function<void(ConstMapVec,MapVec)> const& multiply_in,
                           std::function<void(MapVec)> const& fill_diagonal_in,
                           rokko::skel::mapping_1d const& map)
    : rokko::distributed_mfree_holder(multiply_in, map), fill_diagonal_([this, fill_diagonal_in](double *const x) {
        const int num_local_rows = get_num_local_rows();
        MapVec  X(x, num_local_rows);
        fill_diagonal_in(X);
      }) {}

  ~distributed_mfree_holder() = default;

  void fill_diagonal(double *const x) const override {
    fill_diagonal_(x);
  }

private:
  std::function<void(double *const)> fill_diagonal_;
};

} // end namespace slepc

} // end namespace rokko

#endif // ROKKO_SLEPC_DISTRIBUTED_MFREE_HPP
