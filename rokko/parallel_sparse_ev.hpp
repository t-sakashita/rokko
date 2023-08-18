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

#pragma once

#include <rokko/factory.hpp>
#include <rokko/mapping_1d.hpp>
#include <rokko/distributed_crs_matrix.hpp>
#include <rokko/distributed_mfree.hpp>
#include <rokko/distributed_vector.hpp>
#include <rokko/utility/macro_join.hpp>
#include <rokko/utility/timer.hpp>
#include <rokko/parameters.hpp>

namespace rokko {

namespace detail {
    
class ps_ev_base {
public:
  virtual ~ps_ev_base() = default;
  virtual void initialize(int& argc, char**& argv) = 0;
  virtual void finalize() = 0;
  virtual parameters diagonalize(const rokko::distributed_crs_matrix& mat, rokko::parameters const& params) = 0;
  virtual parameters diagonalize(const rokko::distributed_mfree& mat, rokko::parameters const& params) = 0;
  virtual double eigenvalue(int k) const = 0;
  virtual void eigenvector(int k, std::vector<double>& vec) const = 0;
  virtual void eigenvector(int k, double* vec) const = 0;
  virtual void eigenvector(int k, distributed_vector<double>& vec) const = 0;
  virtual int get_num_conv() const = 0;
};

template<typename SOLVER>
class ps_ev_wrapper : public ps_ev_base {
  using solver_type = SOLVER;
public:
  ps_ev_wrapper() : solver_impl_() {}
  virtual ~ps_ev_wrapper() = default;
  void initialize(int& argc, char**& argv) override { solver_impl_.initialize(argc, argv); }
  void finalize() override { solver_impl_.finalize(); }
  parameters diagonalize(const rokko::distributed_crs_matrix& mat, rokko::parameters const& params) override {
    return solver_impl_.diagonalize(mat, params);
  }
  parameters diagonalize(const rokko::distributed_mfree& mat, rokko::parameters const& params) override {
    return solver_impl_.diagonalize(mat, params);
  }
  double eigenvalue(int k) const override { return solver_impl_.eigenvalue(k); }
  void eigenvector(int k, std::vector<double>& vec) const override { solver_impl_.eigenvector(k, vec); }
  void eigenvector(int k, double* vec) const override { solver_impl_.eigenvector(k, vec); }
  void eigenvector(int k, distributed_vector<double>& vec) const override { solver_impl_.eigenvector(k, vec); }
  int get_num_conv() const override { return solver_impl_.get_num_conv(); }

private:
  solver_type solver_impl_;
};

using ps_solver_factory = factory<ps_ev_base>;
  
} // end namespace detail

class parallel_sparse_ev {
public:
  parallel_sparse_ev(std::string const& solver_name)
    : solver_impl_(detail::ps_solver_factory::instance().make_product(solver_name)), solver_name_(solver_name) {}

  parallel_sparse_ev() : parallel_sparse_ev(default_solver()) {}

  void initialize(int& argc, char**& argv) {
    solver_impl_->initialize(argc, argv);
  }
  void finalize() {
    solver_impl_->finalize();
  }
  rokko::mapping_1d default_mapping(int dim, mpi_comm const& mpi_comm_in) const {
    return rokko::mapping_1d(dim, mpi_comm_in, get_solver_name());
  }
  rokko::mapping_1d default_mapping(int dim, int num_local_rows, mpi_comm const& mpi_comm_in) const {
    return rokko::mapping_1d(dim, num_local_rows, mpi_comm_in, get_solver_name());
  }

  parameters diagonalize(const rokko::distributed_crs_matrix& mat, rokko::parameters const& params) {
    if(get_solver_name() != mat.get_solver_name())
      throw std::invalid_argument(get_solver_name() + "'s diagonalize() : " + mat.get_solver_name() + "'s distributed_crs_matrix is given.");

    return solver_impl_->diagonalize(mat, params);
  }
  parameters diagonalize(const rokko::distributed_mfree& mat, rokko::parameters const& params) {
    return solver_impl_->diagonalize(mat, params);
  }
  double eigenvalue(int k) const { return solver_impl_->eigenvalue(k); }
  void eigenvector(int k, std::vector<double>& vec) const { solver_impl_->eigenvector(k, vec); }
  void eigenvector(int k, double* vec) const { solver_impl_->eigenvector(k, vec); }
  void eigenvector(int k, distributed_vector<double>& vec) const { solver_impl_->eigenvector(k, vec); }
  int get_num_conv() const { return solver_impl_->get_num_conv(); }
  std::string get_solver_name() const { return solver_name_; }
  static std::vector<std::string> solvers() { return detail::ps_solver_factory::product_names(); }
  static const std::string& default_solver() { return detail::ps_solver_factory::default_product_name(); }
private:
  detail::ps_solver_factory::product_pointer_type solver_impl_;
  std::string solver_name_;
};

} // end namespace rokko

#define ROKKO_REGISTER_PARALLEL_SPARSE_SOLVER(solver, name, priority) \
namespace { namespace ROKKO_JOIN(register, __LINE__) { \
struct register_caller { \
  using factory = rokko::detail::ps_solver_factory;     \
  using product = rokko::detail::ps_ev_wrapper<solver>; \
  register_caller() { factory::instance().register_creator<product>(name, priority); } \
} caller; \
} }
