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

#ifndef ROKKO_MAPPING_1D_HPP
#define ROKKO_MAPPING_1D_HPP

#include <rokko/factory.hpp>
#include <rokko/mpi_communicator.hpp>

namespace rokko {

namespace detail {

class ps_mapping_1d_base {
public:
  explicit ps_mapping_1d_base(int dim) {}
  explicit ps_mapping_1d_base(int dim, mpi_comm const& mpi_comm_in) : dim_(dim), mpi_comm_(mpi_comm_in) {}
  ps_mapping_1d_base() {}
  virtual ~ps_mapping_1d_base() {}
  virtual void init(int dim, mpi_comm const& mpi_comm_in) = 0;
  int get_dim() const { return dim_; }
  virtual int get_num_local_rows() const = 0;
  virtual int start_row() const = 0;
  virtual int end_row() const = 0;
  virtual const ps_mapping_1d_base* get_impl() const = 0;
  const rokko::mpi_comm& get_mpi_comm() const { return mpi_comm_; }
  void set_dim(int dim) { dim_ = dim; }
  void set_mpi_comm(rokko::mpi_comm const& mpi_comm_in) { mpi_comm_ = mpi_comm_in; }

private:
  int dim_;
  rokko::mpi_comm mpi_comm_;
};


template<typename MAP>
class ps_mapping_1d_wrapper : public ps_mapping_1d_base {
public:
  using mapping_1d_type = MAP;

  explicit ps_mapping_1d_wrapper(int dim) : map_impl_() {}
  explicit ps_mapping_1d_wrapper(int dim, mpi_comm const& mpi_comm_in) : map_impl_() {}
  ps_mapping_1d_wrapper() : map_impl_() {}
  ~ps_mapping_1d_wrapper() {}
  void init(int dim, mpi_comm const& mpi_comm_in) { map_impl_.init(dim, mpi_comm_in); }
  int get_dim() const { return map_impl_.get_dim(); }
  int get_num_local_rows() const { return map_impl_.get_num_local_rows(); }
  int start_row() const { return map_impl_.start_row(); }
  int end_row() const { return map_impl_.end_row(); }

  const ps_mapping_1d_base* get_impl() const { return &map_impl_; }

private:
  mapping_1d_type map_impl_;
};

using ps_mapping_1d_factory = factory<ps_mapping_1d_base>;

} // end namespace detail

class mapping_1d {
public:
  template<typename SOLVER>
  mapping_1d(int dim, mpi_comm const& mpi_comm_in, SOLVER& solver_in)
    : mapping_1d(dim, mpi_comm_in, solver_in.get_solver_name()) {
  }

  mapping_1d(int dim, mpi_comm const& mpi_comm_in, std::string const& solver_name)
    : solver_name_(solver_name), map_impl_(detail::ps_mapping_1d_factory::instance()->make_product(solver_name)) {
    map_impl_->init(dim, mpi_comm_in);
  }

  mapping_1d(int dim, mpi_comm const& mpi_comm_in)
    : mapping_1d(dim, mpi_comm_in, detail::ps_mapping_1d_factory::instance()->default_product_name()) {}

  int get_dim() const { return map_impl_->get_dim(); }
  int num_local_rows() const { return map_impl_->get_num_local_rows(); }
  int start_row() const { return map_impl_->start_row(); }
  int end_row() const { return map_impl_->end_row(); }
  const rokko::mpi_comm& get_mpi_comm() const { return map_impl_->get_mpi_comm(); }

  std::string get_solver_name() const { return solver_name_; }
  void set_solver_name(std::string const& solver_name) { solver_name_ = solver_name; }

  const detail::ps_mapping_1d_factory::product_pointer_type get_ptr() const { return map_impl_; }

private:
  std::string solver_name_;
  detail::ps_mapping_1d_factory::product_pointer_type map_impl_;
};

} // end namespace rokko

#define ROKKO_REGISTER_PARALLEL_SPARSE_MAPPING_1D(map, name, priority) \
namespace { namespace BOOST_JOIN(register, __LINE__) { \
struct register_caller { \
  using factory = rokko::factory<rokko::detail::ps_mapping_1d_base>;  \
  using product = rokko::detail::ps_mapping_1d_wrapper<map>; \
  register_caller() { factory::instance()->register_creator<product>(name, priority); } \
} caller; \
} }

#endif // ROKKO_MAPPING_1D_HPP
