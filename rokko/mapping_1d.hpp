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

#pragma once

#include <rokko/factory.hpp>
#include <rokko/mpi_communicator.hpp>
#include <rokko/mapping_1d_common.hpp>
#include <rokko/utility/macro_join.hpp>

namespace rokko {

namespace detail {

class ps_mapping_1d_base : virtual public mapping_1d_common {
public:
  explicit ps_mapping_1d_base(int dim) : ps_mapping_1d_base(dim, mpi_comm(MPI_COMM_WORLD)) {}
  explicit ps_mapping_1d_base(int dim, mpi_comm const& mpi_comm_in) : dim_(dim), mpi_comm_(mpi_comm_in) {}
  ps_mapping_1d_base() = default;
  virtual ~ps_mapping_1d_base() = default;
  int get_dim() const override { return dim_; }
  virtual int start_row() const = 0;
  virtual int end_row() const = 0;
  const rokko::mpi_comm& get_mpi_comm() const { return mpi_comm_; }

private:
  int dim_;
  rokko::mpi_comm mpi_comm_;
};

using ps_mapping_1d_factory = factory<ps_mapping_1d_base,int,mpi_comm const&>;
using ps_mapping_1d_factory_num = factory<ps_mapping_1d_base,int,int,mpi_comm const&>;

} // end namespace detail

class mapping_1d {
public:
  mapping_1d(int dim, mpi_comm const& mpi_comm_in, std::string const& solver_name)
    : solver_name_(solver_name), map_impl_(detail::ps_mapping_1d_factory::instance()->make_product(solver_name, dim, mpi_comm_in)) {}

  mapping_1d(int dim, mpi_comm const& mpi_comm_in)
    : mapping_1d(dim, mpi_comm_in, detail::ps_mapping_1d_factory::instance()->default_product_name()) {}

  mapping_1d(int dim, int num_local_rows, mpi_comm const& mpi_comm_in, std::string const& solver_name)
    : solver_name_(solver_name), map_impl_(detail::ps_mapping_1d_factory_num::instance()->make_product(solver_name, dim, num_local_rows, mpi_comm_in)) {}

  mapping_1d(int dim, int num_local_rows, mpi_comm const& mpi_comm_in)
  : mapping_1d(dim, num_local_rows, mpi_comm_in, detail::ps_mapping_1d_factory::instance()->default_product_name()) {}

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
namespace { namespace ROKKO_JOIN(register, __LINE__) { \
struct register_caller { \
  using factory = rokko::detail::ps_mapping_1d_factory; \
  using product = map; \
  register_caller() { factory::instance()->register_creator<product>(name, priority); } \
} caller; \
} } \
 \
namespace { namespace ROKKO_JOIN(register2, __LINE__) { \
struct register_caller { \
  using factory = rokko::detail::ps_mapping_1d_factory_num; \
  using product = map; \
  register_caller() { factory::instance()->register_creator<product>(name, priority); } \
} caller; \
} }
