/*****************************************************************************
*
* Rokko: Integrated Interface for libraries of eigenvalue decomposition
*
* Copyright (C) 2012-2016 by Rokko Developers https://github.com/t-sakashita/rokko
*
* Distributed under the Boost Software License, Version 1.0. (See accompanying
* file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
*
*****************************************************************************/

#ifndef ROKKO_DISTRIBUTED_CRS_MATRIX_HPP
#define ROKKO_DISTRIBUTED_CRS_MATRIX_HPP

#include <rokko/mapping_1d.hpp>
#include <vector>

namespace rokko {
namespace detail {

class ps_crs_base {
public:
  virtual ~ps_crs_base() {}
  virtual void initialize(rokko::mapping_1d const& map, int num_entries_per_row) = 0;
  virtual void insert(int row, std::vector<int> const& cols, std::vector<double> const& values) = 0;
  virtual void insert(int row, int col_size, int* cols, double* const values) = 0;
  virtual void complete() = 0;
  virtual int get_dim() const = 0;
  virtual int num_local_rows() const = 0;
  virtual int start_row() const = 0;
  virtual int end_row() const = 0;
  virtual int get_nnz() const = 0;
  virtual void print() const = 0;
  virtual void output_matrix_market() const = 0;
  virtual ps_crs_base* get_impl() = 0;
  virtual const ps_crs_base* get_impl() const = 0;
};


template<typename CRS>
class ps_crs_wrapper : public ps_crs_base {
  using crs_type = CRS;
public:
  ps_crs_wrapper() : crs_impl_() {}
  virtual ~ps_crs_wrapper() {}
  void initialize(rokko::mapping_1d const& map, int num_entries_per_row) { crs_impl_.initialize(map, num_entries_per_row); }
  void insert(int row, std::vector<int> const& cols, std::vector<double> const& values) {
    crs_impl_.insert(row, cols, values); }
  void insert(int row, int col_size, int* cols, double* const values) { crs_impl_.insert(row, col_size, cols, values); }
  void complete() { crs_impl_.complete(); }
  int get_dim() const { return crs_impl_.get_dim(); }
  int num_local_rows() const { return crs_impl_.num_local_rows(); }
  int start_row() const { return crs_impl_.start_row(); }
  int end_row() const { return crs_impl_.end_row(); }
  int get_nnz() const { return crs_impl_.get_nnz(); }
  void print() const { crs_impl_.print(); }
  void output_matrix_market() const { crs_impl_.print(); }
  MPI_Comm get_comm() const { return crs_impl_.get_comm(); }
  const ps_crs_base* get_impl() const { return &crs_impl_; }
  ps_crs_base* get_impl() { return &crs_impl_; }

private:
  crs_type crs_impl_;
};


using ps_crs_factory = factory<ps_crs_base>;

} // end namespace detail

class distributed_crs_matrix {
public:
  explicit distributed_crs_matrix(rokko::mapping_1d const& map, int num_entries_per_row)
    : solver_name_(map.get_solver_name()), crs_impl_(detail::ps_crs_factory::instance()->make_product(solver_name_)) {
    crs_impl_->initialize(map, num_entries_per_row);
  }
  template<typename SOLVER>
  explicit distributed_crs_matrix(int row_dim, int col_dim, SOLVER& solver_in) : solver_name_(solver_in.get_solver_name()) {
    crs_impl_ = solver_in.create_distributed_crs_matrix(row_dim, col_dim);
  }
  template<typename SOLVER>
  explicit distributed_crs_matrix(int row_dim, int col_dim, int num_entries_per_row, SOLVER& solver_in) : solver_name_(solver_in.get_solver_name()) {
    crs_impl_ = solver_in.create_distributed_crs_matrix(row_dim, col_dim, num_entries_per_row);
  }
  std::string get_solver_name() const { return solver_name_; }
  void insert(int row, std::vector<int> const& cols, std::vector<double> const& values) {
    crs_impl_->insert(row, cols, values);
  }
  void insert(int row, int col_size, int* cols, double* const values) {
    crs_impl_->insert(row, col_size, cols, values);
  }
  void complete() const {
    crs_impl_->complete();
  }
  int get_dim() const {
    return crs_impl_->get_dim();
  }
  int num_local_rows() const {
    return crs_impl_->num_local_rows();
  }
  int start_row() const {
    return crs_impl_->start_row();
  }
  int end_row() const {
    return crs_impl_->end_row();
  }
  int get_nnz() const {
    return crs_impl_->get_nnz();
  }
  void print() const {
    crs_impl_->print();
  }
  void output_matrix_market() const {
    crs_impl_->output_matrix_market();
  }
  const detail::ps_crs_factory::product_pointer_type get_ptr() const { return crs_impl_; }

private:
  std::string solver_name_;
  detail::ps_crs_factory::product_pointer_type crs_impl_;
};

} // end namespace rokko

#define ROKKO_REGISTER_PARALLEL_SPARSE_CRS(crs, name, priority) \
namespace { namespace BOOST_JOIN(register, __LINE__) { \
struct register_caller { \
  using factory = rokko::factory<rokko::detail::ps_crs_base>;  \
  using product = rokko::detail::ps_crs_wrapper<crs>; \
  register_caller() { factory::instance()->register_creator<product>(name, priority); } \
} caller; \
} }

#endif // ROKKO_DISTRIBUTED_CRS_MATRIX_HPP
