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
#include <rokko/utility/tuple_to_array.hpp>
#include <vector>

namespace rokko {
namespace detail {

class ps_crs_base {
public:
  virtual ~ps_crs_base() = default;
  virtual void initialize(rokko::mapping_1d const& map, int num_entries_per_row) = 0;
  virtual void insert(int row, std::vector<int> const& cols, std::vector<double> const& values) = 0;
  virtual void insert(int row, int col_size, int const*const cols, double const*const values) = 0;
  virtual void complete() = 0;
  virtual void extract(int row, std::vector<int>& cols, std::vector<double>& values) const = 0;
  virtual int get_dim() const = 0;
  virtual int num_local_rows() const = 0;
  virtual int start_row() const = 0;
  virtual int end_row() const = 0;
  virtual int get_nnz() const = 0;
  virtual void print() const = 0;
  virtual void output_matrix_market(std::ostream& os = std::cout) const = 0;
  virtual ps_crs_base* get_impl() = 0;
  virtual const ps_crs_base* get_impl() const = 0;
  virtual const ps_mapping_1d_base* get_map() const = 0;
};


template<typename CRS>
class ps_crs_wrapper : public ps_crs_base {
  using crs_type = CRS;
public:
  ps_crs_wrapper() : crs_impl_() {}
  virtual ~ps_crs_wrapper() = default;
  void initialize(rokko::mapping_1d const& map, int num_entries_per_row) { crs_impl_.initialize(map, num_entries_per_row); }
  void insert(int row, std::vector<int> const& cols, std::vector<double> const& values) {
    crs_impl_.insert(row, cols, values);
  }
  void insert(int row, int col_size, int const*const cols, double const*const values) {
    crs_impl_.insert(row, col_size, cols, values);
  }
  void complete() { crs_impl_.complete(); }
  void extract(int row, std::vector<int>& cols, std::vector<double>& values) const {
    crs_impl_.extract(row, cols, values);
  }
  int get_dim() const { return crs_impl_.get_dim(); }
  int num_local_rows() const { return crs_impl_.num_local_rows(); }
  int start_row() const { return crs_impl_.start_row(); }
  int end_row() const { return crs_impl_.end_row(); }
  int get_nnz() const { return crs_impl_.get_nnz(); }
  void print() const { crs_impl_.print(); }
  void output_matrix_market(std::ostream& os = std::cout) const { crs_impl_.output_matrix_market(os); }
  MPI_Comm get_comm() const { return crs_impl_.get_comm(); }
  const ps_crs_base* get_impl() const { return &crs_impl_; }
  ps_crs_base* get_impl() { return &crs_impl_; }
  const ps_mapping_1d_base* get_map() const { return crs_impl_.get_map(); }

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
  std::string get_solver_name() const { return solver_name_; }
  void insert(int row, std::vector<int> const& cols, std::vector<double> const& values) const {
    crs_impl_->insert(row, cols, values);
  }
  template <std::size_t N>
  void insert(int row, std::array<int,N> const& cols, std::array<double,N> const& values) const {
    assert(cols.size() == values.size());
    crs_impl_->insert(row, cols.size(), cols.data(), values.data());
  }
  void insert(int row, int col_size, int const*const cols, double const*const values) const {
    crs_impl_->insert(row, col_size, cols, values);
  }
  void complete() const {
    crs_impl_->complete();
  }
  void extract(int row, std::vector<int>& cols, std::vector<double>& values) const {
    crs_impl_->extract(row, cols, values);
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
  const detail::ps_mapping_1d_base* get_map() const {
    return crs_impl_->get_map();
  }
  void output_matrix_market(std::ostream& os = std::cout) const {
    const auto& comm = get_map()->get_mpi_comm();
    constexpr int root_proc = 0;
    std::vector<int> cols;
    std::vector<double> values;

    const int nnz = get_nnz();
    if (comm.get_myrank() == root_proc) {
      os << "%%MatrixMarket matrix coordinate real general" << std::endl
         << get_dim() << " " << get_dim() << " " << nnz << std::endl;
    }
    comm.barrier();
    for (int global_row=0; global_row<get_dim(); ++global_row) {
      if ((global_row >= start_row()) && (global_row < end_row())) {
        extract(global_row, cols, values);
        for (int i=0; i<cols.size(); ++i) {
          os << global_row + 1 << " " << cols[i] + 1 << " " << values[i] << std::endl;
        }
      }
      comm.barrier();
    }
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
