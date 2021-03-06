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

#ifndef ROKKO_PARALLEL_DENSE_SOLVER_HPP
#define ROKKO_PARALLEL_DENSE_SOLVER_HPP

#include <rokko/factory.hpp>
#include <rokko/matrix_major.hpp>
#include <rokko/mapping_bc.hpp>
#include <rokko/distributed_matrix.hpp>
#include <rokko/eigen3.hpp>
#include <rokko/utility/macro_join.hpp>
#include <rokko/utility/timer.hpp>
#include <rokko/parameters.hpp>

namespace rokko {

namespace detail {
    
class pd_ev_base {
public:
  virtual ~pd_ev_base() = default;
  virtual bool is_available_grid_major(grid_row_major_t const& grid_major) = 0;
  virtual bool is_available_grid_major(grid_col_major_t const& grid_major) = 0;
  virtual void initialize(int& argc, char**& argv) = 0;
  virtual void finalize() = 0;

  // float # eigenvalues/eigenvectors
  virtual parameters diagonalize(distributed_matrix<float, matrix_row_major>& mat,
				 Eigen::VectorXf& eigvals, distributed_matrix<float, matrix_row_major>& eigvecs,
				 parameters const& params) = 0;
  virtual parameters diagonalize(distributed_matrix<float, matrix_col_major>& mat,
				 Eigen::VectorXf& eigvals, distributed_matrix<float, matrix_col_major>& eigvecs,
				 parameters const& params) = 0;
  virtual parameters diagonalize(distributed_matrix<float, matrix_row_major>& mat,
				 Eigen::RefVec<float>& eigvals, distributed_matrix<float, matrix_row_major>& eigvecs,
				 parameters const& params) = 0;
  virtual parameters diagonalize(distributed_matrix<float, matrix_col_major>& mat,
				 Eigen::RefVec<float>& eigvals, distributed_matrix<float, matrix_col_major>& eigvecs,
				 parameters const& params) = 0;
  // float # only eigenvalues
  virtual parameters diagonalize(distributed_matrix<float, matrix_row_major>& mat,
				 Eigen::VectorXf& eigvals,
				 parameters const& params) = 0;
  virtual parameters diagonalize(distributed_matrix<float, matrix_col_major>& mat,
				 Eigen::VectorXf& eigvals,
				 parameters const& params) = 0;
  virtual parameters diagonalize(distributed_matrix<float, matrix_row_major>& mat,
				 Eigen::RefVec<float>& eigvals,
				 parameters const& params) = 0;
  virtual parameters diagonalize(distributed_matrix<float, matrix_col_major>& mat,
				 Eigen::RefVec<float>& eigvals,
				 parameters const& params) = 0;
  // complex float # eigenvalues/eigenvectors
  virtual parameters diagonalize(distributed_matrix<std::complex<float>, matrix_row_major>& mat,
				 Eigen::VectorXf& eigvals, distributed_matrix<std::complex<float>, matrix_row_major>& eigvecs,
				 parameters const& params) = 0;
  virtual parameters diagonalize(distributed_matrix<std::complex<float>, matrix_col_major>& mat,
				 Eigen::VectorXf& eigvals, distributed_matrix<std::complex<float>, matrix_col_major>& eigvecs,
				 parameters const& params) = 0;
  virtual parameters diagonalize(distributed_matrix<std::complex<float>, matrix_row_major>& mat,
				 Eigen::RefVec<float>& eigvals, distributed_matrix<std::complex<float>, matrix_row_major>& eigvecs,
				 parameters const& params) = 0;
  virtual parameters diagonalize(distributed_matrix<std::complex<float>, matrix_col_major>& mat,
				 Eigen::RefVec<float>& eigvals, distributed_matrix<std::complex<float>, matrix_col_major>& eigvecs,
				 parameters const& params) = 0;
  // complex float # only eigenvalues
  virtual parameters diagonalize(distributed_matrix<std::complex<float>, matrix_row_major>& mat,
				 Eigen::VectorXf& eigvals,
				 parameters const& params) = 0;
  virtual parameters diagonalize(distributed_matrix<std::complex<float>, matrix_col_major>& mat,
				 Eigen::VectorXf& eigvals,
				 parameters const& params) = 0;
  virtual parameters diagonalize(distributed_matrix<std::complex<float>, matrix_row_major>& mat,
				 Eigen::RefVec<float>& eigvals,
				 parameters const& params) = 0;
  virtual parameters diagonalize(distributed_matrix<std::complex<float>, matrix_col_major>& mat,
				 Eigen::RefVec<float>& eigvals,
                 parameters const& params) = 0;

  // double # eigenvalues/eigenvectors
  virtual parameters diagonalize(distributed_matrix<double, matrix_row_major>& mat,
				 Eigen::VectorXd& eigvals, distributed_matrix<double, matrix_row_major>& eigvecs,
				 parameters const& params) = 0;
  virtual parameters diagonalize(distributed_matrix<double, matrix_col_major>& mat,
				 Eigen::VectorXd& eigvals, distributed_matrix<double, matrix_col_major>& eigvecs,
				 parameters const& params) = 0;
  virtual parameters diagonalize(distributed_matrix<double, matrix_row_major>& mat,
				 Eigen::RefVec<double>& eigvals, distributed_matrix<double, matrix_row_major>& eigvecs,
				 parameters const& params) = 0;
  virtual parameters diagonalize(distributed_matrix<double, matrix_col_major>& mat,
				 Eigen::RefVec<double>& eigvals, distributed_matrix<double, matrix_col_major>& eigvecs,
				 parameters const& params) = 0;
  // double # only eigenvalues
  virtual parameters diagonalize(distributed_matrix<double, matrix_row_major>& mat,
				 Eigen::VectorXd& eigvals,
				 parameters const& params) = 0;
  virtual parameters diagonalize(distributed_matrix<double, matrix_col_major>& mat,
				 Eigen::VectorXd& eigvals,
				 parameters const& params) = 0;
  virtual parameters diagonalize(distributed_matrix<double, matrix_row_major>& mat,
				 Eigen::RefVec<double>& eigvals,
				 parameters const& params) = 0;
  virtual parameters diagonalize(distributed_matrix<double, matrix_col_major>& mat,
				 Eigen::RefVec<double>& eigvals,
				 parameters const& params) = 0;
  // complex double # eigenvalues/eigenvectors
  virtual parameters diagonalize(distributed_matrix<std::complex<double>, matrix_row_major>& mat,
				 Eigen::VectorXd& eigvals, distributed_matrix<std::complex<double>, matrix_row_major>& eigvecs,
				 parameters const& params) = 0;
  virtual parameters diagonalize(distributed_matrix<std::complex<double>, matrix_col_major>& mat,
				 Eigen::VectorXd& eigvals, distributed_matrix<std::complex<double>, matrix_col_major>& eigvecs,
				 parameters const& params) = 0;
  virtual parameters diagonalize(distributed_matrix<std::complex<double>, matrix_row_major>& mat,
				 Eigen::RefVec<double>& eigvals, distributed_matrix<std::complex<double>, matrix_row_major>& eigvecs,
				 parameters const& params) = 0;
  virtual parameters diagonalize(distributed_matrix<std::complex<double>, matrix_col_major>& mat,
				 Eigen::RefVec<double>& eigvals, distributed_matrix<std::complex<double>, matrix_col_major>& eigvecs,
				 parameters const& params) = 0;

  // complex double # only eigenvalues
  virtual parameters diagonalize(distributed_matrix<std::complex<double>, matrix_row_major>& mat,
				 Eigen::VectorXd& eigvals,
				 parameters const& params) = 0;
  virtual parameters diagonalize(distributed_matrix<std::complex<double>, matrix_col_major>& mat,
				 Eigen::VectorXd& eigvals,
				 parameters const& params) = 0;
  virtual parameters diagonalize(distributed_matrix<std::complex<double>, matrix_row_major>& mat,
				 Eigen::RefVec<double>& eigvals,
				 parameters const& params) = 0;
  virtual parameters diagonalize(distributed_matrix<std::complex<double>, matrix_col_major>& mat,
				 Eigen::RefVec<double>& eigvals,
                 parameters const& params) = 0;

  virtual mapping_bc<matrix_col_major> default_mapping(int dim, grid const& g) const = 0;
};
  
template<typename SOLVER>
class pd_ev_wrapper : public pd_ev_base {
  using solver_type = SOLVER;
public:
  pd_ev_wrapper() : solver_impl_() {}
  virtual ~pd_ev_wrapper() = default;
  bool is_available_grid_major(grid_row_major_t const& grid_major) override {
    return solver_impl_.is_available_grid_major(grid_major);
  }
  bool is_available_grid_major(grid_col_major_t const& grid_major) override {
    return solver_impl_.is_available_grid_major(grid_major);
  }
  void initialize(int& argc, char**& argv) override {
    solver_impl_.initialize(argc, argv);
  }
  void finalize() override { solver_impl_.finalize(); }
  // float # eigenvalues/eigenvectors
  parameters diagonalize(distributed_matrix<float, matrix_row_major>& mat,
			 Eigen::VectorXf& eigvals, distributed_matrix<float, matrix_row_major>& eigvecs,
			 parameters const& params) override {
    return solver_impl_.diagonalize(mat, eigvals, eigvecs, params);
  }
  parameters diagonalize(distributed_matrix<float, matrix_col_major>& mat,
			 Eigen::VectorXf& eigvals, distributed_matrix<float, matrix_col_major>& eigvecs,
			 parameters const& params) override {
    return solver_impl_.diagonalize(mat, eigvals, eigvecs, params);
  }
  parameters diagonalize(distributed_matrix<float, matrix_row_major>& mat,
			 Eigen::RefVec<float>& eigvals, distributed_matrix<float, matrix_row_major>& eigvecs,
			 parameters const& params) override {
    return solver_impl_.diagonalize(mat, eigvals, eigvecs, params);
  }
  parameters diagonalize(distributed_matrix<float, matrix_col_major>& mat,
			 Eigen::RefVec<float>& eigvals, distributed_matrix<float, matrix_col_major>& eigvecs,
			 parameters const& params) override {
    return solver_impl_.diagonalize(mat, eigvals, eigvecs, params);
  }
  // float # only eigenvalues
  parameters diagonalize(distributed_matrix<float, matrix_row_major>& mat,
			 Eigen::VectorXf& eigvals,
			 parameters const& params) override {
    return solver_impl_.diagonalize(mat, eigvals, params);
  }
  parameters diagonalize(distributed_matrix<float, matrix_col_major>& mat,
			 Eigen::VectorXf& eigvals,
			 parameters const& params) override {
    return solver_impl_.diagonalize(mat, eigvals, params);
  }
  parameters diagonalize(distributed_matrix<float, matrix_row_major>& mat,
			 Eigen::RefVec<float>& eigvals,
			 parameters const& params) override {
    return solver_impl_.diagonalize(mat, eigvals, params);
  }
  parameters diagonalize(distributed_matrix<float, matrix_col_major>& mat,
			 Eigen::RefVec<float>& eigvals,
			 parameters const& params) override {
    return solver_impl_.diagonalize(mat, eigvals, params);
  }
  // complex float # eigenvalues/eigenvectors
  parameters diagonalize(distributed_matrix<std::complex<float>, matrix_row_major>& mat,
			 Eigen::VectorXf& eigvals, distributed_matrix<std::complex<float>, matrix_row_major>& eigvecs,
			 parameters const& params) override {
    return solver_impl_.diagonalize(mat, eigvals, eigvecs, params);
  }
  parameters diagonalize(distributed_matrix<std::complex<float>, matrix_col_major>& mat,
			 Eigen::VectorXf& eigvals, distributed_matrix<std::complex<float>, matrix_col_major>& eigvecs,
			 parameters const& params) override {
    return solver_impl_.diagonalize(mat, eigvals, eigvecs, params);
  }
  parameters diagonalize(distributed_matrix<std::complex<float>, matrix_row_major>& mat,
			 Eigen::RefVec<float>& eigvals, distributed_matrix<std::complex<float>, matrix_row_major>& eigvecs,
			 parameters const& params) override {
    return solver_impl_.diagonalize(mat, eigvals, eigvecs, params);
  }
  parameters diagonalize(distributed_matrix<std::complex<float>, matrix_col_major>& mat,
			 Eigen::RefVec<float>& eigvals, distributed_matrix<std::complex<float>, matrix_col_major>& eigvecs,
			 parameters const& params) override {
    return solver_impl_.diagonalize(mat, eigvals, eigvecs, params);
  }
  // complex float # only eigenvalues
  parameters diagonalize(distributed_matrix<std::complex<float>, matrix_row_major>& mat,
			 Eigen::VectorXf& eigvals,
			 parameters const& params) override {
    return solver_impl_.diagonalize(mat, eigvals, params);
  }
  parameters diagonalize(distributed_matrix<std::complex<float>, matrix_col_major>& mat,
			 Eigen::VectorXf& eigvals,
			 parameters const& params) override {
    return solver_impl_.diagonalize(mat, eigvals, params);
  }
  parameters diagonalize(distributed_matrix<std::complex<float>, matrix_row_major>& mat,
			 Eigen::RefVec<float>& eigvals,
			 parameters const& params) override {
    return solver_impl_.diagonalize(mat, eigvals, params);
  }
  parameters diagonalize(distributed_matrix<std::complex<float>, matrix_col_major>& mat,
			 Eigen::RefVec<float>& eigvals,
			 parameters const& params) override {
    return solver_impl_.diagonalize(mat, eigvals, params);
  }
  // double # eigenvalues/eigenvectors
  parameters diagonalize(distributed_matrix<double, matrix_row_major>& mat,
			 Eigen::VectorXd& eigvals, distributed_matrix<double, matrix_row_major>& eigvecs,
			 parameters const& params) override {
    return solver_impl_.diagonalize(mat, eigvals, eigvecs, params);
  }
  parameters diagonalize(distributed_matrix<double, matrix_col_major>& mat,
			 Eigen::VectorXd& eigvals, distributed_matrix<double, matrix_col_major>& eigvecs,
			 parameters const& params) override {
    return solver_impl_.diagonalize(mat, eigvals, eigvecs, params);
  }
  parameters diagonalize(distributed_matrix<double, matrix_row_major>& mat,
			 Eigen::RefVec<double>& eigvals, distributed_matrix<double, matrix_row_major>& eigvecs,
			 parameters const& params) override {
    return solver_impl_.diagonalize(mat, eigvals, eigvecs, params);
  }
  parameters diagonalize(distributed_matrix<double, matrix_col_major>& mat,
			 Eigen::RefVec<double>& eigvals, distributed_matrix<double, matrix_col_major>& eigvecs,
			 parameters const& params) override {
    return solver_impl_.diagonalize(mat, eigvals, eigvecs, params);
  }
  // double # only eigenvalues
  parameters diagonalize(distributed_matrix<double, matrix_row_major>& mat,
			 Eigen::VectorXd& eigvals,
			 parameters const& params) override {
    return solver_impl_.diagonalize(mat, eigvals, params);
  }
  parameters diagonalize(distributed_matrix<double, matrix_col_major>& mat,
			 Eigen::VectorXd& eigvals,
			 parameters const& params) override {
    return solver_impl_.diagonalize(mat, eigvals, params);
  }
  parameters diagonalize(distributed_matrix<double, matrix_row_major>& mat,
			 Eigen::RefVec<double>& eigvals,
			 parameters const& params) override {
    return solver_impl_.diagonalize(mat, eigvals, params);
  }
  parameters diagonalize(distributed_matrix<double, matrix_col_major>& mat,
			 Eigen::RefVec<double>& eigvals,
			 parameters const& params) override {
    return solver_impl_.diagonalize(mat, eigvals, params);
  }
  // complex double # eigenvalues/eigenvectors
  parameters diagonalize(distributed_matrix<std::complex<double>, matrix_row_major>& mat,
			 Eigen::VectorXd& eigvals, distributed_matrix<std::complex<double>, matrix_row_major>& eigvecs,
			 parameters const& params) override {
    return solver_impl_.diagonalize(mat, eigvals, eigvecs, params);
  }
  parameters diagonalize(distributed_matrix<std::complex<double>, matrix_col_major>& mat,
			 Eigen::VectorXd& eigvals, distributed_matrix<std::complex<double>, matrix_col_major>& eigvecs,
			 parameters const& params) override {
    return solver_impl_.diagonalize(mat, eigvals, eigvecs, params);
  }
  parameters diagonalize(distributed_matrix<std::complex<double>, matrix_row_major>& mat,
			 Eigen::RefVec<double>& eigvals, distributed_matrix<std::complex<double>, matrix_row_major>& eigvecs,
			 parameters const& params) override {
    return solver_impl_.diagonalize(mat, eigvals, eigvecs, params);
  }
  parameters diagonalize(distributed_matrix<std::complex<double>, matrix_col_major>& mat,
			 Eigen::RefVec<double>& eigvals, distributed_matrix<std::complex<double>, matrix_col_major>& eigvecs,
			 parameters const& params) override {
    return solver_impl_.diagonalize(mat, eigvals, eigvecs, params);
  }
  // complex double # only eigenvalues
  parameters diagonalize(distributed_matrix<std::complex<double>, matrix_row_major>& mat,
			 Eigen::VectorXd& eigvals,
			 parameters const& params) override {
    return solver_impl_.diagonalize(mat, eigvals, params);
  }
  parameters diagonalize(distributed_matrix<std::complex<double>, matrix_col_major>& mat,
			 Eigen::VectorXd& eigvals,
			 parameters const& params) override {
    return solver_impl_.diagonalize(mat, eigvals, params);
  }
  parameters diagonalize(distributed_matrix<std::complex<double>, matrix_row_major>& mat,
			 Eigen::RefVec<double>& eigvals,
			 parameters const& params) override {
    return solver_impl_.diagonalize(mat, eigvals, params);
  }
  parameters diagonalize(distributed_matrix<std::complex<double>, matrix_col_major>& mat,
			 Eigen::RefVec<double>& eigvals,
			 parameters const& params) override {
    return solver_impl_.diagonalize(mat, eigvals, params);
  }
  // default_mapping
  mapping_bc<matrix_col_major> default_mapping(int dim, grid const& g) const override {
    return solver_impl_.default_mapping(dim, g);
  }
private:
  solver_type solver_impl_;
};

using pd_solver_factory = factory<pd_ev_base>;
  
} // end namespace detail
  
class parallel_dense_ev {
public:
  parallel_dense_ev(std::string const& solver_name) : solver_impl_(detail::pd_solver_factory::instance()->make_product(solver_name)), null_params() {};

  parallel_dense_ev() : parallel_dense_ev(default_solver()) {}

  template <typename GRID_MAJOR>
  bool is_available_grid_major(GRID_MAJOR const& grid_major) {
    return solver_impl_->is_available_grid_major(grid_major);
  }
  void initialize(int& argc, char**& argv) {
    solver_impl_->initialize(argc, argv);
  }
  void finalize() {
    solver_impl_->finalize();
  }
  // with parameters, eigenvalues/eigenvectors
  template<typename T, typename MATRIX_MAJOR, typename VEC>
  parameters diagonalize(distributed_matrix<T, MATRIX_MAJOR>& mat,
			 VEC& eigvals, rokko::distributed_matrix<T, MATRIX_MAJOR>& eigvecs,
			 parameters const& params) {
    return solver_impl_->diagonalize(mat, eigvals, eigvecs, params);
  }
  // with parameters, only eigenvalues
  template<typename T, typename MATRIX_MAJOR, typename VEC>
  parameters diagonalize(distributed_matrix<T, MATRIX_MAJOR>& mat, VEC& eigvals,
			 parameters const& params) {
    return solver_impl_->diagonalize(mat, eigvals, params);
  }
  // no parameters, eigenvalues/eigenvectors
  template<typename T, typename MATRIX_MAJOR, typename VEC>
  parameters diagonalize(distributed_matrix<T, MATRIX_MAJOR>& mat,
			 VEC& eigvals, rokko::distributed_matrix<T, MATRIX_MAJOR>& eigvecs) {
    return solver_impl_->diagonalize(mat, eigvals, eigvecs, null_params);
  }
  // no parameters, only eigenvalues
  template<typename T, typename MATRIX_MAJOR, typename VEC>
  parameters diagonalize(distributed_matrix<T, MATRIX_MAJOR>& mat,
			 VEC& eigvals) {
    return diagonalize(mat, eigvals, null_params);
  }
  mapping_bc<matrix_col_major> default_mapping(int dim, grid const& g) const {
    return solver_impl_->default_mapping(dim, g);
  }
  static std::vector<std::string> solvers() {
    return detail::pd_solver_factory::product_names();
  }
  static const std::string& default_solver() {
    return detail::pd_solver_factory::default_product_name();
  }
private:
  detail::pd_solver_factory::product_pointer_type solver_impl_;
  parameters null_params;
  std::string routine_;
};

} // end namespace rokko

#define ROKKO_REGISTER_PARALLEL_DENSE_SOLVER(solver, name, priority) \
namespace { namespace ROKKO_JOIN(register, __LINE__) { \
struct register_caller { \
  using factory = rokko::detail::pd_solver_factory;     \
  using product = rokko::detail::pd_ev_wrapper<solver>; \
  register_caller() { factory::instance()->register_creator<product>(name, priority); } \
} caller; \
} }

#endif // ROKKO_PARALLEL_DENSE_SOLVER_HPP
