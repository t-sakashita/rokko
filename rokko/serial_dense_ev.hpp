/*****************************************************************************
*
* Rokko: Integrated Interface for libraries of eigenvalue decomposition
*
* Copyright (C) 2012-2015 Rokko Developers https://github.com/t-sakashita/rokko
*
* Distributed under the Boost Software License, Version 1.0. (See accompanying
* file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
*
*****************************************************************************/

#ifndef ROKKO_SERIAL_DENSE_SOLVER_HPP
#define ROKKO_SERIAL_DENSE_SOLVER_HPP

#include <string>
#include <rokko/factory.hpp>
#include <rokko/eigen3.hpp>
#include <rokko/utility/timer.hpp>
#include <rokko/parameters.hpp>

namespace rokko {

namespace detail {

class sd_ev_base {
public:
  virtual ~sd_ev_base() {}
  virtual void initialize(int& argc, char**& argv) = 0;
  virtual void finalize() = 0;
  // -------------- standard eigenvalue probelm ---------------
  // with parameters, eigenvalues/eigenvectors
  virtual parameters diagonalize(Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>& mat,
				 Eigen::VectorXd& eigvals,
				 Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>& eigvecs,
				 parameters const& params) = 0;
  virtual parameters diagonalize(Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::ColMajor>& mat,
				 Eigen::VectorXd& eigvals,
				 Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::ColMajor>& eigvecs,
				 parameters const& params) = 0;
  virtual parameters diagonalize(Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>& mat,
				 std::vector<double>& eigvals,
				 Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>& eigvecs,
				 parameters const& params) = 0;
  virtual parameters diagonalize(Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::ColMajor>& mat,
				 std::vector<double>& eigvals,
				 Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::ColMajor>& eigvecs,
				 parameters const& params) = 0;
  // with parameters, only eigenvalues
  virtual parameters diagonalize(Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>& mat,
				 Eigen::VectorXd& eigvals,
				 parameters const& params) = 0;
  virtual parameters diagonalize(Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::ColMajor>& mat,
				 Eigen::VectorXd& eigvals,
				 parameters const& params) = 0;
  virtual parameters diagonalize(Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>& mat,
				 std::vector<double>& eigvals,
				 parameters const& params) = 0;
  virtual parameters diagonalize(Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::ColMajor>& mat,
				 std::vector<double>& eigvals,
				 parameters const& params) = 0;
};
  
template<typename SOLVER>
class sd_ev_wrapper : public sd_ev_base {
  using solver_type = SOLVER;
public:
  sd_ev_wrapper() : solver_impl_() {}
  virtual ~sd_ev_wrapper() {}
  void initialize(int& argc, char**& argv) {
    solver_impl_.initialize(argc, argv);
  }
  void finalize() { solver_impl_.finalize(); }
  // -------------- standard eigenvalue probelm ---------------
  // with parameters, eigenvalues/eigenvectors
  parameters diagonalize(Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>& mat,
			 Eigen::VectorXd& eigvals, Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>& eigvecs,
			 parameters const& params) {
    return solver_impl_.diagonalize(mat, eigvals, eigvecs, params);
  }
  parameters diagonalize(Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::ColMajor>& mat,
			 Eigen::VectorXd& eigvals, Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::ColMajor>& eigvecs,
			 parameters const& params) {
    return solver_impl_.diagonalize(mat, eigvals, eigvecs, params);
  }
  parameters diagonalize(Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>& mat,
			 std::vector<double>& eigvals, Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>& eigvecs,
			 parameters const& params) {
    return solver_impl_.diagonalize(mat, eigvals, eigvecs, params);
  }
  parameters diagonalize(Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::ColMajor>& mat,
			 std::vector<double>& eigvals, Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::ColMajor>& eigvecs,
			 parameters const& params) {
    return solver_impl_.diagonalize(mat, eigvals, eigvecs, params);
  }
  // with parameters, only eigenvalues
  parameters diagonalize(Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>& mat,
			 Eigen::VectorXd& eigvals,
			 parameters const& params) {
    return solver_impl_.diagonalize(mat, eigvals, params);
  }
  parameters diagonalize(Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::ColMajor>& mat,
			 Eigen::VectorXd& eigvals,
			 parameters const& params) {
    return solver_impl_.diagonalize(mat, eigvals, params);
  }
  parameters diagonalize(Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>& mat,
			 std::vector<double>& eigvals,
			 parameters const& params) {
    return solver_impl_.diagonalize(mat, eigvals, params);
  }
  parameters diagonalize(Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::ColMajor>& mat,
			 std::vector<double>& eigvals,
			 parameters const& params) {
    return solver_impl_.diagonalize(mat, eigvals, params);
  }

private:
  solver_type solver_impl_;
};
    
using sd_solver_factory = factory<sd_ev_base>;
  
} // end namespace detail
  
class serial_dense_ev {
public:
  serial_dense_ev(std::string const& solver_name)
    : solver_impl_(detail::sd_solver_factory::instance()->make_product(solver_name)) {}

  serial_dense_ev() : serial_dense_ev(default_solver()) {}

  void initialize(int& argc, char**& argv) {
    solver_impl_->initialize(argc, argv);
  }
  void finalize() {
    solver_impl_->finalize();
  }
  // -------------- standard eigenvalue probelm ---------------
  // with parameters, eigenvalues/eigenvectors
  template<int MATRIX_MAJOR, typename VEC>
  parameters diagonalize(Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,MATRIX_MAJOR>& mat, VEC& eigvals,
			 Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,MATRIX_MAJOR>& eigvecs,
			 parameters const& params) {
    return solver_impl_->diagonalize(mat, eigvals, eigvecs, params);
  }
  // with parameters, only eigenvalues
  template<int MATRIX_MAJOR, typename VEC>
  parameters diagonalize(Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,MATRIX_MAJOR>& mat, VEC& eigvals,
			 parameters const& params) {
    return solver_impl_->diagonalize(mat, eigvals, params);
  }
  // no parameters, eigenvalues/eigenvectors
  template<int MATRIX_MAJOR, typename VEC>
  parameters diagonalize(Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,MATRIX_MAJOR>& mat, VEC& eigvals,
			 Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,MATRIX_MAJOR>& eigvecs) {
    return solver_impl_->diagonalize(mat, eigvals, eigvecs, null_params);
  }
  // no parameters, only eigenvalues
  template<int MATRIX_MAJOR, typename VEC>
  parameters diagonalize(Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,MATRIX_MAJOR>& mat, VEC& eigvals) {
    return solver_impl_->diagonalize(mat, eigvals, null_params);
  }  
  static std::vector<std::string> solvers() {
    return detail::sd_solver_factory::product_names();
  }
  static std::string default_solver() {
    return detail::sd_solver_factory::default_product_name();
  }
private:
  detail::sd_solver_factory::product_pointer_type solver_impl_;
  parameters null_params;
  //std::string routine_;
};

} // end namespace rokko

#define ROKKO_REGISTER_SERIAL_DENSE_SOLVER(solver, name, priority) \
namespace { namespace BOOST_JOIN(register, __LINE__) { \
struct register_caller { \
  using factory = rokko::factory<rokko::detail::sd_ev_base>; \
  using product = rokko::detail::sd_ev_wrapper<solver>; \
  register_caller() { factory::instance()->register_creator<product>(name, priority); } \
} caller; \
} }

#endif // ROKKO_SERIAL_DENSE_SOLVER_HPP
