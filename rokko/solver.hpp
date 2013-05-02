#ifndef ROKKO_SOLVER_HPP
#define ROKKO_SOLVER_HPP

#include <rokko/solver_factory.hpp>
#include <rokko/distributed_matrix.hpp>
#include <boost/shared_ptr.hpp>

namespace rokko {

class solver {
public:
  solver(std::string const& solver_name) {
    solver_impl_ = solver_factory::instance()->make_solver(solver_name);
  }
  void initialize(int& argc, char**& argv) { solver_impl_->initialize(argc, argv); }
  void finalize() { solver_impl_->finalize(); }
  template<typename MATRIX_MAJOR>
  void diagonalize(rokko::distributed_matrix<MATRIX_MAJOR>& mat, Eigen::VectorXd& eigvals,
    rokko::distributed_matrix<MATRIX_MAJOR>& eigvecs) {
    solver_impl_->diagonalize(mat, eigvals, eigvecs);
  }
private:
  solver_factory::solver_pointer_type solver_impl_;
};

} // end namespace rokko

#endif // ROKKO_SOLVER_HPP
