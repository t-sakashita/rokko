#ifndef ROKKO_SOLVER_FACTORY_H
#define ROKKO_SOLVER_FACTORY_H

#include <rokko/distributed_matrix.hpp>
#include <Eigen/Dense>
#include <boost/noncopyable.hpp>
#include <boost/shared_ptr.hpp>
#include <map>
#include <string>

namespace rokko {

class solver_factory : private boost::noncopyable {
private:
  class solver_base {
  public:
    virtual ~solver_base() {}
    virtual void initialize(int& argc, char**& argv) = 0;
    virtual void finalize() = 0;
    virtual void diagonalize(rokko::distributed_matrix<matrix_row_major>& mat,
      Eigen::VectorXd& eigvals, rokko::distributed_matrix<matrix_row_major>& eigvecs) = 0;
    virtual void diagonalize(rokko::distributed_matrix<matrix_col_major>& mat,
      Eigen::VectorXd& eigvals, rokko::distributed_matrix<matrix_col_major>& eigvecs) = 0;
  };

  template<typename SOLVER>
  class solver_wrapper : public solver_base {
    typedef SOLVER solver_type;
  public:
    solver_wrapper() : solver_impl_() {}
    virtual ~solver_wrapper() {}
    void initialize(int& argc, char**& argv) { solver_impl_.initialize(argc, argv); }
    void finalize() { solver_impl_.finalize(); }
    void diagonalize(rokko::distributed_matrix<matrix_row_major>& mat,
      Eigen::VectorXd& eigvals, rokko::distributed_matrix<matrix_row_major>& eigvecs) {
      solver_impl_.diagonalize(mat, eigvals, eigvecs);
    }
    void diagonalize(rokko::distributed_matrix<matrix_col_major>& mat,
      Eigen::VectorXd& eigvals, rokko::distributed_matrix<matrix_col_major>& eigvecs) {
      solver_impl_.diagonalize(mat, eigvals, eigvecs);
    }
  private:
    solver_type solver_impl_;
  };

  class abstract_solver_creator {
  public:
    virtual ~abstract_solver_creator() {}
    virtual boost::shared_ptr<solver_base> create() const = 0;
  };

  template <typename SOLVER>
  class solver_creator : public abstract_solver_creator {
  public:
    virtual ~solver_creator() {}
    typedef SOLVER solver_type;
    boost::shared_ptr<solver_base> create() const {
      return boost::shared_ptr<solver_base>(new solver_type());
    }
  };

public:  
  typedef boost::shared_ptr<solver_base> solver_pointer_type;

private:
  typedef boost::shared_ptr<abstract_solver_creator> creator_pointer_type;
  typedef std::map<std::string, creator_pointer_type> creator_map_type;
public:
  static solver_pointer_type make_solver(std::string const& name);
  template<typename SOLVER>
  bool register_creator(std::string const& name) {
    bool isnew = (creators_.find(name) == creators_.end());
    creators_[name] = creator_pointer_type(new solver_creator<solver_wrapper<SOLVER> >());
    return isnew;
  }
  bool unregister_creator(std::string const& name);
  static solver_factory* instance();
protected:
  creator_pointer_type make_creator(std::string const& name) const;
private:
  static solver_factory* instance_;
  static int initialize_;
  creator_map_type creators_;
};

} // end namespace rokko

#define ROKKO_REGISTER_SOLVER(solver, name) \
namespace {                               \
  const bool BOOST_JOIN(solver_, __LINE__) = \
    rokko::solver_factory::instance()->register_creator<solver>(name); \
}

#endif // ROKKO_SOLVER_FACTORY_H
