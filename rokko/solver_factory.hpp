/*****************************************************************************
*
* Rokko: Integrated Interface for libraries of eigenvalue decomposition
*
* Copyright (C) 2012-2013 by Tatsuya Sakashita <t-sakashita@issp.u-tokyo.ac.jp>,
*                            Synge Todo <wistaria@comp-phys.org>
*
* Distributed under the Boost Software License, Version 1.0. (See accompanying
* file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
*
*****************************************************************************/

#ifndef ROKKO_SOLVER_FACTORY_H
#define ROKKO_SOLVER_FACTORY_H

#include <rokko/distributed_matrix.hpp>
#include <rokko/localized_vector.hpp>
#include <rokko/utility/timer.hpp>
#include <boost/noncopyable.hpp>
#include <boost/shared_ptr.hpp>
#include <map>
#include <string>
#include <vector>

namespace rokko {

class solver_factory : private boost::noncopyable {
private:
  class solver_base {
  public:
    virtual ~solver_base() {}
    virtual void initialize(int& argc, char**& argv) = 0;
    virtual void finalize() = 0;
    virtual void diagonalize(distributed_matrix<matrix_row_major>& mat,
                             localized_vector& eigvals, distributed_matrix<matrix_row_major>& eigvecs, timer& timer_in) = 0;
    virtual void diagonalize(distributed_matrix<matrix_col_major>& mat,
                             localized_vector& eigvals, distributed_matrix<matrix_col_major>& eigvecs, timer& timer_in) = 0;
    //virtual void optimized_matrix_size(int dim, int nprow, int npcol, int& mb, int& nb, int& lld, int& len_array) = 0;
    virtual void optimized_matrix_size(distributed_matrix<matrix_row_major>& mat) = 0;
    virtual void optimized_matrix_size(distributed_matrix<matrix_col_major>& mat) = 0;
  };

  template<typename SOLVER>
  class solver_wrapper : public solver_base {
    typedef SOLVER solver_type;
  public:
    solver_wrapper() : solver_impl_() {}
    virtual ~solver_wrapper() {}
    void initialize(int& argc, char**& argv) { solver_impl_.initialize(argc, argv); }
    void finalize() { solver_impl_.finalize(); }
    void diagonalize(distributed_matrix<matrix_row_major>& mat,
                     localized_vector& eigvals, distributed_matrix<matrix_row_major>& eigvecs, timer& timer_in) {
      solver_impl_.diagonalize(mat, eigvals, eigvecs, timer_in);
    }
    void diagonalize(distributed_matrix<matrix_col_major>& mat,
                     localized_vector& eigvals, distributed_matrix<matrix_col_major>& eigvecs, timer& timer_in) {
      solver_impl_.diagonalize(mat, eigvals, eigvecs, timer_in);
    }
    void optimized_matrix_size(distributed_matrix<matrix_row_major>& mat) {
      solver_impl_.optimized_matrix_size(mat);
    }
    void optimized_matrix_size(distributed_matrix<matrix_col_major>& mat) {
      solver_impl_.optimized_matrix_size(mat);
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
  std::vector<std::string> solver_names() const;
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
namespace { struct register_caller { register_caller() { rokko::solver_factory::instance()->register_creator<solver>(name); } } caller; }

#endif // ROKKO_SOLVER_FACTORY_H
