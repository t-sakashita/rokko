/*****************************************************************************
*
* Rokko: Integrated Interface for libraries of eigenvalue decomposition
*
* Copyright (C) 2012-2013 by Tatsuya Sakashita <t-sakashita@issp.u-tokyo.ac.jp>,
*                            Synge Todo <wistaria@comp-phys.org>,
*               2013-2013    Ryo IGARASHI <rigarash@issp.u-tokyo.ac.jp>
*
* Distributed under the Boost Software License, Version 1.0. (See accompanying
* file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
*
*****************************************************************************/

#ifndef ROKKO_SERIAL_SOLVER_FACTORY_H
#define ROKKO_SERIAL_SOLVER_FACTORY_H

#include <rokko/localized_matrix.hpp>
#include <rokko/localized_vector.hpp>
#include <rokko/utility/timer.hpp>
#include <boost/noncopyable.hpp>
#include <boost/shared_ptr.hpp>
#include <map>
#include <string>
#include <vector>

namespace rokko {

class serial_solver_factory : private boost::noncopyable {
private:
  class serial_solver_base {
  public:
    virtual ~serial_solver_base() {}
    virtual void initialize(int& argc, char**& argv) = 0;
    virtual void finalize() = 0;
    virtual void diagonalize(localized_matrix<matrix_row_major>& mat,
                             localized_vector& eigvals, localized_matrix<matrix_row_major>& eigvecs, timer& timer_in) = 0;
    virtual void diagonalize(localized_matrix<matrix_col_major>& mat,
                             localized_vector& eigvals, localized_matrix<matrix_col_major>& eigvecs, timer& timer_in) = 0;
  };

  template<typename SOLVER>
  class serial_solver_wrapper : public serial_solver_base {
    typedef SOLVER solver_type;
  public:
    serial_solver_wrapper() : solver_impl_() {}
    virtual ~serial_solver_wrapper() {}
    void initialize(int& argc, char**& argv) { solver_impl_.initialize(argc, argv); }
    void finalize() { solver_impl_.finalize(); }
    void diagonalize(localized_matrix<matrix_row_major>& mat,
                     localized_vector& eigvals, localized_matrix<matrix_row_major>& eigvecs, timer& timer_in) {
      solver_impl_.diagonalize(mat, eigvals, eigvecs, timer_in);
    }
    void diagonalize(localized_matrix<matrix_col_major>& mat,
                     localized_vector& eigvals, localized_matrix<matrix_col_major>& eigvecs, timer& timer_in) {
      solver_impl_.diagonalize(mat, eigvals, eigvecs, timer_in);
    }
  private:
    solver_type solver_impl_;
  };

  class abstract_serial_solver_creator {
  public:
    virtual ~abstract_serial_solver_creator() {}
    virtual boost::shared_ptr<serial_solver_base> create() const = 0;
  };

  template <typename SOLVER>
  class serial_solver_creator : public abstract_serial_solver_creator {
  public:
    virtual ~serial_solver_creator() {}
    typedef SOLVER solver_type;
    boost::shared_ptr<serial_solver_base> create() const {
      return boost::shared_ptr<serial_solver_base>(new solver_type());
    }
  };

public:
  typedef boost::shared_ptr<serial_solver_base> solver_pointer_type;

private:
  typedef boost::shared_ptr<abstract_serial_solver_creator> creator_pointer_type;
  typedef std::map<std::string, creator_pointer_type> creator_map_type;
public:
  static solver_pointer_type make_solver(std::string const& name);
  template<typename SOLVER>
  bool register_creator(std::string const& name) {
    bool isnew = (creators_.find(name) == creators_.end());
    creators_[name] = creator_pointer_type(new serial_solver_creator<serial_solver_wrapper<SOLVER> >());
    return isnew;
  }
  std::vector<std::string> solver_names() const;
  bool unregister_creator(std::string const& name);
  static serial_solver_factory* instance();
protected:
  creator_pointer_type make_creator(std::string const& name) const;
private:
  static serial_solver_factory* instance_;
  static int initialize_;
  creator_map_type creators_;
};

} // end namespace rokko

#define ROKKO_REGISTER_SERIAL_SOLVER(solver, name) \
namespace { struct register_caller { register_caller() { rokko::serial_solver_factory::instance()->register_creator<solver>(name); } } caller; }

#endif // ROKKO_SERIAL_SOLVER_FACTORY_H
