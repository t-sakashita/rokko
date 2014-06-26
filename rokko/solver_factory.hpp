/*****************************************************************************
*
* Rokko: Integrated Interface for libraries of eigenvalue decomposition
*
* Copyright (C) 2012-2014 by Tatsuya Sakashita <t-sakashita@issp.u-tokyo.ac.jp>,
*                            Synge Todo <wistaria@comp-phys.org>,
*               2013-2013    Ryo IGARASHI <rigarash@issp.u-tokyo.ac.jp>
*
* Distributed under the Boost Software License, Version 1.0. (See accompanying
* file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
*
*****************************************************************************/

#ifndef ROKKO_SOLVER_FACTORY_H
#define ROKKO_SOLVER_FACTORY_H

#include <rokko/distributed_matrix.hpp>
#include <rokko/localized_matrix.hpp>
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
  class sd_solver_base {
  public:
    virtual ~sd_solver_base() {}
    virtual void initialize(int& argc, char**& argv) = 0;
    virtual void finalize() = 0;
    virtual void diagonalize(localized_matrix<matrix_row_major>& mat, localized_vector& eigvals,
      localized_matrix<matrix_row_major>& eigvecs, timer& timer_in) = 0;
    virtual void diagonalize(localized_matrix<matrix_col_major>& mat, localized_vector& eigvals,
      localized_matrix<matrix_col_major>& eigvecs, timer& timer_in) = 0;
  };

  template<typename SOLVER>
  class sd_solver_wrapper : public sd_solver_base {
    typedef SOLVER solver_type;
  public:
    sd_solver_wrapper() : solver_impl_() {}
    virtual ~sd_solver_wrapper() {}
    void initialize(int& argc, char**& argv) { solver_impl_.initialize(argc, argv); }
    void finalize() { solver_impl_.finalize(); }
    void diagonalize(localized_matrix<matrix_row_major>& mat, localized_vector& eigvals,
      localized_matrix<matrix_row_major>& eigvecs, timer& timer_in) {
      solver_impl_.diagonalize(mat, eigvals, eigvecs, timer_in);
    }
    void diagonalize(localized_matrix<matrix_col_major>& mat, localized_vector& eigvals,
      localized_matrix<matrix_col_major>& eigvecs, timer& timer_in) {
      solver_impl_.diagonalize(mat, eigvals, eigvecs, timer_in);
    }
  private:
    solver_type solver_impl_;
  };

  class pd_solver_base {
  public:
    virtual ~pd_solver_base() {}
    virtual bool is_available_grid_major(grid_row_major_t const& grid_major) = 0;
    virtual bool is_available_grid_major(grid_col_major_t const& grid_major) = 0;
    virtual void initialize(int& argc, char**& argv) = 0;
    virtual void finalize() = 0;
    virtual void diagonalize(distributed_matrix<matrix_row_major>& mat, localized_vector& eigvals,
      distributed_matrix<matrix_row_major>& eigvecs, timer& timer_in) = 0;
    virtual void diagonalize(distributed_matrix<matrix_col_major>& mat, localized_vector& eigvals,
      distributed_matrix<matrix_col_major>& eigvecs, timer& timer_in) = 0;
    virtual void optimized_matrix_size(distributed_matrix<matrix_row_major>& mat) = 0;
    virtual void optimized_matrix_size(distributed_matrix<matrix_col_major>& mat) = 0;
  };

  template<typename SOLVER>
  class pd_solver_wrapper : public pd_solver_base {
    typedef SOLVER solver_type;
  public:
    pd_solver_wrapper() : solver_impl_() {}
    virtual ~pd_solver_wrapper() {}
    bool is_available_grid_major(grid_row_major_t const& grid_major) {
      return solver_impl_.is_available_grid_major(grid_major);
    }
    bool is_available_grid_major(grid_col_major_t const& grid_major) {
      return solver_impl_.is_available_grid_major(grid_major);
    }
    void initialize(int& argc, char**& argv) { solver_impl_.initialize(argc, argv); }
    void finalize() { solver_impl_.finalize(); }
    void diagonalize(distributed_matrix<matrix_row_major>& mat, localized_vector& eigvals,
      distributed_matrix<matrix_row_major>& eigvecs, timer& timer_in) {
      solver_impl_.diagonalize(mat, eigvals, eigvecs, timer_in);
    }
    void diagonalize(distributed_matrix<matrix_col_major>& mat, localized_vector& eigvals,
      distributed_matrix<matrix_col_major>& eigvecs, timer& timer_in) {
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

  template<typename BASE>
  class abstract_creator {
  public:
    virtual ~abstract_creator() {}
    virtual boost::shared_ptr<BASE> create() const = 0;
  };

  template <typename BASE, typename SOLVER>
  class creator : public abstract_creator<BASE> {
  public:
    virtual ~creator() {}
    boost::shared_ptr<BASE> create() const {
      return boost::shared_ptr<BASE>(new SOLVER());
    }
  };

public:
  typedef boost::shared_ptr<sd_solver_base> serial_dense_solver_pointer_type;
  typedef boost::shared_ptr<pd_solver_base> parallel_dense_solver_pointer_type;

private:
  typedef abstract_creator<sd_solver_base> abstract_sd_solver_creator;
  typedef boost::shared_ptr<abstract_sd_solver_creator> sd_creator_pointer_type;
  typedef std::map<std::string, sd_creator_pointer_type> sd_creator_map_type;
  typedef abstract_creator<pd_solver_base> abstract_pd_solver_creator;
  typedef boost::shared_ptr<abstract_pd_solver_creator> pd_creator_pointer_type;
  typedef std::map<std::string, pd_creator_pointer_type> pd_creator_map_type;
public:
  solver_factory() : sd_largest_priority_(0), pd_largest_priority_(0) {}
  static serial_dense_solver_pointer_type make_serial_dense_solver(std::string const& name);
  static serial_dense_solver_pointer_type make_serial_dense_solver();
  static parallel_dense_solver_pointer_type make_parallel_dense_solver(std::string const& name);
  static parallel_dense_solver_pointer_type make_parallel_dense_solver();
  template<typename SOLVER>
  bool register_serial_dense_creator(std::string const& name, int priority = 0) {
    bool isnew = (sd_creators_.find(name) == sd_creators_.end());
    sd_creators_[name] = sd_creator_pointer_type(
      new creator<sd_solver_base, sd_solver_wrapper<SOLVER> >());
    if (priority >= sd_largest_priority_) {
      sd_largest_priority_ = priority;
      default_sd_solver_ = name;
    }
    return isnew;
  }
  template<typename SOLVER>
  bool register_parallel_dense_creator(std::string const& name, int priority = 0) {
    bool isnew = (pd_creators_.find(name) == pd_creators_.end());
    pd_creators_[name] = pd_creator_pointer_type(
      new creator<pd_solver_base, pd_solver_wrapper<SOLVER> >());
    if (priority >= pd_largest_priority_) {
      pd_largest_priority_ = priority;
      default_pd_solver_ = name;
    }
    return isnew;
  }
  bool unregister_serial_dense_creator(std::string const& name);
  bool unregister_parallel_dense_creator(std::string const& name);
  static std::vector<std::string> serial_dense_solver_names();
  static std::string default_serial_dense_solver_name();
  static std::vector<std::string> parallel_dense_solver_names();
  static std::string default_parallel_dense_solver_name();
  static solver_factory* instance();
protected:
  sd_creator_pointer_type make_sd_creator(std::string const& name) const;
  pd_creator_pointer_type make_pd_creator(std::string const& name) const;
private:
  static solver_factory* instance_;
  sd_creator_map_type sd_creators_;
  pd_creator_map_type pd_creators_;
  int sd_largest_priority_;
  int pd_largest_priority_;
  std::string default_sd_solver_;
  std::string default_pd_solver_;
};

} // end namespace rokko

#define ROKKO_REGISTER_SERIAL_DENSE_SOLVER(solver, name, priority) \
  namespace { namespace BOOST_JOIN(solver_register, __LINE__) { struct register_caller { register_caller() { rokko::solver_factory::instance()->register_serial_dense_creator<solver>(name, priority); } } caller; } }

#define ROKKO_REGISTER_PARALLEL_DENSE_SOLVER(solver, name, priority) \
  namespace { namespace BOOST_JOIN(solver_register, __LINE__) { struct register_caller { register_caller() { rokko::solver_factory::instance()->register_parallel_dense_creator<solver>(name, priority); } } caller; } }

#endif // ROKKO_SOLVER_FACTORY_H
