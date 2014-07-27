/*****************************************************************************
*
* Rokko: Integrated Interface for libraries of eigenvalue decomposition
*
* Copyright (C) 2012-2014 by Tatsuya Sakashita <t-sakashita@issp.u-tokyo.ac.jp>,
*                            Synge Todo <wistaria@comp-phys.org>
*
* Distributed under the Boost Software License, Version 1.0. (See accompanying
* file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
*
*****************************************************************************/

#ifndef ROKKO_DISTRIBUTED_MFREE_HPP
#define ROKKO_DISTRIBUTED_MFREE_HPP

#include <rokko/factory.hpp>
#include <rokko/mapping_1d.hpp>

namespace rokko {


class distributed_operator {
public:
  distributed_operator() {}
  ~distributed_operator() {}
  void multiply(const double* x, double* y) const {};
};


namespace detail {

class df_matrix_base {
public:
  df_matrix_base() {}
  ~df_matrix_base() {}
  virtual void initialize(mapping_1d const& map) = 0;
  virtual void define_operator(rokko::distributed_operator& op) = 0;
};
    
typedef factory<df_matrix_base> df_matrix_factory;
  
} // end namespace detail

class distributed_mfree {
public:
  distributed_mfree(std::string const& solver_name, mapping_1d const& map) {
    matrix_impl_ = detail::df_matrix_factory::instance()->make_product(solver_name);
    matrix_impl_->initialize(map);
  }
  distributed_mfree() {
    matrix_impl_ = detail::df_matrix_factory::instance()->make_product();
  }
  void define_operator(rokko::distributed_operator& op) {
    matrix_impl_->define_operator(op);
  }
// private:
  detail::df_matrix_factory::product_pointer_type matrix_impl_;
};

} // end namespace rokko

#endif // ROKKO_DISTRIBUTED_MFREE_HPP

#define ROKKO_REGISTER_DISTRIBUTED_MFREE(matrix, name, priority)    \
namespace { namespace BOOST_JOIN(register, __LINE__) { \
struct register_caller { \
  typedef rokko::factory<rokko::detail::df_matrix_base> factory; \
  register_caller() { factory::instance()->register_creator<matrix>(name, priority); } \
} caller; \
} }
