/*****************************************************************************
*
* Rokko: Integrated Interface for libraries of eigenvalue decomposition
*
* Copyright (C) 2017 by Rokko Developers https://github.com/t-sakashita/rokko
*
* Distributed under the Boost Software License, Version 1.0. (See accompanying
* file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
*
*****************************************************************************/

#ifndef ROKKO_TRAITS_NORM_T_HPP
#define ROKKO_TRAITS_NORM_T_HPP

#include "value_t.hpp"
#include <complex>
#include <boost/numeric/ublas/traits.hpp>

namespace rokko {

template<typename T>
struct norm_t {
  using value_type = typename value_t<T>::type;
  using type = typename boost::numeric::ublas::type_traits<value_type>::real_type;
};

} // namespace rokko

#endif // ROKKO_TRAITS_NORM_T_HPP
