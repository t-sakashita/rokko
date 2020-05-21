/*****************************************************************************
*
* Rokko: Integrated Interface for libraries of eigenvalue decomposition
*
* Copyright (C) 2012-2020 by Rokko Developers https://github.com/t-sakashita/rokko
*
* Distributed under the Boost Software License, Version 1.0. (See accompanying
* file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
*
*****************************************************************************/

#ifndef ROKKO_TRAITS_REAL_T_HPP
#define ROKKO_TRAITS_REAL_T_HPP

#include "value_t.hpp"
#include <complex>

namespace rokko {

namespace detail {

template<typename T>
struct norm_scalar_type_traits {
};

template<>
struct norm_scalar_type_traits<float> {
  using type = float;
};

template<>
struct norm_scalar_type_traits<std::complex<float>> {
  using type = float;
};

template<>
struct norm_scalar_type_traits<double> {
  using type = double;
};

template<>
struct norm_scalar_type_traits<std::complex<double>> {
  using type = double;
};

template<typename T>
struct real_type_traits {
  using value_type = rokko::value_t<T>;
  using type = typename norm_scalar_type_traits<value_type>::type;
};

} // end namespace detail

template<typename T>
using real_t = typename detail::real_type_traits<T>::type;

} // namespace rokko

#endif // ROKKO_TRAITS_REAL_T_HPP
