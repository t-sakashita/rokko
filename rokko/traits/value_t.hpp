/*****************************************************************************
*
* Rokko: Integrated Interface for libraries of eigenvalue decomposition
*
* Copyright (C) 2017-2019 by Rokko Developers https://github.com/t-sakashita/rokko
*
* Distributed under the Boost Software License, Version 1.0. (See accompanying
* file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
*
*****************************************************************************/

#ifndef ROKKO_TRAITS_VALUE_T_HPP
#define ROKKO_TRAITS_VALUE_T_HPP

#include <complex>

namespace rokko {

namespace detail {

template<typename T>
struct value_type_traits {
  using type = typename T::value_type;
};

template<typename T>
struct value_type_traits<T*> {
  using type = T;
};

template<typename T>
struct value_type_traits<T**> {
  using type = T;
};

template<>
struct value_type_traits<float> {
  using type = float;
};

template<>
struct value_type_traits<std::complex<float>> {
  using type = std::complex<float>;
};

template<>
struct value_type_traits<double> {
  using type = double;
};

template<>
struct value_type_traits<std::complex<double>> {
  using type = std::complex<double>;
};

} // end namespace detail

template<typename T>
using value_t = typename detail::value_type_traits<T>::type;

} // namespace rokko

#endif // ROKKO_TRAITS_VALUE_T_HPP
