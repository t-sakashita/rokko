/*****************************************************************************
*
* Rokko: Integrated Interface for libraries of eigenvalue decomposition
*
* Copyright (C) 2010-2019 Rokko Developers https://github.com/t-sakashita/rokko
*
* Distributed under the Boost Software License, Version 1.0. (See accompanying
* file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
*
*****************************************************************************/

#ifndef ROKKO_UTIL_TUPLE_TO_ARRAY_HPP
#define ROKKO_UTIL_TUPLE_TO_ARRAY_HPP

#include <array>
#include <tuple>
#include <utility>

namespace rokko {

template<typename T, std::size_t N, typename TUPLE, std::size_t... I>
constexpr decltype(auto) to_array_impl(const TUPLE& a, std::index_sequence<I...>) {
  return std::array<T,N>{std::get<I>(a)...};
}

template<typename HEAD, typename... T>
constexpr decltype(auto) to_array(const std::tuple<HEAD, T...>& a) {
  using TUPLE = std::tuple<HEAD, T...>;
  constexpr auto N = sizeof...(T) + 1;
  return to_array_impl<HEAD, N, TUPLE>(a, std::make_index_sequence<N>());
}

} // namespace rokko

#endif // ROKKO_UTIL_TUPLE_TO_ARRAY_HPP
