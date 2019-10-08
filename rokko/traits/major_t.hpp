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

#ifndef ROKKO_TRAITS_MAJOR_T_HPP
#define ROKKO_TRAITS_MAJOR_T_HPP

#include <complex>

namespace rokko {

template<typename T>
struct major_t {
  using type = typename T::major_type;
};

} // namespace rokko

#endif // ROKKO_TRAITS_MAJOR_T_HPP
