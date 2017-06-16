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

#ifndef ROKKO_MATRIX_TRAITS_HPP
#define ROKKO_MATRIX_TRAITS_HPP

namespace rokko {

template<typename MATRIX>
struct matrix_traits {
  typedef typename MATRIX::value_type value_type;
  typedef typename MATRIX::major_type major_type;
};

} // namespace rokko

#endif // ROKKO_MATRIX_TRAITS_HPP
