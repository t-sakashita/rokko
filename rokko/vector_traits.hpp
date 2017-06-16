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

#ifndef ROKKO_VECTOR_TRAITS_HPP
#define ROKKO_VECTOR_TRAITS_HPP

namespace rokko {

template<typename VECTOR>
struct vector_traits {
  typedef typename VECTOR::value_type value_type;
};

} // namespace rokko

#endif // ROKKO_VECTOR_TRAITS_HPP
