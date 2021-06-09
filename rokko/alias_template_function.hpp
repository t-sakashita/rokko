/*****************************************************************************
*
* Rokko: Integrated Interface for libraries of eigenvalue decomposition
*
* Copyright (C) 2012-2020 Rokko Developers https://github.com/t-sakashita/rokko
*
* Distributed under the Boost Software License, Version 1.0. (See accompanying
* file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
*
*****************************************************************************/

#pragma once

#define ALIAS_TEMPLATE_FUNCTION(target, source)  \
template<typename... Types> \
inline auto target(Types&&... args) -> decltype(source(std::forward<Types>(args)...)) { \
  return source(std::forward<Types>(args)...); \
}
