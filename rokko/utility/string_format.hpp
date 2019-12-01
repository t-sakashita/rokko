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

#ifndef ROKKO_UTIL_STRING_FORMAT_HPP
#define ROKKO_UTIL_STRING_FORMAT_HPP

#include <cstdio>
#include <string>
#include <vector>

namespace rokko {

template <typename ... Args>
std::string format(const std::string& fmt, Args ... args) {
  size_t len = std::snprintf( nullptr, 0, fmt.c_str(), args ... );
  std::vector<char> buf(len + 1);
  std::snprintf(buf.data(), len + 1, fmt.c_str(), args ... );
  return std::string(buf.data(), buf.data() + len);
}

} // namespace rokko

#endif // ROKKO_UTIL_STRING_FORMAT_HPP
