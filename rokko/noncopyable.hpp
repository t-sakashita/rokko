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

#ifndef ROKKO_NONCOPYABLE_HPP
#define ROKKO_NONCOPYABLE_HPP

namespace rokko {

template <class T>
class noncopyable {
protected:
  noncopyable() = default;
  ~noncopyable() = default;

public:
  noncopyable(const noncopyable&) = delete;
  noncopyable& operator=(const noncopyable&) = delete;
};

} // end namespace rokko

#endif // ROKKO_NONCOPYABLE_HPP
