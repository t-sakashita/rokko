/*****************************************************************************
*
* Rokko: Integrated Interface for libraries of eigenvalue decomposition
*
* Copyright (C) 2012-2015 Rokko Developers https://github.com/t-sakashita/rokko
*
* Distributed under the Boost Software License, Version 1.0. (See accompanying
* file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
*
*****************************************************************************/

#ifndef ROKKO_UTILITY_SOLVER_NAME_HPP
#define ROKKO_UTILITY_SOLVER_NAME_HPP

#include <iostream>
#include <vector>
#include <boost/throw_exception.hpp>
#include <boost/algorithm/string.hpp>

namespace rokko {

void split_solver_name(std::string const& str, std::string& library, std::string& routine) {
  std::vector<std::string> v;
  boost::algorithm::split(v, str, boost::is_any_of(":"));
  if (v.size() > 2) {
    BOOST_THROW_EXCEPTION(std::invalid_argument("split_solver_name() : More than colon were detected from given string."));
  }
  if (v.size() >= 1)
    library = v[0];
  if (v.size() == 2)
    routine = v[1];
}
    
} // namespace rokko

#endif // ROKKO_UTILITY_SOLVER_NAME_HPP
