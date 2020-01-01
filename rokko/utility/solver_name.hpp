/*****************************************************************************
*
* Rokko: Integrated Interface for libraries of eigenvalue decomposition
*
* Copyright (C) 2012-2019 Rokko Developers https://github.com/t-sakashita/rokko
*
* Distributed under the Boost Software License, Version 1.0. (See accompanying
* file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
*
*****************************************************************************/

#ifndef ROKKO_UTILITY_SOLVER_NAME_HPP
#define ROKKO_UTILITY_SOLVER_NAME_HPP

#include <regex>
#include <vector>

namespace rokko {

void split_solver_name(std::string const& str, std::string& library, std::string& routine) {
  const std::regex separator{":"};
  auto it = std::sregex_token_iterator{str.cbegin(), str.cend(), separator, -1};
  const auto end_it = std::sregex_token_iterator{};
  if (it != end_it)
    library = *it++;

  if (it != end_it)
    routine = *it++;

  if (it != end_it)
    throw std::invalid_argument("split_solver_name() : More than one colon were detected from given string.");
}

} // namespace rokko

#endif // ROKKO_UTILITY_SOLVER_NAME_HPP
