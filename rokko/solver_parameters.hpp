/*****************************************************************************
*
* Rokko: Integrated Interface for libraries of eigenvalue decomposition
*
* Copyright (C) 2012-2016 by Rokko Developers https://github.com/t-sakashita/rokko
*
* Distributed under the Boost Software License, Version 1.0. (See accompanying
* file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
*
*****************************************************************************/

#ifndef ROKKO_SOLVER_PARAMETERS_HPP
#define ROKKO_SOLVER_PARAMETERS_HPP

namespace rokko {

#define ARRAY_SIZE(array) (sizeof *ARRAY_SIZE_(&(array)))

#define ARRAY_END(array) &array[(sizeof *ARRAY_SIZE_(&(array)))]

template <typename T, size_t N>
char (*ARRAY_SIZE_(T (*)[N]))[N];

static const std::vector<std::string> rokko_solver_keys{ "num_eigenvalues", "routine" };

bool is_rokko_solver_key(std::string const& key) {
  if (std::find(rokko_solver_keys.begin(), rokko_solver_keys.end(), key) != rokko_solver_keys.end()) {
    //std::cout << key << " is rokko key" << std::endl;
    return true;
  } else {
    //std::cout << key << " is NOT rokko key" << std::endl;
    return false;
  }
}

} // end namespace rokko

#endif // ROKKO_SOLVER_PARAMETERS_HPP
