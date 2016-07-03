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

static const char* const rokko_solver_keys[] =  { "num_eigenvalues", "routine" };

bool is_rokko_solver_key(std::string const& key) {
  std::vector<std::string> vec(&rokko_solver_keys[0], ARRAY_END(rokko_solver_keys));
  if (std::find(vec.begin(), vec.end(), key) != vec.end()) {
    //std::cout << key << " is rokko key" << std::endl;
    return true;
  } else {
    //std::cout << key << " is NOT rokko key" << std::endl;
    return false;
  }
}

} // end namespace rokko

#endif // ROKKO_SOLVER_PARAMETERS_HPP
