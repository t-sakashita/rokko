/*****************************************************************************
*
* Rokko: Integrated Interface for libraries of eigenvalue decomposition
*
* Copyright (C) 2012-2019 by Rokko Developers https://github.com/t-sakashita/rokko
*
* Distributed under the Boost Software License, Version 1.0. (See accompanying
* file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
*
*****************************************************************************/

#ifndef ROKKO_SOLVER_PARAMETERS_HPP
#define ROKKO_SOLVER_PARAMETERS_HPP

namespace rokko {

static const std::vector<std::string> rokko_solver_keys{ "num_eigenvalues", "max_block_size", "routine", "block_size", "conv_tol", "max_iters", "verbose" };

bool is_rokko_solver_key(std::string const& key) {
  return std::find(rokko_solver_keys.begin(), rokko_solver_keys.end(), key) != rokko_solver_keys.end();
}

} // end namespace rokko

#endif // ROKKO_SOLVER_PARAMETERS_HPP
