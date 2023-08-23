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

#include <rokko/parameters.hpp>
#include <rokko/elpa/elpa.h>

namespace rokko {
namespace elpa {

auto get_nev(parameters const& params, int default_nev) {
  if (params.defined("num_eigvals"))
    return params.get<int>("num_eigvals");
  if (params.defined("nev"))
    return params.get<int>("nev");
  else
    return default_nev;
}

int get_kernel(parameters const& params) {
  return params.defined("kernel") ? params.get<int>("kernel") : ELPA_2STAGE_REAL_GENERIC_SIMPLE;
}

int get_blocked_qr(parameters const& params) {
  return params.defined("blocked_qr") ? static_cast<int>(params.get<bool>("blocked_qr")) : 0;
}

int get_solver(parameters const& params) {
  const std::string routine = params.defined("routine") ? params.get_string("routine") : "";

  if (routine=="elpa1") {
    return ELPA_SOLVER_1STAGE;
  } else if (routine=="elpa2") {
    return ELPA_SOLVER_2STAGE;
  } else if (routine=="") {  // default
    return ELPA_SOLVER_2STAGE;
  } else {
    throw std::invalid_argument("elpa::get_solver() : " + routine + " is invalid routine name");
  }
}

} // namespace elpa
} // namespace rokko
