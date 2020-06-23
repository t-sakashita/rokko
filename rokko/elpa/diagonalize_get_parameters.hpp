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

#ifndef ROKKO_ELPA_DIAGONALIZE_GET_PARAMETERS_HPP
#define ROKKO_ELPA_DIAGONALIZE_GET_PARAMETERS_HPP

#include <rokko/parameters.hpp>
#include <rokko/elpa/elpa.h>

namespace rokko {
namespace elpa {

void get_nev(parameters const& params, int& nev) {
  if (params.defined("num_eigvals"))
    nev = params.get<int>("num_eigvals");
  if (params.defined("nev"))
    nev = params.get<int>("nev");
}

int get_kernel(parameters const& params) {
  return params.defined("kernel") ? params.get<int>("kernel") : ELPA_2STAGE_REAL_GENERIC_SIMPLE;
}

int get_blocked_qr(parameters const& params) {
  return params.defined("blocked_qr") ? static_cast<int>(params.get<bool>("blocked_qr")) : 0;
}

} // namespace elpa
} // namespace rokko

#endif // ROKKO_ELPA_DIAGONALIZE_GET_PARAMETERS_HPP

