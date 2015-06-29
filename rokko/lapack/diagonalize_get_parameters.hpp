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

#ifndef ROKKO_LAPACK_DIAGONALIZE_GET_PARAMMETERS_HPP
#define ROKKO_LAPACK_DIAGONALIZE_GET_PARAMMETERS_HPP

namespace rokko {
namespace lapack {

void get_matrix_part(rokko::parameters const& params, std::string& matrix_part, char& uplow) {
  if (params.defined("uplow"))
    matrix_part = params.get_string("uplow");
  if (params.defined("matrix_part"))
    matrix_part = params.get_string("matrix_part");
  if ((matrix_part[0] == 'u') || (matrix_part[0] == 'U'))
    matrix_part = "upper";  uplow = 'U';
  if ((matrix_part[0] == 'l') || (matrix_part[0] == 'L'))
    matrix_part = "lower";  uplow = 'L';
}


void get_eigenvalues_range(rokko::parameters const& params, std::string& matrix_part, char& range, double vu, double vl, int iu, int il, bool& is_upper_value, bool& is_lower_value, bool& is_upper_index, bool& is_lower_index) {
  is_upper_value = get_key(params, "upper_value", vu);
  is_upper_index = get_key(params, "upper_value", iu);
  is_lower_value = get_key(params, "lower_value", vl);
  is_lower_index = get_key(params, "lower_value", il);
  if (is_upper_index && is_lower_index)   range = 'I';
  if (is_upper_value && is_lower_value)   range = 'V';
  if (is_upper_index && is_lower_value) {
    std::cerr << "error: sepcify either of a pair of upper_value and lower_value or a pair of upper_index and lower_index";
    throw;
  }
}

} // namespace lapack
} // namespace rokko

#endif // ROKKO_LAPACK_DIAGONALIZE_GET_PARAMETERS_HPP

