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

#include <rokko/parameters.hpp>

namespace rokko {
namespace lapack {

char get_matrix_part(parameters const& params) {
  std::string matrix_part;

  if (params.defined("matrix_part"))
    matrix_part = params.get_string("matrix_part");
  else if (params.defined("uplow"))
    matrix_part = params.get_string("uplow");

  if (!matrix_part.empty()) {
    char matrix_part_letter = matrix_part[0];
    if ((matrix_part_letter == 'u') || (matrix_part_letter == 'U'))
      return 'U';
    if ((matrix_part_letter == 'l') || (matrix_part_letter == 'L'))
      return 'L';
  }
  return 'U';  // default
}
  
std::string get_matrix_part(char const& uplow) {
  if (uplow == 'U')
    return "upper";
  else if (uplow == 'L')
    return "lower";
  else
    return "none";
}

template<typename T>
char get_eigenvalues_range(parameters const& params, T& vl, T& vu, int& il, int& iu) {
  bool is_lower_value = get_key(params, "lower_value", vl);
  bool is_lower_index = get_key(params, "lower_index", il);
  bool is_upper_value = get_key(params, "upper_value", vu);
  bool is_upper_index = get_key(params, "upper_index", iu);

  if (is_lower_index && is_upper_index)   return 'I';
  if (is_lower_value && is_upper_value)   return 'V';
  else if (!(is_lower_index && is_lower_value && is_upper_index && is_upper_value))
    return 'A';
  else {
    throw std::invalid_argument("error: sepcify either of a pair of upper_value and lower_value or a pair of upper_index and lower_index");
  }
}

void print_verbose(std::string const& routine, char const& jobz, char const& uplow) {
  std::cout << "specified solver: " << routine << std::endl;
  if (jobz == 'V') {
    std::cout << "eigenvalues / eigenvectors were requested" << std::endl;
  } else if (jobz == 'N') {
    std::cout << "only eigenvalues were requested" << std::endl;
  } else {
    std::cout << "jobz '" << jobz << "' is not valid." << std::endl;
  }

  std::cout << "The " << get_matrix_part(uplow) << " part of the matrix was used" << std::endl;
}

template<typename T>
void print_verbose(std::string const& routine, char const& jobz, char const& range, char const& uplow,
		   T vl, T vu, int il, int iu) {
  std::cout << "specified solver: " << routine << std::endl;
  if (jobz == 'V') {
    std::cout << "eigenvalues / eigenvectors were requested" << std::endl;
  } else if (jobz == 'N') {
    std::cout << "only eigenvalues were requested" << std::endl;
  } else {
    std::cout << "jobz '" << jobz << "' is not valid." << std::endl;
  }

  if (range == 'A') {
    std::cout << "All eigenvalues were requested" << std::endl;
  } else if (range == 'V') {
    std::cout << "Eigenvalues contained in the interval [" << vl << ", " << vu << "]" << " were requested" << std::endl;
  } else if (range == 'I') {
    std::cout << "Eigenvalues from " << il << "th" << " to " << iu << "th" << " were requested" << std::endl;
  } else {
    std::cout << "range '" << range << "' is not valid." << std::endl;
  }

  std::cout << "The " << get_matrix_part(uplow) << " part of the matrix was used" << std::endl;
}

template<typename T>
void print_verbose(std::string const& routine, char const& jobz, char const& range, char const& uplow,
		   T vl, T vu, int il, int iu,
		   parameters const& params) {
  print_verbose(routine, jobz, range, uplow, vl, vu, il, iu);
  if (routine == "dsyevx") {
    std::cout << "The number of found eigenvalues = " << params.get_string("m") << std::endl;
    //std::cout << "ifail = " << params.get_string("ifail") << std::endl;
  }
  if (!params.defined("abstol")) {
    //std::cout << "abstol was not specified, so used optimal value for bisection method: 2 * LAPACKE_dlamch('S')" << std::endl;
    std::cout << "specified value: abstol=" << params.get_string("abstol") << std::endl;
  }
}

bool is_interval(parameters const& params) {
  if (params.defined("lower_value") || params.defined("upper_value") || params.defined("lower_index") || params.defined("upper_index")) {
    return true;
  }
  else if (params.defined("vl") || params.defined("vu") || params.defined("il") || params.defined("iu")) {
    return true;
  } else {
    return false;
  }
}

template
char get_eigenvalues_range(parameters const& params, double& vl, double& vu, int& il, int& iu);

void print_verbose(std::string const& routine, char const& jobz, char const& uplow);

template
void print_verbose(std::string const& routine, char const& jobz, char const& range, char const& uplow,
		   double vl, double vu, int il, int iu);

template void print_verbose<double>(std::string const& routine, char const& jobz, char const& range, char const& uplow,
		   double vl, double vu, int il, int iu,
		   parameters const& params);

} // namespace lapack
} // namespace rokko
