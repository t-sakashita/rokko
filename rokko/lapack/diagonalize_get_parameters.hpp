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

#pragma once

namespace rokko {
namespace lapack {

char get_matrix_part(parameters const& params);

std::string get_matrix_part(char const& uplow);

template<typename T>
char get_eigenvalues_range(parameters const& params, T& vl, T& vu, int& il, int& iu);

void print_verbose(std::string const& routine, char const& jobz, char const& uplow);

template<typename T>
void print_verbose(std::string const& routine, char const& jobz, char const& range, char const& uplow,
		   T vl, T vu, int il, int iu);

template<typename T>
void print_verbose(std::string const& routine, char const& jobz, char const& range, char const& uplow,
		   T vl, T vu, int il, int iu,
		   parameters const& params);

bool is_interval(parameters const& params);

} // namespace lapack
} // namespace rokko
