/*****************************************************************************
*
* Rokko: Integrated Interface for libraries of eigenvalue decomposition
*
* Copyright (C) 2013-2019 Rokko Developers https://github.com/t-sakashita/rokko
*
* Distributed under the Boost Software License, Version 1.0. (See accompanying
* file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
*
*****************************************************************************/

#ifndef ROKKO_UTILITY_XYZ_LATTICE_HPP
#define ROKKO_UTILITY_XYZ_LATTICE_HPP

#include <string>
#include <list>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <tuple>
#include <boost/algorithm/string.hpp>
#include <boost/tokenizer.hpp>

#include <rokko/utility/lattice.hpp>

namespace rokko {

void read_lattice_stream(std::ifstream& ifs, int& num_sites, std::vector<std::pair<int, int>>& lattice, std::vector<std::tuple<double, double, double>>& coupling) {
  std::size_t num_bonds;
  std::istringstream is;
  if (detail::read_line_with_comment(ifs, is)) {
    is >> num_sites >> num_bonds;
  }
  std::cout << "num_sites=" << num_sites << " num_bonds=" << num_bonds << std::endl;

  bool offset1 = detail::read_offset_info(ifs);

  do {
    int j, k;
    if (detail::read_line_with_comment(ifs, is)) {
      is >> j >> k;
      std::cout << "j=" << j << " k=" << k << std::endl;
      if (offset1)  lattice.push_back(std::make_pair(j-1, k-1));
      else  lattice.push_back(std::make_pair(j, k));
      //std::cout << "back()=" << lattice.back().first << ", " << lattice.back().second << std::endl;
      if ((lattice.back().first < 0) || (lattice.back().first >= num_sites)) {
        throw std::invalid_argument("read_lattice_stream() : first index of " + std::to_string(lattice.size() - 1) + "-th bond \"" + std::to_string(lattice.back().first) + "\" is out of range");
      } else if ((lattice.back().second < 0) || (lattice.back().second >= num_sites)) {
        throw std::invalid_argument("read_lattice_stream() : second index of " + std::to_string(lattice.size() - 1) + "-th bond \"" + std::to_string(lattice.back().second) + "\" is out of range");
      }
    }
  } while (lattice.size() < num_bonds);

  do {
    double jx, jy, jz;
    if (detail::read_line_with_comment(ifs, is)) {
      is >> jx >> jy >> jz;
      std::cout << "jx=" << jx << " jy=" << jy << " jz=" << jz << std::endl;
      coupling.push_back(std::make_tuple(jx, jy, jz));
    }
  } while (coupling.size() < num_bonds);
}

void read_lattice_file(std::string const& filename, int& num_sites, std::vector<std::pair<int, int>>& lattice, std::vector<std::tuple<double, double, double>>& coupling) {
  std::ifstream ifs(filename);
  if (!ifs) {
    throw std::runtime_error("read_lattice_file() : can't open file \"" + filename + "\"");
  }
  return read_lattice_stream(ifs, num_sites, lattice, coupling);
}

} // namespace rokko

#endif // ROKKO_UTILITY_XYZ_LATTICE_HPP
