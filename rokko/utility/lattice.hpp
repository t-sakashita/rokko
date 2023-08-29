/*****************************************************************************
*
* Rokko: Integrated Interface for libraries of eigenvalue decomposition
*
* Copyright (C) 2013-2023 Rokko Developers https://github.com/t-sakashita/rokko
*
* Distributed under the Boost Software License, Version 1.0. (See accompanying
* file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
*
*****************************************************************************/

#pragma once

#include <sstream>      // for std::istringstream
#include <tuple>
#include <rokko/utility/lattice_file_common.hpp>

namespace rokko {

auto read_lattice_stream(std::ifstream& ifs) {
  const auto [num_sites, num_bonds] = detail::read_num_sites_bonds(ifs);
  std::cout << "num_sites=" << num_sites << " num_bonds=" << num_bonds << std::endl;

  const auto offset1 = detail::read_offset_info(ifs);

  std::vector<std::pair<int,int>> lattice;
  std::istringstream is;
  do {
    int j, k;
    if (detail::read_line_with_comment(ifs, is)) {
      is >> j >> k;
      std::cout << "j=" << j << " k=" << k << std::endl;
      if (offset1)  lattice.emplace_back(std::make_pair(j-1, k-1));
      else  lattice.emplace_back(std::make_pair(j, k));
      //std::cout << "back()=" << lattice.back().first << ", " << lattice.back().second << std::endl;
      if ((lattice.back().first < 0) || (lattice.back().first >= num_sites)) {
        throw std::invalid_argument("read_lattice_stream() : first index of " + std::to_string(lattice.size() - 1) + "-th bond \"" + std::to_string(lattice.back().first) + "\" is out of range");
      } else if ((lattice.back().second < 0) || (lattice.back().second >= num_sites)) {
        throw std::invalid_argument("read_lattice_stream() : second index of " + std::to_string(lattice.size() - 1) + "-th bond \"" + std::to_string(lattice.back().second) + "\" is out of range");
      }
    }
  } while (lattice.size() < num_bonds);

  return std::tuple(num_sites, lattice);
}

auto read_lattice_file(std::string const& filename) {
  std::ifstream ifs(filename);
  if (!ifs) {
    throw std::runtime_error("read_lattice_file() : can't open file \"" + filename + "\"");
  }
  return read_lattice_stream(ifs);
}

auto create_ladder_lattice_1dim(int len_ladder) {
  std::vector<std::pair<int, int>> lattice;

  const auto L = 2 * len_ladder;
  for (std::size_t i = 0; i < (len_ladder-1); ++i) {
    lattice.emplace_back(std::make_pair(i, i+1));
  }
  for (std::size_t i = len_ladder; i < (L-1); ++i) {
    lattice.emplace_back(std::make_pair(i, i+1));
  }
  for (std::size_t i = 0; i < len_ladder; ++i) {
    lattice.emplace_back(std::make_pair(i, i+len_ladder));
  }

  return lattice;
}

void output_lattice(std::ostream& os, std::vector<std::pair<int, int>> const& lattice) {
  for (std::size_t i=0; i<lattice.size(); ++i) {
    os << "no=" << i << " <" << lattice[i].first << ", " << lattice[i].second << ">" << std::endl;
  }
}

void print_lattice(std::vector<std::pair<int, int>> const& lattice) {
  output_lattice(std::cout, lattice);
}

} // namespace rokko
