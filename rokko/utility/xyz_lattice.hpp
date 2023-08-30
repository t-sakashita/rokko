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

#include <rokko/utility/lattice_file_common.hpp>

namespace rokko {

auto read_lattice_stream(std::ifstream& ifs) {
  const auto [num_sites, num_bonds] = detail::read_num_sites_bonds(ifs);
  std::cout << "num_sites=" << num_sites << " num_bonds=" << num_bonds << std::endl;

  const auto offset1 = detail::read_offset_info(ifs);

  std::vector<std::pair<int, int>> lattice;
  std::vector<std::tuple<double, double, double>> coupling;
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

  do {
    double jx, jy, jz;
    if (detail::read_line_with_comment(ifs, is)) {
      is >> jx >> jy >> jz;
      std::cout << "jx=" << jx << " jy=" << jy << " jz=" << jz << std::endl;
      coupling.emplace_back(std::make_tuple(jx, jy, jz));
    }
  } while (coupling.size() < num_bonds);

  return std::tuple(num_sites, lattice, coupling);
}

auto read_lattice_file(std::string const& filename) {
  std::ifstream ifs(filename);
  if (!ifs) {
    throw std::runtime_error("read_lattice_file() : can't open file \"" + filename + "\"");
  }
  return read_lattice_stream(ifs);
}

auto print_lattice_coupling(std::size_t num_sites, std::vector<std::pair<int, int>> const& lattice, std::vector<std::tuple<double, double, double>> const& coupling) {
  std::cout << "num_sites=" << num_sites << " num_bonds=" << lattice.size() << std::endl;

  for (auto i=0; i<lattice.size(); ++i) {
    std::cout << "No=" << i << "  bond=<" << lattice[i].first << "," << lattice[i].second << ">   coupling=(" << std::get<0>(coupling[i]) << "," << std::get<1>(coupling[i]) << "," << std::get<2>(coupling[i]) << ")" << std::endl;
  }
}

} // namespace rokko
