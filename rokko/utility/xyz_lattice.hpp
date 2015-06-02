/*****************************************************************************
*
* Rokko: Integrated Interface for libraries of eigenvalue decomposition
*
* Copyright (C) 2013-2015 Rokko Developers https://github.com/t-sakashita/rokko
*
* Distributed under the Boost Software License, Version 1.0. (See accompanying
* file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
*
*****************************************************************************/

#ifndef ROKKO_UTILITY_XYZ_LATTICE_HPP
#define ROKKO_UTILITY_XYZ_LATTICE_HPP

#include <string>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <boost/tuple/tuple.hpp>
#include <boost/algorithm/string.hpp>

namespace rokko {

int read_lattice_stream(std::ifstream& ifs, int& num_sites, std::vector<std::pair<int, int> >& lattice, std::vector<boost::tuple<double, double, double> >& coupling) {
  int num_bonds;

  std::string str_line;
  std::ifstream::pos_type file_pos;
  getline(ifs, str_line);
  std::istringstream is(str_line);
  if ((str_line.find("#") != 0) && (str_line.size() != 0)) { // not comment
    is >> num_sites >> num_bonds;
  }
  std::cout << "num_sites=" << num_sites << " num_bonds=" << num_bonds << std::endl;
  bool start_index1 = false;
  do {
    file_pos = ifs.tellg();
    getline(ifs, str_line);
  } while (!str_line.empty());

  if (str_line.find("start_index = 1") == 0) {
    start_index1 = true;
    std::cout << "start_index = 1" << std::endl;
  } else if (str_line.find("start_index = 0") == 0) {
    start_index1 = false;
    std::cout << "start_index = 0" << std::endl;
  } else {
    std::cout << "else file_pos" << std::endl;
    ifs.seekg(file_pos);
  }

  do {
    int j, k;
    getline(ifs, str_line);
    boost::trim(str_line);
    std::istringstream is(str_line);
    if ((str_line.find("#") != 0) && (str_line.size() != 0)) { // not comment
      is >> j >> k;
      std::cout << "j=" << j << " k=" << k << std::endl;
      if (start_index1)  lattice.push_back(std::make_pair(j-1, k-1));
      else  lattice.push_back(std::make_pair(j, k));
      if (lattice.back().first >= num_sites) {
	std::cerr << "error: first index of"  << lattice.size() - 1 << "-th bond \"" << lattice.back().first << "\" is out of range" << std::endl;
	return 1;
      } else if (lattice.back().second >= num_sites) {
	std::cerr << "error: second index of " << lattice.size() - 1 << "-th bond \"" << lattice.back().second << "\" is out of range" << std::endl;
	return 2;
      }
    } else {
      std::cout << "comment:" << str_line << std::endl;
    }
  } while (lattice.size() < num_bonds);
  do {
    double jx, jy, jz;
    getline(ifs, str_line);
    boost::trim(str_line);
    std::istringstream is(str_line);
    if ((str_line.find("#") != 0) && (str_line.size() != 0)) { // not comment
      is >> jx >> jy >> jz;
      std::cout << "jx=" << jx << " jy=" << jy << " jz=" << jz << std::endl;
      coupling.push_back(boost::make_tuple(jx, jy, jz));
    } else {      
      std::cout << "comment:" << str_line << std::endl;
    }
  } while (coupling.size() < num_bonds);
  return 0;
}

int read_lattice_file(std::string const& filename, int& num_sites, std::vector<std::pair<int, int> >& lattice, std::vector<boost::tuple<double, double, double> >& coupling) {
  std::ifstream ifs(filename);
  if (!ifs) {
    std::cout << "can't open file" << std::endl;
    return 1;
  }
  return read_lattice_stream(ifs, num_sites, lattice, coupling);
}

} // namespace rokko

#endif // ROKKO_UTILITY_XYZ_LATTICE_HPP
