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
#include <list>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <boost/tuple/tuple.hpp>
#include <boost/algorithm/string.hpp>

namespace rokko {

namespace detail {

bool read_line_with_comment(std::ifstream& ifs, std::istringstream& is) {
  std::string str_line;
  getline(ifs, str_line);
  std::list<std::string> list_string;
  boost::split(list_string, str_line, boost::is_any_of("#"));
  std::string trimed_str = list_string.front();
  boost::trim(trimed_str);
  is.clear();
  is.str(trimed_str);
  //std::cout << "string:" << trimed_str << std::endl;
  //std::cout << "comment:" << list_string.back() << std::endl;
  if (!trimed_str.empty()) { // the sentece is not just comment
    return true;
  } else { // the sentence is just comment
    return false;
  }
}

} // namespace detail

void read_lattice_stream(std::ifstream& ifs, int& num_sites, std::vector<std::pair<int, int> >& lattice, std::vector<boost::tuple<double, double, double> >& coupling) {
  std::string str_line;
  int num_bonds;
  std::istringstream is;
  std::ifstream::pos_type file_pos;
  if (detail::read_line_with_comment(ifs, is)) {
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
    //std::cout << "else file_pos" << std::endl;
    ifs.seekg(file_pos);
  }

  do {
    int j, k;
    if (detail::read_line_with_comment(ifs, is)) {
      is >> j >> k;
      std::cout << "j=" << j << " k=" << k << std::endl;
      if (start_index1)  lattice.push_back(std::make_pair(j-1, k-1));
      else  lattice.push_back(std::make_pair(j, k));
      if (lattice.back().first >= num_sites) {
	std::cerr << "error: first index of"  << lattice.size() - 1 << "-th bond \"" << lattice.back().first << "\" is out of range" << std::endl;
	throw 1;
      } else if (lattice.back().second >= num_sites) {
	std::cerr << "error: second index of " << lattice.size() - 1 << "-th bond \"" << lattice.back().second << "\" is out of range" << std::endl;
	throw 2;
      }
    }
  } while (lattice.size() < num_bonds);
  do {
    double jx, jy, jz;
    if (detail::read_line_with_comment(ifs, is)) {
      is >> jx >> jy >> jz;
      std::cout << "jx=" << jx << " jy=" << jy << " jz=" << jz << std::endl;
      coupling.push_back(boost::make_tuple(jx, jy, jz));
    }
  } while (coupling.size() < num_bonds);
}

void read_lattice_file(std::string const& filename, int& num_sites, std::vector<std::pair<int, int> >& lattice, std::vector<boost::tuple<double, double, double> >& coupling) {
  std::ifstream ifs(filename.c_str());
  if (!ifs) {
    std::cerr << "can't open file" << std::endl;
    throw 1;
  }
  return read_lattice_stream(ifs, num_sites, lattice, coupling);
}

} // namespace rokko

#endif // ROKKO_UTILITY_XYZ_LATTICE_HPP
