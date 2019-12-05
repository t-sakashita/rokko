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

#ifndef ROKKO_UTILITY_LATTICE_HPP
#define ROKKO_UTILITY_LATTICE_HPP

#include <string>
#include <sstream>      // for std::istringstream
#include <list>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <tuple>
#include <boost/algorithm/string.hpp>
#include <boost/tokenizer.hpp>

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
  return !trimed_str.empty(); // empty means the sentence is just comment
}


bool detect_offset_info(std::string const& str_line, bool& offset1) {
  boost::char_separator<char>  sep(" ", "=", boost::drop_empty_tokens);
  boost::tokenizer<boost::char_separator<char>>  tokens0(str_line, sep);
  std::vector<std::string> tokens;
  std::copy(tokens0.begin(), tokens0.end(), std::back_inserter<std::vector<std::string>>(tokens));
  if (tokens.size() < 3) {
    return false;
  }

  if (tokens[0]=="offset" && tokens[1]=="=") {
    if (tokens[2]=="1") {
      offset1 = true;
      std::cout << "offset = 1" << std::endl;
    } else if (tokens[2]=="0") {
      offset1 = false;
      std::cout << "offset = 0" << std::endl;
    } else
      throw std::invalid_argument("detail::detect_offset_info() : give 0 or 1 after 'offset='");
    return true;
  } else {
    return false;
  }
}

void read_offset_info(std::ifstream& ifs, bool& offset1) {
  std::string str_line;
  std::ifstream::pos_type file_pos;
  do {
    file_pos = ifs.tellg();
    getline(ifs, str_line);
    if (detect_offset_info(str_line, offset1)) {
      return;
    }
  } while (str_line.empty());
  ifs.seekg(file_pos);  // resotre file position
}

} // namespace detail

void read_lattice_stream(std::ifstream& ifs, int& num_sites, std::vector<std::pair<int, int>>& lattice) {
  std::size_t num_bonds;
  std::istringstream is;
  if (detail::read_line_with_comment(ifs, is)) {
    is >> num_sites >> num_bonds;
  }
  std::cout << "num_sites=" << num_sites << " num_bonds=" << num_bonds << std::endl;

  bool offset1 = false;
  detail::read_offset_info(ifs, offset1);

  do {
    int j, k;
    if (detail::read_line_with_comment(ifs, is)) {
      is >> j >> k;
      std::cout << "j=" << j << " k=" << k << std::endl;
      if (offset1)  lattice.push_back(std::make_pair(j-1, k-1));
      else  lattice.push_back(std::make_pair(j, k));
      //std::cout << "back()=" << lattice.back().first << ", " << lattice.back().second << std::endl;
      if ((lattice.back().first < 0) || (lattice.back().first >= num_sites)) {
        std::stringstream msg;
        msg << "read_lattice_stream() : first index of " << lattice.size() - 1 << "-th bond \"" << lattice.back().first << "\" is out of range";
        throw std::invalid_argument(msg.str());
      } else if ((lattice.back().second < 0) || (lattice.back().second >= num_sites)) {
        std::stringstream msg;
        msg << "read_lattice_stream() : second index of " << lattice.size() - 1 << "-th bond \"" << lattice.back().first << "\" is out of range";
        throw std::invalid_argument(msg.str());
      }
    }
  } while (lattice.size() < num_bonds);
}

void read_lattice_file(std::string const& filename, int& num_sites, std::vector<std::pair<int, int>>& lattice) {
  std::ifstream ifs(filename.c_str());
  if (!ifs) {
    throw std::runtime_error("read_lattice_file() : can't open file \"" + filename + "\"");
  }
  return read_lattice_stream(ifs, num_sites, lattice);
}

void ladder_lattice_1dim(int len_ladder, std::vector<std::pair<int, int>>& lattice) {
  int L = 2 * len_ladder;
  int dim = 1 << L;
  for (std::size_t i = 0; i < (len_ladder-1); ++i) {
    lattice.push_back(std::make_pair(i, i+1));
  }
  for (std::size_t i = len_ladder; i < (L-1); ++i) {
    lattice.push_back(std::make_pair(i, i+1));
  }
  for (std::size_t i = 0; i < len_ladder; ++i) {
    lattice.push_back(std::make_pair(i, i+len_ladder));
  }
}

void print_lattice(std::vector<std::pair<int, int>>& lattice) {
  for (std::size_t i=0; i<lattice.size(); ++i) {
    std::cout << "no=" << i << " <" << lattice[i].first << ", " << lattice[i].second << ">" << std::endl;
  }
}

void output_lattice(std::ostream& os, std::vector<std::pair<int, int>>& lattice) {
  for (std::size_t i=0; i<lattice.size(); ++i) {
    os << "no=" << i << " <" << lattice[i].first << ", " << lattice[i].second << ">" << std::endl;
  }
}

} // namespace rokko

#endif // ROKKO_UTILITY_LATTICE_HPP
