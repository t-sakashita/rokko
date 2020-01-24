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

#include <regex>
#include <string>
#include <sstream>      // for std::istringstream
#include <list>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <tuple>
#include <rokko/utility/string_trim.hpp>

namespace rokko {

namespace detail {

std::string retrieve_before_comment(std::string const& str) {
  std::regex separator{"#"};
  auto it = std::sregex_token_iterator{str.cbegin(), str.cend(), separator, -1};
  const auto end_it = std::sregex_token_iterator{};
  if (it != end_it)
    return *it;
  else
    return str;
}

bool read_line_with_comment(std::ifstream& ifs, std::istringstream& is) {
  std::string str_line;
  std::getline(ifs, str_line);
  std::string trimed_str = trim_copy(retrieve_before_comment(str_line));
  is.clear();
  is.str(trimed_str);
  //std::cout << "no_comment:" << trimed_str << std::endl;
  //std::cout << "str_line:" << str_line << std::endl;
  return !trimed_str.empty(); // empty means the sentence is just comment
}

std::vector<std::string> split_by_symbol(std::string const& str_line) {
  const std::regex separator{"="};
  std::vector<std::string> vec;
  for (std::sregex_token_iterator it{str_line.cbegin(), str_line.cend(), separator, -1}, end; it != end; ++it) {
    vec.emplace_back(trim_copy(*it));
  }

  return vec;
}

bool detect_offset_info(std::string const& str_line, bool& offset1) {
  auto tokens = split_by_symbol(str_line);
  if (tokens.size() < 2) {
    return false;
  }

  if (tokens[0]=="offset") {
    if (tokens[1]=="1") {
      offset1 = true;
      std::cout << "offset = 1" << std::endl;
    } else if (tokens[1]=="0") {
      offset1 = false;
      std::cout << "offset = 0" << std::endl;
    } else
      throw std::invalid_argument("detail::detect_offset_info() : give 0 or 1 after 'offset='");
    return true;
  } else {
    return false;
  }
}

bool read_offset_info(std::ifstream& ifs) {
  bool offset1 = false;
  std::string str_line;
  while(true) {
    std::ifstream::pos_type file_pos = ifs.tellg();  // save file position
    std::getline(ifs, str_line);
    if (detect_offset_info(str_line, offset1)) {
      break;
    } else if (!str_line.empty()) {  // if str_line is first sentence
      ifs.seekg(file_pos);  // resotre file position
      break;
    }
  }
  return offset1;
}

} // namespace detail

void read_lattice_stream(std::ifstream& ifs, int& num_sites, std::vector<std::pair<int, int>>& lattice) {
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
}

void read_lattice_file(std::string const& filename, int& num_sites, std::vector<std::pair<int, int>>& lattice) {
  std::ifstream ifs(filename);
  if (!ifs) {
    throw std::runtime_error("read_lattice_file() : can't open file \"" + filename + "\"");
  }
  return read_lattice_stream(ifs, num_sites, lattice);
}

void create_ladder_lattice_1dim(int len_ladder, std::vector<std::pair<int, int>>& lattice) {
  int L = 2 * len_ladder;
  for (std::size_t i = 0; i < (len_ladder-1); ++i) {
    lattice.emplace_back(std::make_pair(i, i+1));
  }
  for (std::size_t i = len_ladder; i < (L-1); ++i) {
    lattice.emplace_back(std::make_pair(i, i+1));
  }
  for (std::size_t i = 0; i < len_ladder; ++i) {
    lattice.emplace_back(std::make_pair(i, i+len_ladder));
  }
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

#endif // ROKKO_UTILITY_LATTICE_HPP
