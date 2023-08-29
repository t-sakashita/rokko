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

#include <regex>
#include <string>
#include <sstream>      // for std::istringstream
#include <list>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <rokko/utility/string_trim.hpp>

namespace rokko {

namespace detail {

std::string retrieve_before_comment(std::string const& str) {
  const std::regex separator{"#"};
  const auto it = std::sregex_token_iterator{str.cbegin(), str.cend(), separator, -1};
  const auto end_it = std::sregex_token_iterator{};
  if (it != end_it)
    return *it;
  else
    return str;
}

bool read_line_with_comment(std::ifstream& ifs, std::istringstream& is) {
  std::string str_line;
  std::getline(ifs, str_line);
  const std::string trimed_str = trim_copy(retrieve_before_comment(str_line));
  is.clear();
  is.str(trimed_str);
  //std::cout << "no_comment:" << trimed_str << std::endl;
  //std::cout << "str_line:" << str_line << std::endl;
  return !trimed_str.empty(); // empty means the sentence is just comment
}

auto split_by_symbol(std::string const& str_line) {
  const std::regex separator{"="};
  std::vector<std::string> vec;
  for (std::sregex_token_iterator it{str_line.cbegin(), str_line.cend(), separator, -1}, end; it != end; ++it) {
    vec.emplace_back(trim_copy(*it));
  }

  return vec;
}

bool detect_offset_info(std::string const& str_line, bool& offset1) {
  const auto tokens = split_by_symbol(str_line);
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
    const std::ifstream::pos_type file_pos = ifs.tellg();  // save file position
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

} // namespace rokko
