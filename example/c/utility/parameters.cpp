/*****************************************************************************
*
* Rokko: Integrated Interface for libraries of eigenvalue decomposition
*
* Copyright (C) 2015 Rokko Developers https://github.com/t-sakashita/rokko
*
* Distributed under the Boost Software License, Version 1.0. (See accompanying
* file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
*
*****************************************************************************/

#include <rokko/parameters.hpp>
#include <iostream>

int main(int argc, char **argv) {
  rokko_parameters params;
  char p = 'V';
  // set some parameters
  set_double(params, "T", 0.95);
  set_double(params, "Pi", 3.14);
  set_char_string(params, "A", "test");
  //params.set("C", 'A');
  set_char(params, "C", p);

  // get double
  double t = params.get<double>("T");

  // get string
  std::string a = params.get<std::string>("A");

  // get char
  char c = params.get<char>("C");

  // is "T" defined?
  bool t_defined = params.defined("T");
  
  // is "S" defined?
  bool s_defined = params.defined("S");

  // clear paramter "Pi"
  params.clear("Pi");
  
  // output list of parameters as string
  std::list<std::string> keys = params.keys();
  std::cout << "number of parameters = " << keys.size() << std::endl;
  std::cout << "key-value pairs:\n";
  BOOST_FOREACH(std::string const& key, keys) {
    std::cout << "  " << key << " = " << params.get_string(key) << std::endl;
  }
}
