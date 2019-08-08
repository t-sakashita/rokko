/*****************************************************************************
*
* Rokko: Integrated Interface for libraries of eigenvalue decomposition
*
* Copyright (C) 2012-2019 Rokko Developers https://github.com/t-sakashita/rokko
*
* Distributed under the Boost Software License, Version 1.0. (See accompanying
* file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
*
*****************************************************************************/

#include <rokko/parameters.hpp>
#include <iostream>

int main() {
  rokko::parameters params;
  char p = 'V';
  // set some parameters
  params.set("T", 0.95);
  params.set("Pi", 3.14);
  params.set("A", "test");
  //params.set("C", 'A');
  params.set("C", p);

  // get double
  /* double t = */ params.get<double>("T");

  // get string
  /* std::string a = */ params.get<std::string>("A");

  // get char
  /* char c = */ params.get<char>("C");

  // is "T" defined?
  /* bool t_defined = */ params.defined("T");
  
  // is "S" defined?
  /* bool s_defined = */ params.defined("S");

  // clear paramter "Pi"
  params.clear("Pi");
  
  // output list of parameters as string
  std::list<std::string> keys = params.keys();
  std::cout << "number of parameters = " << keys.size() << std::endl;
  std::cout << "key-value pairs:\n";
  for(auto const& key : keys) {
    std::cout << "  " << key << " = " << params.get_string(key) << std::endl;
  }
}
