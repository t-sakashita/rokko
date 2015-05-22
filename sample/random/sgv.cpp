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

#include <rokko/rokko.hpp>
#include <boost/random.hpp>

typedef rokko::localized_vector<double> vector_t;

int main(int argc, char *argv[]) {
  int n = 6;
  vector_t v(n);
  boost::mt19937 engine(12345l);
  boost::variate_generator<boost::mt19937&, boost::normal_distribution<> >
    gauss(engine, boost::normal_distribution<>());
  for (int i = 0; i < n; ++i) v(i) = gauss();
  std::cout << "v: " << v << std::endl;
}
