/*****************************************************************************
*
* Rokko: Integrated Interface for libraries of eigenvalue decomposition
*
* Copyright (C) 2014 by Synge Todo <wistaria@comp-phys.org>
*
* Distributed under the Boost Software License, Version 1.0. (See accompanying
* file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
*
*****************************************************************************/

#ifndef OPTIONS_HPP
#define OPTIONS_HPP

#include <iosfwd>

struct options {
  unsigned int N;
  std::string solver;
  bool valid;

  options(unsigned int argc, char *argv[], unsigned int N_def, std::string const& solver_def, bool print = true) :
    N(N_def), solver(solver_def), valid(true) {
    for (unsigned int i = 1; i < argc; ++i) {
      switch (argv[i][0]) {
      case '-' :
        switch (argv[i][1]) {
        case 'N' :
          if (++i == argc) { usage(); return; }
          N = std::stoi(argv[i]); break;
        case 's' :
          if (++i == argc) { usage(); return; }
          solver = argv[i]; break;
        case 'h' :
          usage(std::cout); return;
        default :
          usage(); return;
        }
        break;
      default :
        usage(); return;
      }
    }
    if (N <= 0) {
      if (print) {
        std::cerr << "invalid parameter\n"; usage();
      }
      return;
    }
    if (print) {
      std::cout << "System Size = " << N << '\n';
      std::cout << "Solver Name = " << solver << '\n';
    }
  }

  void usage(std::ostream& os = std::cerr) {
    os << "[command line options]\n"
       << "  -N int    System Size\n"
       << "  -s solver Solver Name\n"
       << "  -h        this help\n";
    valid = false;
  }
};

#endif
