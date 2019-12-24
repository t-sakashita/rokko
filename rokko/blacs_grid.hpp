/*****************************************************************************
*
* Rokko: Integrated Interface for libraries of eigenvalue decomposition
*
* Copyright (C) 2012-2019 by Rokko Developers https://github.com/t-sakashita/rokko
*
* Distributed under the Boost Software License, Version 1.0. (See accompanying
* file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
*
*****************************************************************************/

#ifndef ROKKO_BLACS_GRID_HPP
#define ROKKO_BLACS_GRID_HPP

#include <mpi.h>
#include <array>
#include <rokko/blacs.hpp>

namespace rokko {

class blacs_grid {
public:
  blacs_grid() {}

  void set_blacs_grid(MPI_Comm comm, bool is_row, std::array<int,2> const& grid_size) {
    blacs_context = blacs::sys2blacs_handle(comm);
    char char_grid_major = is_row ? 'R' : 'C';
    blacs::gridinit(blacs_context, char_grid_major, grid_size[0], grid_size[1]);
  }

  int get_blacs_context() const { return blacs_context; }

private:
  int blacs_context;
};

} // namespace rokko

#endif // ROKKO_BLACS_GRID_HPP
