#
# Rokko: Integrated Interface for libraries of eigenvalue decomposition
#
# Copyright (C) 2015-2019 by Rokko Developers https://github.com/t-sakashita/rokko
#
# Distributed under the Boost Software License, Version 1.0. (See accompanying
# file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
#

import mpi4py
import numpy
import pyrokko

# This program must be run with 4 MPI processes.
grid_dim = 2

# row-major
g = pyrokko.grid(pyrokko.grid_row_major)
assert( g.nprocs == 4 )
grid_mat = numpy.arange(grid_dim**2, dtype='int').reshape(grid_dim, grid_dim)
assert( g.myrank == grid_mat[g.mine] )

# Defalut: row-major
g = pyrokko.grid()
assert( g.myrank == grid_mat[g.mine] )

# col-major
g = pyrokko.grid(pyrokko.grid_col_major)
grid_mat = numpy.arange(grid_dim**2, dtype='int').reshape(grid_dim, grid_dim, order='F')
assert( g.myrank == grid_mat[g.mine] )

mpi4py.MPI.Finalize()
