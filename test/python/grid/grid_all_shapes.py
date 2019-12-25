#
# Rokko: Integrated Interface for libraries of eigenvalue decomposition
#
# Copyright (C) 2015-2019 by Rokko Developers https://github.com/t-sakashita/rokko
#
# Distributed under the Boost Software License, Version 1.0. (See accompanying
# file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
#

import mpi4py
import mpi4py.MPI
import pyrokko
import numpy

nprocs = mpi4py.MPI.COMM_WORLD.size
dict_major = {pyrokko.grid_row_major:'C', pyrokko.grid_col_major:'F'}

for proc in range(1,nprocs+1):
    if (nprocs % proc == 0):
        nprow = nprocs // proc
        npcol = nprocs // nprow
        for major in [pyrokko.grid_row_major, pyrokko.grid_col_major]:
            grid_mat = numpy.arange(nprow*npcol, dtype='int').reshape(nprow, npcol, order=dict_major[major])
            g = pyrokko.grid(mpi4py.MPI.COMM_WORLD, (nprow,npcol), major)
            assert( g.myrank == grid_mat[g.mine] )

mpi4py.MPI.Finalize()
