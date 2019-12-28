#
# Rokko: Integrated Interface for libraries of eigenvalue decomposition
#
# Copyright (C) 2015-2019 by Rokko Developers https://github.com/t-sakashita/rokko
#
# Distributed under the Boost Software License, Version 1.0. (See accompanying
# file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
#

import mpi4py
import pyrokko
import numpy

dim = 10

g = pyrokko.grid(pyrokko.grid_row_major)
map = pyrokko.mapping_bc(dim, 2, g, pyrokko.matrix_major.col)
mat = pyrokko.distributed_matrix(map)
pyrokko.matrix012.generate(mat)

mat_loc = numpy.ndarray((dim, dim), order='F')
pyrokko.gather(mat, mat_loc, 0)

if (g.myrank == 0):
    mat_ref = numpy.ndarray((dim, dim), order='F')
    pyrokko.matrix012.generate(mat_ref)
    assert((mat_loc == mat_ref).all())

mpi4py.MPI.Finalize()
