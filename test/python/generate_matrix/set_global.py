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

dim = 8

g = pyrokko.grid(pyrokko.grid_row_major)
map = pyrokko.mapping_bc(dim, 2, 10, g, pyrokko.matrix_major.col)  # where lld=10
mat = pyrokko.distributed_matrix(map)
for i in range(map.global_shape[0]):
    for j in range(map.global_shape[1]):
        mat.set_global(i, j, min(i, j) + 1)

dim_proc = dim if g.myrank == 0 else 2  # Want to fix 2 to 0
mat_loc = numpy.ndarray((dim_proc, dim_proc), order='F')
pyrokko.gather(mat, mat_loc, 0)

if g.myrank == 0:
    mat_ref = numpy.ndarray((dim, dim), order='F')
    pyrokko.minij_matrix.generate(mat_ref)
    assert(mat_loc.all() == mat_ref.all())

mpi4py.MPI.Finalize()
