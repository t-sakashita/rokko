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
for local_i in range(map.local_shape[0]):
    global_i = map.translate_l2g_row(local_i)
    for local_j in range(map.local_shape[1]):
        global_j = map.translate_l2g_col(local_j)
        mat.set_local(local_i, local_j, min(global_i, global_j) + 1)

dim_proc = dim if g.myrank == 0 else 2  # Want to fix 2 to 0
mat_loc = numpy.ndarray((dim_proc, dim_proc), order='F')
pyrokko.gather(mat, mat_loc, 0)

if g.myrank == 0:
    mat_ref = numpy.ndarray((dim, dim), order='F')
    pyrokko.minij_matrix.generate(mat_ref)
    assert((mat_loc == mat_ref).all())

mpi4py.MPI.Finalize()
