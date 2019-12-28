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

mat_loc = numpy.ndarray((dim, dim), order='F')
pyrokko.matrix012.generate(mat_loc)

pyrokko.scatter(mat_loc, mat, 0)

mat_ref = pyrokko.distributed_matrix(map)
pyrokko.matrix012.generate(mat_ref)

if g.myrank == 0:
    for local_i in range(map.local_shape[0]):
        for local_j in range(map.local_shape[1]):
            assert(mat.get_local(local_i, local_j) == mat_ref.get_local(local_i, local_j))

mpi4py.MPI.Finalize()
