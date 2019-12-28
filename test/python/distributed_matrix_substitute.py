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
import numpy as np

dim = 8

g = pyrokko.grid(pyrokko.grid_row_major)
map = pyrokko.mapping_bc(dim, 2, 10, g, pyrokko.matrix_major.row)  # where lld=10
mat = pyrokko.distributed_matrix(map)
flattened_size = map.local_shape[0] * map.local_shape[1]

# substitute row-major matrix to row-major distributed_matrix
mat_r = np.arange(flattened_size, dtype='float').reshape(map.local_shape)
mat.ndarray = mat_r
assert((mat.ndarray == mat_r).all())
mat.print()
if(g.myrank == 0): print("")

# substitute col-major matrix to row-major distributed_matrix
mat_c = np.arange(flattened_size, dtype='float').reshape(map.local_shape, order='F')
mat.ndarray = mat_c
assert((mat.ndarray == mat_c).all())
mpi4py.MPI.COMM_WORLD.barrier()
mat.print();
if(g.myrank == 0): print("")

mat.ndarray = mat_r.T
assert((mat.ndarray == mat_c).all())
mpi4py.MPI.COMM_WORLD.barrier()
mat.print()

mpi4py.MPI.Finalize()
