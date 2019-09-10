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

dim = 4

g = pyrokko.grid(pyrokko.grid_row_major)
map = pyrokko.mapping_bc(dim, 2, g, pyrokko.matrix_major.col)
mat = pyrokko.distributed_matrix(map)
pyrokko.matrix012.generate(mat)
mat_loc = pyrokko.localized_matrix(dim, dim, pyrokko.matrix_major.col)

pyrokko.gather(mat, mat_loc, 0)

if (mpi4py.MPI.COMM_WORLD.Get_rank() == 0):
    mat_ref = np.arange(dim**2, dtype='float').reshape(mat_loc.local_shape)
    assert(mat_loc.ndarray.all() == mat_ref.all())

mpi4py.MPI.Finalize()
