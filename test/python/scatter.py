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

mat_loc = pyrokko.localized_matrix(dim, dim, pyrokko.matrix_major.col)
mat_loc.ndarray = np.arange(dim**2, dtype='float').reshape(mat_loc.local_shape)

g = pyrokko.grid(pyrokko.grid_row_major)
map = pyrokko.mapping_bc(dim, 2, g, pyrokko.matrix_major.col)
mat = pyrokko.distributed_matrix(map)

pyrokko.scatter(mat_loc, mat, 0)

mat_ref = pyrokko.distributed_matrix(map)
pyrokko.matrix012.generate(mat_ref)

assert(mat.ndarray.all() == mat_ref.ndarray.all())

mpi4py.MPI.Finalize()
