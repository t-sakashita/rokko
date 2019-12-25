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

mat_loc = np.arange(dim**2, dtype='float').reshape((dim,dim), order='F')
#mat_loc = np.arange(dim**2, dtype='float').reshape((dim,dim), order='C')

g = pyrokko.grid(pyrokko.grid_row_major)
map = pyrokko.mapping_bc(dim, 2, g, pyrokko.matrix_major.col)
#map = pyrokko.mapping_bc(dim, 2, g, pyrokko.matrix_major.row)
mat = pyrokko.distributed_matrix(map)

pyrokko.scatter(mat_loc, mat, 0)

if g.myrank == 0:
    print("mat_loc=")
    print(mat_loc)
    print("mat=")
mat.print()

mpi4py.MPI.Finalize()
