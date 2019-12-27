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

myrank = mpi4py.MPI.COMM_WORLD.rank

dim = 10

g = pyrokko.grid(pyrokko.grid_row_major)
map = pyrokko.mapping_bc(dim, 2, g, pyrokko.matrix_major.col)
matA = pyrokko.distributed_matrix(map)
matC = pyrokko.distributed_matrix(map)
pyrokko.minij_matrix.generate(matA)

pyrokko.product(1., matA, False, matA, False, 0., matC);

lmatC = numpy.ndarray((dim, dim), order='F')
pyrokko.gather(matC, lmatC, 0);

if myrank == 0:
    mat_ref = numpy.ndarray((dim, dim), order='F')
    pyrokko.minij_matrix.generate(mat_ref)
    assert((lmatC == numpy.dot(mat_ref, mat_ref)).all())
