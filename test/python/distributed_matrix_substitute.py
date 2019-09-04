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

dim = 4
g = pyrokko.grid()
#map = pyrokko.mapping_bc(dim, 2, g)
map = pyrokko.mapping_bc(dim, 2, g, pyrokko.matrix_major.row)
mat = pyrokko.distributed_matrix(map)
mat.set_local(1, 0, 5.)
mat.print()

A = mat.ndarray

print(A)
print(A.flags)
print(A.strides)
A[0,0] = 9
A[0,1] = 2
print(A)
mat.print()


B = numpy.full_like(A, 7.)
mat.ndarray = B
mat.print()

mpi4py.MPI.Finalize()
