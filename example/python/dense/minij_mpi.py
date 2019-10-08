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
g = pyrokko.grid()

solver = pyrokko.parallel_dense_ev()
map = solver.default_mapping(dim, g)
mat = pyrokko.distributed_matrix(map)
pyrokko.minij_matrix.generate(mat)

eigval = numpy.ndarray(dim)
eigvec = pyrokko.distributed_matrix(map)
params = pyrokko.parameters()

solver.initialize()
solver.diagonalize(mat, eigval, eigvec, params)
eigvec_loc = numpy.ndarray((dim,dim), order='F')
pyrokko.gather(eigvec, eigvec_loc, 0)

if (mpi4py.MPI.COMM_WORLD.Get_rank() == 0):
    print("eigenvalues:")
    print(eigval)
    print("eigenvectors:")
    print(eigvec_loc)

solver.finalize()
mpi4py.MPI.Finalize()
