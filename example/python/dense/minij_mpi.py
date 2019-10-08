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

dim = 10
g = pyrokko.grid()

solver = pyrokko.parallel_dense_ev()
map = solver.default_mapping(dim, g)
mat = pyrokko.distributed_matrix(map)
pyrokko.minij_matrix.generate(mat)

eigval = pyrokko.localized_vector(dim);
eigvec = pyrokko.distributed_matrix(map)
params = pyrokko.parameters()

solver.initialize()
solver.diagonalize(mat, eigval, eigvec, params)
eigvec_loc = pyrokko.localized_matrix(dim,dim)
pyrokko.gather(eigvec, eigvec_loc, 0)

if (mpi4py.MPI.COMM_WORLD.Get_rank() == 0):
    print("eigenvalues:")
    eigval.print()
    print("eigenvectors:")
    eigvec_loc.print()

solver.finalize()
mpi4py.MPI.Finalize()