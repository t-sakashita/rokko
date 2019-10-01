#
# Rokko: Integrated Interface for libraries of eigenvalue decomposition
#
# Copyright (C) 2015-2019 by Rokko Developers https://github.com/t-sakashita/rokko
#
# Distributed under the Boost Software License, Version 1.0. (See accompanying
# file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
#

import mpi4py.MPI
from pyrokko import *
import heisenberg


# Main program
L = 10
lattice = [(i, (i+1) % L) for i in range(0, L)]
mat = heisenberg.mfree(L, lattice)

params = parameters()
params.set("Block Size", 5);
params.set("Maximum Iterations", 500);
params.set("Convergence Tolerance", 1.0e-8);
params.set("num_eigenvalues", 10)
params.set("routine", "LOBPCG");
params.set("verbose", True)

solver_name = "anasazi"
solver = parallel_sparse_ev(solver_name)
params_out = solver.diagonalize(mat, params)
num_conv = params_out.get("num_conv")

eig_val = solver.eigenvalue(0)
eig_vec = solver.eigenvector(0)

if (mpi4py.MPI.COMM_WORLD.Get_rank() == 0):
    print("number of converged eigenpairs = {}".format(num_conv))
    print("Computed Eigenvalue :")
    print("{:30.20f}".format(eig_val))
    print("Computed Eigenvector :")
    print(eig_vec)

solver.finalize()
mpi4py.MPI.Finalize()
