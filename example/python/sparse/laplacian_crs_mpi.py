#
# Rokko: Integrated Interface for libraries of eigenvalue decomposition
#
# Copyright (C) 2015-2020 by Rokko Developers https://github.com/t-sakashita/rokko
#
# Distributed under the Boost Software License, Version 1.0. (See accompanying
# file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
#

import mpi4py.MPI
from pyrokko import *

dim = 100
solver_name = "anasazi"

if (mpi4py.MPI.COMM_WORLD.Get_rank() == 0):
	print("solver name = {}".format(solver_name))
	print("dimension = {}".format(dim))

solver = parallel_sparse_ev(solver_name)
map = solver.default_mapping(dim, mpi4py.MPI.COMM_WORLD)
mat = distributed_crs_matrix(map, 3)

print("row_start={}, row_end={}".format(mat.start_row, mat.start_row))

if mat.start_row == 0:
    mat.insert(0, [0, 1], [1., -1.])

for row in range(max(1,mat.start_row), min(mat.end_row,dim-1)):
    mat.insert(row, [row-1, row, row+1], [-1., 2., -1.]);

if mat.end_row == dim:
    mat.insert(dim-1, [dim-2, dim-1], [-1., 2.]);

mat.complete()
mat.print()

params = parameters()
params.set("Block Size", 5);
params.set("Maximum Iterations", 500);
params.set("Convergence Tolerance", 1.0e-8);
params.set("num_eigenvalues", 10)

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
