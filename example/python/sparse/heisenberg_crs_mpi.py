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

L = 10
dim = 1 << L
solver_name = "anasazi"

lattice = [(i, (i+1) % L) for i in range(0, L)]

if (mpi4py.MPI.COMM_WORLD.Get_rank() == 0):
	print("solver name = {}".format(solver_name))
	print("Eigenvalue decomposition of antiferromagnetic Heisenberg chain")
	print("L = {}".format(L))
	print("dimension = {}".format(dim))

solver = parallel_sparse_ev(solver_name)

mat = distributed_crs_matrix(dim, dim, solver)

print("row_start={}, row_end={}".format(mat.start_row, mat.start_row))

for row in range(mat.start_row, mat.end_row):
    cols = []
    values = []
    diag = 0.
    for l in range(0, L):
        i = lattice[l][0]
        j = lattice[l][1]
        m1 = 1 << i
        m2 = 1 << j
        m3 = m1 + m2
        if ((row & m3) == m1) or ((row & m3) == m2):
            cols.append(row ^ m3)
            values.append(0.5)
            diag += -0.25
        else:
            diag += 0.25
    if diag != 0.:
        cols.append(row)
        values.append(diag)
    mat.insert(row, cols, values)

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
