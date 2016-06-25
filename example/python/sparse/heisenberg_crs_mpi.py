#
# Rokko: Integrated Interface for libraries of eigenvalue decomposition
#
# Copyright (C) 2015-2016 by Rokko Developers https://github.com/t-sakashita/rokko
#
# Distributed under the Boost Software License, Version 1.0. (See accompanying
# file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
#

from mpi4py import MPI
from rokko import *

L = 8
dim = 1 << L
solver_name = "anasazi"

lattice_first  = [x           for x in range(0, L)]
lattice_second = [(x + 1) % L for x in range(0, L)]

if (MPI.COMM_WORLD.Get_rank() == 0):
	print("solver name = %s" % solver_name)
	print("Eigenvalue decomposition of antiferromagnetic Heisenberg chain")
	print("L = %d" % L)
	print("dimension = %d" % dim)

solver = rokko_parallel_sparse_ev(solver_name, 0, None)

mat = rokko_distributed_crs_matrix(dim, dim, solver)

Z = rokko_distributed_crs_matrix(dim, dim, solver)

row_start = mat.start_row()
row_end = mat.end_row()
print("row_start=%d row_end=%d" % (row_start, row_end))

cols =   [0   for x in range(0, dim)]
values = [0.0 for x in range(0, dim)]

for row in range(row_start, row_end):
	count = 0
	diag = 0
	for l in range(0, L):
		i = lattice_first[l]
		j = lattice_second[l]
		m1 = 1 << i
		m2 = 1 << j
		m3 = m1 + m2
		if (row & m3) == m1 or (row & m3) == m2:
			cols[count] = row ^ m3
			values[count] = 0.5
			count += 1
			diag += -0.25
		else:
			diag += 0.25

	cols[count] = row
	values[count] = diag
	count += 1
	mat.insert(row, count, cols, values)

mat.complete()

params = rokko_parameters()
#params.set("routine", routine);
params.set("Block Size", 5);
params.set("Maximum Iterations", 500);
params.set("Convergence Tolerance", 1.0e-8);
params.set("num_eigenvalues", 10)

params_out = solver.diagonalize_distributed_crs_matrix(mat, params)
print("num_conv=%d" % params_out.get("num_conv"))
keys = params_out.keys()
num_conv = solver.num_conv()

eig_val = solver.eigenvalue(0)
num_local_rows = mat.num_local_rows()
print("num_local_rows=%d" % num_local_rows)

eig_vec = [0.0 for x in range(0, num_local_rows)]

solver.eigenvector(0, eig_vec)

if (MPI.COMM_WORLD.Get_rank() == 0):
	print("number of converged eigenpairs = %d" % num_conv)
	print("Computed Eigenvalue = ")
	print("%30.20f" % eig_val)
	print("Computed Eigenvector =")
	for i in range(0, num_local_rows):
		print("%30.20f " % eig_vec[i]),

#print "rank = ", MPI.COMM_WORLD.Get_rank()

#MPI.Finalize()
