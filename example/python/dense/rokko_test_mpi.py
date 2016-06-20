#
# Rokko: Integrated Interface for libraries of eigenvalue decomposition
#
# Copyright (C) 2015 by Rokko Developers https://github.com/t-sakashita/rokko
#
# Distributed under the Boost Software License, Version 1.0. (See accompanying
# file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
#

from mpi4py import MPI
from rokko import *

solver = rokko_parallel_dense_ev("scalapack", 0, None)
#solver = rokko_parallel_dense_ev("elemental", 0, None)

dim = 50

rokko_grid_row_major = rokko.grid_row_major
rokko_grid_col_major = rokko.grid_col_major
rokko_matrix_col_major = rokko.matrix_col_major

grid = rokko_grid(MPI.COMM_WORLD, rokko_grid_col_major)
#map = solver.default_mapping(dim, grid)
map = rokko_mapping_bc(dim, 1, grid)
mat = rokko_distributed_matrix(map)
Z = rokko_distributed_matrix(map)
w = rokko_localized_vector(dim)

rokko_frank_matrix_generate_distributed_matrix(mat)
mat.show()

solver.diagonalize_distributed_matrix(mat, w, Z)

if (MPI.COMM_WORLD.Get_rank() == 0):
	print("Computed Eigenvalues =\n");

	for i in range(0, dim):
		print(w.get(i))

# print "rank = ", MPI.COMM_WORLD.Get_rank()

MPI.Finalize()
