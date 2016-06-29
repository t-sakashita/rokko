#
# Rokko: Integrated Interface for libraries of eigenvalue decomposition
#
# Copyright (C) 2015-2016 by Rokko Developers https://github.com/t-sakashita/rokko
#
# Distributed under the Boost Software License, Version 1.0. (See accompanying
# file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
#

from rokko import *

solver = rokko_serial_dense_ev("lapack", 0, None)

dim = 10
rokko_matrix_col_major = rokko.matrix_col_major

mat = rokko_localized_matrix(dim, dim, rokko_matrix_col_major)
Z = rokko_localized_matrix(dim, dim, rokko_matrix_col_major)
w = rokko_localized_vector(dim) 

rokko_generate_frank_matrix(mat)
mat.show()

params = rokko_parameters()
params.set("routine", "dsyevr")
b = True
params.set("verbose", 1)
solver.diagonalize(mat, w, Z, params)

print("Computed Eigenvalues =\n");

for i in range(0, dim):
	print(w.get(i))