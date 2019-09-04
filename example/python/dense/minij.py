#
# Rokko: Integrated Interface for libraries of eigenvalue decomposition
#
# Copyright (C) 2015-2019 by Rokko Developers https://github.com/t-sakashita/rokko
#
# Distributed under the Boost Software License, Version 1.0. (See accompanying
# file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
#

from pyrokko import *

dim = 10

mat = localized_matrix(dim, dim, matrix_major.col)
eigvec = localized_matrix(dim, dim, matrix_major.col)
eigval = localized_vector(dim)

minij_matrix.generate(mat)
mat.print()

params = parameters()
params.set("routine", "dsyevr")
params.set("verbose", True)
solver = serial_dense_ev("lapack")
solver.diagonalize(mat, eigval, eigvec, params)

print("Computed eigenvalues:")
eigval.print()
print("eigenvectors:")
eigvec.print()
