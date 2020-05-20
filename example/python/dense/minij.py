#
# Rokko: Integrated Interface for libraries of eigenvalue decomposition
#
# Copyright (C) 2015-2019 by Rokko Developers https://github.com/t-sakashita/rokko
#
# Distributed under the Boost Software License, Version 1.0. (See accompanying
# file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
#

from pyrokko import *
import numpy

dim = 5

mat = numpy.ndarray((dim, dim))#, order='F')
eigval = numpy.ndarray(dim)
eigvec = numpy.ndarray((dim, dim))#, order='F')

minij_matrix.generate(mat)
print(mat)

params = parameters()
params.set("routine", "syev")
params.set("verbose", True)
solver = serial_dense_ev("lapack")
solver.diagonalize(mat, eigval, eigvec, params)
#solver.diagonalize(mat, eigval, params)

print("Computed eigenvalues:")
print(eigval)
print("eigenvectors:")
print(eigvec)
