#
# Rokko: Integrated Interface for libraries of eigenvalue decomposition
#
# Copyright (C) 2015-2019 by Rokko Developers https://github.com/t-sakashita/rokko
#
# Distributed under the Boost Software License, Version 1.0. (See accompanying
# file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
#

from pyrokko import *
import numpy as np

dim = 2
mat = localized_matrix(dim, dim, matrix_major.row)
assert(mat.major == "row")

# substitute row-major matrix to row-major localized_matrix
mat_r = np.arange(dim**2, dtype='float').reshape(dim, dim)
mat.ndarray = mat_r
assert(mat.ndarray.all() == mat_r.all())

# substitute col-major matrix to row-major localized_matrix
mat_c = np.arange(dim**2, dtype='float').reshape(dim, dim, order='F')
mat.ndarray = mat_c
assert(mat.ndarray.all() == mat_c.all())

mat.ndarray = mat_r.T
assert(mat.ndarray.all() == mat_c.all())
