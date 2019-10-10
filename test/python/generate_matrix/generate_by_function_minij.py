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

f = lambda i,j : min(i,j) + 1

dim = 5
mat = numpy.ndarray((dim, dim))
generate(mat, f)
print(mat)

mat_ref = numpy.ndarray((dim, dim))
minij_matrix.generate(mat_ref)

assert((mat == mat_ref).all())
