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

mat = numpy.ndarray((dim, dim))
matrix012.generate(mat)
print(mat)
mat_ref = numpy.arange(dim**2, dtype='float').reshape(dim, dim)
assert(mat.all() == mat_ref.all())

mat_col = numpy.ndarray((dim, dim), order='F')
matrix012.generate(mat_col)
print(mat_col)
mat_col_ref = numpy.arange(dim**2, dtype='float').reshape(dim, dim, order='F')
assert(mat_col.all() == mat_col_ref.all())
