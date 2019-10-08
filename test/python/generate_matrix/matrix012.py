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

mat_c = numpy.ndarray((dim, dim), order='F')
matrix012.generate(mat_c)
print(mat_c)
