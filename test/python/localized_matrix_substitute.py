#
# Rokko: Integrated Interface for libraries of eigenvalue decomposition
#
# Copyright (C) 2015-2019 by Rokko Developers https://github.com/t-sakashita/rokko
#
# Distributed under the Boost Software License, Version 1.0. (See accompanying
# file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
#

from pyrokko import *
from numpy import *

#mat = localized_matrix(2,2, matrix_major.row)
mat = localized_matrix(2,2, matrix_major.col)

A = mat.ndarray

B = full_like(A, 7.)
print(B)

mat.ndarray = B
mat.print()
