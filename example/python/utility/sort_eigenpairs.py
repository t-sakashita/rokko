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
import random

dim = 5

eigval = numpy.arange(dim, dtype='float')
random.shuffle(eigval)
eigvec = numpy.eye(dim)#, order='F')

eigval_sorted = numpy.ndarray(dim)
eigvec_sorted = numpy.ndarray((dim, dim))#, order='F')

sort_eigenpairs(eigval, eigvec, eigval_sorted, eigvec_sorted, True)
print("eigval=", eigval)
print("eigval_sorted=", eigval_sorted)
print("eigvec_sorted:")
print(eigvec_sorted)

