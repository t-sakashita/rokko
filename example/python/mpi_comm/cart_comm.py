#
# Rokko: Integrated Interface for libraries of eigenvalue decomposition
#
# Copyright (C) 2015-2020 by Rokko Developers https://github.com/t-sakashita/rokko
#
# Distributed under the Boost Software License, Version 1.0. (See accompanying
# file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
#

from mpi4py import MPI

def find_square_root_like_divisor(n):
    import math
    for i in range(int(math.sqrt(n)), 1, -1):
        if ((n % i) == 0):  break
    return i

def get_square_like_dims(nprocs):
    dim0 = find_square_root_like_divisor(nprocs)
    return [dim0, nprocs // dim0]


def create(comm_in = MPI.COMM_WORLD):
    dims = get_square_like_dims(comm_in.size)
    periods = [False, False]

    return comm_in.Create_cart(dims, periods=periods)


if __name__ == "__main__":
    comm = create()
    print (comm.rank)
