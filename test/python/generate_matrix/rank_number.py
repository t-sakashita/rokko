#
# Rokko: Integrated Interface for libraries of eigenvalue decomposition
#
# Copyright (C) 2015-2019 by Rokko Developers https://github.com/t-sakashita/rokko
#
# Distributed under the Boost Software License, Version 1.0. (See accompanying
# file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
#

import mpi4py
import mpi4py.MPI
import pyrokko
import numpy

def test(global_size, block_size, grid_major):
    myrank = mpi4py.MPI.COMM_WORLD.rank
    g = pyrokko.grid(grid_major)
    grid_size = g.shape

    map = pyrokko.mapping_bc(global_size, block_size, g, pyrokko.matrix_major.col)
    mat = pyrokko.distributed_matrix(map)
    for local_i in range(map.local_shape[0]):
        for local_j in range(map.local_shape[1]):
            mat.set_local(local_i, local_j, myrank)

    dim_proc = global_size if g.myrank == 0 else (2,2)  # Want to fix 2 to 0
    mat_loc = numpy.ndarray(dim_proc, order='F')
    pyrokko.gather(mat, mat_loc, 0)

    if myrank == 0:
        for i in range(0, global_size[0], block_size[0]):
            i_block = i // block_size[0]
            ip = i_block % grid_size[0]
            block0 = min(global_size[0] - i, block_size[0])
            for j in range(0, global_size[1], block_size[1]):
                j_block = j // block_size[1]
                jp = j_block % grid_size[1]
                block1 = min(global_size[1] - j, block_size[1])
                rank = g.calculate_rank_form_coords(ip, jp)
                assert( (mat_loc[i:i+block0, j:j+block1] == rank).all() )


def test_both_grid_major(global_size, block_size):
    for grid_major in [pyrokko.grid_row_major, pyrokko.grid_col_major]:
       test(global_size, block_size, grid_major)

## main program
test_both_grid_major((20,8), (4,3))
test_both_grid_major((8,20), (4,3))
test_both_grid_major((5,7), (9,12))
test_both_grid_major((7,9), (9,5))

mpi4py.MPI.Finalize()
