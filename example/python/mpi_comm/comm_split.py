#
# Rokko: Integrated Interface for libraries of eigenvalue decomposition
#
# Copyright (C) 2015-2020 by Rokko Developers https://github.com/t-sakashita/rokko
#
# Distributed under the Boost Software License, Version 1.0. (See accompanying
# file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
#

from mpi4py import MPI

def create_even_odd_comm(comm_in = MPI.COMM_WORLD):
    world_rank = comm_in.rank

    color = world_rank % 2
    return MPI.COMM_WORLD.Split(color, world_rank)


if __name__ == "__main__":
    comm = create_even_odd_comm()

    for i in range(MPI.COMM_WORLD.size):
        if (MPI.COMM_WORLD.rank == i) :
            print ("MPI_COMM : rank {} of {}, new comm : rank {} of {}".format(MPI.COMM_WORLD.rank, MPI.COMM_WORLD.size, comm.rank, comm.size))
        MPI.COMM_WORLD.Barrier()

    comm.Free()
