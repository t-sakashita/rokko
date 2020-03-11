#
# Rokko: Integrated Interface for libraries of eigenvalue decomposition
#
# Copyright (C) 2015-2020 by Rokko Developers https://github.com/t-sakashita/rokko
#
# Distributed under the Boost Software License, Version 1.0. (See accompanying
# file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
#

from mpi4py import MPI

def create_even_comm(comm_in = MPI.COMM_WORLD):
    nprocs = MPI.COMM_WORLD.size

    Neven = (nprocs + 1) // 2
    even_members = [2 * i for i in range(Neven)]

    world_group = comm_in.Get_group()
    even_group = world_group.Incl(even_members)
    even_comm = comm_in.Create(even_group)
    world_group.Free()
    even_group.Free()
    return even_comm


def create_odd_comm(comm_in = MPI.COMM_WORLD):
    nprocs = MPI.COMM_WORLD.size

    Nodd = nprocs // 2
    odd_members = [2 * i + 1 for i in range(Nodd)]

    world_group = comm_in.Get_group()
    odd_group = world_group.Incl(odd_members)
    odd_comm = comm_in.Create(odd_group)
    world_group.Free()
    odd_group.Free()
    return odd_comm


def create_even_odd_comm(comm_in = MPI.COMM_WORLD):
    return create_even_comm(comm_in) if ((comm_in.rank % 2) == 0) else create_odd_comm(comm_in)


if __name__ == "__main__":
    comm = create_even_comm()
    if (comm != MPI.COMM_NULL):
        print (comm.rank)
