#
# Rokko: Integrated Interface for libraries of eigenvalue decomposition
#
# Copyright (C) 2015-2020 by Rokko Developers https://github.com/t-sakashita/rokko
#
# Distributed under the Boost Software License, Version 1.0. (See accompanying
# file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
#

import mpi4py.MPI
from pyrokko import *

class heisenberg_op(distributed_mfree):
    def __init__(self, L, lattice):
        self.L = L
        self.lattice = lattice
        self.comm = mpi4py.MPI.COMM_WORLD
        self.nprocs = self.comm.Get_size()
        self.myrank = self.comm.Get_rank()

        p = self.get_pow_of_2(self.nprocs)
        if self.nprocs != (1 << p):
            raise ValueError("This program can be run only for powers of 2")
        self.local_pow = self.L - p
        
        dim = 1 << L
        distributed_mfree.__init__(self, self.multiply, dim, self.comm)


    def get_pow_of_2(self, n):
        assert(n >= 1)
        p = 0
        while n > 1:
            n //= 2
            p += 1
        return p

    def multiply(self, x, y):
        local_pow = self.local_pow
        myrank = self.myrank
        
        for (i,j) in lattice:

            if i < local_pow:
                if j < local_pow:
                    m1 = 1 << i
                    m2 = 1 << j
                    m3 = m1 + m2

                    for k in range(self.num_local_rows):
                        if (k & m3) == m1:  # when (bit i == 1, bit j == 0) or (bit i == 0, bit j == 1)
                            y[k] += 0.5 * x[k^m3] - 0.25 * x[k]
                        elif (k & m3) == m2:
                            y[k] += 0.5 * x[k^m3] - 0.25 * x[k]
                        else:
                            y[k] += 0.25 * x[k]
                else:
                    m = 1 << (j - local_pow)
                    buffer = self.comm.sendrecv(x, dest=myrank^m, source=myrank^m)
                    m1 = 1 << i
                    if (myrank & m) == m:
                        for k in range(self.num_local_rows):
                            if (k & m1) == m1:
                                y[k] += 0.25 * x[k]
                            else:
                                y[k] += 0.5 * buffer[k^m1] - 0.25 * x[k]
                    else:
                        for k in range(self.num_local_rows):
                            if (k & m1) == m1:
                                y[k] += 0.5 * buffer[k^m1] - 0.25 * x[k]
                            else:
                                y[k] += 0.25 * x[k]
            else:
                if j < local_pow:
                    m = 1 << (i - local_pow)
                    buffer = self.comm.sendrecv(x, dest=myrank^m, source=myrank^m)
                    m1 = 1 << j
                    if (myrank & m) == m:
                        for k in range(self.num_local_rows):
                            if (k & m1) == m1:
                                y[k] += 0.25 * x[k]
                            else:
                                y[k] += 0.5 * buffer[k^m1] - 0.25 * x[k]
                    else:
                        for k in range(self.num_local_rows):
                            if (k & m1) == m1:
                                y[k] += 0.5 * buffer[k^m1] - 0.25 * x[k]
                            else:
                                y[k] += 0.25 * x[k]
                else:
                    m = (1 << (i - local_pow)) + (1 << (j - local_pow))
                    if ((myrank & m) != m) and ((myrank & m) != 0):
                        buffer = self.comm.sendrecv(x, dest=myrank^m, source=myrank^m)
                        for k in range(self.num_local_rows):
                            y[k] += 0.5 * buffer[k] - 0.25 * x[k]
                    else:
                        for k in range(self.num_local_rows):
                            y[k] += 0.25 * x[k]


# Main program
L = 10
lattice = [(i, (i+1) % L) for i in range(0, L)]
mat = heisenberg_op(L, lattice)

params = parameters()
params.set("Block Size", 5);
params.set("Maximum Iterations", 500);
params.set("Convergence Tolerance", 1.0e-8);
params.set("num_eigenvalues", 10)
params.set("routine", "LOBPCG");
params.set("verbose", True)

solver_name = "anasazi"
solver = parallel_sparse_ev(solver_name)
params_out = solver.diagonalize(mat, params)
num_conv = params_out.get("num_conv")

eig_val = solver.eigenvalue(0)
eig_vec = solver.eigenvector(0)

if (mpi4py.MPI.COMM_WORLD.Get_rank() == 0):
    print("number of converged eigenpairs = {}".format(num_conv))
    print("Computed Eigenvalue :")
    print("{:30.20f}".format(eig_val))
    print("Computed Eigenvector :")
    print(eig_vec)

solver.finalize()
mpi4py.MPI.Finalize()
