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
import numpy

class laplacian_op(distributed_mfree):
    def __init__(self, dim):
        self.comm = mpi4py.MPI.COMM_WORLD
        self.nprocs = self.comm.Get_size()
        self.myrank = self.comm.Get_rank()
        distributed_mfree.__init__(self, self.multiply, dim, self.comm)

        self.is_first_proc = self.start_row == 0
        self.is_last_proc = self.end_row == dim
        self.end_k = self.num_local_rows - 1

    def multiply(self, x, y):
        if self.num_local_rows == 0:
            return

        if (not self.is_first_proc) and (self.nprocs != 1):
            self.comm.send(x[0], dest=self.myrank-1)
            self.buf_m = self.comm.recv(source=self.myrank-1)

        if (not self.is_last_proc) and (self.nprocs != 1):
            self.buf_p = self.comm.recv(source=self.myrank+1)
            self.comm.send(x[self.end_k], dest=self.myrank+1)

        if self.is_first_proc:
            if self.num_local_rows != 1:
                y[0] = x[0] - x[1]
                if self.nprocs != 1:
                    y[self.end_k] = - x[self.end_k - 1] + 2 * x[self.end_k] - self.buf_p
            else:
                y[0] = x[0] - self.buf_p

        if self.is_last_proc:
            if self.num_local_rows != 1:
                if self.nprocs != 1:
                    y[0] = - self.buf_m + 2 * x[0] - x[1];
                y[self.end_k] = 2 * x[self.end_k] - x[self.end_k - 1]
            else:
                y[self.end_k] = 2 * x[self.end_k] - self.buf_m

        if not(self.is_first_proc or self.is_last_proc):  # neither first or last process
            if self.num_local_rows != 1:
                y[0] = - self.buf_m + 2 * x[0] - x[1]
                y[self.end_k] = - x[self.end_k - 1] + 2 * x[self.end_k] - self.buf_p
            else:
                y[0] = - self.buf_m + 2 * x[0] - self.buf_p

        # from 1 to end-1
        for k in range(1, self.end_k):
            y[k] = - x[k-1] + 2 * x[k] - x[k+1]


# Main program
dim = 100
mat = laplacian_op(dim)

solver_name = "anasazi"
solver = parallel_sparse_ev(solver_name)

params = parameters()
params.set("Block Size", 5);
params.set("Maximum Iterations", 500);
params.set("Convergence Tolerance", 1.0e-8);
params.set("num_eigenvalues", 10)
params.set("routine", "lobpcg");
params.set("verbose", True)

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
