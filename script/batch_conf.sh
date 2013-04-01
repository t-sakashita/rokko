#!/bin/bash -x
#PJM -L "rscunit=unit1"
#PJM -L "rscgrp=large"
#PJM --rsc-list "node=8"
#PJM --rsc-list "elapse=02:00:00"
#PJM --mpi "proc=8"
#PJM -s

#. /work/system/Env_base

MPI_HOME=/opt/FJSVfxlang/1.2.0
TOFU_LIB=/opt/FJSVpxtof/sparc64fx/lib64
export PATH="${MPI_HOME}/bin:${PATH}"
export LD_LIBRARY_PATH=${MPI_HOME}/lib64:${TOFU_LIB}:${LD_LIBRARY_PATH}
export THREAD_STACK_SIZE=65536
export FLIB_FASTOMP=TRUE

export OMP_NUM_THREADS=16
#export PARALLEL=2

cd /k/home/users/zs4/wget_try/petsc-3.3-p6

mpiexec ./conftest-arch-fujitsu-sparc64fx-opt


