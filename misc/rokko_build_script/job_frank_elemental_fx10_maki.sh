#!/bin/bash
#------ pjsub option --------#
#PJM -L "rscgrp=debug"
#PJM -L "node=4"
#PJM --mpi "proc=4"
#PJM -L "elapse=10:00"
###PJM -o "${HOME}/sample_elemental.o"
#PJM -j

#------- Program execution -------#
export OMP_NUM_THREADS=8
cd $HOME/build/rokko/sample_elemental

mpiexec ./diagonalize_elemental

