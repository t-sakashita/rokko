#!/bin/bash

##wget http://ftp.mcs.anl.gov/pub/petsc/release-snapshots/petsc-3.3-p6.tar.gz
##tar xvf petsc-3.3-p6.tar.gz
cd petsc-3.3-p6

export PETSC_DIR=$PWD
export PETSC_ARCH=arch-fujitsu-sparc64fx-opt
echo petsc_dir=$PETSC_DIR

./configure \
--prefix="../install_build/" \
--with-cc="mpifccpx" --CFLAGS="-mt -Xg" --COPTFLAGS="-Kfast" \
--with-cxx="mpiFCCpx" --CXXFLAGS="-mt" --CXXOPTFLAGS="-Kfast" \
--with-fc="mpifrtpx" --FFLAGS="-Kthreadsafe" --FOPTFLAGS="-Kfast" \
--with-blas-lapack-lib="-SSL2" \
--with-bx=0 --bwith-c++-support --with-info=1 --with-debugging=0 --known-mpi-shared-libraries=0 --with-valgrind=0 \
--with-batch

cd ..
pjsub batch_conf.sh


