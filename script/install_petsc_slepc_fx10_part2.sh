#!/bin/bash

# run this script after install_petsc_slepc_fx10_part1.sh
# assuming arch-fujitsu-sparc64fx-opt.py was made by running the batch job
SCRIPT_DIR=${0%/*}
INSTALL_DIR=$WORK/rokko_lib
BUILD_DIR=$WORK/build

cd $BUILD_DIR/petsc-3.4.2

# When installing petsc, PETSC_DIR must be an current dir name
export PETSC_DIR=$PWD
export PETSC_ARCH=arch-fujitsu-sparc64fx-opt

python ./reconfigure-arch-fujitsu-sparc64fx-opt.py

make all test
make install

# finished to install PETSc

# Now, let us install SLEPc
#wget http://www.grycap.upv.es/slepc/download/download.php?filename=slepc-3.4.2.tar.gz
tar xvf slepc-3.4.2.tar.gz
cd $BUILD_DIR/slepc-3.4.2

export SLEPC_DIR=$PWD
export PETSC_DIR=$INSTALL_DIR
# we need to reset PETSC_ARCH to specify a install prefix before configuring SLEPc.
unset PETSC_ARCH
./configure --prefix=$INSTALL_DIR

export PETSC_ARCH=arch-installed-petsc
make
make test
make install
