#!/bin/bash
# Reference: http://www.openpetascale.org/source/PETSc%20Installation%20Guide%20(June%202013).pdf

SCRIPT_DIR=${0%/*}
SOURCE_DIR=$WORK/source
INSTALL_DIR=$WORK/rokko_lib/

#wget -O - http://ftp.mcs.anl.gov/pub/petsc/release-snapshots/petsc-3.4.2.tar.gz | tar zxf -

tar xvf $SOURCE_DIR/petsc-3.4.2.tar.gz

cd $WORK/build/petsc-3.4.2
#wget -O - http://www.openpetascale.org/source/petsc-3.4.1_fujitsu_patches.tar.gz | tar zxf -
tar xvf petsc-3.4.1_fujitsu_patches.tar.gz

patch -p0 < petsc-3.4.1_setCompilers.patch

export PETSC_DIR=$PWD
export PETSC_ARCH=arch-fujitsu-sparc64fx-opt
echo petsc_dir=$PETSC_DIR

./configure \
--prefix=$INSTALL_DIR \
--with-cc="mpifccpx" --CFLAGS="-mt -Xg" --COPTFLAGS="-Kfast" \
--with-cxx="mpiFCCpx" --CXXFLAGS="-mt" --CXXOPTFLAGS="-Kfast" \
--with-fc="mpifrtpx" --FFLAGS="-Kthreadsafe" --FOPTFLAGS="-Kfast" \
--LDFLAGS="-lmpi_f77 -lmpi_f90" \
--with-blas-lapack-lib="-SSL2" \
--with-x=0 --bwith-c++-support --with-info=1 --with-debugging=0 --known-mpi-shared-libraries=0 --with-valgrind=0 \
--with-batch=1 --with-info=1

pjsub $SCRIPT_DIR/batch_conf.sh

#cd ..


