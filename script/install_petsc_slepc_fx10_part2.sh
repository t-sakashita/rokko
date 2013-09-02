SCRIPT_DIR=${0%/*}
INSTALL_DIR=/work/n0004/n000402/rokko_lib

cd petsc-3.4.2

# When installing petsc, PETSC_DIR must be an current dir name
export PETSC_DIR=$PWD
export PETSC_ARCH=arch-fujitsu-sparc64fx-opt

#python ./reconfigure-arch-fujitsu-sparc64fx-opt.py

#make all test
#make install

cd ..


#$INSTALL_DIR

# finished to install PETSc

# installing SLEPc
##wget http://www.grycap.upv.es/slepc/download/download.php?filename=slepc-3.4.2.tar.gz
tar xvf slepc-3.4.2.tar.gz
cd slepc-3.4.2

export SLEPC_DIR=$PWD
export PETSC_DIR=$INSTALL_DIR
unset PETSC_ARCH
#export PETSC_ARCH=
#arch-fujitsu-sparc64fx-opt

./configure --prefix=$INSTALL_DIR

export PETSC_ARCH=arch-installed-petsc
make
make test
make install
