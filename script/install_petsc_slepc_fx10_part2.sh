SCRIPT_DIR=${0%/*}
INSTALL_DIR=/work/n0004/n000402/rokko_lib

cd petsc-3.4.2

# When installing petsc, PETSC_DIR must be an current dir name
export PETSC_DIR=$PWD
export PETSC_ARCH=arch-fujitsu-sparc64fx-opt

#export PETSC_ARCH=  #arch-installed-petsc

python ./reconfigure-arch-fujitsu-sparc64fx-opt.py

make all test
make install

cd ..


#$INSTALL_DIR

# finished to install PETSc

# start to install SLEPc
##wget http://www.grycap.upv.es/slepc/download/download.php?filename=slepc-3.4.2.tar.gz
##tar xvf slepc-3.4.2.tar.gz
cd slepc-3.4.2

export PETSC_DIR=$INSTALL_DIR
export SLEPC_DIR=$PWD
export SLEPC_ARCH=  #"arch-installed-petsc"
#arch-fujitsu-sparc64fx-opt

./configure
make
#make test
