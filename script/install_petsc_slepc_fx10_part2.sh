cd petsc-3.3-p6

#./reconfigure-arch-fujitsu-sparc64fx-opt.py

#export PETSC_ARCH=arch-fujitsu-sparc64fx-opt

export PETSC_ARCH=                                                                                                                                               
export PETSC_DIR=$PWD


#make all test                                                                                                                                                                             
#make install                                                                                                                                                                              

cd ..

#export PETSC_DIR=../install_build

# PETScのインストールおわり

# SLEPcのインストール
##wget http://www.grycap.upv.es/slepc/download/download.php?filename=slepc-3.3-p3.tar.gz
##tar xvf slepc-3.3-p3.tar.gz
cd slepc-3.3-p3
export SLEPC_DIR=$PWD
#export SLEPC_ARCH=arch-fujitsu-sparc64fx-opt

./configure
make
#make test 
