#!/bin/bash -x

PREFIX="$1"
if test -z "$PREFIX"; then
    if test -d "/opt/nano/rokko"; then
	PREFIX=/opt/nano/rokko
    elif test -d "/opt/rokko"; then
	PREFIX=/opt/nano/rokko
    fi
fi
echo "PREFIX = $PREFIX"

SCRIPT_DIR=`dirname $0`
mkdir -p $HOME/build
mkdir -p $PREFIX/include $PREFIX/lib

# eigen_s

cd $HOME/build
rm -rf eigen_s
wget http://ccse.jaea.go.jp/ja/download/eigenk_files/eigen_s.zip
unzip eigen_s.zip
rm -f eigen_s.zip

cd eigen_s
make FF=mpif90 FC=mpif90 CC=mpicc CCFLAG="-O3 -openmp -fpp" SCALAPACK="-lmkl_scalapack_lp64 -lmkl_blacs_intelmpi_lp64 -mkl=parallel"
make FF=mpif90 FC=mpif90 CC=mpicc CCFLAG="-O3 -openmp -fpp" SCALAPACK="-lmkl_scalapack_lp64 -lmkl_blacs_intelmpi_lp64 -mkl=parallel" a.out

cp -p communication_s.mod $PREFIX/include
cp -p libEigen_s.a $PREFIX/lib

# eigen_sx

cd $HOME/build
rm -rf eigen_sx
wget http://ccse.jaea.go.jp/ja/download/eigenk_files/eigen_sx.zip 
unzip eigen_sx.zip
rm -rf eigen_sx.zip

cd eigen_sx
patch -p1 < $SCRIPT_DIR/eigen_sx.patch

make FF=mpif90 FC=mpif90 CC=mpicc CCFLAG="-O3 -openmp -fpp" SCALAPACK="-lmkl_scalapack_lp64 -lmkl_blacs_intelmpi_lp64 -mkl=parallel"
make FF=mpif90 FC=mpif90 CC=mpicc CCFLAG="-O3 -openmp -fpp" SCALAPACK="-lmkl_scalapack_lp64 -lmkl_blacs_intelmpi_lp64 -mkl=parallel" a.out

cp -p communication_sx.mod $PREFIX/include
cp -p libEigen_sx.a $PREFIX/lib
