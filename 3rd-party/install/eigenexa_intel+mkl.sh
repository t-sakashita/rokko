#!/bin/bash -x

PREFIX="$1"
test -z "$PREFIX" && PREFIX=$HOME/opt/rokko
echo "PREFIX = $PREFIX"
DEBUG="$2"
test -z "$DEBUG" && DEBUG=0
echo "DEBUG = $DEBUG"
dir=`dirname $0`
SCRIPT_DIR=`cd $dir && pwd`

mkdir -p $HOME/build
cd $HOME/build
rm -rf EigenExa-1.3 EigenExa-1.3-build
if test -f $HOME/source/EigenExa-1.3.tgz; then
  tar zxf $HOME/source/EigenExa-1.3.tgz
else
  wget -O - http://www.aics.riken.jp/labs/lpnctrt/EigenExa-1.3.tgz | tar zxf -
fi
cd $HOME/build/EigenExa-1.3
patch -p1 < $SCRIPT_DIR/EigenExa-1.3.patch

cd $HOME/build
mkdir -p EigenExa-1.3-build && cd EigenExa-1.3-build
if [ "$DEBUG" = "0" ]; then
    OPTFLAGS="-O3 -xSSE3"
else
    OPTFLAGS="-O0 -g"
fi
if [ `which mpicxx > /dev/null 2>&1; echo $?` = 0 ]; then
    cmake -DCMAKE_INSTALL_PREFIX=$PREFIX \
	-DCMAKE_C_COMPILER=mpicc -DCMAKE_Fortran_COMPILER=mpif90 \
	-DCMAKE_C_FLAGS="$OPTFLAGS" -DCMAKE_Fortran_FLAGS="$OPTFLAGS" \
	-DSCALAPACK_LIB="-lmkl_scalapack_lp64 -lmkl_blacs_intelmpi_lp64 -mkl=parallel" \
	$HOME/build/EigenExa-1.3
else
    cmake -DCMAKE_INSTALL_PREFIX=$PREFIX \
	-DCMAKE_C_COMPILER=icc -DCMAKE_Fortran_COMPILER=ifort \
	-DCMAKE_C_FLAGS="$OPTFLAGS" -DCMAKE_Fortran_FLAGS="$OPTFLAGS" \
	-DMPI_C_INCLUDE_PATH="/usr/include" \
	-DMPI_C_LIBRARIES="-lmpi" -DMPI_Fortran_LIBRARIES="-lmpi" \
	-DSCALAPACK_LIB="-lmkl_scalapack_lp64 -lmkl_blacs_intelmpi_lp64 -mkl=parallel" \
	$HOME/build/EigenExa-1.3
fi
make -j4 VERBOSE=1
make install
