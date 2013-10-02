#!/bin/bash -x

PREFIX="$1"
test -z "$PREFIX" && PREFIX=$HOME/opt/rokko
echo "PREFIX = $PREFIX"
DEBUG="$2"
test -z "$DEBUG" && DEBUG=0
echo "DEBUG = $DEGUG"
dir=`dirname $0`
SCRIPT_DIR=`cd $dir && pwd`

mkdir -p $HOME/build
cd $HOME/build
rm -rf elemental-0.80 elemental-0.80-build
if test -f $HOME/source/elemental-0.80.tgz; then
  tar zxf $HOME/source/elemental-0.80.tgz
else
  wget -O - http://elemental.googlecode.com/files/elemental-0.80.tgz | tar zxf -
fi
cd $HOME/build/elemental-0.80
patch -p1 < $SCRIPT_DIR/elemental-0.80.patch

cd $HOME/build
mkdir elemental-0.80-build && cd elemental-0.80-build
if [ "$DEBUG" = "0" ]; then
    OPTFLAGS="-O3 -xSSE3"
else
    OPTFLAGS="-O0 -g"
fi
if [ `which mpicxx > /dev/null 2>&1; echo $?` = 0 ]; then
  cmake -DCMAKE_INSTALL_PREFIX=$PREFIX \
      -DSHARED_LIBRARIES=ON \
      -DCMAKE_CXX_COMPILER=mpicxx -DCMAKE_C_COMPILER=mpicc -DCMAKE_Fortran_COMPILER=mpif90 \
      -DCMAKE_CXX_FLAGS="$OPTFLAGS" -DCMAKE_C_FLAGS="$OPTFLAGS" -DCMAKE_Fortran_FLAGS="$OPTFLAGS" \
      -DMATH_LIBS="-mkl=parallel;-lifcore" \
      $HOME/build/elemental-0.80
else
  cmake -DCMAKE_INSTALL_PREFIX=$PREFIX \
      -DSHARED_LIBRARIES=ON \
      -DCMAKE_CXX_COMPILER=icpc -DCMAKE_C_COMPILER=icc -DCMAKE_Fortran_COMPILER=ifort \
      -DCMAKE_CXX_FLAGS="$OPTFLAGS" -DCMAKE_C_FLAGS="$OPTFLAGS" -DCMAKE_Fortran_FLAGS="$OPTFLAGS" \
      -DMPI_C_INCLUDE_PATH="/usr/include" -DMPI_CXX_INCLUDE_PATH="/usr/include" \
      -DMPI_C_LIBRARIES="-lmpi" -DMPI_CXX_LIBRARIES="-lmpi++;-lmpi" -DMPI_Fortran_LIBRARIES="-lmpi" \
      -DMATH_LIBS="-mkl=parallel;-lifcore" \
      $HOME/build/elemental-0.80
fi
make -j4 all VERBOSE=1
make install
