#!/bin/bash

SCRIPT_DIR=$(cd "$(dirname $0)"; pwd)
. $SCRIPT_DIR/../util.sh
. $SCRIPT_DIR/version.sh
set_prefix

sh $SCRIPT_DIR/setup.sh

MKL=$(echo $MKLROOT | cut -d: -f1)
unset PETSC_DIR
BUILD_TYPES="Release Debug"
for build_type in $BUILD_TYPES; do
  PREFIX=$PREFIX_ROKKO/petsc-$PETSC_VERSION-$PETSC_RK_REVISION/$build_type
  cd $BUILD_DIR
  cp -rp petsc-$PETSC_VERSION petsc-$PETSC_VERSION-build-$build_type
  cd petsc-$PETSC_VERSION-build-$build_type
  if [ $build_type = "Release" ]; then
    if [ `which mpicxx > /dev/null 2>&1; echo $?` = 0 ]; then
      ./configure --prefix=$PREFIX \
        --with-cxx=mpicxx --with-cc=mpicc --with-fc=mpif90 \
	--CFLAGS="-mt" --CXXFLAGS="-mt" --FFLAGS="-mt" \
        --COPTFLAGS="-O3" --CXXOPTFLAGS="-O3" --FOPTFLAGS="-O3" \
        --with-mpiexec="mpiexec" \
        --with-blas-lapack-dir="$MKL/bin" \
        --with-make-np=1 \
        --with-c++-support=1 --with-debugging=0
    else
      ./configure --prefix=$PREFIX \
        --with-cxx=icpc --with-cc=icc --with-fc=ifort \
        --COPTFLAGS="-O3" --CXXOPTFLAGS="-O3" --FOPTFLAGS="-O3" \
        --with-mpi-include="/usr/include" --with-mpi-lib="-lmpi_mt" \
        --with-blas-lapack-dir="$MKL/bin" \
        --with-make-np=1 \
        --with-c++-support=1 --with-debugging=0
    fi
  else
    if [ `which mpicxx > /dev/null 2>&1; echo $?` = 0 ]; then
      ./configure --prefix=$PREFIX \
        --with-cxx=mpicxx --with-cc=mpicc --with-fc=mpif90 \
	--CFLAGS="-mt" --CXXFLAGS="-mt" --FFLAGS="-mt" \
        --COPTFLAGS="-O0 -g" --CXXOPTFLAGS="-O0 -g" --FOPTFLAGS="-O0 -g" \
        --with-mpiexec="mpiexec" \
        --with-blas-lapack-dir="$MKL/bin" \
	--with-make-np=1 \
        --with-c++-support=1 --with-debugging=1
    else
      ./configure --prefix=$PREFIX \
        --with-cxx=icpc --with-cc=icc --with-fc=ifort \
        --COPTFLAGS="-O0 -g" --CXXOPTFLAGS="-O0 -g" --FOPTFLAGS="-O0 -g" \
        --with-mpi-include="/usr/include" --with-mpi-lib="-lmpi_mt" \
        --with-blas-lapack-dir="$MKL/bin" \
        --with-make-np=1 \
        --with-c++-support=1 --with-debugging=1
    fi
  fi
  check make
  $SUDO env LD_LIBRARY_PATH=$LD_LIBRARY_PATH make install
done

DATE=$(date +%Y%m%d-%H%M%S)
for build_type in $BUILD_TYPES; do
  PREFIX=$PREFIX_ROKKO/petsc-$PETSC_VERSION-$PETSC_RK_REVISION/$build_type
  cat << EOF > $BUILD_DIR/petscvars.sh
# petsc $(basename $0 .sh) $PETSC_VERSION $PETSC_RK_REVISION $DATE
export PETSC_ROOT=$PREFIX_ROKKO/petsc-$PETSC_VERSION-$PETSC_RK_REVISION
export LD_LIBRARY_PATH=$PREFIX/lib:\$LD_LIBRARY_PATH
EOF
  $SUDO cp -f $BUILD_DIR/petscvars.sh $PREFIX
done

cat << EOF > $BUILD_DIR/petscvars.sh
# petsc $(basename $0 .sh) $PETSC_VERSION $PETSC_RK_REVISION $DATE
export PETSC_ROOT=$PREFIX_ROKKO/petsc-$PETSC_VERSION-$PETSC_RK_REVISION
EOF
$SUDO cp -f $BUILD_DIR/petscvars.sh $PREFIX_ROKKO/petsc-$PETSC_VERSION-$PETSC_RK_REVISION
