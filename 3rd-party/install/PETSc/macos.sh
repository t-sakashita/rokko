#!/bin/bash

SCRIPT_DIR=$(cd "$(dirname $0)"; pwd)
. $SCRIPT_DIR/../util.sh
. $SCRIPT_DIR/version.sh
set_prefix

sh $SCRIPT_DIR/setup.sh

unset PETSC_DIR
BUILD_TYPES="Release Debug"
for build_type in $BUILD_TYPES; do
  PREFIX=$PREFIX_ROKKO/petsc-$PETSC_VERSION-$PETSC_PATCH_VERSION/$build_type
  cd $BUILD_DIR
  cp -rp petsc-$PETSC_VERSION petsc-$PETSC_VERSION-build-$build_type
  cd petsc-$PETSC_VERSION-build-$build_type
  if [ $build_type == "Release" ]; then
    ./configure --prefix=$PREFIX \
      --with-cc=mpicc --COPTFLAGS="-O3" --with-cxx=mpicxx --CXXOPTFLAGS="-O3" --with-fc=mpif77 --FOPTFLAGS="-O3" --with-mpiexec=mpiexec \
    --with-blas-lapack-lib="-framework Accelerate" \
    --with-blacs-dir=$PREFIX_ROKKO/scalapack-2.0.2-0/$build_type --with-scalapack-dir=$PREFIX_ROKKO/scalapack-2.0.2-0/$build_type --with-c++-support=1 --with-debugging=0
    check make PETSC_DIR=$BUILD_DIR/petsc-$PETSC_VERSION-build-$build_type PETSC_ARCH=arch-darwin-c-debug
    $SUDO make install PETSC_DIR=$BUILD_DIR/petsc-$PETSC_VERSION-build-$build_type PETSC_ARCH=arch-darwin-c-debug
  else
    ./configure --prefix=$PREFIX \
      --with-cc=mpicc --COPTFLAGS="-g -O0" --with-cxx=mpicxx --CXXOPTFLAGS="-g -O0" --with-fc=mpif77 --FOPTFLAGS="-g -O0" --with-mpiexec=mpiexec \
    --with-blas-lapack-lib="-framework Accelerate" \
    --with-blacs-dir=$PREFIX_ROKKO/scalapack-2.0.2-0/$build_type --with-scalapack-dir=$PREFIX_ROKKO/scalapack-2.0.2-0/$build_type --with-c++-support=1 --with-debugging=1
    check make PETSC_DIR=$BUILD_DIR/petsc-$PETSC_VERSION-build-$build_type PETSC_ARCH=arch-darwin-c-debug
    $SUDO make install PETSC_DIR=$BUILD_DIR/petsc-$PETSC_VERSION-build-$build_type PETSC_ARCH=arch-darwin-c-debug
  fi
  cat << EOF > $BUILD_DIR/petscvars.sh
export PETSC_ROOT=$PREFIX_ROKKO/petsc-$PETSC_VERSION-$PETSC_PATCH_VERSION
export LD_LIBRARY_PATH=$PREFIX/lib:\$LD_LIBRARY_PATH
EOF
  $SUDO cp -f $BUILD_DIR/petscvars.sh $PREFIX
done

cat << EOF > $BUILD_DIR/petscvars.sh
export PETSC_ROOT=$PREFIX_ROKKO/petsc-$PETSC_VERSION-$PETSC_PATCH_VERSION
EOF
$SUDO cp -f $BUILD_DIR/petscvars.sh $PREFIX_ROKKO/petsc-$PETSC_VERSION-$PETSC_PATCH_VERSION
