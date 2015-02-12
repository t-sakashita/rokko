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
      --COPTFLAGS="-O3" --CXXOPTFLAGS="-O3" --FOPTFLAGS="-O3" \
      --with-scalapack-lib="/usr/lib64/openmpi/lib/libscalapack.so /usr/lib64/openmpi/lib/libmpiblacs.so /usr/lib64/openmpi/lib/libmpiblacsF77init.so" \
      --with-c++-support=1 --with-debugging=0
  else
    ./configure --prefix=$PREFIX \
      --COPTFLAGS="-g -O0" --CXXOPTFLAGS="-g -O0" --FOPTFLAGS="-g -O0" \
      --with-scalapack-lib="/usr/lib64/openmpi/lib/libscalapack.so /usr/lib64/openmpi/lib/libmpiblacs.so /usr/lib64/openmpi/lib/libmpiblacsF77init.so" \
      --with-c++-support=1 --with-debugging=1
  fi
  check make
  $SUDO env LD_LIBRARY_PATH=$LD_LIBRARY_PATH make install
  cat << EOF > $BUILD_DIR/petscvars.sh
export PETSC_ROOT=$PREFIX_ROKKO/petsc-$PETSC_VERSION-$PETSC_PATCH_VERSION
export LD_LIBRARY_PATH=$PREFIX/lib:$LD_LIBRARY_PATH
EOF
  $SUDO cp -f $BUILD_DIR/petscvars.sh $PREFIX
done

cat << EOF > $BUILD_DIR/petscvars.sh
export PETSC_ROOT=$PREFIX_ROKKO/petsc-$PETSC_VERSION-$PETSC_PATCH_VERSION
EOF
$SUDO cp -f $BUILD_DIR/petscvars.sh $PREFIX_ROKKO/petsc-$PETSC_VERSION-$PETSC_PATCH_VERSION
