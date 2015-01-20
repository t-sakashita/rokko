#!/bin/bash

SCRIPT_DIR=$(cd "$(dirname $0)"; pwd)
. $SCRIPT_DIR/../util.sh
. $SCRIPT_DIR/version.sh
set_prefix

sh $SCRIPT_DIR/setup.sh

source $PREFIX_ROKKO/rokkoenv.d/petscvars.sh
BUILD_TYPES="Release Debug"
for build_type in $BUILD_TYPES; do
  PREFIX=$PREFIX_ROKKO/slpec-$SLEPC_VERSION-$SLEPC_PATCH_VERSION/$build_type
  cd $BUILD_DIR
  cp -rp slepc-$SLEPC_VERSION slepc-$SLEPC_VERSION-build-$build_type
  cd slepc-$SLEPC_VERSION-build-$build_type
  export SLEPC_DIR=$PWD
  export PETSC_DIR=$PETSC_ROOT/$build_type
  ./configure --prefix=$PREFIX
  make
  $SUDO make install SLEPC_DIR=$PWD PETSC_DIR=$PETSC_ROOT/$build_type
  cat << EOF > $BUILD_DIR/slepcvars.sh
export SLEPC_ROOT=$PREFIX_ROKKO/slepc-$SLEPC_VERSION-$SLEPC_PATCH_VERSION
export LD_LIBRARY_PATH=$PREFIX/lib:$LD_LIBRARY_PATH
EOF
  $SUDO cp -f $BUILD_DIR/slepcvars.sh $PREFIX
done

cat << EOF > $BUILD_DIR/slepcvars.sh
export SLEPC_ROOT=$PREFIX_ROKKO/slepc-$SLEPC_VERSION-$SLEPC_PATCH_VERSION
EOF
$SUDO cp -f $BUILD_DIR/slepcvars.sh $PREFIX_ROKKO/slepc-$SLEPC_VERSION-$SLEPC_PATCH_VERSION
