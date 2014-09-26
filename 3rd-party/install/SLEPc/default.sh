#!/bin/bash

SCRIPT_DIR=$(cd "$(dirname $0)"; pwd)
. $SCRIPT_DIR/../util.sh
. $SCRIPT_DIR/version.sh
set_prefix
set_build_dir

. $PREFIX_OPT/env.sh

sh $SCRIPT_DIR/setup.sh

BUILD_TYPES="Release Debug"

for build_type in $BUILD_TYPES; do
  cd $BUILD_DIR
  cp -rp slepc-$SLEPC_VERSION slepc-$SLEPC_VERSION-build-$build_type
  cd slepc-$SLEPC_VERSION-build-$build_type
  ./configure --prefix=$PREFIX_ROKKO/$build_type
  make SLEPC_DIR=$PWD PETSC_DIR=$PREFIX_ROKKO/$build_type
  $SUDO LD_LIBRARY_PATH=$LD_LIBRARY_PATH make SLEPC_DIR=$PWD PETSC_DIR=$PREFIX_ROKKO/$build_type install
done
