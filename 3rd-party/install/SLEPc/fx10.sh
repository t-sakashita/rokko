#!/bin/bash

SCRIPT_DIR=$(cd "$(dirname $0)"; pwd)
. $SCRIPT_DIR/../util.sh
. $SCRIPT_DIR/version.sh
set_prefix

sh $SCRIPT_DIR/setup.sh

. $PREFIX_ROKKO/rokkoenv-Linux-s64fx.d/petscvars.sh
BUILD_TYPES="Release Debug"
for build_type in $BUILD_TYPES; do
  PREFIX=$PREFIX_ROKKO/slepc-$SLEPC_VERSION-$SLEPC_RK_REVISION/Linux-s64fx/$build_type
  cd $BUILD_DIR
  cp -rp slepc-$SLEPC_VERSION slepc-$SLEPC_VERSION-build-Linux-s64fx-$build_type
  cd slepc-$SLEPC_VERSION-build-Linux-s64fx-$build_type
  export SLEPC_DIR=$PWD
  export PETSC_DIR=$PETSC_ROOT/$build_type
  ./configure --prefix=$PREFIX
  make
  $SUDO env LD_LIBRARY_PATH=$LD_LIBRARY_PATH make install SLEPC_DIR=$PWD PETSC_DIR=$PETSC_ROOT/$build_type
done

DATE=$(date +%Y%m%d-%H%M%S)
for build_type in $BUILD_TYPES; do
  PREFIX=$PREFIX_ROKKO/slepc-$SLEPC_VERSION-$SLEPC_RK_REVISION/Linux-s64fx/$build_type
  cat << EOF > $BUILD_DIR/slepcvars.sh
# slepc $(basename $0 .sh) $SLEPC_VERSION $SLEPC_RK_REVISION $DATE
export SLEPC_ROOT=$PREFIX_ROKKO/slepc-$SLEPC_VERSION-$SLEPC_RK_REVISION/Linux-s64fx
export LD_LIBRARY_PATH=$PREFIX/lib:\$LD_LIBRARY_PATH
EOF
  $SUDO cp -f $BUILD_DIR/slepcvars.sh $PREFIX
done

cat << EOF > $BUILD_DIR/slepcvars.sh
# slepc $(basename $0 .sh) $SLEPC_VERSION $SLEPC_RK_REVISION $DATE
export SLEPC_ROOT=$PREFIX_ROKKO/slepc-$SLEPC_VERSION-$SLEPC_RK_REVISION/Linux-s64fx
EOF
$SUDO cp -f $BUILD_DIR/slepcvars.sh $PREFIX_ROKKO/slepc-$SLEPC_VERSION-$SLEPC_RK_REVISION/Linux-s64fx
