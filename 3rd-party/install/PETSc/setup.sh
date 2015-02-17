#!/bin/bash

SCRIPT_DIR=$(cd "$(dirname $0)"; pwd)
. $SCRIPT_DIR/../util.sh
. $SCRIPT_DIR/version.sh
set_prefix

cd $BUILD_DIR
rm -rf petsc-$PETSC_VERSION*
if [ -f $HOME/source/petsc-$PETSC_VERSION.tar.gz ]; then
  check tar zxf $HOME/source/petsc-$PETSC_VERSION.tar.gz
else
  check wget http://ftp.mcs.anl.gov/pub/petsc/release-snapshots/petsc-$PETSC_VERSION.tar.gz
  check tar zxf petsc-$PETSC_VERSION.tar.gz
fi
cd petsc-$PETSC_VERSION
patch -p1 < $SCRIPT_DIR/petsc-3.5.2.patch
