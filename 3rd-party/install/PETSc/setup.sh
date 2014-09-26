#!/bin/bash

SCRIPT_DIR=$(cd "$(dirname $0)"; pwd)
. $SCRIPT_DIR/../util.sh
. $SCRIPT_DIR/version.sh
set_build_dir

cd $BUILD_DIR
rm -rf petsc-$PETSC_VERSION*
if [ -f $HOME/source/petsc-$PETSC_VERSION.tar.gz ]; then
  check tar zxf $HOME/source/petsc-$PETSC_VERSION.tar.gz
else
  check wget -O - http://ftp.mcs.anl.gov/pub/petsc/release-snapshots/petsc-$PETSC_VERSION.tar.gz | tar zxf -
fi
