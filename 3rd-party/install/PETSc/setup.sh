#!/bin/bash

SCRIPT_DIR=$(cd "$(dirname $0)"; pwd)
. $SCRIPT_DIR/../util.sh
. $SCRIPT_DIR/version.sh
set_prefix

sh $SCRIPT_DIR/download.sh

cd $BUILD_DIR
if [ -f $SOURCE_DIR/petsc-$PETSC_VERSION.tar.gz ]; then
  check mkdir -p petsc-$PETSC_VERSION
  tar zxf $SOURCE_DIR/petsc-$PETSC_VERSION.tar.gz -C petsc-$PETSC_VERSION --strip-components=1
  cd petsc-$PETSC_VERSION
  if [ -f $SCRIPT_DIR/petsc-$PETSC_VERSION.patch ]; then
    echo "applying petsc-$PETSC_VERSION.patch"
    check patch -p1 < $SCRIPT_DIR/petsc-$PETSC_VERSION.patch
  fi
fi
