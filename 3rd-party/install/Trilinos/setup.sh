#!/bin/bash

SCRIPT_DIR=$(cd "$(dirname $0)"; pwd)
. $SCRIPT_DIR/../util.sh
. $SCRIPT_DIR/version.sh
set_prefix

cd $BUILD_DIR

rm -rf trilinos-$TRILINOS_VERSION-Source*
if [ -f $HOME/source/trilinos-$TRILINOS_VERSION-Source.tar.bz2 ]; then
  check tar jxf $HOME/source/trilinos-$TRILINOS_VERSION-Source.tar.bz2
else
  echo "Error: trilinos-$TRILINOS_VERSION-Source.tar.bz2 not found"
  exit 127
fi

if [ -f $SCRIPT_DIR/trilinos-$TRILINOS_VERSION.patch ]; then
  cd trilinos-$TRILINOS_VERSION-Source
  check patch -p1 < $SCRIPT_DIR/trilinos-$TRILINOS_VERSION.patch
fi
