#!/bin/bash

SCRIPT_DIR=$(cd "$(dirname $0)"; pwd)
. $SCRIPT_DIR/../util.sh
. $SCRIPT_DIR/version.sh
set_prefix

cd $BUILD_DIR
rm -rf trilinos-$TRILINOS_VERSION-Source*
if [ -f $SOURCE_DIR/trilinos-$TRILINOS_VERSION-Source.tar.bz2 ]; then
  check tar jxf $SOURCE_DIR/trilinos-$TRILINOS_VERSION-Source.tar.bz2
elif [ -f $SOURCE_DIR/trilinos-release-$TRILINOS_VERSION_DASHED.tar.gz ]; then
  check tar zxf $SOURCE_DIR/trilinos-release-$TRILINOS_VERSION_DASHED.tar.gz
else
  check wget -O - https://github.com/trilinos/Trilinos/archive/trilinos-release-$TRILINOS_VERSION_DASHED.tar.gz | tar zxf -
fi

if [ -d trilinos-$TRILINOS_VERSION-Source ]; then :; else
  mv Trilinos-trilinos-release-$TRILINOS_VERSION_DASHED trilinos-$TRILINOS_VERSION-Source
fi

if [ -f $SCRIPT_DIR/trilinos-$TRILINOS_VERSION.patch ]; then
  cd trilinos-$TRILINOS_VERSION-Source
  check patch -p1 < $SCRIPT_DIR/trilinos-$TRILINOS_VERSION.patch
fi
