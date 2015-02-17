#!/bin/bash

SCRIPT_DIR=$(cd "$(dirname $0)"; pwd)
. $SCRIPT_DIR/../util.sh
. $SCRIPT_DIR/version.sh
set_prefix

cd $BUILD_DIR
rm -rf trilinos-$TRILINOS_VERSION-Source*
check tar jxf $HOME/source/trilinos-$TRILINOS_VERSION-Source.tar.bz2
cd trilinos-$TRILINOS_VERSION-Source
patch -p1 < $SCRIPT_DIR/trilinos-$TRILINOS_VERSION.patch
