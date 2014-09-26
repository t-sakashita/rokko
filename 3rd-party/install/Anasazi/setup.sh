#!/bin/bash

SCRIPT_DIR=$(cd "$(dirname $0)"; pwd)
. $SCRIPT_DIR/../util.sh
. $SCRIPT_DIR/version.sh
set_build_dir

cd $BUILD_DIR
rm -rf trilinos-$TRILINOS_VERSION-Source*
check tar jxf $HOME/source/trilinos-$TRILINOS_VERSION-Source.tar.bz2
