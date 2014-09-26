#!/bin/bash

SCRIPT_DIR=$(cd "$(dirname $0)"; pwd)
. $SCRIPT_DIR/../util.sh
. $SCRIPT_DIR/version.sh
set_build_dir

cd $BUILD_DIR
rm -rf elpa_lib-$ELPA_VERSION*
check tar zxf $HOME/source/elpa_lib-$ELPA_VERSION.tar.gz
cd $BUILD_DIR/elpa_lib-$ELPA_VERSION
patch -p1 < $SCRIPT_DIR/elpa_lib-201305.patch
