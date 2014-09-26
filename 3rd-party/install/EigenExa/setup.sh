#!/bin/bash -x

SCRIPT_DIR=$(cd "$(dirname $0)"; pwd)
. $SCRIPT_DIR/../util.sh
. $SCRIPT_DIR/version.sh
set_build_dir

cd $BUILD_DIR
rm -rf EigenExa-$EIGENEXA_VERSION*
if test -f $HOME/source/EigenExa-$EIGENEXA_VERSION.tgz; then
  tar zxf $HOME/source/EigenExa-$EIGENEXA_VERSION.tgz
else
  wget -O - http://www.aics.riken.jp/labs/lpnctrt/EigenExa-$EIGENEXA_VERSION.tgz | tar zxf -
fi
cd $BUILD_DIR/EigenExa-$EIGENEXA_VERSION
patch -p1 < $SCRIPT_DIR/EigenExa-$EIGENEXA_VERSION.patch
