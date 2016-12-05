#!/bin/bash

SCRIPT_DIR=$(cd "$(dirname $0)"; pwd)
. $SCRIPT_DIR/../util.sh
. $SCRIPT_DIR/version.sh
set_prefix

cd $BUILD_DIR
rm -rf elpa-$ELPA_VERSION*
if test -f $SOURCE_DIR/elpa-$ELPA_VERSION.tar.gz; then
    tar zxf $SOURCE_DIR/elpa-$ELPA_VERSION.tar.gz
else
    wget -O - http://elpa.mpcdf.mpg.de/html/Releases/$ELPA_VERSION/elpa-$ELPA_VERSION.tar.gz | tar zxf -
fi

if [ -f $SCRIPT_DIR/elpa-$ELPA_VERSION.patch ]; then
  cd elpa-$ELPA_VERSION
  check patch -p1 < $SCRIPT_DIR/elpa-$ELPA_VERSION.patch
fi

