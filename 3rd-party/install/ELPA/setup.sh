#!/bin/bash -x

SCRIPT_DIR=$(cd "$(dirname $0)"; pwd)
. $SCRIPT_DIR/../util.sh
. $SCRIPT_DIR/version.sh
set_prefix

cd $BUILD_DIR
rm -rf elpa-$ELPA_VERSION*
if test -f $SOURCE_DIR/elpa-$ELPA_VERSION.tar.gz; then :; else
  wget -O $SOURCE_DIR/elpa-$ELPA_VERSION.tar.gz http://elpa.mpcdf.mpg.de/html/Releases/$ELPA_VERSION/elpa-$ELPA_VERSION.tar.gz
fi
tar zxf $SOURCE_DIR/elpa-$ELPA_VERSION.tar.gz

if [ -f $SCRIPT_DIR/elpa-$ELPA_VERSION.patch ]; then
  cd elpa-$ELPA_VERSION
  check patch -p1 < $SCRIPT_DIR/elpa-$ELPA_VERSION.patch
fi
