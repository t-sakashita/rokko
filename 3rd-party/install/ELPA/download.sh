#!/bin/bash -x

SCRIPT_DIR=$(cd "$(dirname $0)"; pwd)
. $SCRIPT_DIR/../util.sh
. $SCRIPT_DIR/version.sh
set_prefix

FILE=elpa-$ELPA_VERSION.tar.gz
if [ -f $SOURCE_DIR/$FILE ]; then :; else
  wget -O $SOURCE_DIR/$FILE https://elpa.mpcdf.mpg.de/software/tarball-archive/Releases/$ELPA_VERSION/$FILE
fi
