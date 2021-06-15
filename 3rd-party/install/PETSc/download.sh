#!/bin/bash -x

SCRIPT_DIR=$(cd "$(dirname $0)"; pwd)
. $SCRIPT_DIR/../util.sh
. $SCRIPT_DIR/version.sh
set_prefix

FILE="petsc-$PETSC_VERSION.tar.gz"
if [ -f $SOURCE_DIR/$FILE ]; then :; else
  wget -O $SOURCE_DIR/$FILE http://ftp.mcs.anl.gov/pub/petsc/release-snapshots/$FILE
fi
