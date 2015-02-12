#!/bin/sh

SCRIPT_DIR=$(cd "$(dirname $0)"; pwd)
. $SCRIPT_DIR/../util.sh
set_prefix

BUILD_TYPES="Release Debug"
for build_type in $BUILD_TYPES; do
  cat << EOF > $BUILD_DIR/rokkoenv-$build_type.sh
export ROKKO_SOLVER_ROOT=$PREFIX_ROKKO
for i in $PREFIX_ROKKO/rokkoenv-$build_type.d/*.sh ; do
  if [ -r "\$i" ]; then
    if [ "\${-#*i}" != "\$-" ]; then
      . "\$i"
    else
      . "\$i" >/dev/null 2>&1
    fi
  fi
done
unset i
EOF
  $SUDO mkdir -p $PREFIX_ROKKO/rokkoenv-$build_type.d
  $SUDO cp -f $BUILD_DIR/rokkoenv-$build_type.sh $PREFIX_ROKKO
done

cat << EOF > $BUILD_DIR/rokkoenv.sh
export ROKKO_SOLVER_ROOT=$PREFIX_ROKKO
for i in $PREFIX_ROKKO/rokkoenv.d/*.sh ; do
  if [ -r "\$i" ]; then
    if [ "\${-#*i}" != "\$-" ]; then
      . "\$i"
    else
      . "\$i" >/dev/null 2>&1
    fi
  fi
done
unset i
EOF
$SUDO mkdir -p $PREFIX_ROKKO/rokkoenv.d
$SUDO cp -f $BUILD_DIR/rokkoenv.sh $PREFIX_ROKKO
