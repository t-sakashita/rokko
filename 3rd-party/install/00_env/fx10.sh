#!/bin/sh

SCRIPT_DIR=$(cd "$(dirname $0)"; pwd)
. $SCRIPT_DIR/../util.sh
set_prefix

ARCH_TYPES="Linux-s64fx Linux-x86_64"
for arch_type in $ARCH_TYPES; do
  BUILD_TYPES="Release Debug"
  for build_type in $BUILD_TYPES; do
    cat << EOF > $BUILD_DIR/rokkoenv-$arch_type-$build_type.sh
export ROKKO_SOLVER_ROOT=$PREFIX_ROKKO/$arch_type
for i in $PREFIX_ROKKO/rokkoenv-$arch_type-$build_type.d/*.sh ; do
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
    $SUDO mkdir -p $PREFIX_ROKKO/rokkoenv-$arch_type-$build_type.d
    $SUDO cp -f $BUILD_DIR/rokkoenv-$arch_type-$build_type.sh $PREFIX_ROKKO
  done

  cat << EOF > $BUILD_DIR/rokkoenv-$arch_type.sh
export ROKKO_SOLVER_ROOT=$PREFIX_ROKKO/$arch_type
for i in $PREFIX_ROKKO/rokkoenv-$arch_type.d/*.sh ; do
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
  $SUDO mkdir -p $PREFIX_ROKKO/rokkoenv-$arch_type.d
  $SUDO cp -f $BUILD_DIR/rokkoenv-$arch_type.sh $PREFIX_ROKKO
done
