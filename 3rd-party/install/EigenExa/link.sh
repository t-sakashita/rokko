#!/bin/sh

SCRIPT_DIR=$(cd "$(dirname $0)"; pwd)
. $SCRIPT_DIR/../util.sh
. $SCRIPT_DIR/version.sh
set_prefix

if [ -d $PREFIX_ROKKO/eigenexa-$EIGENEXA_VERSION-$EIGENEXA_RK_REVISION ]; then
  $SUDO rm -f $PREFIX_ROKKO/eigenexa
  $SUDO ln -s $PREFIX_ROKKO/eigenexa-$EIGENEXA_VERSION-$EIGENEXA_RK_REVISION $PREFIX_ROKKO/eigenexa
fi

BUILD_TYPES="Release Debug"
for build_type in $BUILD_TYPES; do
  if [ -f $PREFIX_ROKKO/eigenexa/$build_type/eigenexavars.sh ]; then
    $SUDO rm -f $PREFIX_ROKKO/rokkoenv-$build_type.d/eigenexavars.sh
    $SUDO ln -s $PREFIX_ROKKO/eigenexa/$build_type/eigenexavars.sh $PREFIX_ROKKO/rokkoenv-$build_type.d/eigenexavars.sh
  fi
done

if [ -f $PREFIX_ROKKO/eigenexa/eigenexavars.sh ]; then
  $SUDO rm -f $PREFIX_ROKKO/rokkoenv.d/eigenexavars.sh
  $SUDO ln -s $PREFIX_ROKKO/eigenexa/eigenexavars.sh $PREFIX_ROKKO/rokkoenv.d/eigenexavars.sh
fi

ARCH_TYPES="Linux-s64fx Linux-x86_64"
for arch_type in $ARCH_TYPES; do
  BUILD_TYPES="Release Debug"
  for build_type in $BUILD_TYPES; do
    if [ -f $PREFIX_ROKKO/eigenexa/$arch_type/$build_type/eigenexavars.sh ]; then
      $SUDO rm -f $PREFIX_ROKKO/rokkoenv-$arch_type-$build_type.d/eigenexavars.sh
      $SUDO ln -s $PREFIX_ROKKO/eigenexa/$arch_type/$build_type/eigenexavars.sh $PREFIX_ROKKO/rokkoenv-$arch_type-$build_type.d/eigenexavars.sh
    fi
  done

  if [ -f $PREFIX_ROKKO/eigenexa/$arch_type/eigenexavars.sh ]; then
    $SUDO rm -f $PREFIX_ROKKO/rokkoenv-$arch_type.d/eigenexavars.sh
    $SUDO ln -s $PREFIX_ROKKO/eigenexa/$arch_type/eigenexavars.sh $PREFIX_ROKKO/rokkoenv-$arch_type.d/eigenexavars.sh
  fi
done
