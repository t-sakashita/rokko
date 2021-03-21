#!/bin/sh

SCRIPT_DIR=$(cd "$(dirname $0)"; pwd)
. $SCRIPT_DIR/../util.sh
. $SCRIPT_DIR/version.sh
set_prefix

__NAME__="scalapack"
__VERSION__=${SCALAPACK_VERSION}
__RK_REVISION__=${SCALAPACK_RK_REVISION}

BUILD_TYPES="Release Debug"
for build_type in $BUILD_TYPES; do
  if [ -f $PREFIX_ROKKO/${__NAME__}/${__NAME__}-${__VERSION__}-${__RK_REVISION__}/$build_type/${__NAME__}vars.sh ]; then
    $SUDO rm -f $PREFIX_ROKKO/rokkoenv-$build_type.d/${__NAME__}vars.sh
    $SUDO ln -s $PREFIX_ROKKO/${__NAME__}/${__NAME__}-${__VERSION__}-${__RK_REVISION__}/$build_type/${__NAME__}vars.sh $PREFIX_ROKKO/rokkoenv-$build_type.d/${__NAME__}vars.sh
  fi
done

if [ -f $PREFIX_ROKKO/${__NAME__}/${__NAME__}-${__VERSION__}-${__RK_REVISION__}/${__NAME__}vars.sh ]; then
  $SUDO rm -f $PREFIX_ROKKO/rokkoenv.d/${__NAME__}vars.sh
  $SUDO ln -s $PREFIX_ROKKO/${__NAME__}/${__NAME__}-${__VERSION__}-${__RK_REVISION__}/${__NAME__}vars.sh $PREFIX_ROKKO/rokkoenv.d/${__NAME__}vars.sh
fi

ARCH_TYPES="Linux-s64fx Linux-x86_64"
for arch_type in $ARCH_TYPES; do
  BUILD_TYPES="Release Debug"
  for build_type in $BUILD_TYPES; do
    if [ -f $PREFIX_ROKKO/${__NAME__}/${__NAME__}-${__VERSION__}-${__RK_REVISION__}/$arch_type/$build_type/${__NAME__}vars.sh ]; then
      $SUDO rm -f $PREFIX_ROKKO/rokkoenv-$arch_type-$build_type.d/${__NAME__}vars.sh
      $SUDO ln -s $PREFIX_ROKKO/${__NAME__}/${__NAME__}-${__VERSION__}-${__RK_REVISION__}/$arch_type/$build_type/${__NAME__}vars.sh $PREFIX_ROKKO/rokkoenv-$arch_type-$build_type.d/${__NAME__}vars.sh
    fi
  done

  if [ -f $PREFIX_ROKKO/${__NAME__}/${__NAME__}-${__VERSION__}-${__RK_REVISION__}/$arch_type/${__NAME__}vars.sh ]; then
    $SUDO rm -f $PREFIX_ROKKO/rokkoenv-$arch_type.d/${__NAME__}vars.sh
    $SUDO ln -s $PREFIX_ROKKO/${__NAME__}/${__NAME__}-${__VERSION__}-${__RK_REVISION__}/$arch_type/${__NAME__}vars.sh $PREFIX_ROKKO/rokkoenv-$arch_type.d/${__NAME__}vars.sh
  fi
done
