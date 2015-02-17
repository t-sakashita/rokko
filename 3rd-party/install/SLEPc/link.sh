#!/bin/sh

SCRIPT_DIR=$(cd "$(dirname $0)"; pwd)
. $SCRIPT_DIR/../util.sh
. $SCRIPT_DIR/version.sh
set_prefix

if [ -d $PREFIX_ROKKO/slepc-$SLEPC_VERSION-$SLEPC_PATCH_VERSION ]; then
  $SUDO rm -f $PREFIX_ROKKO/slepc
  $SUDO ln -s $PREFIX_ROKKO/slepc-$SLEPC_VERSION-$SLEPC_PATCH_VERSION $PREFIX_ROKKO/slepc
fi

BUILD_TYPES="Release Debug"
for build_type in $BUILD_TYPES; do
  if [ -f $PREFIX_ROKKO/slepc/$build_type/slepcvars.sh ]; then
    $SUDO rm -f $PREFIX_ROKKO/rokkoenv-$build_type.d/slepcvars.sh
    $SUDO ln -s $PREFIX_ROKKO/slepc/$build_type/slepcvars.sh $PREFIX_ROKKO/rokkoenv-$build_type.d/slepcvars.sh
  fi
done

if [ -f $PREFIX_ROKKO/slepc/slepcvars.sh ]; then
  $SUDO rm -f $PREFIX_ROKKO/rokkoenv.d/slepcvars.sh
  $SUDO ln -s $PREFIX_ROKKO/slepc/slepcvars.sh $PREFIX_ROKKO/rokkoenv.d/slepcvars.sh
fi

ARCH_TYPES="Linux-s64fx Linux-x86_64"
for arch_type in $ARCH_TYPES; do
  BUILD_TYPES="Release Debug"
  for build_type in $BUILD_TYPES; do
    if [ -f $PREFIX_ROKKO/slepc/$arch_type/$build_type/slepcvars.sh ]; then
      $SUDO rm -f $PREFIX_ROKKO/rokkoenv-$arch_type-$build_type.d/slepcvars.sh
      $SUDO ln -s $PREFIX_ROKKO/slepc/$arch_type/$build_type/slepcvars.sh $PREFIX_ROKKO/rokkoenv-$arch_type-$build_type.d/slepcvars.sh
    fi
  done

  if [ -f $PREFIX_ROKKO/slepc/$arch_type/slepcvars.sh ]; then
    $SUDO rm -f $PREFIX_ROKKO/rokkoenv-$arch_type.d/slepcvars.sh
    $SUDO ln -s $PREFIX_ROKKO/slepc/$arch_type/slepcvars.sh $PREFIX_ROKKO/rokkoenv-$arch_type.d/slepcvars.sh
  fi
done
