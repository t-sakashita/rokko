#!/bin/sh

SCRIPT_DIR=$(cd "$(dirname $0)"; pwd)
. $SCRIPT_DIR/../util.sh
. $SCRIPT_DIR/version.sh
set_prefix

if [ -d $PREFIX_ROKKO/petsc-$PETSC_VERSION-$PETSC_PATCH_VERSION ]; then
  $SUDO rm -f $PREFIX_ROKKO/petsc
  $SUDO ln -s $PREFIX_ROKKO/petsc-$PETSC_VERSION-$PETSC_PATCH_VERSION $PREFIX_ROKKO/petsc
fi

BUILD_TYPES="Release Debug"
for build_type in $BUILD_TYPES; do
  if [ -f $PREFIX_ROKKO/petsc/$build_type/petscvars.sh ]; then
    $SUDO rm -f $PREFIX_ROKKO/rokkoenv-$build_type.d/petscvars.sh
    $SUDO ln -s $PREFIX_ROKKO/petsc/$build_type/petscvars.sh $PREFIX_ROKKO/rokkoenv-$build_type.d/petscvars.sh
  fi
done

if [ -f $PREFIX_ROKKO/petsc/petscvars.sh ]; then
  $SUDO rm -f $PREFIX_ROKKO/rokkoenv.d/petscvars.sh
  $SUDO ln -s $PREFIX_ROKKO/petsc/petscvars.sh $PREFIX_ROKKO/rokkoenv.d/petscvars.sh
fi

ARCH_TYPES="Linux-s64fx Linux-x86_64"
for arch_type in $ARCH_TYPES; do
  BUILD_TYPES="Release Debug"
  for build_type in $BUILD_TYPES; do
    if [ -f $PREFIX_ROKKO/petsc/$arch_type/$build_type/petscvars.sh ]; then
      $SUDO rm -f $PREFIX_ROKKO/rokkoenv-$arch_type-$build_type.d/petscvars.sh
      $SUDO ln -s $PREFIX_ROKKO/petsc/$arch_type/$build_type/petscvars.sh $PREFIX_ROKKO/rokkoenv-$arch_type-$build_type.d/petscvars.sh
    fi
  done

  if [ -f $PREFIX_ROKKO/petsc/$arch_type/petscvars.sh ]; then
    $SUDO rm -f $PREFIX_ROKKO/rokkoenv-$arch_type.d/petscvars.sh
    $SUDO ln -s $PREFIX_ROKKO/petsc/$arch_type/petscvars.sh $PREFIX_ROKKO/rokkoenv-$arch_type.d/petscvars.sh
  fi
done
