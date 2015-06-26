#!/bin/sh

SCRIPT_DIR=$(cd "$(dirname $0)"; pwd)
. $SCRIPT_DIR/../util.sh
. $SCRIPT_DIR/version.sh
set_prefix

if [ -d $PREFIX_ROKKO/trilinos-$TRILINOS_VERSION-$TRILINOS_RK_REVISION ]; then
  $SUDO rm -f $PREFIX_ROKKO/trilinos
  $SUDO ln -s $PREFIX_ROKKO/trilinos-$TRILINOS_VERSION-$TRILINOS_RK_REVISION $PREFIX_ROKKO/trilinos
fi

BUILD_TYPES="Release Debug"
for build_type in $BUILD_TYPES; do
  if [ -f $PREFIX_ROKKO/trilinos/$build_type/trilinosvars.sh ]; then
    $SUDO rm -f $PREFIX_ROKKO/rokkoenv-$build_type.d/trilinosvars.sh
    $SUDO ln -s $PREFIX_ROKKO/trilinos/$build_type/trilinosvars.sh $PREFIX_ROKKO/rokkoenv-$build_type.d/trilinosvars.sh
  fi
done

if [ -f $PREFIX_ROKKO/trilinos/trilinosvars.sh ]; then
  $SUDO rm -f $PREFIX_ROKKO/rokkoenv.d/trilinosvars.sh
  $SUDO ln -s $PREFIX_ROKKO/trilinos/trilinosvars.sh $PREFIX_ROKKO/rokkoenv.d/trilinosvars.sh
fi

ARCH_TYPES="Linux-s64fx Linux-x86_64"
for arch_type in $ARCH_TYPES; do
  BUILD_TYPES="Release Debug"
  for build_type in $BUILD_TYPES; do
    if [ -f $PREFIX_ROKKO/trilinos/$arch_type/$build_type/trilinosvars.sh ]; then
      $SUDO rm -f $PREFIX_ROKKO/rokkoenv-$arch_type-$build_type.d/trilinosvars.sh
      $SUDO ln -s $PREFIX_ROKKO/trilinos/$arch_type/$build_type/trilinosvars.sh $PREFIX_ROKKO/rokkoenv-$arch_type-$build_type.d/trilinosvars.sh
    fi
  done

  if [ -f $PREFIX_ROKKO/trilinos/$arch_type/trilinosvars.sh ]; then
    $SUDO rm -f $PREFIX_ROKKO/rokkoenv-$arch_type.d/trilinosvars.sh
    $SUDO ln -s $PREFIX_ROKKO/trilinos/$arch_type/trilinosvars.sh $PREFIX_ROKKO/rokkoenv-$arch_type.d/trilinosvars.sh
  fi
done
