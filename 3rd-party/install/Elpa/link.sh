#!/bin/sh

SCRIPT_DIR=$(cd "$(dirname $0)"; pwd)
. $SCRIPT_DIR/../util.sh
. $SCRIPT_DIR/version.sh
set_prefix

if [ -d $PREFIX_ROKKO/elpa-$ELPA_VERSION-$ELPA_RK_REVISION ]; then
  $SUDO rm -f $PREFIX_ROKKO/elpa
  $SUDO ln -s $PREFIX_ROKKO/elpa-$ELPA_VERSION-$ELPA_RK_REVISION $PREFIX_ROKKO/elpa
fi

BUILD_TYPES="Release Debug"
for build_type in $BUILD_TYPES; do
  if [ -f $PREFIX_ROKKO/elpa/$build_type/elpavars.sh ]; then
    $SUDO rm -f $PREFIX_ROKKO/rokkoenv-$build_type.d/elpavars.sh
    $SUDO ln -s $PREFIX_ROKKO/elpa/$build_type/elpavars.sh $PREFIX_ROKKO/rokkoenv-$build_type.d/elpavars.sh
  fi
done

if [ -f $PREFIX_ROKKO/elpa/elpavars.sh ]; then
  $SUDO rm -f $PREFIX_ROKKO/rokkoenv.d/elpavars.sh
  $SUDO ln -s $PREFIX_ROKKO/elpa/elpavars.sh $PREFIX_ROKKO/rokkoenv.d/elpavars.sh
fi

ARCH_TYPES="Linux-s64fx Linux-x86_64"
for arch_type in $ARCH_TYPES; do
  BUILD_TYPES="Release Debug"
  for build_type in $BUILD_TYPES; do
    if [ -f $PREFIX_ROKKO/elpa/$arch_type/$build_type/elpavars.sh ]; then
      $SUDO rm -f $PREFIX_ROKKO/rokkoenv-$arch_type-$build_type.d/elpavars.sh
      $SUDO ln -s $PREFIX_ROKKO/elpa/$arch_type/$build_type/elpavars.sh $PREFIX_ROKKO/rokkoenv-$arch_type-$build_type.d/elpavars.sh
    fi
  done

  if [ -f $PREFIX_ROKKO/elpa/$arch_type/elpavars.sh ]; then
    $SUDO rm -f $PREFIX_ROKKO/rokkoenv-$arch_type.d/elpavars.sh
    $SUDO ln -s $PREFIX_ROKKO/elpa/$arch_type/elpavars.sh $PREFIX_ROKKO/rokkoenv-$arch_type.d/elpavars.sh
  fi
done
