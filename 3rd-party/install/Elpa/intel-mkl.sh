#!/bin/bash

SCRIPT_DIR=$(cd "$(dirname $0)"; pwd)
. $SCRIPT_DIR/../util.sh
. $SCRIPT_DIR/version.sh
set_prefix
set_build_dir

. $PREFIX_OPT/env.sh

sh $SCRIPT_DIR/setup.sh

BUILD_TYPES="Release Debug"
for build_type in $BUILD_TYPES; do
  cd $BUILD_DIR
  rm -rf elpa_lib-$ELPA_VERSION-build-$build_type
  mkdir -p elpa_lib-$ELPA_VERSION-build-$build_type
  cd elpa_lib-$ELPA_VERSION-build-$build_type
  if [ `which mpicxx > /dev/null 2>&1; echo $?` = 0 ]; then
    check cmake -DCMAKE_BUILD_TYPE=$build_type -DCMAKE_INSTALL_PREFIX=$PREFIX_ROKKO/$build_type \
      -DCMAKE_C_COMPILER=mpicc -DCMAKE_Fortran_COMPILER=mpif90 \
      -DSCALAPACK_LIB="-lmkl_scalapack_lp64 -lmkl_blacs_intelmpi_lp64 -mkl=parallel" \
      $BUILD_DIR/elpa_lib-$ELPA_VERSION
  else
    check cmake -DCMAKE_BUILD_TYPE=$build_type -DCMAKE_INSTALL_PREFIX=$PREFIX_ROKKO/$build_type \
      -DCMAKE_C_COMPILER=icc -DCMAKE_Fortran_COMPILER=ifort \
      -DSCALAPACK_LIB="-lmkl_scalapack_lp64 -lmkl_blacs_intelmpi_lp64 -mkl=parallel" \
      $BUILD_DIR/elpa_lib-$ELPA_VERSION
  fi
  check make VERBOSE=1 -j4
  $SUDO make install
done
