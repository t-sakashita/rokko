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
  mkdir -p trilinos-$TRILINOS_VERSION-Source-build-$build_type
  cd trilinos-$TRILINOS_VERSION-Source-build-$build_type
  check cmake -DCMAKE_BUILD_TYPE="$build_type" -DCMAKE_INSTALL_PREFIX=$PREFIX_ROKKO/$build_type \
    -DBUILD_SHARED_LIBS=ON -DBUILD_STATIC_LIBS=OFF \
    -DTPL_ENABLE_MPI=ON -DMPI_USE_COMPILER_WRAPPERS=ON \
    -DCMAKE_Fortran_COMPILER=mpif90 -DCMAKE_CXX_COMPILER=mpicxx -DCMAKE_C_COMPILER=mpicc \
    -DTPL_BLAS_LIBRARIES="-framework veclib" -DTPL_LAPACK_LIBRARIES="-framework veclib" \
    -DTrilinos_ENABLE_Anasazi=ON -DTrilinos_ENABLE_Didasko=ON -DTrilinos_ENABLE_EXAMPLES=ON -DTrilinos_ENABLE_TESTS=ON \
    -DCMAKE_INSTALL_RPATH="$PREFIX_ROKKO/$build_type/lib" -DCMAKE_SKIP_BUILD_RPATH=ON -DCMAKE_BUILD_WITH_INSTALL_RPATH=OFF -DCMAKE_INSTALL_RPATH_USE_LINK_PATH=ON -DCMAKE_MACOSX_RPATH=1 \
    $BUILD_DIR/trilinos-$TRILINOS_VERSION-Source
  check make VERBOSE=1 -j4
  $SUDO make install
done
