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
  if [ `which mpicxx > /dev/null 2>&1; echo $?` = 0 ]; then
    check cmake -DCMAKE_BUILD_TYPE="$build_type" -DCMAKE_INSTALL_PREFIX=$PREFIX_ROKKO/$build_type \
      -DBUILD_SHARED_LIBS=ON -DBUILD_STATIC_LIBS=OFF \
      -DTPL_ENABLE_MPI=ON -DMPI_USE_COMPILER_WRAPPERS=ON \
      -DCMAKE_CXX_COMPILER=mpicxx -DCMAKE_C_COMPILER=mpicc -DCMAKE_Fortran_COMPILER=mpif90 \
      -DBLAS_LIBRARY_DIRS="$MKLROOT/lib/intel64" -DTPL_BLAS_LIBRARIES="-mkl=parallel" \
      -DLAPACK_LIBRARY_DIRS="$MKLROOT/lib/intel64" -DTPL_LAPACK_LIBRARIES="-mkl=parallel" \
      -DTrilinos_ENABLE_Anasazi=ON -DTrilinos_ENABLE_Didasko=ON -DTrilinos_ENABLE_EXAMPLES=ON -DTrilinos_ENABLE_TESTS=ON \
      $BUILD_DIR/trilinos-$TRILINOS_VERSION-Source
  else
    check cmake -DCMAKE_BUILD_TYPE="$build_type" -DCMAKE_INSTALL_PREFIX=$PREFIX_ROKKO/$build_type \
      -DBUILD_SHARED_LIBS=ON -DBUILD_STATIC_LIBS=OFF \
      -DTPL_ENABLE_MPI=ON -DMPI_USE_COMPILER_WRAPPERS=ON \
      -DCMAKE_CXX_COMPILER=icpc -DCMAKE_C_COMPILER=icc -DCMAKE_Fortran_COMPILER=ifort \
      -DMPI_C_INCLUDE_PATH="/usr/include" -DMPI_CXX_INCLUDE_PATH="/usr/include" \
      -DMPI_C_LIBRARIES="-lmpi" -DMPI_CXX_LIBRARIES="-lmpi++;-lmpi" -DMPI_Fortran_LIBRARIES="-lmpi" \
      -DBLAS_LIBRARY_DIRS="$MKLROOT/lib/intel64" -DTPL_BLAS_LIBRARIES="-mkl=parallel" \
      -DLAPACK_LIBRARY_DIRS="$MKLROOT/lib/intel64" -DTPL_LAPACK_LIBRARIES="-mkl=parallel" \
      -DTrilinos_ENABLE_Anasazi=ON -DTrilinos_ENABLE_Didasko=ON \
      $BUILD_DIR/trilinos-$TRILINOS_VERSION-Source
  fi
  check make VERBOSE=1 -j4
  $SUDO make install
done
