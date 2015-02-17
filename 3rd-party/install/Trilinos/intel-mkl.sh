#!/bin/bash

SCRIPT_DIR=$(cd "$(dirname $0)"; pwd)
. $SCRIPT_DIR/../util.sh
. $SCRIPT_DIR/version.sh
set_prefix

sh $SCRIPT_DIR/setup.sh

BUILD_TYPES="Release Debug"
for build_type in $BUILD_TYPES; do
  PREFIX=$PREFIX_ROKKO/trilinos-$TRILINOS_VERSION-$TRILINOS_PATCH_VERSION/$build_type
  cd $BUILD_DIR
  mkdir -p trilinos-$TRILINOS_VERSION-Source-build-$build_type
  cd trilinos-$TRILINOS_VERSION-Source-build-$build_type
  if [ `which mpicxx > /dev/null 2>&1; echo $?` = 0 ]; then
    check cmake -DCMAKE_BUILD_TYPE="$build_type" -DCMAKE_INSTALL_PREFIX=$PREFIX \
      -DBUILD_SHARED_LIBS=ON -DBUILD_STATIC_LIBS=OFF \
      -DTPL_ENABLE_MPI=ON -DMPI_USE_COMPILER_WRAPPERS=ON \
      -DCMAKE_CXX_COMPILER=mpicxx -DCMAKE_C_COMPILER=mpicc -DCMAKE_Fortran_COMPILER=mpif90 \
      -DBLAS_LIBRARY_DIRS="$MKLROOT/lib/intel64" -DTPL_BLAS_LIBRARIES="-mkl=parallel" \
      -DLAPACK_LIBRARY_DIRS="$MKLROOT/lib/intel64" -DTPL_LAPACK_LIBRARIES="-mkl=parallel" \
      -DTrilinos_ENABLE_Anasazi=ON -DTrilinos_ENABLE_Didasko=ON -DTrilinos_ENABLE_EXAMPLES=ON -DTrilinos_ENABLE_TESTS=ON \
      $BUILD_DIR/trilinos-$TRILINOS_VERSION-Source
  else
    check cmake -DCMAKE_BUILD_TYPE="$build_type" -DCMAKE_INSTALL_PREFIX=$PREFIX \
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
  cat << EOF > $BUILD_DIR/trilinosvars.sh
export TRILINOS_ROOT=$PREFIX_ROKKO/trilinos-$TRILINOS_VERSION-$TRILINOS_PATCH_VERSION
export LD_LIBRARY_PATH=$PREFIX/lib:\$LD_LIBRARY_PATH
EOF
  $SUDO cp -f $BUILD_DIR/trilinosvars.sh $PREFIX
done

cat << EOF > $BUILD_DIR/trilinosvars.sh
export TRILINOS_ROOT=$PREFIX_ROKKO/trilinos-$TRILINOS_VERSION-$TRILINOS_PATCH_VERSION
EOF
$SUDO cp -f $BUILD_DIR/trilinosvars.sh $PREFIX_ROKKO/trilinos-$TRILINOS_VERSION-$TRILINOS_PATCH_VERSION
