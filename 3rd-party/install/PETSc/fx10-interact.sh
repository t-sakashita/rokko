#!/bin/bash

SCRIPT_DIR=$(cd "$(dirname $0)"; pwd)
. $SCRIPT_DIR/../util.sh
. $SCRIPT_DIR/version.sh
set_prefix

sh $SCRIPT_DIR/setup.sh

unset PETSC_DIR
BUILD_TYPES="Release Debug"
for build_type in $BUILD_TYPES; do
  PREFIX_BACKEND=$PREFIX_ROKKO/petsc-$PETSC_VERSION-$PETSC_PATCH_VERSION/Linux-s64fx/$build_type
  cd $BUILD_DIR
  cp -rp petsc-$PETSC_VERSION petsc-$PETSC_VERSION-build-Linux-s64fx-$build_type
  cd petsc-$PETSC_VERSION-build-Linux-s64fx-$build_type
  if [ $build_type == "Release" ]; then
    check ./configure --prefix=$PREFIX_BACKEND \
      --with-cc="mpifcc" --CFLAGS="-mt -Xg" --COPTFLAGS="-Kfast" \
      --with-cxx="mpiFCC" --CXXFLAGS="-mt -Xg" --CXXOPTFLAGS="-Kfast" \
      --with-fc="mpifrt" --FFLAGS="-Kthreadsafe" --FOPTFLAGS="-Kfast" \
      --LDFLAGS="-lmpi_f77 -lmpi_f90" \
      --with-blas-lapack-lib="-SSL2" \
      --with-x=0 --with-c++-support=1 --with-info=1 --with-debugging=0 --known-mpi-shared-libraries=0 --with-valgrind=0
  else
    check ./configure --prefix=$PREFIX_BACKEND \
      --with-cc="mpifcc" --CFLAGS="-mt -Xg" --COPTFLAGS="-g -O0" \
      --with-cxx="mpiFCC" --CXXFLAGS="-mt -Xg" --CXXOPTFLAGS="-g -O0" \
      --with-fc="mpifrt" --FFLAGS="-Kthreadsafe" --FOPTFLAGS="-g -O0" \
      --LDFLAGS="-lmpi_f77 -lmpi_f90" \
      --with-blas-lapack-lib="-SSL2" \
      --with-x=0 --with-c++-support=1 --with-info=1 --with-debugging=1 --known-mpi-shared-libraries=0 --with-valgrind=0
  fi
  check make
  $SUDO make install
  cat << EOF > $BUILD_DIR/petscvars.sh
export PETSC_ROOT=$PREFIX_ROKKO/petsc-$PETSC_VERSION-$PETSC_PATCH_VERSION/Linux-s64fx
export LD_LIBRARY_PATH=$PREFIX_BACKEND/lib:\$LD_LIBRARY_PATH
EOF
  $SUDO cp -f $BUILD_DIR/petscvars.sh $PREFIX_BACKEND
done

cat << EOF > $BUILD_DIR/petscvars.sh
export PETSC_ROOT=$PREFIX_ROKKO/petsc-$PETSC_VERSION-$PETSC_PATCH_VERSION/Linux-s64fx
EOF
$SUDO cp -f $BUILD_DIR/petscvars.sh $PREFIX_ROKKO/petsc-$PETSC_VERSION-$PETSC_PATCH_VERSION/Linux-s64fx
