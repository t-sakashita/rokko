#!/bin/bash

SCRIPT_DIR=$(cd "$(dirname $0)"; pwd)
. $SCRIPT_DIR/../util.sh
. $SCRIPT_DIR/version.sh
set_prefix

sh $SCRIPT_DIR/setup.sh

BUILD_TYPES="Release Debug"
for build_type in $BUILD_TYPES; do
  PREFIX=$PREFIX_ROKKO/elemental-$ELEMENTAL_VERSION-$ELEMENTAL_PATCH_VERSION/$build_type
  cd $BUILD_DIR
  rm -rf Elemental-$ELEMENTAL_VERSION-build-$build_type
  mkdir -p Elemental-$ELEMENTAL_VERSION-build-$build_type && cd Elemental-$ELEMENTAL_VERSION-build-$build_type
  if [ `which mpicxx > /dev/null 2>&1; echo $?` = 0 ]; then
    check cmake -DCMAKE_BUILD_TYPE="Hybrid$build_type" -DCMAKE_INSTALL_PREFIX=$PREFIX \
      -DELEM_SHARED_LIBRARIES=ON \
      -DCMAKE_CXX_COMPILER=mpicxx -DCMAKE_C_COMPILER=mpicc -DCMAKE_Fortran_COMPILER=mpif90 \
      -DMATH_LIBS="-mkl=parallel;-lifcore" \
    $BUILD_DIR/Elemental-$ELEMENTAL_VERSION
  else
    check cmake -DCMAKE_BUILD_TYPE="Hybrid$build_type" -DCMAKE_INSTALL_PREFIX=$PREFIX_ROKKO \
      -DELEM_SHARED_LIBRARIES=ON \
      -DCMAKE_CXX_COMPILER=icpc -DCMAKE_C_COMPILER=icc -DCMAKE_Fortran_COMPILER=ifort \
      -DMPI_C_INCLUDE_PATH="/usr/include" -DMPI_CXX_INCLUDE_PATH="/usr/include" \
      -DMPI_C_LIBRARIES="-lmpi" -DMPI_CXX_LIBRARIES="-lmpi++;-lmpi" -DMPI_Fortran_LIBRARIES="-lmpi" \
      -DMATH_LIBS="-mkl=parallel;-lifcore" \
    $BUILD_DIR/Elemental-$ELEMENTAL_VERSION
  fi
  check make VERBOSE=1 -j4
  $SUDO make install
  cat << EOF > $BUILD_DIR/elementalvars.sh
export ELEMENTAL_ROOT=$PREFIX_ROKKO/elemental-$ELEMENTAL_VERSION-$ELEMENTAL_PATCH_VERSION
export LD_LIBRARY_PATH=$PREFIX/lib:$LD_LIBRARY_PATH
EOF
  $SUDO cp -f $BUILD_DIR/elementalvars.sh $PREFIX
done

cat << EOF > $BUILD_DIR/elementalvars.sh
export ELEMENTAL_ROOT=$PREFIX_ROKKO/elemental-$ELEMENTAL_VERSION-$ELEMENTAL_PATCH_VERSION
EOF
$SUDO cp -f $BUILD_DIR/elementalvars.sh $PREFIX_ROKKO/elemental-$ELEMENTAL_VERSION-$ELEMENTAL_PATCH_VERSION
