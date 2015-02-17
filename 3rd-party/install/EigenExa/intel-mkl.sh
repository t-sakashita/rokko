#!/bin/bash -x

SCRIPT_DIR=$(cd "$(dirname $0)"; pwd)
. $SCRIPT_DIR/../util.sh
. $SCRIPT_DIR/version.sh
set_prefix

sh $SCRIPT_DIR/setup.sh

BUILD_TYPES="Release Debug"
for build_type in $BUILD_TYPES; do
  PREFIX=$PREFIX_ROKKO/eigenexa-$EIGENEXA_VERSION-$EIGENEXA_PATCH_VERSION/$build_type
  cd $BUILD_DIR
  mkdir -p EigenExa-$EIGENEXA_VERSION-build-$build_type
  cd EigenExa-$EIGENEXA_VERSION-build-$build_type
  if [ `which mpicxx > /dev/null 2>&1; echo $?` = 0 ]; then
    check cmake -DCMAKE_BUILD_TYPE=$build_type -DCMAKE_INSTALL_PREFIX=$PREFIX \
      -DCMAKE_C_COMPILER=mpicc -DCMAKE_Fortran_COMPILER=mpif90 \
      -DSCALAPACK_LIB="-lmkl_scalapack_lp64 -lmkl_blacs_intelmpi_lp64 -mkl=parallel" \
      $BUILD_DIR/EigenExa-$EIGENEXA_VERSION
  else
    check cmake -DCMAKE_BUILD_TYPE=$build_type -DCMAKE_INSTALL_PREFIX=$PREFIX \
      -DCMAKE_C_COMPILER=icc -DCMAKE_Fortran_COMPILER=ifort \
      -DMPI_C_INCLUDE_PATH="/usr/include" \
      -DMPI_C_LIBRARIES="-lmpi" -DMPI_Fortran_LIBRARIES="-lmpi" \
      -DSCALAPACK_LIB="-lmkl_scalapack_lp64 -lmkl_blacs_intelmpi_lp64 -mkl=parallel" \
      $BUILD_DIR/EigenExa-$EIGENEXA_VERSION
  fi
  check make VERBOSE=1 -j4
  $SUDO make install
  cat << EOF > $BUILD_DIR/eigenexavars.sh
export EIGENEXA_ROOT=$PREFIX_ROKKO/eigenexa-$EIGENEXA_VERSION-$EIGENEXA_PATCH_VERSION
export LD_LIBRARY_PATH=$PREFIX/lib:\$LD_LIBRARY_PATH
EOF
  $SUDO cp -f $BUILD_DIR/eigenexavars.sh $PREFIX
done

cat << EOF > $BUILD_DIR/eigenexavars.sh
export EIGENEXA_ROOT=$PREFIX_ROKKO/eigenexa-$EIGENEXA_VERSION-$EIGENEXA_PATCH_VERSION
EOF
$SUDO cp -f $BUILD_DIR/eigenexavars.sh $PREFIX_ROKKO/eigenexa-$EIGENEXA_VERSION-$EIGENEXA_PATCH_VERSION
