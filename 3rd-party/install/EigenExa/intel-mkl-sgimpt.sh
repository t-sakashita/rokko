#!/bin/bash -x

SCRIPT_DIR=$(cd "$(dirname $0)"; pwd)
. $SCRIPT_DIR/../util.sh
. $SCRIPT_DIR/version.sh
set_prefix

sh $SCRIPT_DIR/setup.sh

cd $BUILD_DIR/EigenExa-$EIGENEXA_VERSION

BUILD_TYPES="Release Debug"
for build_type in $BUILD_TYPES; do
  PREFIX=$PREFIX_ROKKO/eigenexa/eigenexa-$EIGENEXA_VERSION-$EIGENEXA_RK_REVISION/$build_type
  cd $BUILD_DIR
  mkdir -p EigenExa-$EIGENEXA_VERSION-build-$build_type
  cd EigenExa-$EIGENEXA_VERSION-build-$build_type
  check cmake -DCMAKE_BUILD_TYPE=$build_type -DCMAKE_INSTALL_PREFIX=$PREFIX \
      -DCMAKE_C_COMPILER=mpicc -DCMAKE_Fortran_COMPILER=mpif90 \
      -DCMAKE_C_FLAGS="-mt" -DCMAKE_Fortran_FLAGS="-mt" \
      -DMPI_C_COMPILER=mpicc -DMPI_Fortran_COMPILER=mpif90 \
      -DMPI_C_INCLUDE_PATH="/home/app/mpt/mpt-2.14-p11333/include" -DMPI_Fortran_INCLUDE_PATH="/home/app/mpt/mpt-2.14-p11333/include" \
      -DMPI_C_LIBRARIES="-mt" -DMPI_Fortran_LIBRARIES="-mt" \
      -DSCALAPACK_LIB="-lmkl_scalapack_lp64 -lmkl_blacs_sgimpt_lp64 -mkl=parallel" \
      $BUILD_DIR/EigenExa-$EIGENEXA_VERSION
  check make -j2
  $SUDO make install
done

DATE=$(date +%Y%m%d-%H%M%S)
for build_type in $BUILD_TYPES; do
  PREFIX=$PREFIX_ROKKO/eigenexa/eigenexa-$EIGENEXA_VERSION-$EIGENEXA_RK_REVISION/$build_type
  cat << EOF > $BUILD_DIR/eigenexavars.sh
# eigenexa $(basename $0 .sh) $EIGENEXA_VERSION $EIGENEXA_RK_REVISION $DATE
export EIGENEXA_ROOT=$PREFIX_ROKKO/eigenexa/eigenexa-$EIGENEXA_VERSION-$EIGENEXA_RK_REVISION
export LD_LIBRARY_PATH=$PREFIX/lib:\$LD_LIBRARY_PATH
EOF
  $SUDO cp -f $BUILD_DIR/eigenexavars.sh $PREFIX
done

cat << EOF > $BUILD_DIR/eigenexavars.sh
# eigenexa $(basename $0 .sh) $EIGENEXA_VERSION $EIGENEXA_RK_REVISION $DATE
export EIGENEXA_ROOT=$PREFIX_ROKKO/eigenexa/eigenexa-$EIGENEXA_VERSION-$EIGENEXA_RK_REVISION
EOF
$SUDO cp -f $BUILD_DIR/eigenexavars.sh $PREFIX_ROKKO/eigenexa/eigenexa-$EIGENEXA_VERSION-$EIGENEXA_RK_REVISION
