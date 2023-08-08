#!/bin/bash -x

SCRIPT_DIR=$(cd "$(dirname $0)"; pwd)
. $SCRIPT_DIR/../util.sh
. $SCRIPT_DIR/version.sh
set_prefix

sh $SCRIPT_DIR/setup.sh

BUILD_TYPES="Release Debug"
for build_type in $BUILD_TYPES; do
  PREFIX=$PREFIX_ROKKO/eigenexa/eigenexa-$EIGENEXA_VERSION-$EIGENEXA_RK_REVISION/$build_type
  cd $BUILD_DIR
  mkdir -p eigenexa-$EIGENEXA_VERSION-build-$build_type
  cd eigenexa-$EIGENEXA_VERSION-build-$build_type
  check cmake -DCMAKE_BUILD_TYPE=$build_type -DCMAKE_INSTALL_PREFIX=$PREFIX \
    -DCMAKE_C_COMPILER=mpicc -DCMAKE_Fortran_COMPILER=mpif90 \
    -DSCALAPACK_LIB="-L/usr/lib/x86_64-linux-gnu  /usr/lib/x86_64-linux-gnu/libscalapack-openmpi.so -llapack -lblas" \
    $BUILD_DIR/eigenexa-$EIGENEXA_VERSION
  check make -j4 VERBOSE=1
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
