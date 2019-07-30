#!/bin/bash -x

SCRIPT_DIR=$(cd "$(dirname $0)"; pwd)
. $SCRIPT_DIR/../util.sh
. $SCRIPT_DIR/version.sh
set_prefix

sh $SCRIPT_DIR/setup.sh

cd $BUILD_DIR/EigenExa-$EIGENEXA_VERSION

BUILD_TYPES="Release Debug"
for build_type in $BUILD_TYPES; do
  PREFIX=$PREFIX_ROKKO/eigenexa-$EIGENEXA_VERSION-$EIGENEXA_RK_REVISION/$build_type
  cd $BUILD_DIR
  mkdir -p EigenExa-$EIGENEXA_VERSION-build-$build_type
  cd EigenExa-$EIGENEXA_VERSION-build-$build_type
  check cmake -DCMAKE_BUILD_TYPE=$build_type -DCMAKE_INSTALL_PREFIX=$PREFIX \
    -DCMAKE_C_COMPILER=mpicc -DCMAKE_Fortran_COMPILER=mpif90 \
    -DSCALAPACK_LIB="-L$PREFIX_ROKKO/scalapack/$build_type/lib -L/opt/local/lib -lscalapack -Wl,-framework -Wl,Accelerate" \
    -DCMAKE_INSTALL_RPATH="$PREFIX/lib" -DCMAKE_SKIP_BUILD_RPATH=OFF -DCMAKE_BUILD_WITH_INSTALL_RPATH=OFF -DCMAKE_INSTALL_RPATH_USE_LINK_PATH=ON -DCMAKE_MACOSX_RPATH=1 \
    $BUILD_DIR/EigenExa-$EIGENEXA_VERSION
  check make VERBOSE=1 -j4
  $SUDO make install
done

DATE=$(date +%Y%m%d-%H%M%S)
for build_type in $BUILD_TYPES; do
  PREFIX=$PREFIX_ROKKO/eigenexa-$EIGENEXA_VERSION-$EIGENEXA_RK_REVISION/$build_type
  cat << EOF > $BUILD_DIR/eigenexavars.sh
# eigenexa $(basename $0 .sh) $EIGENEXA_VERSION $EIGENEXA_RK_REVISION $DATE
export EIGENEXA_ROOT=$PREFIX_ROKKO/eigenexa-$EIGENEXA_VERSION-$EIGENEXA_RK_REVISION
export LD_LIBRARY_PATH=$PREFIX/lib:\$LD_LIBRARY_PATH
EOF
  $SUDO cp -f $BUILD_DIR/eigenexavars.sh $PREFIX
done

cat << EOF > $BUILD_DIR/eigenexavars.sh
# eigenexa $(basename $0 .sh) $EIGENEXA_VERSION $EIGENEXA_RK_REVISION $DATE
export EIGENEXA_ROOT=$PREFIX_ROKKO/eigenexa-$EIGENEXA_VERSION-$EIGENEXA_RK_REVISION
EOF
$SUDO cp -f $BUILD_DIR/eigenexavars.sh $PREFIX_ROKKO/eigenexa-$EIGENEXA_VERSION-$EIGENEXA_RK_REVISION
