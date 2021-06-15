#!/bin/bash -x

SCRIPT_DIR=$(cd "$(dirname $0)"; pwd)
. $SCRIPT_DIR/../util.sh
. $SCRIPT_DIR/version.sh
set_prefix

sh $SCRIPT_DIR/setup.sh

BUILD_TYPES="Release Debug"
for build_type in $BUILD_TYPES; do
  PREFIX_BACKEND=$PREFIX_ROKKO/eigenexa/eigenexa-$EIGENEXA_VERSION-$EIGENEXA_RK_REVISION/Linux-s64fx/$build_type
  cd $BUILD_DIR
  mkdir -p eigenexa-$EIGENEXA_VERSION-build-Linux-s64fx-$build_type
  cd eigenexa-$EIGENEXA_VERSION-build-Linux-s64fx-$build_type
  check cmake -DCMAKE_BUILD_TYPE=$build_type -DCMAKE_INSTALL_PREFIX=$PREFIX_BACKEND \
    -DCMAKE_C_COMPILER=mpifccpx -DCMAKE_Fortran_COMPILER=mpifrtpx \
    -DCMAKE_C_FLAGS="-Kfast -Xg -KPIC" -DCMAKE_Fortran_FLAGS="-Kfast -KPIC -Kocl -Ksimd -KXFILL -Cpp" -DOpenMP_C_FLAGS="-Kopenmp" \
    -DSCALAPACK_LIB="-SCALAPACK -SSL2BLAMP" \
    -DUSE_C_LINKER=ON \
    $BUILD_DIR/eigenexa-$EIGENEXA_VERSION
  check make -j4
  $SUDO make install
done

DATE=$(date +%Y%m%d-%H%M%S)
for build_type in $BUILD_TYPES; do
  PREFIX_BACKEND=$PREFIX_ROKKO/eigenexa/eigenexa-$EIGENEXA_VERSION-$EIGENEXA_RK_REVISION/Linux-s64fx/$build_type
  $SUDO make install
  cat << EOF > $BUILD_DIR/eigenexavars.sh
# eigenexa $(basename $0 .sh) $EIGENEXA_VERSION $EIGENEXA_RK_REVISION $DATE
export EIGENEXA_ROOT=$PREFIX_ROKKO/eigenexa/eigenexa-$EIGENEXA_VERSION-$EIGENEXA_RK_REVISION/Linux-s64fx
export LD_LIBRARY_PATH=$PREFIX_BACKEND/lib:\$LD_LIBRARY_PATH
EOF
  $SUDO cp -f $BUILD_DIR/eigenexavars.sh $PREFIX_BACKEND
done

cat << EOF > $BUILD_DIR/eigenexavars.sh
# eigenexa $(basename $0 .sh) $EIGENEXA_VERSION $EIGENEXA_RK_REVISION $DATE
export EIGENEXA_ROOT=$PREFIX_ROKKO/eigenexa/eigenexa-$EIGENEXA_VERSION-$EIGENEXA_RK_REVISION/Linux-s64fx
EOF
$SUDO cp -f $BUILD_DIR/eigenexavars.sh $PREFIX_ROKKO/eigenexa/eigenexa-$EIGENEXA_VERSION-$EIGENEXA_RK_REVISION/Linux-s64fx
