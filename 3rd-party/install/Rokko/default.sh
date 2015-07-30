SCRIPT_DIR=$(cd "$(dirname $0)"; pwd)
. $SCRIPT_DIR/../util.sh
. $SCRIPT_DIR/version.sh
set_prefix

ROKKO_SOURCE_DIR=$(cd "$SCRIPT_DIR/../../.."; pwd)

BUILD_TYPES="Release Debug"
for build_type in $BUILD_TYPES; do
  PREFIX=$PREFIX_ROKKO/rokko-$ROKKO_VERSION-$ROKKO_RK_REVISION/$build_type
  cd $BUILD_DIR
  rm -rf Rokko-$ROKKO_VERSION-build-$build_type
  mkdir -p Rokko-$ROKKO_VERSION-build-$build_type && cd Rokko-$ROKKO_VERSION-build-$build_type
  check cmake -DCMAKE_BUILD_TYPE="$build_type" -DCMAKE_INSTALL_PREFIX=$PREFIX \
    -DROKKO_SOLVER_ROOT=$PREFIX_ROKKO \
    $ROKKO_SOURCE_DIR
  check make VERBOSE=1 -j2
  $SUDO make install
done

DATE=$(date +%Y%m%d-%H%M%S)
for build_type in $BUILD_TYPES; do
  PREFIX=$PREFIX_ROKKO/rokko-$ROKKO_VERSION-$ROKKO_RK_REVISION/$build_type
  cat << EOF > $BUILD_DIR/rokkovars.sh
# rokko $(basename $0 .sh) $ROKKO_VERSION $ROKKO_RK_REVISION $DATE
export ROKKO_ROOT=$PREFIX_ROKKO/rokko-$ROKKO_VERSION-$ROKKO_RK_REVISION
export LD_LIBRARY_PATH=$PREFIX/lib:\$LD_LIBRARY_PATH
EOF
  $SUDO cp -f $BUILD_DIR/rokkovars.sh $PREFIX
done

cat << EOF > $BUILD_DIR/rokkovars.sh
# rokko $(basename $0 .sh) $ROKKO_VERSION $ROKKO_RK_REVISION $DATE
export ROKKO_ROOT=$PREFIX_ROKKO/rokko-$ROKKO_VERSION-$ROKKO_RK_REVISION
EOF
$SUDO cp -f $BUILD_DIR/rokkovars.sh $PREFIX_ROKKO/rokko-$ROKKO_VERSION-$ROKKO_RK_REVISION
