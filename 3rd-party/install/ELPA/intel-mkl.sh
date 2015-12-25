#!/bin/bash

SCRIPT_DIR=$(cd "$(dirname $0)"; pwd)
. $SCRIPT_DIR/../util.sh
. $SCRIPT_DIR/version.sh
set_prefix

sh $SCRIPT_DIR/setup.sh

BUILD_TYPES="Release Debug"
for build_type in $BUILD_TYPES; do
  PREFIX=$PREFIX_ROKKO/elpa-$ELPA_VERSION-$ELPA_RK_REVISION/$build_type
  cd $BUILD_DIR
  cp -rp elpa-$ELPA_VERSION elpa-$ELPA_VERSION-build-$build_type
  cd elpa-$ELPA_VERSION-build-$build_type
  if [ $build_type == "Release" ]; then
      FLAGS="-O3"
  else
      FLAGS="-g"
  fi
  ./configure SCALAPACK_LDFLAGS="-lmkl_scalapack_lp64 -lmkl_blacs_intelmpi_lp64 -mkl=parallel" \
	      FCFLAGS=$FLAGS CFLAGS=$FLAGS CXXFLAGS=$FLAGS \
	      --enable-openmp  --prefix=$PREFIX
  check make -j4
  $SUDO make install
done

DATE=$(date +%Y%m%d-%H%M%S)
for build_type in $BUILD_TYPES; do
  PREFIX=$PREFIX_ROKKO/elpa-$ELPA_VERSION-$ELPA_RK_REVISION/$build_type
  cat << EOF > $BUILD_DIR/elpavars.sh
# elpa $(basename $0 .sh) $ELPA_VERSION $ELPA_RK_REVISION $DATE
export ELPA_ROOT=$PREFIX_ROKKO/elpa-$ELPA_VERSION-$ELPA_RK_REVISION
export LD_LIBRARY_PATH=$PREFIX/lib:\$LD_LIBRARY_PATH
EOF
  $SUDO cp -f $BUILD_DIR/elpavars.sh $PREFIX
done

cat << EOF > $BUILD_DIR/elpavars.sh
# elpa $(basename $0 .sh) $ELPA_VERSION $ELPA_RK_REVISION $DATE
export ELPA_ROOT=$PREFIX_ROKKO/elpa-$ELPA_VERSION-$ELPA_RK_REVISION
EOF
$SUDO cp -f $BUILD_DIR/elpavars.sh $PREFIX_ROKKO/elpa-$ELPA_VERSION-$ELPA_RK_REVISION
