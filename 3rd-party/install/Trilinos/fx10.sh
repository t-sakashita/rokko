#!/bin/bash

SCRIPT_DIR=$(cd "$(dirname $0)"; pwd)
. $SCRIPT_DIR/../util.sh
. $SCRIPT_DIR/version.sh
set_prefix

sh $SCRIPT_DIR/setup.sh

BUILD_TYPES="Release Debug"
for build_type in $BUILD_TYPES; do
  PREFIX_BACKEND=$PREFIX_ROKKO/trilinos-$TRILINOS_VERSION-$TRILINOS_RK_REVISION/Linux-s64fx/$build_type
  cd $BUILD_DIR
  mkdir -p trilinos-$TRILINOS_VERSION-Source-build-Linux-s64fx-$build_type
  cd trilinos-$TRILINOS_VERSION-Source-build-Linux-s64fx-$build_type
  check cmake -DCMAKE_BUILD_TYPE="$build_type" -DCMAKE_INSTALL_PREFIX=$PREFIX_BACKEND \
    -DBUILD_SHARED_LIBS=ON -DBUILD_STATIC_LIBS=OFF \
    -DTPL_ENABLE_MPI=ON \
    -DCMAKE_Fortran_COMPILER=mpifrtpx -DCMAKE_CXX_COMPILER=mpiFCCpx -DCMAKE_C_COMPILER=mpifccpx \
    -DCMAKE_CXX_FLAGS="-Kfast -KPIC -Kopenmp -Xg -mt" -DCMAKE_C_FLAGS="-Kfast -KPIC -Kopenmp -Xg -mt" \
    -DCMAKE_Fortran_FLAGS="-Kfast -KPIC -Cpp"\
    -DTPL_BLAS_LIBRARIES="-SSL2BLAMP --linkfortran" \
    -DTPL_LAPACK_LIBRARIES="-SSL2BLAMP --linkfortran" \
    -DCMAKE_EXE_LINKER_FLAGS="--linkfortran" \
    -DTrilinos_EXTRA_LINK_FLAGS="--linkfortran" \
    -DTrilinos_SKIP_FORTRANCINTERFACE_VERIFY_TEST=ON \
    -DTrilinos_ENABLE_Anasazi=ON -DTrilinos_ENABLE_Didasko=ON -DTrilinos_ENABLE_EXAMPLES=ON -DTrilinos_ENABLE_TESTS=ON \
    $BUILD_DIR/trilinos-$TRILINOS_VERSION-Source
  check make VERBOSE=1 -j4
  $SUDO make install
done

DATE=$(date +%Y%m%d-%H%M%S)
for build_type in $BUILD_TYPES; do
  PREFIX_BACKEND=$PREFIX_ROKKO/trilinos-$TRILINOS_VERSION-$TRILINOS_RK_REVISION/Linux-s64fx/$build_type
  cat << EOF > $BUILD_DIR/trilinosvars.sh
# trilinos $(basename $0 .sh) $TRILINOS_VERSION $TRILINOS_RK_REVISION $DATE
export TRILINOS_ROOT=$PREFIX_ROKKO/trilinos-$TRILINOS_VERSION-$TRILINOS_RK_REVISION/Linux-s64fx
export LD_LIBRARY_PATH=$PREFIX_BACKEND/lib:\$LD_LIBRARY_PATH
EOF
  $SUDO cp -f $BUILD_DIR/trilinosvars.sh $PREFIX_BACKEND
done

cat << EOF > $BUILD_DIR/trilinosvars.sh
# trilinos $(basename $0 .sh) $TRILINOS_VERSION $TRILINOS_RK_REVISION $DATE
export TRILINOS_ROOT=$PREFIX_ROKKO/trilinos-$TRILINOS_VERSION-$TRILINOS_RK_REVISION/Linux-s64fx
EOF
$SUDO cp -f $BUILD_DIR/trilinosvars.sh $PREFIX_ROKKO/trilinos-$TRILINOS_VERSION-$TRILINOS_RK_REVISION/Linux-s64fx
