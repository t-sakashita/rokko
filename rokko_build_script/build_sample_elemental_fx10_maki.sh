
rm -rf $WORK/build/rokko

cmake ~/development/rokko/ \
-DCMAKE_CXX_COMPILER=mpiFCCpx \
-DCMAKE_C_COMPILER=mpifccpx -DCMAKE_Fortran_COMPILER=mpifrtpx \
-DCMAKE_CXX_FLAGS="-Kfast -Xg -mt" -DCMAKE_C_FLAGS="-Kfast -Xg -mt" \
-DCMAKE_BUILD_TYPE=Debug \
-DELEMENTAL_DIR="$HOME/lib/rokko" \
-DEIGEN3_INCLUDE_DIR="$HOME/lib/eigen-eigen-2249f9c22fe8" \
-DBOOST_INCLUDE_DIR=/global/nano/alps/boost_1_52_0-r2491 \
-DSCALAPACK_LIB="-SCALAPACK -SSL2 --linkfortran"


make diagonalize_elemental VERBOSE=1

