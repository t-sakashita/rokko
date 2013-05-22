
rm $WORK/build/rokko/*

#cmake ~/development/rokko/ \
#-DCMAKE_CXX_COMPILER=mpiFCCpx \
#-DCMAKE_C_COMPILER=mpifccpx -DCMAKE_Fortran_COMPILER=mpifrtpx \
#-DCMAKE_CXX_FLAGS="-Kfast -Xg -mt" -DCMAKE_C_FLAGS="-Kfast -Xg -mt" \
#-DCMAKE_BUILD_TYPE=Debug \
#-DELEMENTAL_DIR="$HOME/lib/rokko" \
#-DEIGEN3_INCLUDE_DIR="$HOME/lib/eigen-eigen-2249f9c22fe8" \
#-DBOOST_INCLUDE_DIR=/global/nano/alps/boost_1_52_0-r2491 \
#-DSCALAPACK_LIB="-SCALAPACK -SSL2 --linkfortran"

cmake /global/home/n0004/n000402/development/rokko -DBOOST_INCLUDE_DIR=/global/nano/alps/boost_1_52_0-r2491 -DCMAKE_CXX_COMPILER=mpiFCCpx -DCMAKE_CXX_FLAGS="-Kfast -Xg -mt" -DCMAKE_C_COMPILER=mpifccpx -DCMAKE_C_FLAGS="-Kfast -Xg -mt" -DCMAKE_Fortran_COMPILER=mpifrtpx -DCMAKE_Fortran_FLAGS="-Kfast -KPIC -mt" -DELEMENTAL_DIR="/work/n0004/n000402/rokko_lib" -DOpenMP_CXX_FLAGS=-Kopenmp -DOpenMP_C_FLAGS=-Kopenmp -DSCALAPACK_LIB="-SCALAPACK -SSL2"


make diagonalize_elemental VERBOSE=1

