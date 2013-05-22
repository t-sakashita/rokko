mkdir -p rokko
rm -f rokko/CMakeCache.txt
rm -f rokko/CPackConfig.cmake
rm -f rokko/CPackSourceConfig.cmake
rm -rf rokko/CMakeFiles/2.*
cd rokko
cmake /global/home/n0004/n000402/development/rokko -DBOOST_INCLUDE_DIR=/global/nano/alps/boost_1_52_0-r2491 -DCMAKE_CXX_COMPILER=mpiFCCpx -DCMAKE_CXX_FLAGS="-Kfast -Xg -mt" -DCMAKE_C_COMPILER=mpifccpx -DCMAKE_C_FLAGS="-Kfast -Xg -mt" -DCMAKE_Fortran_COMPILER=mpifrtpx -DCMAKE_Fortran_FLAGS="-Kfast -KPIC -mt" -DCMAKE_INSTALL_PREFIX=/global/work/n0004/n000402/build -DELEMENTAL_DIR=$WORK/rokko_lib -DOpenMP_CXX_FLAGS=-Kopenmp -DOpenMP_C_FLAGS=-Kopenmp -DSCALAPACK_LIB="-SCALAPACK -SSL2"
cp -fp CMakeCache.txt ../cmake-rokko.cache
