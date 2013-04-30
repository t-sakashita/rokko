#!/bin/bash -x

cmake -D CMAKE_CXX_COMPILER=mpiFCCpx \
-D CMAKE_C_COMPILER=mpifccpx -D
-D CMAKE_CXX_FLAGS="-Xg -mt" -D CMAKE_C_FLAGS="-Xg -mt" -D MATH_LIBS="-SSL2 --linkfortran" \
-D CMAKE_INSTALL_PREFIX="$HOME/lib/rokko" \
..

#CMAKE_Fortran_COMPILER=mpifrtpx \

make all
make install

