#!/bin/bash

tar xvf trilinos-11.0.3-Source.tar.bz2 
patch -u -p1 -d ./trilinos-11.0.3-Source/ < ./trilinos_anasazi_fx10.patch

cd trilinos-11.0.3-Source
rm -Rf ./build/*
mkdir build
cd build

cmake \
-D TPL_ENABLE_MPI:BOOL=ON  \
-D CMAKE_Fortran_COMPILER:FILEPATH="mpifrtpx" -D CMAKE_CXX_COMPILER:FILEPATH="mpiFCCpx" -D CMAKE_C_COMPILER:FILEPATH="mpifccpx" \
-D CMAKE_CXX_FLAGS="-Kopenmp -Xg -mt" \
-D CMAKE_C_FLAGS="-Kopenmp -Xg -mt" \
-D TPL_BLAS_LIBRARIES:STRING="-SSL2BLAMP --linkfortran" \
-D TPL_LAPACK_LIBRARIES:STRING="-SSL2BLAMP --linkfortran" \
-D CMAKE_INSTALL_PREFIX:PATH="../install_build/" \
-D CMAKE_EXE_LINKER_FLAGS="--linkfortran" \
-D Trilinos_EXTRA_LINK_FLAGS:STRING="--linkfortran" \
-D Trilinos_ENABLE_Anasazi:BOOL=ON \
..

make -j2 2>&1 | tee make.log

