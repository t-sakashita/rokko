#!/bin/bash -x

wget http://elemental.googlecode.com/files/elemental-0.78-p1.tgz
tar xvf elemental-0.78-p1.tgz
#mkdir -p $HOME/build
#cd $HOME/build
cmake -D CMAKE_CXX_COMPILER=mpicxx -D CMAKE_C_COMPILER=mpicc -D CMAKE_Fortran_COMPILER=mpif90 \
-D MATH_LIBS="-mkl=parallel" \
-D CMAKE_INSTALL_PREFIX=$HOME/installed_elemental_debug \
-D ELEM_EXAMPLES=ON -D ELEM_TESTS=ON \
-D CMAKE_BUILD_TYPE=PureDebug \
..

make all
make install

