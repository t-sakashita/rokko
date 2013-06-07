#!/bin/bash -x

mkdir -p $WORK/build
cd $WORK/build
tar xvf $WORK/source/elemental-0.79-p1.tar
rm -rf build_elemental-0.79-p1
mkdir -p build_elemental-0.79-p1
cd build_elemental-0.79-p1

cmake \
-D CMAKE_CXX_COMPILER=icpc -D CMAKE_C_COMPILER=icc -D CMAKE_Fortran_COMPILER=ifort \
-D MATH_LIBS="-mkl=parallel;-lifcore" \
-D IFCORE_LIB="-lifcore" \
-D CMAKE_INSTALL_PREFIX="$WORK/rokko_lib" \
-D ELEM_EXAMPLES=ON -D ELEM_TESTS=ON \
-DMPI_C_INCLUDE_PATH="/usr/include" -DMPI_CXX_INCLUDE_PATH="/usr/include" \
-DMPI_C_LIBRARIES="-lmpi" -DMPI_CXX_LIBRARIES="-lmpi++;-lmpi" -DMPI_Fortran_LIBRARIES="-lmpi" \
-DSHARED_LIBRARIES=ON \
$WORK/build/elemental-0.79-p1

# The following option is not needed, because "-O3" is automatically specified by make of Elemental.
#-D CMAKE_CXX_FLAGS="-O" -D CMAKE_C_FLAGS="-O" \

# The following option is not needed:
# -D CMAKE_EXE_LINKER_FLAGS="-lmpi" \
# This is because "-lmpi" "-lmpi++" for is automatically specified by make of Elemental.
 #-- Found MPI_C: /usr/lib64/libmpi.so 
 #-- Found MPI_CXX: /usr/lib64/libmpi.so /usr/lib64/libmpi++.so 
 #-- Found MPI_Fortran: /usr/lib64/libmpi.so 

make all
make install

