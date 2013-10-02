#!/bin/bash -x

PREFIX="$1"
test -z "$PREFIX" && PREFIX=$HOME/opt/rokko
echo "PREFIX = $PREFIX"
DEBUG="$2"
test -z "$DEBUG" && DEBUG=0
echo "DEBUG = $DEGUG"

mkdir -p $HOME/build
cd $HOME/build
rm -rf trilinos-11.4.1-Source trilinos-11.4.1-Source-build
tar jxf $HOME/source/trilinos-11.4.1-Source.tar.bz2

cd $HOME/build
mkdir -p trilinos-11.4.1-Source-build && cd trilinos-11.4.1-Source-build
if [ "$DEBUG" = "0" ]; then
    OPTFLAGS="-O3 -xSSE3"
else
    OPTFLAGS="-O0 -g"
fi
if [ `which mpicxx > /dev/null 2>&1; echo $?` = 0 ]; then
    cmake -DCMAKE_INSTALL_PREFIX=$PREFIX \
	-DBUILD_SHARED_LIBS=ON -DBUILD_STATIC_LIBS=OFF \
	-DTPL_ENABLE_MPI=ON -DMPI_USE_COMPILER_WRAPPERS=ON \
	-DCMAKE_CXX_COMPILER=mpicxx -DCMAKE_C_COMPILER=mpicc -DCMAKE_Fortran_COMPILER=mpif90 \
	-DCMAKE_CXX_FLAGS="$OPTFLAGS" -DCMAKE_C_FLAGS="$OPTFLAGS" -DCMAKE_Fortran_FLAGS="$OPTFLAGS" \
	-DBLAS_LIBRARY_DIRS="$MKLROOT/lib/intel64" -DTPL_BLAS_LIBRARIES="-mkl=parallel" \
	-DLAPACK_LIBRARY_DIRS="$MKLROOT/lib/intel64" -DTPL_LAPACK_LIBRARIES="-mkl=parallel" \
	-DTrilinos_ENABLE_Anasazi=ON \
	-DTrilinos_ENABLE_Didasko=ON \
	$HOME/build/trilinos-11.4.1-Source
else
    cmake -DCMAKE_INSTALL_PREFIX=$PREFIX \
	-DBUILD_SHARED_LIBS=ON -DBUILD_STATIC_LIBS=OFF \
	-DTPL_ENABLE_MPI=ON -DMPI_USE_COMPILER_WRAPPERS=ON \
	-DCMAKE_CXX_COMPILER=icpc -DCMAKE_C_COMPILER=icc -DMPI_Fortran_COMPILER=ifort \
	-DCMAKE_CXX_FLAGS="$OPTFLAGS" -DCMAKE_C_FLAGS="$OPTFLAGS" -DCMAKE_Fortran_FLAGS="$OPTFLAGS" \
      -DMPI_C_INCLUDE_PATH="/usr/include" -DMPI_CXX_INCLUDE_PATH="/usr/include" \
      -DMPI_C_LIBRARIES="-lmpi" -DMPI_CXX_LIBRARIES="-lmpi++;-lmpi" -DMPI_Fortran_LIBRARIES="-lmpi" \
	-DBLAS_LIBRARY_DIRS="$MKLROOT/lib/intel64" -DTPL_BLAS_LIBRARIES="-mkl=parallel" \
	-DLAPACK_LIBRARY_DIRS="$MKLROOT/lib/intel64" -DTPL_LAPACK_LIBRARIES="-mkl=parallel" \
	-DTrilinos_ENABLE_Anasazi=ON \
	-DTrilinos_ENABLE_Didasko=ON \
	$HOME/build/trilinos-11.4.1-Source
fi

make -j4 VERBOSE=1
make install
