#!/bin/bash -x

PREFIX="$1"
test -z "$PREFIX" && PREFIX=$HOME/opt/rokko
echo "PREFIX = $PREFIX"

mkdir -p $HOME/build
cd $HOME/build
rm -rf scalapack_installer_1.0.2
if test -f $HOME/source/scalapack_installer_1.0.2.tgz; then
  tar zxf $HOME/source/scalapack_installer_1.0.2.tgz
else
  wget -O - http://www.netlib.org/scalapack/scalapack_installer.tgz | tar zxf -
fi

cd scalapack_installer_1.0.2
patch -p1 << EOF
*** scalapack_installer_1.0.2.orig/netlib.py    2012-05-02 14:44:38.000000000 +0900
--- scalapack_installer_1.0.2/netlib.py 2013-02-15 11:58:47.000000000 +0900
***************
*** 10,17 ****
    lapacklib   = ""                # the Lapack library
    lapclib     = ""                # the Lapack C interface
    lapackinc   = ""                # the Lapack headers
!   cc          = "cc"              # the C compiler for plasma
!   fc          = "gfortran"        # the Fortran compiler for core_lapack
    ranlib      = ""                # Ranlib
    arflags     = "rc"              # ar flags
    ldflags_c   = ""                # loader flags when main program is in C
--- 10,17 ----
    lapacklib   = ""                # the Lapack library
    lapclib     = ""                # the Lapack C interface
    lapackinc   = ""                # the Lapack headers
!   cc          = "gcc-mp-4.7"              # the C compiler for plasma
!   fc          = "gfortran-mp-4.7"        # the Fortran compiler for core_lapack
    ranlib      = ""                # Ranlib
    arflags     = "rc"              # ar flags
    ldflags_c   = ""                # loader flags when main program is in C
EOF

./setup.py --mpicc=openmpicc --mpif90=openmpif90 --mpiincdir=/opt/local/include/openmpi --lapacklib="-L/opt/local/lib -llapack -lptcblas -lptf77blas -latlas" --ldflags_c="-lgfortran" --mpirun=openmpiexec

mkdir -p $PREFIX/lib
cp -p build/scalapack-2.0.2/libscalapack.a $PREFIX/lib
ln -s libscalapack.a $PREFIX/lib/libblacs.a
