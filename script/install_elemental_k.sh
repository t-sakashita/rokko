#!/bin/bash -x

PREFIX="$1"
test -z "$PREFIX" && PREFIX=$HOME/rokko
echo "PREFIX = $PREFIX"

mkdir -p $PREFIX/build
cd $PREFIX/build
rm -rf elemental-0.80
wget -O - http://elemental.googlecode.com/files/elemental-0.80.tgz | tar zxf -
#tar xvf $PREFIX/source/elemental-0.80.tgz
cd elemental-0.80
patch -p1 << EOF
--- elemental-0.80.orig/CMakeLists.txt
+++ elemental-0.80/CMakeLists.txt
@@ -185,7 +185,7 @@
   set(HAVE_F90_INTERFACE FALSE)
   if(CMAKE_Fortran_COMPILER_SUPPORTS_F90)
     include(FortranCInterface)
-    FortranCInterface_VERIFY(CXX)
+# :   FortranCInterface_VERIFY(CXX)
     if(FortranCInterface_VERIFIED_CXX)
       set(HAVE_F90_INTERFACE TRUE)
       FortranCInterface_HEADER(
EOF

cd $PREFIX/build
rm -rf elemental-0.80-build
mkdir elemental-0.80-build
cd elemental-0.80-build

cmake -DCMAKE_C_COMPILER=mpifccpx -DCMAKE_C_FLAGS="-Kfast -Xg -mt" -DCMAKE_CXX_COMPILER=mpiFCCpx -DCMAKE_CXX_FLAGS="-Kfast -Xg -mt" -DCMAKE_Fortran_COMPILER=mpifrtpx -DCMAKE_Fortran_FLAGS="-Kfast -mt" -DMATH_LIBS="-SSL2 --linkfortran" -DSHARED_LIBRARIES=ON -DCMAKE_INSTALL_PREFIX=$PREFIX $PREFIX/build/elemental-0.80

make -j8 2>&1 | tee make.log
make install
