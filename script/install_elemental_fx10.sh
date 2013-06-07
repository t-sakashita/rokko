#!/bin/bash -x

PREFIX="$1"
test -z "$PREFIX" && PREFIX=$HOME/rokko
echo "PREFIX = $PREFIX"

mkdir -p $HOME/build
cd $HOME/build
rm -rf elemental-0.79-p1
tar xvf $HOME/source/elemental-0.79-p1.tgz
cd elemental-0.79-p1
patch -p1 << EOF
--- elemental-0.79-p1.orig/CMakeLists.txt
+++ elemental-0.79-p1/CMakeLists.txt
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

cd $HOME/build
rm -rf elemental-0.79-p1-build
mkdir elemental-0.79-p1-build
cd elemental-0.79-p1-build

cmake -DCMAKE_C_COMPILER=mpifccpx -DCMAKE_C_FLAGS="-Kfast -Xg -mt" -DCMAKE_CXX_COMPILER=mpiFCCpx -DCMAKE_CXX_FLAGS="-Kfast -Xg -mt" -DCMAKE_Fortran_COMPILER=mpifrtpx -DCMAKE_Fortran_FLAGS="-Kfast -mt" -DMATH_LIBS="-SSL2 --linkfortran" -DSHARED_LIBRARIES=ON -DCMAKE_INSTALL_PREFIX=$PREFIX $HOME/build/elemental-0.79-p1

make -j8 2>&1 | tee make.log
make install
