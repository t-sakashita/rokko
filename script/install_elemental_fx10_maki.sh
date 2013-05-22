#!/bin/bash -x

cd $WORK/build
tar xvf $WORK/source/elemental-0.79-p1.tar

cd $WORK/build/elemental-0.79-p1
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

rm -rf $WORK/build/elemental-0.79-p1-build
mkdir $WORK/build/elemental-0.79-p1-build
cd $WORK/build/elemental-0.79-p1-build

cmake -D CMAKE_CXX_COMPILER=mpiFCCpx \
-D CMAKE_C_COMPILER=mpifccpx \
-D CMAKE_CXX_FLAGS="-Xg -mt" -D CMAKE_C_FLAGS="-Xg -mt" -D MATH_LIBS="-SSL2 --linkfortran" \
-D CMAKE_INSTALL_PREFIX="$WORK/rokko_lib" \
$WORK/build/elemental-0.79-p1

#CMAKE_Fortran_COMPILER=mpifrtpx \

make all
make install

