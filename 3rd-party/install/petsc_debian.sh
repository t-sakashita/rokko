#!/bin/bash -x

PREFIX="$1"
test -z "$PREFIX" && PREFIX=$HOME/opt/rokko
echo "PREFIX = $PREFIX"

mkdir -p $HOME/build
cd $HOME/build
rm -rf petsc-3.4.2
if test -f $HOME/source/petsc-3.4.2.tar.gz; then
  tar zxf $HOME/source/petsc-3.4.2.tar.gz
else
  wget -O - http://ftp.mcs.anl.gov/pub/petsc/release-snapshots/petsc-3.4.2.tar.gz | tar zxf -
fi

cd $HOME/build/petsc-3.4.2
unset PETSC_DIR
./configure --prefix=$PREFIX --COPTFLAGS="-O3" --CXXOPTFLAGS="-O3" --FOPTFLAGS="-O3" --with-c++-support=1 --with-debugging=0

make # -j option can not be specified
make install
