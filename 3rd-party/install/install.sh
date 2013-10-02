#!/bin/sh

ARCH="$1"
if test -z "$ARCH"; then
    echo "Error: $0 arch [prefix]"
    exit 127
fi

PREFIX="$2"
test -z "$PREFIX" && PREFIX=$HOME/opt

dir=`dirname $0`
SCRIPT_DIR=`cd $dir && pwd`

PACKAGES="scalapack eigenexa elemental elpa anasazi"
for p in $PACKAGES; do
    if test -f "$SCRIPT_DIR/${p}_${ARCH}.sh"; then
	echo "Running: $SCRIPT_DIR/${p}_${ARCH}.sh"
	time sh $SCRIPT_DIR/${p}_${ARCH}.sh $PREFIX
	echo "Done: $SCRIPT_DIR/${p}_${ARCH}.sh"
    fi
done

if test -f "$SCRIPT_DIR/petsc_${ARCH}.sh"; then
    echo "Running: $SCRIPT_DIR/petsc_${ARCH}.sh"
    time sh $SCRIPT_DIR/petsc_${ARCH}.sh $PREFIX
    echo "Done: $SCRIPT_DIR/petsc_${ARCH}.sh"
    echo "Running: $SCRIPT_DIR/slepc.sh"
    time sh $SCRIPT_DIR/slepc.sh $PREFIX
    echo "Done: $SCRIPT_DIR/slepc.sh"
fi
