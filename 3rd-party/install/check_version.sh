#!/bin/sh

SCRIPT_DIR=`dirname $0`
. $SCRIPT_DIR/util.sh
set_prefix

SOLVERS="eigenexa elemental elpa petsc scalapack slepc trilinos"

echo "ROKKO_SOLVER_ROOT: $PREFIX_ROKKO"
awk '$1=="#" && $2=="env" {print}' "$PREFIX_ROKKO/rokkoenv.sh"
for s in $SOLVERS; do
  if [ -f "$PREFIX_ROKKO/rokkoenv.d/${s}vars.sh" ]; then
    awk '$1=="#" && $2==s {print}' s="$s" "$PREFIX_ROKKO/rokkoenv.d/${s}vars.sh"
  fi
done
