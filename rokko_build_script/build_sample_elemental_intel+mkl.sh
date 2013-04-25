
cmake ~/development/rokko/  -DCMAKE_CXX_COMPILER="mpicxx" -DCMAKE_Fortran_COMPILER="mpif90" -DCMAKE_BUILD_TYPE=Debug -DELEMENTAL_DIR="/opt/nano/rokko"

make make diagonalize_elemental VERBOSE=1

