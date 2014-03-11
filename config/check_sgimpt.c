#include <mpi.h>

#ifdef OMPI_MPI_H // OpenMPI
  nonexistent();
#endif

#ifdef MPICH2 // MPICH2 or Intel MPI
  nonexistent();
#endif

#ifdef MVAPICH2_VERSION //MVAPICH2
  nonexistent();
#endif

int main(int argc, char** argv) {
  return 0;
}

