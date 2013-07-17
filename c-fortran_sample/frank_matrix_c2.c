#include <mpi.h>
#include <stdio.h>
#include "wrapper_c2.h"
#include "frank_matrix_wrapper_c2.h"
#include "timer_wrapper_c2.h"
int main(int argc, char *argv[]) {
  MPI_Init(&argc, &argv);

  if (argc <= 2) {
    fprintf(stderr, "error: %s solver_name matrix_size\n", argv[0]);
    MPI_Abort(MPI_COMM_WORLD, 34);
  }
  void* timer;
  timer = initialize_timer();
  timer_registrate(timer, 1, "diagonalize");
    

  //bin/  rokko::timer timer;
  //timer.registrate( 1, "diagonalize");
  //mkl_set_num_threads(4);
    
  char *solver_name;
  solver_name = argv[1];
  printf("solver_name= %s \n",solver_name);

    
  int dim = atoi(argv[2]);
  printf("dim= %d \n",dim);

  void* solver_;
  solver_ = initialize_solver(solver_name, argc, argv);
  

  MPI_Comm comm = MPI_COMM_WORLD;
  void* g;
  g = initialize_grid_col_major(comm);
  int myrank = grid_get_myrank(g);
  int nprocs = grid_get_nprocs(g);

  const int root = 0;

  void* mat;
  void* w;
  void* Z;
  mat = initialize_distributed_matrix_col_major(dim, dim, g, solver_);
  w = initialize_localized_vector(dim);
  Z = initialize_distributed_matrix_col_major(dim, dim, g, solver_);

  printf("finished matrix generation\n");
  int count;
  for (count=0; count<3; ++count) {
    generate_distributed_matrix_col_major(mat);

    MPI_Barrier(MPI_COMM_WORLD);
    //timer.start(1);
    solver_diagonalize_matrix_col_major(solver_, mat, w, Z, timer);
    MPI_Barrier(MPI_COMM_WORLD);
    //timer.stop(1);
  }
  printf("finished matrix generation frank\n");

#ifndef NDEBUG
  if (myrank == root) {
    printf("Computed Eigenvalues=\n");
    int i;
    for (i = 0; i < dim; i++){
      printf( "%.20f\n", localized_vector_get_element(w,i));
    }
  }
#endif

  if (myrank == 0) {
    printf("num_procs = %i\n", nprocs);
    printf("num_threads = %i\n", omp_get_max_threads());
    printf("solver_name = %s\n", solver_name);
    printf("matrix = frank\n");
    printf("dim = %i\n",dim);
    printf("time = %d\n", timer_get_average(timer, 1));
    //    printf("rokko_version = %s\n", ROKKO_VERSION);
    //    std::cout << "hostname = " << boost::asio::ip::host_name() << std::endl;
    //std::time_t now = std::time(0);
    //std::cout << "date = " << ctime(&now)<< std::endl;
  }

  delete_distributed_matrix_col_major(mat);
  delete_localized_vector(w);
  delete_distributed_matrix_col_major(Z);
  delete_timer(timer);
  delete_solver(solver_);
  MPI_Finalize();
  return 0;
}
