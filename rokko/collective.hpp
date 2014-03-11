#ifndef ROKKO_COLLECTIVE_HPP
#define ROKKO_COLLECTIVE_HPP

#include <mpi.h>

#include <rokko/distributed_matrix.hpp>
#include <rokko/localized_matrix.hpp>

namespace rokko {

// calculation of the size and the block number of distributed_matrix for a given (proc_row, proc_col)
template<typename MATRIX_MAJOR>
void calculate_local_matrix_size(const rokko::distributed_matrix<MATRIX_MAJOR>& mat, int proc_row, int proc_col, int& local_num_block_rows, int& local_num_block_cols, int& local_matrix_rows, int& local_matrix_cols, int& local_rest_block_rows, int& local_rest_block_cols)
{
  int m_global = mat.get_m_global();  int n_global = mat.get_n_global();  int mb = mat.get_mb();  int nb = mat.get_nb();
  int m_local = mat.calculate_row_size(proc_row);
  int n_local = mat.calculate_col_size(proc_col);
  int nprow = mat.get_nprow();  int npcol = mat.get_npcol();

  int tmp = m_global / mb;
  local_num_block_rows = (tmp + nprow-1 - proc_row) / nprow;

  tmp = n_global / nb;
  local_num_block_cols = (tmp + npcol-1 - proc_col) / npcol;

  int rest_block_row = ((m_global / mb) % nprow) % nprow; // 最後のブロックを持つプロセスの次のプロセス
  if (proc_row == rest_block_row)
    local_rest_block_rows = m_global % mb;
  else
    local_rest_block_rows = 0;

  int rest_block_col = ((n_global / nb) % npcol) % npcol; // 最後のブロックを持つプロセス//の次のプロセス
  if (proc_col == rest_block_col)
    local_rest_block_cols = n_global % nb;
  else
    local_rest_block_cols = 0;

  local_matrix_rows = local_num_block_rows * mb + local_rest_block_rows;
  local_matrix_cols = local_num_block_cols * nb + local_rest_block_cols;
}


template<typename MATRIX_MAJOR>
void create_struct_local(const rokko::distributed_matrix<MATRIX_MAJOR>& mat, MPI_Datatype& local_array_type)
{
  int m_global = mat.get_m_global();  int n_global = mat.get_n_global();  int mb = mat.get_mb();  int nb = mat.get_nb();
  int m_local = mat.get_m_local();  int n_local = mat.get_n_local();
  int myrank = mat.get_myrank();
  int myrow = mat.get_myrow();  int mycol = mat.get_mycol(); int nprow = mat.get_nprow();  int npcol = mat.get_npcol();
  int lld = mat.get_lld();

  int count_max = n_local;
  int*          array_of_blocklengths = new int[count_max];
  MPI_Aint*     array_of_displacements = new MPI_Aint[count_max];
  MPI_Datatype* array_of_types = new MPI_Datatype[count_max];

  int count = 0;
  for (int i=0; i<n_local; ++i) {
    array_of_blocklengths[count] = m_local;
    array_of_displacements[count] = i * lld * sizeof(double);
    array_of_types[count] = MPI_DOUBLE;
    //if (myrank == 0)
    //  std::cout << "verify: count=" << count << "  length=" << array_of_blocklengths[count] << "  disp="  << (int)array_of_displacements[count] << std::endl;
    ++count;
  }

#ifndef NDEBUG
  // print out struct of local matrix
  //  if (myrank == 0) {
    for (int i=0; i<count; ++i) {
      //cout << "proc=" << proc << "  count=" << count << " lengthhhhh=" << array_of_blocklengths[i] << std::endl;
      // printf("local Type proc=%d count=%d:  length=%3d  disp=%3d\n", myrank, i, array_of_blocklengths[i], (int)array_of_displacements[i]/8);
    }
    //  }
#endif

  MPI_Type_create_struct(count, array_of_blocklengths, array_of_displacements, array_of_types, &local_array_type);
  MPI_Type_commit(&local_array_type);

  delete[] array_of_blocklengths;
  array_of_blocklengths = 0;
  delete[] array_of_displacements;
  array_of_displacements = 0;
  delete[] array_of_types;
  array_of_types = 0;
}

// struct of global matrix
template<typename MATRIX_MAJOR>
void create_struct_global(const rokko::distributed_matrix<MATRIX_MAJOR>& mat, MPI_Datatype& global_array_type, int proc)
{
  int m_global = mat.get_m_global();  int n_global = mat.get_n_global();  int mb = mat.get_mb();  int nb = mat.get_nb();
  int nprow = mat.get_nprow();  int npcol = mat.get_npcol();

  int myrank = mat.get_myrank();
  int myrow, mycol;
  myrow = mat.get_grid().calculate_grid_row(proc);
  mycol = mat.get_grid().calculate_grid_col(proc);
  int m_local = mat.calculate_row_size(myrow);
  int n_local = mat.calculate_col_size(mycol);


  int num_block_rows, num_block_cols, local_matrix_rows, local_matrix_cols, rest_num_block_rows, rest_num_block_cols;
  calculate_local_matrix_size(mat, myrow, mycol, num_block_rows, num_block_cols, local_matrix_rows, local_matrix_cols, rest_num_block_rows, rest_num_block_cols);

  int type_block_rows = (m_local + mb - 1) / mb;   // 切り上げ
  // std::cout << "myrank=" << myrank << " type_block_rows=" << type_block_rows << " kn_local=" << n_local << std::endl;
  int count_max = type_block_rows * n_local;

  //std::cout << "proc=" << myrank << "count_max=" << count_max << std::endl;
  //cout << "proc=" << myrank << std::endl;
  //cout << "count_max=" << count_max << std::endl;
  //cout << "type_block_rows=" << type_block_rows << std::endl;
  //cout << "num_block_rows=" << num_block_rows << std::endl;
  //cout << "local_matrix_rows=" << local_matrix_rows << std::endl;
  //cout << "rest_num_block_rows=" << rest_num_block_rows << std::endl;
  //cout << "type_block_cols=" << type_block_cols << std::endl;
  //cout << "num_block_cols=" << num_block_cols << std::endl;
  //cout << "local_matrix_cols=" << local_matrix_cols << std::endl;
  //cout << "rest_num_block_cols=" << rest_num_block_cols << std::endl;

  int*          array_of_blocklengths = new int[count_max];
  MPI_Aint*     array_of_displacements = new MPI_Aint[count_max];
  MPI_Datatype* array_of_types = new MPI_Datatype[count_max];

  int count = 0;

  for (int i=0; i<num_block_cols; ++i) {
    for (int k=0; k<nb; ++k) {
      for (int j=0; j<num_block_rows; ++j) {
	array_of_blocklengths[count] = mb;
	  array_of_displacements[count] = ( ((i*npcol + mycol)*nb+k) * m_global + (j * nprow + myrow) * mb ) * sizeof(double);
	  array_of_types[count] = MPI_DOUBLE;
#ifndef NDEBUG
          //	  if (myrank == 0)
	    std::cout << "verify: count=" << count << "  length=" << array_of_blocklengths[count] << "  disp="  << (int)array_of_displacements[count] << std::endl;
#endif          
	  ++count;
      }
      if (rest_num_block_rows != 0) {
	array_of_blocklengths[count] = rest_num_block_rows;
	array_of_displacements[count] = ( ((i*npcol + mycol)*nb+k) * m_global + (num_block_rows * nprow + myrow) * mb ) * sizeof(double);
#ifndef NDEBUG
        //	if (myrank == 0)
	    std::cout << "amari: count=" << count << "  length=" << array_of_blocklengths[count] << "  disp="  << array_of_displacements[count] << std::endl;
#endif
	array_of_types[count] = MPI_DOUBLE;
	++count;
      }
    }
  }
  //cout << "before amari: count=" << count << std::endl;
  // reminder of column
  for (int k=0; k<rest_num_block_cols; ++k) {
    for (int j=0; j<num_block_rows; ++j) {
      array_of_blocklengths[count] = mb;
      array_of_displacements[count] = ( ((num_block_cols * npcol + mycol)*nb+k) * m_global + (j * nprow + myrow) * mb ) * sizeof(double);
      array_of_types[count] = MPI_DOUBLE;
      ++count;
    }
    if (rest_num_block_rows != 0) {
#ifndef NDEBUG
      std::cout << "rest: count=" << count << "  disp="  << array_of_displacements[count-1] << std::endl;
#endif
      array_of_blocklengths[count] = rest_num_block_rows;
      array_of_displacements[count] = ( ((num_block_cols * npcol + mycol)*nb+k) * m_global + (num_block_rows * nprow + myrow) * mb ) * sizeof(double);
#ifndef NDEBUG
      //      if (myrank == 0)
	std::cout << "rest_amari: count=" << count << "  disp="  << array_of_displacements[count] << std::endl;
#endif
      array_of_types[count] = MPI_DOUBLE;
      ++count;
    }
  }

#ifndef NDEBUG
  // print out struct of global matrix
  //  if (myrank == 0) {
    for (int i=0; i<count; ++i) {
      printf("global Type proc=%d count=%d:  length=%3d  disp=%3d\n", proc, i, array_of_blocklengths[i], (int)array_of_displacements[i]/8);
    }
    std::cout << "num_block_rows=" << num_block_rows << std::endl;
    //  }
#endif
    
    MPI_Type_create_struct(count, array_of_blocklengths, array_of_displacements, array_of_types, &global_array_type);
    MPI_Type_commit(&global_array_type);

  // std::cout << "count_max=" << count_max << std::endl;
  delete[] array_of_blocklengths;
  array_of_blocklengths = 0;
  delete[] array_of_displacements;
  array_of_displacements = 0;
  delete[] array_of_types;
  array_of_types = 0;
}

template<typename LOC_MATRIX_MAJOR, typename DIST_MATRIX_MAJOR>
void copy_g2l_root(localized_matrix<LOC_MATRIX_MAJOR>& mat_global, const rokko::distributed_matrix<DIST_MATRIX_MAJOR>& mat)
{
  int m_global = mat.m_global;  int n_global = mat.n_global;  int mb = mat.mb;  int nb = mat.nb;
  int m_local = mat.m_local;  int n_local = mat.n_local;
  int myrow = mat.myrow;  int mycol = mat.mycol; int nprow = mat.nprow;  int npcol = mat.npcol;
  int lld = mat.lld;

  int num_block_rows, num_block_cols, local_matrix_rows, local_matrix_cols, rest_num_block_rows, rest_num_block_cols;
  calculate_local_matrix_size(mat, myrow, mycol, num_block_rows, num_block_cols, local_matrix_rows, local_matrix_cols, rest_num_block_rows, rest_num_block_cols);

  int count = 0;
  for (int i=0; i<num_block_cols; ++i) {
    for (int k=0; k<nb; ++k) {
      for (int j=0; j<num_block_rows; ++j) {
	//array_of_blocklengths[count] = mb;
 	//array_of_displacements[count] = ( (i*nb+k) * lld + j * mb ) * sizeof(double);
	//array_of_types[count] = MPI_DOUBLE;
	for (int c=0; c<mb; ++c) {
	  mat.array[(i*nb+k) * lld + j * mb + c] = mat_global( (j*nprow+myrow)*mb+c, (i*npcol+mycol)*nb+k );
	}
	++count;
      }
      if (rest_num_block_rows != 0) {
	//array_of_blocklengths[count] = rest_num_block_rows;
	//array_of_displacements[count] = ( (i*nb+k) * lld + num_block_rows * mb ) * sizeof(double);
	//array_of_types[count] = MPI_DOUBLE;
	for (int c=0; c<rest_num_block_rows; ++c) {
	  mat.array[(i*nb+k) * lld + num_block_rows * mb + c] = mat_global( (num_block_rows*nprow+myrow)*mb+c, (i*npcol+mycol)*nb+k );
	}
	++count;
      }
    }
  }
  // reminder of column
  for (int k=0; k<rest_num_block_cols; ++k) {
    for (int j=0; j<num_block_rows; ++j) {
      //array_of_blocklengths[count] = mb;
      //array_of_displacements[count] = ( (num_block_cols * nb+k) * lld + j * mb ) * sizeof(double);
      //array_of_types[count] = MPI_DOUBLE;
      for (int c=0; c<mb; ++c) {
	mat.array[(num_block_cols*nb+k) * lld + j*mb+c] = mat_global( (j*nprow+myrow)*mb+c, (num_block_cols*npcol+mycol)*nb+k );
      }
      ++count;
    }
    if (rest_num_block_rows != 0) {
      //array_of_blocklengths[count] = rest_num_block_rows;
      //array_of_displacements[count] = ( (num_block_cols * nb+k) * lld + num_block_rows * mb ) * sizeof(double);
      //array_of_types[count] = MPI_DOUBLE;
      for (int c=0; c<rest_num_block_rows; ++c) {
	mat.array[(num_block_cols*nb+k)*lld + num_block_rows*mb+c] = mat_global( (num_block_rows*nprow+myrow)*mb+c, (num_block_cols*npcol+mycol)*nb+k );
      }
      ++count;
    }
  }
}

template<typename DIST_MATRIX_MAJOR, typename LOC_MATRIX_MAJOR>
void copy_l2g_root(const rokko::distributed_matrix<DIST_MATRIX_MAJOR>& mat, localized_matrix<LOC_MATRIX_MAJOR>& mat_global)
{
  int m_global = mat.m_global;  int n_global = mat.n_global;  int mb = mat.mb;  int nb = mat.nb;
  int m_local = mat.m_local;  int n_local = mat.n_local;
  int myrow = mat.myrow;  int mycol = mat.mycol; int nprow = mat.nprow;  int npcol = mat.npcol;
  int lld = mat.lld;

  int num_block_rows, num_block_cols, local_matrix_rows, local_matrix_cols, rest_num_block_rows, rest_num_block_cols;
  calculate_local_matrix_size(mat, myrow, mycol, num_block_rows, num_block_cols, local_matrix_rows, local_matrix_cols, rest_num_block_rows, rest_num_block_cols);

  int count = 0;
  for (int i=0; i<num_block_cols; ++i) {
    for (int k=0; k<nb; ++k) {
      for (int j=0; j<num_block_rows; ++j) {
	//array_of_blocklengths[count] = mb;
 	//array_of_displacements[count] = ( (i*nb+k) * lld + j * mb ) * sizeof(double);
	//array_of_types[count] = MPI_DOUBLE;
	for (int c=0; c<mb; ++c) {
	  mat_global( (j*nprow+myrow)*mb+c, (i*npcol+mycol)*nb+k ) = mat.array[(i*nb+k) * lld + j * mb + c];
	}
	++count;
      }
      if (rest_num_block_rows != 0) {
	//array_of_blocklengths[count] = rest_num_block_rows;
	//array_of_displacements[count] = ( (i*nb+k) * lld + num_block_rows * mb ) * sizeof(double);
	//array_of_types[count] = MPI_DOUBLE;
	for (int c=0; c<rest_num_block_rows; ++c) {
#ifndef NDEBUG
	  std::cout << "ROW=" << (num_block_rows*nprow+myrow)*mb+c << "  COL=" << (i*npcol+mycol)*nb+k << "rest_num_block_rows=" << rest_num_block_rows << std::endl;
#endif
	  mat_global( (num_block_rows*nprow+myrow)*mb+c, (i*npcol+mycol)*nb+k ) = mat.array[(i*nb+k) * lld + num_block_rows * mb + c];
	}
	++count;
      }
    }
  }
  // reminder of column
  for (int k=0; k<rest_num_block_cols; ++k) {
    for (int j=0; j<num_block_rows; ++j) {
      //array_of_blocklengths[count] = mb;
      //array_of_displacements[count] = ( (num_block_cols * nb+k) * lld + j * mb ) * sizeof(double);
      //array_of_types[count] = MPI_DOUBLE;
      for (int c=0; c<mb; ++c) {
	 mat_global( (j*nprow+myrow)*mb+c, (num_block_cols*npcol+mycol)*nb+k ) = mat.array[(num_block_cols*nb+k) * lld + j*mb+c];
      }
      ++count;
    }
    if (rest_num_block_rows != 0) {
      //array_of_blocklengths[count] = rest_num_block_rows;
      //array_of_displacements[count] = ( (num_block_cols * nb+k) * lld + num_block_rows * mb ) * sizeof(double);
      //array_of_types[count] = MPI_DOUBLE;
      for (int c=0; c<rest_num_block_rows; ++c) {
	 mat_global( (num_block_rows*nprow+myrow)*mb+c, (num_block_cols*npcol+mycol)*nb+k ) = mat.array[(num_block_cols*nb+k)*lld + num_block_rows*mb+c];
      }
      ++count;
    }
  }
}

// gather to pointer specified array
template<typename DIST_MATRIX_MAJOR>
int gather(rokko::distributed_matrix<DIST_MATRIX_MAJOR> const& mat, double* global_array, int root) {
  MPI_Status  status;
  int ierr;

  int local_matrix_rows, local_matrix_cols;
  int local_num_block_rows, local_num_block_cols;

  int local_rest_block_rows, local_rest_block_cols;
  int count_max;

  MPI_Comm cart_comm = MPI_COMM_WORLD;

  int m_global = mat.get_m_global();  int n_global = mat.get_n_global();  int mb = mat.get_mb();  int nb = mat.get_nb();
  int m_local = mat.get_m_local();  int n_local = mat.get_n_local();
  int myrow = mat.get_myrow();  int mycol = mat.get_mycol(); int nprow = mat.get_nprow();  int npcol = mat.get_npcol();
  int myrank = mat.get_myrank(); int nprocs = mat.get_nprocs();

  const double* local_array = mat.get_array_pointer();

  MPI_Datatype local_array_type, global_array_type;

  create_struct_local(mat, local_array_type);

  int rank_recv = root;  // gather to the process specified by rank_recv
  int sendcount = 1;
  int recvcount = 1;

  for (int proc = 0; proc < nprocs; ++proc) {
    if ((myrank == proc) && (myrank != root)) {
      ierr = MPI_Send(const_cast<double*>(local_array), sendcount, local_array_type, root, 0, cart_comm);
      if (ierr != 0) {
	printf("Error with Recv (Gather). ierr=%d\nExiting\n", ierr);
	MPI_Abort(MPI_COMM_WORLD,78);
	exit(78);
      }
    }
    if ((proc != root) &&  (myrank == root)) {
      create_struct_global(mat, global_array_type, proc);
      ierr = MPI_Recv(global_array, recvcount, global_array_type, proc, 0, cart_comm, &status);
      if (ierr != 0) {
	printf("Error with Recv (Gather). ierr=%d\nExiting\n", ierr);
	MPI_Abort(MPI_COMM_WORLD,78);
	exit(78);
      }
      MPI_Type_free(&global_array_type);
    }

    if ((proc == root) && (myrank == root)) {
      create_struct_global(mat, global_array_type, root);
      ierr = MPI_Sendrecv(const_cast<double*>(local_array), sendcount, local_array_type, root, 0, global_array, recvcount, global_array_type, root, 0, cart_comm, &status);
      if (ierr != 0) {
      	printf("Error with Sendrecv (Gather). ierr=%d\nExiting\n", ierr);
      	MPI_Abort(MPI_COMM_WORLD,78);
      	exit(78);
      }
      MPI_Type_free(&global_array_type);
      //copy_l2g_root(mat, mat_global);
    }

  } // for (int proc = 0; proc < nprocs; ++proc)

  MPI_Type_free(&local_array_type);

  return ierr;
}


template<typename DIST_MATRIX_MAJOR, typename LOC_MATRIX_MAJOR>
int gather(rokko::distributed_matrix<DIST_MATRIX_MAJOR> const& mat, localized_matrix<LOC_MATRIX_MAJOR>& mat_global, int root) {
  double* global_array;

  if (mat.get_myrank() == root) {
    mat_global.resize(mat.get_m_global(), mat.get_n_global());
    global_array = &mat_global(0,0);
  }

  return gather(mat, global_array, root);
}

template<typename DIST_MATRIX_MAJOR>
int scatter(double* global_array, distributed_matrix<DIST_MATRIX_MAJOR>& mat, int root) {
  MPI_Status  status;
  int ierr;

  int local_matrix_rows, local_matrix_cols;
  int local_num_block_rows, local_num_block_cols;

  int local_rest_block_rows, local_rest_block_cols;
  int count_max;

  MPI_Comm cart_comm = MPI_COMM_WORLD;

  int m_global = mat.get_m_global();
  int n_global = mat.get_n_global();
  int mb = mat.get_mb();
  int nb = mat.get_nb();
  int m_local = mat.get_m_local();
  int n_local = mat.get_n_local();
  int myrow = mat.get_myrow();
  int mycol = mat.get_mycol();
  int nprow = mat.get_nprow();
  int npcol = mat.get_npcol();
  int myrank = mat.get_myrank();
  int nprocs = mat.get_nprocs();

  double* local_array = mat.get_array_pointer();

  MPI_Datatype local_array_type, global_array_type;
  create_struct_local(mat, local_array_type);

  int sendcount = 1, recvcount = 1;

  for (int proc = 0; proc < nprocs; ++proc) {
    create_struct_global(mat, global_array_type, proc);
    if ((proc != root) &&  (myrank == root)) {
      ierr = MPI_Send(const_cast<double*>(global_array), sendcount, global_array_type, proc, 0, cart_comm);
      if (ierr != 0) {
	printf("Error with Recv (Scatter). ierr=%d\nExiting\n", ierr);
	MPI_Abort(MPI_COMM_WORLD,68);
	exit(78);
      }
      MPI_Type_free(&global_array_type);
    }

    if ((myrank == proc) && (myrank != root)) {
      ierr = MPI_Recv(local_array, recvcount, local_array_type, root, 0, cart_comm, &status);
      if (ierr != 0) {
	printf("Error with Recv (Scatter). ierr=%d\nExiting\n", ierr);
	MPI_Abort(MPI_COMM_WORLD,69);
	exit(78);
      }
    }

    if ((proc == root) && (myrank == root)) {
      create_struct_global(mat, global_array_type, root);
      ierr = MPI_Sendrecv(const_cast<double*>(global_array), sendcount, global_array_type, root, 0, local_array, recvcount, local_array_type, root, 0, cart_comm, &status);
      if (ierr != 0) {
      printf("Error with Sendrecv (Scatter). ierr=%d\nExiting\n", ierr);
      	MPI_Abort(MPI_COMM_WORLD,70);
	exit(78);
      }
      MPI_Type_free(&global_array_type);
      //copy_g2l_root(mat_global, mat);
    }

  } // for (int proc = 0; proc < nprocs; ++proc)

  MPI_Type_free(&local_array_type);

  return ierr;
}

template<typename LOC_MATRIX_MAJOR, typename DIST_MATRIX_MAJOR>
int scatter(localized_matrix<LOC_MATRIX_MAJOR>& mat_global, distributed_matrix<DIST_MATRIX_MAJOR>& mat, int root) {
  double* global_array;

  if (mat.get_myrank() == root) {
    mat_global.resize(mat.get_m_global(), mat.get_n_global());
    global_array = &mat_global(0,0);
  }

  return scatter(global_array, mat, root);
}

} // namespace rokko

#endif // ROKKO_COLLECTIVE_HPP

