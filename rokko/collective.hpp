#ifndef ROKKO_COLLECTIVE_H
#define ROKKO_COLLECTIVE_H

#include <mpi.h>

// Eigen3に関するヘッダファイル
#include <Eigen/Dense>

namespace rokko {

// ローカル行列のサイズ，ブロック数の計算
void calculate_local_matrix_size(const rokko::distributed_matrix& mat, int proc_row, int proc_col, int& local_num_block_rows, int& local_num_block_cols, int& local_matrix_rows, int& local_matrix_cols, int& local_rest_block_rows, int& local_rest_block_cols)
{
  int m_global = mat.m_global;  int n_global = mat.n_global;  int mb = mat.mb;  int nb = mat.nb;
  int m_local = mat.m_local;  int n_local = mat.n_local;
  int myrow = mat.myrow;  int mycol = mat.mycol; int nprow = mat.nprow;  int npcol = mat.npcol;

  int tmp = m_global / mb;
  local_num_block_rows = (tmp + nprow-1 - proc_row) / nprow;

  tmp = n_global / nb;
  local_num_block_cols = (tmp + npcol-1 - proc_col) / npcol;

  //int rest_block_row = ((m_global / mb) % nprow + 2) % nprow; // 最後のブロックを持つプロセスの次のプロセス
  int rest_block_row = ((m_global / mb) % nprow) % nprow; // 最後のブロックを持つプロセスの次のプロセス
  if (proc_row == rest_block_row)
    local_rest_block_rows = m_global % mb;
  else
    local_rest_block_rows = 0;

  //int rest_block_col = ((n_global / nb) % npcol + 2) % npcol; // 最後のブロックを持つプロセスの次のプロセス
  int rest_block_col = ((n_global / nb) % npcol) % npcol; // 最後のブロックを持つプロセス//の次のプロセス
  if (proc_col == rest_block_col)
    local_rest_block_cols = n_global % nb;
  else
    local_rest_block_cols = 0;

  local_matrix_rows = local_num_block_rows * mb + local_rest_block_rows;
  local_matrix_cols = local_num_block_cols * nb + local_rest_block_cols;
}

// struct of global matrix
  void create_struct_global(const rokko::distributed_matrix& mat, MPI_Datatype& global_array_type, MPI_Comm& cart_comm, int myrank_cart)
{
  int numprocs_cart;
  MPI_Comm_size(cart_comm, &numprocs_cart);

  int coords[2];

  int m_global = mat.m_global;  int n_global = mat.n_global;  int mb = mat.mb;  int nb = mat.nb;
  int m_local = mat.m_local;  int n_local = mat.n_local;
  int myrow = mat.myrow;  int mycol = mat.mycol; int nprow = mat.nprow;  int npcol = mat.npcol;

  MPI_Cart_coords(cart_comm, myrank_cart, 2, coords);
  myrow = coords[0];  mycol = coords[1];
  int num_block_rows, num_block_cols, local_matrix_rows, local_matrix_cols, rest_num_block_rows, rest_num_block_cols;
  calculate_local_matrix_size(mat, myrow, mycol, num_block_rows, num_block_cols, local_matrix_rows, local_matrix_cols, rest_num_block_rows, rest_num_block_cols);

  //const int type_block_rows = ( (m_global + nprow * mb - 1) / (nprow * mb) ) * mb;   // 切り上げ
  const int type_block_rows = (m_local + mb - 1) / mb;   // 切り上げ
  const int type_block_cols = (n_local + nb - 1) / nb;   // 切り上げ
  int count_max = 4 * type_block_rows * type_block_cols;

  //cout << "proc=" << myrank_cart << "count_max=" << count_max << endl;

  cout << "proc=" << myrank_cart << endl;
  cout << "count_max=" << count_max << endl;
  cout << "type_block_rows=" << type_block_rows << endl;
  cout << "num_block_rows=" << num_block_rows << endl;
  cout << "local_matrix_rows=" << local_matrix_rows << endl;
  cout << "rest_num_block_rows=" << rest_num_block_rows << endl;
  cout << "type_block_cols=" << type_block_cols << endl;
  cout << "num_block_cols=" << num_block_cols << endl;
  cout << "local_matrix_cols=" << local_matrix_cols << endl;
  cout << "rest_num_block_cols=" << rest_num_block_cols << endl;

  int*          array_of_blocklengths = new int[count_max];
  MPI_Aint*     array_of_displacements = new MPI_Aint[count_max];
  MPI_Datatype* array_of_types = new MPI_Datatype[count_max];

  int count = 0;

  for (int i=0; i<num_block_cols; ++i) {
    for (int k=0; k<nb; ++k) {
      //for (int j=0; j<num_block_rows; ++j) {
      for (int j=0; j<num_block_rows; ++j) {
	array_of_blocklengths[count] = mb;
	  array_of_displacements[count] = ( ((i*npcol + mycol)*nb+k) * n_global + (j * nprow + myrow) * mb ) * sizeof(double);
	  array_of_types[count] = MPI_DOUBLE;
	  if (myrank_cart == 0)
	    cout << "verify: count=" << count << "  length=" << array_of_blocklengths[count] << "  disp="  << (int)array_of_displacements[count] << endl;
	  ++count;
      }
      if (rest_num_block_rows != 0) {
	array_of_blocklengths[count] = rest_num_block_rows;
	array_of_displacements[count] = ( ((i*npcol + mycol)*nb+k) * n_global + (num_block_rows * nprow + myrow) * mb ) * sizeof(double);
	if (myrank_cart == 0)
	    cout << "amari: count=" << count << "  length=" << array_of_blocklengths[count] << "  disp="  << array_of_displacements[count] << endl;
	array_of_types[count] = MPI_DOUBLE;
	++count;
      }
    }
  }
  //cout << "before amari: count=" << count << endl;
  // 行のあまり
  for (int k=0; k<rest_num_block_cols; ++k) {
    cout << "ddddddddddddddddddddddddddddddd" << endl;
    for (int j=0; j<num_block_rows; ++j) {
      array_of_blocklengths[count] = mb;
      array_of_displacements[count] = ( ((num_block_cols * npcol + mycol)*nb+k) * n_global + (j * nprow + myrow) * mb ) * sizeof(double);
      array_of_types[count] = MPI_DOUBLE;
      ++count;
    }
    if (rest_num_block_rows != 0) {
      cout << "rest: count=" << count << "  disp="  << array_of_displacements[count-1] << endl;
      array_of_blocklengths[count] = rest_num_block_rows;
      array_of_displacements[count] = ( ((num_block_cols * npcol + mycol)*nb+k) * n_global + (num_block_rows * nprow + myrow) * mb ) * sizeof(double);
      if (myrank_cart == 0)
	cout << "rest_amari: count=" << count << "  disp="  << array_of_displacements[count] << endl;
      array_of_types[count] = MPI_DOUBLE;
      ++count;
    }
  }

  //if (myrank_cart == 0)
  //cout << "myrank" << myrank_cart << "  imano_count=" << count << "   num_block_rows_r=" << num_block_rows_r << "  num_block_cols=" << num_block_cols << "  num_block_rows=" << num_block_rows << endl;

  // print out struct of global matrix
  if (myrank_cart == 0) {
    for (int i=0; i<count; ++i) {
	//cout << "proc=" << proc << "  count=" << count << " lengthhhhh=" << array_of_blocklengths[i] << endl;
      printf("global Type proc=%d count=%d:  length=%3d  disp=%3d\n", myrank_cart, i, array_of_blocklengths[i], (int)array_of_displacements[i]/8);
    }
    cout << "num_block_rows=" << num_block_rows << endl;
  }

  MPI_Type_struct(count, array_of_blocklengths, array_of_displacements, array_of_types, &global_array_type);
  MPI_Type_commit(&global_array_type);


  delete[] array_of_blocklengths;
  array_of_blocklengths = NULL;
  delete[] array_of_displacements;
  array_of_displacements = NULL;
  delete[] array_of_types;
  array_of_types = NULL;
}

// struct of localally distributed matrices
void create_struct_local(const rokko::distributed_matrix& mat, MPI_Datatype& local_array_type, MPI_Comm& cart_comm)
{
  int m_global = mat.m_global;  int n_global = mat.n_global;  int mb = mat.mb;  int nb = mat.nb;
  int m_local = mat.m_local;  int n_local = mat.n_local;
  int myrow = mat.myrow;  int mycol = mat.mycol; int nprow = mat.nprow;  int npcol = mat.npcol;

  int numprocs_cart;
  int myrank_cart;
  int coords[2];
  MPI_Comm_rank(cart_comm, &myrank_cart);
  MPI_Comm_size(cart_comm, &numprocs_cart);
  MPI_Cart_coords(cart_comm, myrank_cart, 2, coords);
  myrow = coords[0];  mycol = coords[1];

  int num_block_rows, num_block_cols, local_matrix_rows, local_matrix_cols, rest_num_block_rows, rest_num_block_cols;
  calculate_local_matrix_size(mat, myrow, mycol, num_block_rows, num_block_cols, local_matrix_rows, local_matrix_cols, rest_num_block_rows, rest_num_block_cols);

  const int type_block_rows = ( (m_global + nprow * mb - 1) / (nprow * mb) ) * mb;   // 切り上げ
  const int type_block_cols = (n_local + nb - 1) / nb;   // 切り上げ
  int count_max = 4 * type_block_rows * type_block_cols;

  /*
  if (myrank_cart == 0) {
    cout << "count_max=" << count_max << endl;
    cout << "type_block_rows=" << type_block_rows << endl;
    cout << "num_block_rows=" << num_block_rows << endl;
    cout << "local_matrix_rows=" << local_matrix_rows << endl;
    cout << "rest_num_block_rows=" << rest_num_block_rows << endl;
    cout << "type_block_cols=" << type_block_cols << endl;
    cout << "num_block_cols=" << num_block_cols << endl;
    cout << "local_matrix_cols=" << local_matrix_cols << endl;
    cout << "rest_num_block_cols=" << rest_num_block_cols << endl;
  }
  */

  int*          array_of_blocklengths = new int[count_max];
  MPI_Aint*     array_of_displacements = new MPI_Aint[count_max];
  MPI_Datatype* array_of_types = new MPI_Datatype[count_max];

  int count = 0;
  for (int i=0; i<num_block_cols; ++i) {
    for (int k=0; k<nb; ++k) {
      for (int j=0; j<num_block_rows; ++j) {
	array_of_blocklengths[count] = mb;
 	array_of_displacements[count] = ( (i*nb+k) * local_matrix_rows + j * mb ) * sizeof(double);
	array_of_types[count] = MPI_DOUBLE;
	++count;
      }
      if (rest_num_block_rows != 0) {
	array_of_blocklengths[count] = rest_num_block_rows;
	array_of_displacements[count] = ( (i*nb+k) * local_matrix_rows + num_block_rows * mb ) * sizeof(double);
	array_of_types[count] = MPI_DOUBLE;
	++count;
      }
    }
  }
  // 行のあまり
  for (int k=0; k<rest_num_block_cols; ++k) {
    for (int j=0; j<num_block_rows; ++j) {
      array_of_blocklengths[count] = mb;
      array_of_displacements[count] = ( (num_block_cols * nb+k) * local_matrix_rows + j * mb ) * sizeof(double);
      array_of_types[count] = MPI_DOUBLE;
      ++count;
    }
    if (rest_num_block_rows != 0) {
      array_of_blocklengths[count] = rest_num_block_rows;
      array_of_displacements[count] = ( (num_block_cols * nb+k) * local_matrix_rows + num_block_rows * mb ) * sizeof(double);
      array_of_types[count] = MPI_DOUBLE;
      ++count;
    }
  }
  //cout << "count=" << count << endl;

  MPI_Barrier(MPI_COMM_WORLD);

  // print out struct of localally distributed matrices
  if (myrank_cart == 0) {
    for (int i=0; i<count; ++i) {
      printf("local Type proc=%d count=%d:  length=%3d  disp=%3d\n", myrank_cart, i, array_of_blocklengths[i], (int)array_of_displacements[i]/8);
    }
  }

  MPI_Type_struct(count, array_of_blocklengths, array_of_displacements, array_of_types, &local_array_type);
  MPI_Type_commit(&local_array_type);

  delete[] array_of_blocklengths;
  array_of_blocklengths = NULL;
  delete[] array_of_displacements;
  array_of_displacements = NULL;
  delete[] array_of_types;
  array_of_types = NULL;
}



int gather(const rokko::distributed_matrix& mat, Eigen::MatrixXd& mat_global, int root)
{
  double* global_array;

  if (mat.g.myrank == root) {
    mat_global.resize(mat.m_global, mat.n_global);
    global_array = &mat_global(0,0);  // 本当に、内部では連続な配列になっているか？
    //global_array = mat_global.data();  // 本当に、内部では連続な配列になっているか
    for (int ii=0; ii<mat.m_global * mat.n_global; ++ii) {
      global_array[ii] = ii;
    }
  }

  int numprocs;
  int coords[2];
  int myrank_cart, numprocs_cart;

  MPI_Status  status;
  int ierr;

  int local_matrix_rows, local_matrix_cols;
  int local_num_block_rows, local_num_block_cols;

  int local_rest_block_rows, local_rest_block_cols;
  int count_max;

  MPI_Comm cart_comm;

  int m_global = mat.m_global;  int n_global = mat.n_global;  int mb = mat.mb;  int nb = mat.nb;
  int m_local = mat.m_local;  int n_local = mat.n_local;
  int myrow = mat.myrow;  int mycol = mat.mycol; int nprow = mat.nprow;  int npcol = mat.npcol;
  double* local_array = mat.array;

  int array_of_psizes[2] = {0, 0};  // initial values should be 0
  MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
  MPI_Dims_create(numprocs, 2, array_of_psizes);  // 2-dimensional grid
  int periodic[2] = {0, 0};

  array_of_psizes[0] = nprow;
  array_of_psizes[1] = npcol;
  MPI_Cart_create(MPI_COMM_WORLD, 2, array_of_psizes, periodic, false, &cart_comm);

  MPI_Comm_rank(cart_comm, &myrank_cart);
  MPI_Comm_size(cart_comm, &numprocs_cart);
  MPI_Cart_coords(cart_comm, myrank_cart, 2, coords);
  myrow = coords[0];  mycol = coords[1];

  MPI_Datatype global_array_type;

  int rank_recv = root;  // プロセスrank_recvに集約
  int sendcount = m_local * n_local;
  int recvcount = 1;

  for (int proc = 0; proc < numprocs_cart; ++proc) {
    // Todo: copy routine in root.
    if ((myrank_cart == proc) && (myrank_cart != root)) {
      ierr = MPI_Send(local_array, sendcount, MPI_DOUBLE, root, 0, cart_comm);
      if (ierr != 0) {
	printf("Error with Recv (Gather). ierr=%d\nExiting\n", ierr);
	MPI_Abort(MPI_COMM_WORLD,78);
	exit(78);
      }
    }
    if ((proc != root) &&  (myrank_cart == root)) {
      create_struct_global(mat, global_array_type, cart_comm, proc);
      ierr = MPI_Recv(global_array, recvcount, global_array_type, proc, 0, cart_comm, &status);
      if (ierr != 0) {
	printf("Error with Recv (Gather). ierr=%d\nExiting\n", ierr);
	MPI_Abort(MPI_COMM_WORLD,78);
	exit(78);
      }
      MPI_Type_free(&global_array_type);
      global_array_type = NULL;
    }

    if ((proc == root) && (myrank_cart == root)) {
      create_struct_global(mat, global_array_type, cart_comm, root);
      ierr = MPI_Sendrecv(local_array, sendcount, MPI_DOUBLE, root, 0, global_array, recvcount, global_array_type, root, 0, cart_comm, &status);
      if (ierr != 0) {
	printf("Error with Sendrecv (Gather). ierr=%d\nExiting\n", ierr);
	MPI_Abort(MPI_COMM_WORLD,78);
	exit(78);
      }
      MPI_Type_free(&global_array_type);
      global_array_type = NULL;
    }

  } // for (int proc = 0; proc < numprocs_cart; ++proc)

  MPI_Comm_free(&cart_comm);

  return ierr;
}


int scatter(const rokko::distributed_matrix& mat, Eigen::MatrixXd& mat_global, int root)
{
  double* global_array;


  if (mat.g.myrank == root) {
    // mat_globalは初期化，値の代入後に呼び出されているから，resizeはしない
    //mat_global.resize(mat.m_global, mat.n_global);
    global_array = &mat_global(0,0);  // 本当に、内部では連続な配列になっているか？
    //global_array = mat_global.data();  // 本当に、内部では連続な配列になっているか
    //  for (int ii=0; ii<mat.m_global * mat.n_global; ++ii) {
    //    global_array[ii] = ii;
    //  }
  }

  int numprocs;
  int myrank_cart, numprocs_cart;
  int coords[2];

  MPI_Status  status;
  int ierr;

  int local_matrix_rows, local_matrix_cols;
  int local_num_block_rows, local_num_block_cols;

  int local_rest_block_rows, local_rest_block_cols;
  int count_max;

  MPI_Comm cart_comm;

  int m_global = mat.m_global;  int n_global = mat.n_global;  int mb = mat.mb;  int nb = mat.nb;
  int m_local = mat.m_local;  int n_local = mat.n_local;
  int myrow = mat.myrow;  int mycol = mat.mycol; int nprow = mat.nprow;  int npcol = mat.npcol;
  double* local_array = mat.array;

  int array_of_psizes[2] = {0, 0};  // initial values should be 0
  MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
  MPI_Dims_create(numprocs, 2, array_of_psizes);  // 2-dimensional grid
  int periodic[2] = {0, 0};

  array_of_psizes[0] = nprow;
  array_of_psizes[1] = npcol;
  MPI_Cart_create(MPI_COMM_WORLD, 2, array_of_psizes, periodic, false, &cart_comm);

  MPI_Comm_rank(cart_comm, &myrank_cart);
  MPI_Comm_size(cart_comm, &numprocs_cart);
  MPI_Cart_coords(cart_comm, myrank_cart, 2, coords);
  myrow = coords[0];  mycol = coords[1];

  MPI_Datatype global_array_type;
  int sendcount = 1, recvcount = m_local * n_local;

  for (int proc = 0; proc < numprocs_cart; ++proc) {
    // Todo: copy routine in root.
    if ((proc != root) &&  (myrank_cart == root)) {
      create_struct_global(mat, global_array_type, cart_comm, proc);
      ierr = MPI_Send(global_array, sendcount, global_array_type, proc, 0, cart_comm);
      if (ierr != 0) {
	printf("Error with Recv (Scatter). ierr=%d\nExiting\n", ierr);
	MPI_Abort(MPI_COMM_WORLD,68);
	exit(78);
      }
      MPI_Type_free(&global_array_type);
      global_array_type = NULL;
    }

    if ((myrank_cart == proc) && (myrank_cart != root)) {
      ierr = MPI_Recv(local_array, recvcount, MPI_DOUBLE, root, 0, cart_comm, &status);
      if (ierr != 0) {
	printf("Error with Recv (Scatter). ierr=%d\nExiting\n", ierr);
	MPI_Abort(MPI_COMM_WORLD,69);
	exit(78);
      }
    }

    if ((proc == root) && (myrank_cart == root)) {
      create_struct_global(mat, global_array_type, cart_comm, root);
      ierr = MPI_Sendrecv(global_array, sendcount, global_array_type, root, 0, local_array, recvcount, MPI_DOUBLE, root, 0, cart_comm, &status);
      if (ierr != 0) {
	printf("Error with Sendrecv (Scatter). ierr=%d\nExiting\n", ierr);
	MPI_Abort(MPI_COMM_WORLD,70);
	exit(78);
      }
      MPI_Type_free(&global_array_type);
      global_array_type = NULL;
    }

  } // for (int proc = 0; proc < numprocs_cart; ++proc)

  MPI_Comm_free(&cart_comm);

  return ierr;
}


} // namespace rokko

#endif // ROKKO_COLLECTIVE_H


