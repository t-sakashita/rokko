#ifndef ROKKO_COLLECTIVE_H
#define ROKKO_COLLECTIVE_H

#include <mpi.h>

// Eigen3に関するヘッダファイル
#include <Eigen/Dense>

namespace rokko {

// ローカル行列のサイズ，ブロック数の計算
void calculate_local_matrix_size(const rokko::distributed_matrix& mat, int proc_row, int proc_col, int& local_num_block_rows, int& local_num_block_cols, int& local_matrix_rows, int& local_matrix_cols, int& local_rest_block_rows, int& local_rest_block_cols)
{
  int m_global = mat.m_global;  int n_global = mat.n_global;  int blockrows = mat.mb;  int blockcols = mat.nb;
  int m_local = mat.m_local;  int n_local = mat.n_local;
  int myrow = mat.myrow;  int mycol = mat.mycol; int nprow = mat.nprow;  int npcol = mat.npcol;

  int tmp = m_global / blockrows;
  local_num_block_rows = (tmp + nprow-1 - proc_row) / nprow;

  tmp = n_global / blockcols;
  local_num_block_cols = (tmp + npcol-1 - proc_col) / npcol;
  int rest_block_row = ((m_global / blockrows) % nprow + 2) % nprow; // 最後のブロックを持つプロセスの次のプロセス
  if (proc_row == rest_block_row)
    local_rest_block_rows = m_global % blockrows;
  else
    local_rest_block_rows = 0;

  int rest_block_col = ((n_global / blockcols) % npcol + 2) % npcol; // 最後のブロックを持つプロセスの次のプロセス
  if (proc_col == rest_block_col)
    local_rest_block_cols = n_global % blockcols;
  else
    local_rest_block_cols = 0;

  local_matrix_rows = local_num_block_rows * blockrows + local_rest_block_rows;
  local_matrix_cols = local_num_block_cols * blockcols + local_rest_block_cols;
}

// struct of global matrix
void create_struct_global(const rokko::distributed_matrix& mat, MPI_Datatype& global_array_type, MPI_Comm& cart_comm)
{
  int numprocs_cart;
  int coords[2];

  int myrank_cart;

  int m_global = mat.m_global;  int n_global = mat.n_global;  int blockrows = mat.mb;  int blockcols = mat.nb;
  int m_local = mat.m_local;  int n_local = mat.n_local;
  int myrow = mat.myrow;  int mycol = mat.mycol; int nprow = mat.nprow;  int npcol = mat.npcol;

  MPI_Comm_rank(cart_comm, &myrank_cart);
  MPI_Comm_size(cart_comm, &numprocs_cart);
  MPI_Cart_coords(cart_comm, myrank_cart, 2, coords);
  myrow = coords[0];  mycol = coords[1];
  int num_block_rows, num_block_cols, local_matrix_rows, local_matrix_cols, rest_num_block_rows, rest_num_block_cols;
  calculate_local_matrix_size(mat, myrow, mycol, num_block_rows, num_block_cols, local_matrix_rows, local_matrix_cols, rest_num_block_rows, rest_num_block_cols);

  const int type_block_rows = ( (m_global + nprow * blockrows - 1) / (nprow * blockrows) ) * blockrows;   // 切り上げ
  const int type_block_cols = (n_global + blockcols - 1) / blockcols;   // 切り上げ
  int count_max = type_block_rows * type_block_cols;

  cout << "count_max=" << count_max << endl;

  int*          array_of_blocklengths = new int[count_max];
  MPI_Aint*     array_of_displacements = new MPI_Aint[count_max];
  MPI_Datatype* array_of_types = new MPI_Datatype[count_max];

  int count = 0;

  for (int i=0; i<num_block_rows; ++i) {
    for (int k=0; k<blockrows; ++k) {
      for (int j=0; j<num_block_cols; ++j) {
	array_of_blocklengths[count] = blockcols;
	  array_of_displacements[count] = ( ((i*nprow + myrow)*blockrows+k) * n_global + (j * npcol + mycol) * blockcols ) * sizeof(double);
	  array_of_types[count] = MPI_DOUBLE;
	  if (myrank_cart == 0)
	    cout << "verify: count=" << count << "  length=" << array_of_blocklengths[count] << "  disp="  << (int)array_of_displacements[count] << endl;
	  count++;
      }
      if (rest_num_block_cols != 0) {
	array_of_blocklengths[count] = rest_num_block_cols;
	array_of_displacements[count] = ( ((i*nprow + myrow)*blockrows+k) * n_global + (num_block_cols * npcol + mycol) * blockcols ) * sizeof(double);
	if (myrank_cart == 0)
	    cout << "amari: count=" << count << "  length=" << array_of_blocklengths[count] << "  disp="  << array_of_displacements[count] << endl;
	array_of_types[count] = MPI_DOUBLE;
	count++;
      }
    }
  }
  //cout << "before amari: count=" << count << endl;
  // 行のあまり
  for (int k=0; k<rest_num_block_rows; ++k) {
    cout << "ddddddddddddddddddddddddddddddd" << endl;
    for (int j=0; j<num_block_cols; ++j) {
      array_of_blocklengths[count] = blockcols;
      array_of_displacements[count] = ( ((num_block_rows * nprow + myrow)*blockrows+k) * n_global + (j * npcol + mycol) * blockcols ) * sizeof(double);
      array_of_types[count] = MPI_DOUBLE;
      count++;
    }
    if (rest_num_block_cols != 0) {
      cout << "rest: count=" << count << "  disp="  << array_of_displacements[count] << endl;
      array_of_blocklengths[count] = rest_num_block_cols;
      array_of_displacements[count] = ( ((num_block_rows * nprow + myrow)*blockrows+k) * n_global + (num_block_cols * npcol + mycol) * blockcols ) * sizeof(double);
      if (myrank_cart == 0)
	cout << "rest_amari: count=" << count << "  disp="  << array_of_displacements[count] << endl;
      array_of_types[count] = MPI_DOUBLE;
      count++;
    }
  }

  MPI_Barrier(MPI_COMM_WORLD);
  //if (myrank_cart == 0)
  //cout << "myrank" << myrank_cart << "  imano_count=" << count << "   num_block_cols_r=" << num_block_cols_r << "  num_block_rows=" << num_block_rows << "  num_block_cols=" << num_block_cols << endl;
  MPI_Barrier(MPI_COMM_WORLD);

  //cout << "blockcols=" << blockcols << endl;
  // print out struct of global matrix
  if (myrank_cart == 2) {
    for (int i=0; i<count; ++i) {
	//cout << "proc=" << proc << "  count=" << count << " lengthhhhh=" << array_of_blocklengths[i] << endl;
      printf("global Type proc=%d count=%d:  length=%3d  disp=%3d\n", myrank_cart, i, array_of_blocklengths[i], (int)array_of_displacements[i]);
    }
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
  int numprocs_cart;
  int coords[2];

  int myrank_cart;

  int m_global = mat.m_global;  int n_global = mat.n_global;  int blockrows = mat.mb;  int blockcols = mat.nb;
  int m_local = mat.m_local;  int n_local = mat.n_local;
  int myrow = mat.myrow;  int mycol = mat.mycol; int nprow = mat.nprow;  int npcol = mat.npcol;

  MPI_Comm_rank(cart_comm, &myrank_cart);
  MPI_Comm_size(cart_comm, &numprocs_cart);
  MPI_Cart_coords(cart_comm, myrank_cart, 2, coords);
  myrow = coords[0];  mycol = coords[1];
  int num_block_rows, num_block_cols, local_matrix_rows, local_matrix_cols, rest_num_block_rows, rest_num_block_cols;
  calculate_local_matrix_size(mat, myrow, mycol, num_block_rows, num_block_cols, local_matrix_rows, local_matrix_cols, rest_num_block_rows, rest_num_block_cols);

  const int type_block_rows = ( (m_global + nprow * blockrows - 1) / (nprow * blockrows) ) * blockrows;   // 切り上げ
  const int type_block_cols = (n_global + blockcols - 1) / blockcols;   // 切り上げ
  int count_max = type_block_rows * type_block_cols;

  cout << "count_max=" << count_max << endl;

  int*          array_of_blocklengths = new int[count_max];
  MPI_Aint*     array_of_displacements = new MPI_Aint[count_max];
  MPI_Datatype* array_of_types = new MPI_Datatype[count_max];

  int count = 0;
  for (int i=0; i<num_block_rows; ++i) {
    for (int k=0; k<blockrows; ++k) {
      for (int j=0; j<num_block_cols; ++j) {
	array_of_blocklengths[count] = blockcols;
 	array_of_displacements[count] = ( (i*blockrows+k) * local_matrix_cols + j * blockcols ) * sizeof(double);
	array_of_types[count] = MPI_DOUBLE;
	++count;
      }
      if (rest_num_block_cols != 0) {
	array_of_blocklengths[count] = rest_num_block_cols;
	array_of_displacements[count] = ( (i*blockrows+k) * local_matrix_cols + num_block_cols * blockcols ) * sizeof(double);
	array_of_types[count] = MPI_DOUBLE;
	++count;
      }
    }
  }
  // 行のあまり
  for (int k=0; k<rest_num_block_rows; ++k) {
    for (int j=0; j<num_block_cols; ++j) {
      array_of_blocklengths[count] = blockcols;
      array_of_displacements[count] = ( (num_block_rows * blockrows+k) * local_matrix_cols + j * blockcols ) * sizeof(double);
      array_of_types[count] = MPI_DOUBLE;
      ++count;
    }
    if (rest_num_block_cols != 0) {
      array_of_blocklengths[count] = rest_num_block_cols;
      array_of_displacements[count] = ( (num_block_rows * blockrows+k) * local_matrix_cols + num_block_cols * blockcols ) * sizeof(double);
      array_of_types[count] = MPI_DOUBLE;
      ++count;
    }
  }
  //cout << "count=" << count << endl;

  MPI_Barrier(MPI_COMM_WORLD);

  // print out struct of localally distributed matrices
  if (myrank_cart == 0) {
    for (int i=0; i<count; ++i) {
      printf("local Type proc=%d count=%d:  length=%3d  disp=%3d\n", myrank_cart, i, array_of_blocklengths[i], (int)array_of_displacements[i]);
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
  }

  // variables of class Distributed_Matrix
  int mb, nb;

  // variables of class Grid
  //int myrank, nprocs;
  int ictxt;

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

  int m_global = mat.m_global;  int n_global = mat.n_global;  int blockrows = mat.mb;  int blockcols = mat.nb;
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

  MPI_Datatype global_array_type, local_array_type;

  create_struct_global(mat, global_array_type, cart_comm);
  create_struct_local(mat, local_array_type, cart_comm);

  int rank_recv = root;  // プロセスrank_recvに集約
  int sendcount, recvcount;

  for (int proc = 0; proc < numprocs_cart; ++proc) {
    /*
    int rank_send = proc; // 全プロセス（ルートプロセスも含む）から集約

    int dest, source;
    if (myrank_cart == rank_recv) {
      source = rank_send;
      recvcount = 1;
    }
    else {  // 自プロセスが受信者ではない場合(送信者である場合を含む)
      source = MPI_PROC_NULL;
      recvcount = 0;
    }
    if (myrank_cart == rank_send) {
      dest = rank_recv;
      sendcount = 1;
    }
    else {  // 自プロセスが送信者ではない場合(受信者である場合を含む)
      dest = MPI_PROC_NULL;
      sendcount = 0;
    }
    //cout << "myrank_cart=" << myrank_cart << "  dest=" << dest << "  source=" << source << endl;
    */

    if (myrank_cart == root) {
      ierr = MPI_Recv(global_array, recvcount, global_array_type, proc, proc, cart_comm, &status);
      if (ierr != 0) {
	printf("Error with Recv (Scatter). ierr=%d\nExiting\n", ierr);
	MPI_Abort(MPI_COMM_WORLD,78);
	exit(78);
      }
    }
    if (myrank_cart == proc) {
      ierr = MPI_Send(local_array, sendcount, local_array_type, root, proc, cart_comm);
      if (ierr != 0) {
	printf("Error with Recv (Scatter). ierr=%d\nExiting\n", ierr);
	MPI_Abort(MPI_COMM_WORLD,78);
	exit(78);
      }
    }

    //ierr = MPI_Sendrecv(local_array, sendcount, local_array_type, dest, proc, global_array, recvcount, global_array_type, source, proc, cart_comm, &status);

    //ierr = 0;
    /*
    if (ierr != 0) {
      printf("Error with Sendrecv (Scatter). ierr=%d\nExiting\n", ierr);
      MPI_Abort(MPI_COMM_WORLD,78);
      exit(78);
    }
    */
  } // for (int proc = 0; proc < numprocs_cart; ++proc)

  MPI_Type_free(&local_array_type);
  MPI_Type_free(&global_array_type);

  MPI_Comm_free(&cart_comm);
  MPI_Barrier(MPI_COMM_WORLD);
  return ierr;

}

} // namespace rokko

#endif // ROKKO_COLLECTIVE_H


