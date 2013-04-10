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

// 送信・受信のデータ型の構造体の作成
  void create_struct(const rokko::distributed_matrix& mat, MPI_Datatype* global_array_type, MPI_Datatype* local_array_type, MPI_Comm& cart_comm)
{
  int numprocs_cart;
  int coords[2];

  int myrank_cart;

  MPI_Comm_rank(cart_comm, &myrank_cart);
  MPI_Comm_size(cart_comm, &numprocs_cart);
  //MPI_Cart_coords(cart_comm, myrank_cart, 2, coords);

  int m_global = mat.m_global;  int n_global = mat.n_global;  int blockrows = mat.mb;  int blockcols = mat.nb;
  int m_local = mat.m_local;  int n_local = mat.n_local;
  int myrow = mat.myrow;  int mycol = mat.mycol; int nprow = mat.nprow;  int npcol = mat.npcol;

  const int type_block_rows = ( (m_global + nprow * blockrows - 1) / (nprow * blockrows) ) * blockrows;   // 切り上げ
  const int type_block_cols = (n_global + blockcols - 1) / blockcols;   // 切り上げ
  int count_max = type_block_rows * type_block_cols;

  cout << "count_max=" << count_max << endl;

  int*          array_of_blocklengths = new int[count_max];
  MPI_Aint*     array_of_displacements = new MPI_Aint[count_max];
  MPI_Datatype* array_of_types = new MPI_Datatype[count_max];

  for (int proc = 0; proc < numprocs_cart; ++proc) {
    MPI_Cart_coords(cart_comm, proc, 2, coords);
    int proc_row = coords[0], proc_col = coords[1];
    int proc_num_block_rows, proc_num_block_cols, proc_local_matrix_rows, proc_local_matrix_cols, proc_rest_block_rows, proc_rest_block_cols;
    calculate_local_matrix_size(mat, proc_row, proc_col, proc_num_block_rows, proc_num_block_cols, proc_local_matrix_rows, proc_local_matrix_cols, proc_rest_block_rows, proc_rest_block_cols);

    // 送信データの構造体
    int count = 0;

    for (int i=0; i<proc_num_block_rows; ++i) {
      for (int k=0; k<blockrows; ++k) {
	for (int j=0; j<proc_num_block_cols; ++j) {
	  array_of_blocklengths[count] = blockcols;
	  array_of_displacements[count] = ( ((i*nprow + proc_row)*blockrows+k) * n_global + (j * npcol + proc_col) * blockcols ) * sizeof(double);
	  array_of_types[count] = MPI_DOUBLE;
	  if (myrank_cart == 0)
	    cout << "verify: count=" << count << "  length=" << array_of_blocklengths[count] << "  disp="  << (int)array_of_displacements[count] << endl;
	  count++;
	}
	if (proc_rest_block_cols != 0) {
	  array_of_blocklengths[count] = proc_rest_block_cols;
	  array_of_displacements[count] = ( ((i*nprow + proc_row)*blockrows+k) * n_global + (proc_num_block_cols * npcol + proc_col) * blockcols ) * sizeof(double);
          if (myrank_cart == 0)
	    cout << "amari: count=" << count << "  length=" << array_of_blocklengths[count] << "  disp="  << array_of_displacements[count] << endl;
	  array_of_types[count] = MPI_DOUBLE;
	  count++;
	}
      }
    }
    //cout << "before amari: count=" << count << endl;
    // 行のあまり
    for (int k=0; k<proc_rest_block_rows; ++k) {
      cout << "ddddddddddddddddddddddddddddddd" << endl;
      for (int j=0; j<proc_num_block_cols; ++j) {
	array_of_blocklengths[count] = blockcols;
	array_of_displacements[count] = ( ((proc_num_block_rows * nprow + proc_row)*blockrows+k) * n_global + (j * npcol + proc_col) * blockcols ) * sizeof(double);
	array_of_types[count] = MPI_DOUBLE;
	count++;
      }
      if (proc_rest_block_cols != 0) {
	cout << "rest: count=" << count << "  disp="  << array_of_displacements[count] << endl;
	array_of_blocklengths[count] = proc_rest_block_cols;
	array_of_displacements[count] = ( ((proc_num_block_rows * nprow + proc_row)*blockrows+k) * n_global + (proc_num_block_cols * npcol + proc_col) * blockcols ) * sizeof(double);
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
    // 作成した送信データTypeの表示
    if (myrank_cart == 0) {
      for (int i=0; i<count; ++i) {
	//cout << "proc=" << proc << "  count=" << count << " lengthhhhh=" << array_of_blocklengths[i] << endl;
	printf("send Type proc=%d count=%d:  length=%3d  disp=%3d\n", proc, i, array_of_blocklengths[i], (int)array_of_displacements[i]);
      }
    }

    MPI_Type_struct(count, array_of_blocklengths, array_of_displacements, array_of_types, &global_array_type[proc]);
    MPI_Type_commit(&global_array_type[proc]);

    // 受信データの構造体
    count = 0;
    for (int i=0; i<proc_num_block_rows; ++i) {
      for (int k=0; k<blockrows; ++k) {
	for (int j=0; j<proc_num_block_cols; ++j) {
	  array_of_blocklengths[count] = blockcols;
	  array_of_displacements[count] = ( (i*blockrows+k) * proc_local_matrix_cols + j * blockcols ) * sizeof(double);
	  array_of_types[count] = MPI_DOUBLE;
	  count++;
	}
	if (proc_rest_block_cols != 0) {
	  array_of_blocklengths[count] = proc_rest_block_cols;
	  array_of_displacements[count] = ( (i*blockrows+k) * proc_local_matrix_cols + proc_num_block_cols * blockcols ) * sizeof(double);
	  array_of_types[count] = MPI_DOUBLE;
	  count++;
	}
      }
    }
    // 行のあまり
    for (int k=0; k<proc_rest_block_rows; ++k) {
      for (int j=0; j<proc_num_block_cols; ++j) {
	array_of_blocklengths[count] = blockcols;
	array_of_displacements[count] = ( (proc_num_block_rows * blockrows+k) * proc_local_matrix_cols + j * blockcols ) * sizeof(double);
	array_of_types[count] = MPI_DOUBLE;
	++count;
      }
      if (proc_rest_block_cols != 0) {
	array_of_blocklengths[count] = proc_rest_block_cols;
	array_of_displacements[count] = ( (proc_num_block_rows * blockrows+k) * proc_local_matrix_cols + proc_num_block_cols * blockcols ) * sizeof(double);
	array_of_types[count] = MPI_DOUBLE;
	++count;
      }
    }
    //cout << "count=" << count << endl;

    MPI_Barrier(MPI_COMM_WORLD);

    // 作成した受信データTypeの表示
    if (myrank_cart == 0) {
      for (int i=0; i<count; ++i) {
	printf("recv Type proc=%d count=%d:  length=%3d  disp=%3d\n", proc, i, array_of_blocklengths[i], (int)array_of_displacements[i]);
      }
    }

    MPI_Type_struct(count, array_of_blocklengths, array_of_displacements, array_of_types, &local_array_type[proc]);
    MPI_Type_commit(&local_array_type[proc]);
  } // end of "for (int proc = 0; proc < numprocs_cart; ++proc)"

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

  calculate_local_matrix_size(mat, myrow, mycol, local_num_block_rows, local_num_block_cols, local_matrix_rows, local_matrix_cols, local_rest_block_rows, local_rest_block_cols);
  cout << "local_num_block_rows=" << local_num_block_rows << "  local_num_block_cols=" << local_num_block_cols << endl;

  MPI_Datatype* global_array_type = new MPI_Datatype[numprocs_cart];
  MPI_Datatype* local_array_type = new MPI_Datatype[numprocs_cart];

  create_struct(mat, global_array_type, local_array_type, cart_comm);

  int rank_recv = 0;  // プロセスrank_recvに集約
  int sendcount, recvcount;

  for (int proc = 0; proc < numprocs_cart; ++proc) {
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
    ierr = MPI_Sendrecv(local_array, sendcount, local_array_type[proc], dest, 0, global_array, recvcount, global_array_type[proc], source, 0, cart_comm, &status);
    //ierr = 0;

    if (ierr != 0) {
      printf("Problem with Sendrecv (Scatter).\nExiting\n");
      MPI_Abort(MPI_COMM_WORLD,78);
      exit(78);
    }
  } // for (int proc = 0; proc < numprocs_cart; ++proc)

  MPI_Barrier(MPI_COMM_WORLD);
  return ierr;

}

} // namespace rokko

/*
  ~block_cyclic_adaptor()
  {
  cout << "Destructor ~block_cyclic_adaptor()" << endl;
    delete[] global_array_type;
    global_array_type = NULL;
    delete[] local_array_type;
    local_array_type = NULL;

  }

  */

#endif // ROKKO_COLLECTIVE_H


