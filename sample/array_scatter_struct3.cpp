#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

using namespace std;

#define COLS  7
#define ROWS  7

class Distributed_matrix {

public:
  Distributed_matrix();
  int create_struct();
  int gather();
  int scatter();

  int create_global_matrix();
  int create_local_matrix();
  void print_matrix();

private:
  int numprocs;
  int coords[2];

  int myrank_cart, numprocs_cart;

  MPI_Status  status;
  int ierr;

  char* a;  // グローバル行列
  char* b;  // ローカル行列

  int nprows, npcols;
  int myrow, mycol;
  int blockrows, blockcols;
  int local_matrix_rows, local_matrix_cols;
  int num_block_rows, num_block_cols;
  int my_rest_block_rows, my_rest_block_cols;
  int count_max;

  MPI_Comm cart_comm;
  MPI_Datatype* global_dist_array;
  MPI_Datatype* local_dist_array;
};

Distributed_matrix::Distributed_matrix()
{
  int array_of_psizes[2] = {0, 0};  // 2-dimensional grid
  MPI_Comm_size(MPI_COMM_WORLD, &numprocs);  
  MPI_Dims_create(numprocs, 2, array_of_psizes);  // 2次元グリッド
  int periodic[2] = {0, 0};
  nprows = 2;  npcols = 2;   // 例として
  array_of_psizes[0] = nprows;
  array_of_psizes[1] = npcols;
  MPI_Cart_create(MPI_COMM_WORLD, 2, array_of_psizes, periodic, false, &cart_comm);
  
  MPI_Comm_rank(cart_comm, &myrank_cart);
  MPI_Comm_size(cart_comm, &numprocs_cart);
  MPI_Cart_coords(cart_comm, myrank_cart, 2, coords);
  myrow = coords[0];  mycol = coords[1];
  
  blockrows = 2;  //BLOCKROWS/2;  // ブロック長(行方向)任意の長さ
  blockcols = 3;  //BLOCKCOLS/2;  // ブロック長(列方向)任意の長さ
}


// 送信・受信のデータ型の構造体の作成
int Distributed_matrix::create_struct()
{
  global_dist_array = new MPI_Datatype[numprocs_cart];
  local_dist_array = new MPI_Datatype[numprocs_cart];

  for (int proc = 0; proc < numprocs_cart; ++proc) {
    MPI_Cart_coords(cart_comm, proc, 2, coords);
    int proc_row = coords[0], proc_col = coords[1];

    int tmp = ROWS / blockrows;
    num_block_rows = (tmp + nprows-1 - proc_row) / nprows;

    tmp = COLS / blockcols;
    num_block_cols = (tmp + npcols-1 - proc_col) / npcols;
    cout << "num_block_rows=" << num_block_rows << endl;
    cout << "num_block_cols=" << num_block_cols << endl;

    int amari_proc_row = ((ROWS / blockrows) % nprows + 2) % nprows; // 最後のブロックを持つプロセスの次のプロセス  % rank番号は0から始まることに注意．割り切れている場合は，my_rest_block_rows=0となるので問題ない．
    int rest_rows = ROWS % blockrows;
    my_rest_block_rows = 0;
    if (proc_row == amari_proc_row)
      my_rest_block_rows = rest_rows;

    int amari_proc_col = ((COLS / blockcols) % npcols + 2) % npcols; // 最後のブロックを持つプロセスの次のプロセス
    int rest_cols = COLS % blockcols;
    my_rest_block_cols = 0;
    if (proc_col == amari_proc_col)
      my_rest_block_cols = rest_cols;

    cout << "my_rest_block_cols=" << my_rest_block_cols << "  my_rest_block_rows=" << my_rest_block_rows << endl;
    int num_block_cols_r = (COLS + npcols * blockcols - 1) /  (npcols * blockcols);
    cout << "num_block_cols_r=" << num_block_cols_r << endl;

    local_matrix_rows = num_block_rows * blockrows + my_rest_block_rows;
    local_matrix_cols = num_block_cols * blockcols + my_rest_block_cols;
    cout << "local_matrix_rows=" << local_matrix_rows << " local_matrix_cols=" << local_matrix_cols << endl;

    count_max = local_matrix_rows * num_block_cols_r;
    int* array_of_blocklengths = new int[count_max];
    MPI_Aint* array_of_displacements = new MPI_Aint[count_max];
    MPI_Datatype* array_of_types = new MPI_Datatype[count_max];

    // 送信データの構造体
    int count = 0;

    for (int i=0; i<num_block_rows; ++i) {
      for (int k=0; k<blockrows; ++k) {
	for (int j=0; j<num_block_cols; ++j) {
	  array_of_blocklengths[count] = blockcols;
	  array_of_displacements[count] = ((i*nprows + proc_row)*blockrows+k) * COLS + (j * npcols + proc_col) * blockcols;
	  array_of_types[count] = MPI_CHAR;
	  count++;
	}
	if (my_rest_block_cols != 0) {
	  array_of_blocklengths[count] = my_rest_block_cols;
	  array_of_displacements[count] = ((i*nprows + proc_row)*blockrows+k) * COLS + (num_block_cols * npcols + proc_col) * blockcols;
	  cout << "amari: count=" << count << "  length=" << array_of_blocklengths[count] << "  disp="  << array_of_displacements[count] << endl;
	  array_of_types[count] = MPI_CHAR;
	  count++;
	}
      }
    }
    cout << "before amari: count=" << count << endl;
    // 行のあまり
    for (int k=0; k<my_rest_block_rows; ++k) {
      for (int j=0; j<num_block_cols; ++j) {
	array_of_blocklengths[count] = blockcols;
	array_of_displacements[count] = ((num_block_rows * nprows + proc_row)*blockrows+k) * COLS + (j * npcols + proc_col) * blockcols;
	array_of_types[count] = MPI_CHAR;
	count++;
      }
      if (my_rest_block_cols != 0) {
	array_of_blocklengths[count] = my_rest_block_cols;
	array_of_displacements[count] = ((num_block_rows * nprows + proc_row)*blockrows+k) * COLS + (num_block_cols * npcols + proc_col) * blockcols;
	cout << "amari: count=" << count << "  disp="  << array_of_displacements[count] << endl;
	array_of_types[count] = MPI_CHAR;
	count++;
      }
    }

    MPI_Barrier(MPI_COMM_WORLD);
    cout << "myrank" << myrank_cart << "  imano_count=" << count << "   num_block_cols_r=" << num_block_cols_r << "  num_block_rows=" << num_block_rows << "  num_block_cols=" << num_block_cols << endl;
    MPI_Barrier(MPI_COMM_WORLD);

    // 作成した送信データTypeの表示
    if (myrank_cart == 0) {
      for (int i=0; i<count; ++i) {
	printf("send Type count %d:  length=%3d  disp=%3d\n", i, array_of_blocklengths[i], array_of_displacements[i]);
      }
    }

    MPI_Type_struct(count, array_of_blocklengths, array_of_displacements, array_of_types, &global_dist_array[proc]);
    MPI_Type_commit(&global_dist_array[proc]);
 
    // 受信データの構造体
    count = 0;
    for (int i=0; i<num_block_rows; ++i) {
      for (int k=0; k<blockrows; ++k) {
	for (int j=0; j<num_block_cols; ++j) {
	  array_of_blocklengths[count] = blockcols;
	  array_of_displacements[count] = (i*blockrows+k) * local_matrix_cols + j * blockcols;
	  array_of_types[count] = MPI_CHAR;
	  count++;
	}
	if (my_rest_block_cols != 0) {
	  array_of_blocklengths[count] = my_rest_block_cols;
	  array_of_displacements[count] = (i*blockrows+k) * local_matrix_cols + num_block_cols * blockcols;
	  array_of_types[count] = MPI_CHAR;
	  count++;
	}
      }
    }
    // 行のあまり
    for (int k=0; k<my_rest_block_rows; ++k) {
      for (int j=0; j<num_block_cols; ++j) {
	array_of_blocklengths[count] = blockcols;
	array_of_displacements[count] = (num_block_rows * blockrows+k) * local_matrix_cols + j * blockcols;
	array_of_types[count] = MPI_CHAR;
	count++;
      }
      if (my_rest_block_cols != 0) {
	array_of_blocklengths[count] = my_rest_block_cols;
	array_of_displacements[count] = (num_block_rows * blockrows+k) * local_matrix_cols + num_block_cols * blockcols;
	array_of_types[count] = MPI_CHAR;
	count++;
      }
    }
    cout << "count=" << count << endl;
    
    MPI_Barrier(MPI_COMM_WORLD);
  
    // 作成した受信データTypeの表示
    if (myrank_cart == 0) {
      for (int i=0; i<count; ++i) {
	printf("recv Type count %d:  length=%3d  disp=%3d\n", i, array_of_blocklengths[i], array_of_displacements[i]);
      }
    }
    
    MPI_Type_struct(count, array_of_blocklengths, array_of_displacements, array_of_types, &local_dist_array[proc]);
    MPI_Type_commit(&local_dist_array[proc]);
  } // end of "for (int proc = 0; proc < numprocs_cart; ++proc)"

  MPI_Barrier(MPI_COMM_WORLD);
}

//　送受信(Scatter)関数
int Distributed_matrix::scatter()
{
  int sendcount, recvcount;

  int rank_send = 0;  // プロセスrank_sendから分散

  for (int proc = 0; proc < numprocs_cart; ++proc) {
    //cout << endl << endl << "##myrank_cart=" << myrank_cart << "  proc_row=" << proc_row << "  proc_col=" << proc_col << endl;
    
    int rank_recv = proc;
    
    int dest, source;
    if (myrank_cart == rank_recv) {
      source = rank_send;
      recvcount = 1;
    } 
    else {  // 自プロセスが受信者ではない場合(送信者である場合を含む)
      source = MPI_PROC_NULL;
      recvcount = 0;
    }
    if (myrank_cart == rank_send)  {
      dest = rank_recv;
      sendcount = 1;
    }
    else {  // 自プロセスが送信者ではない場合(受信者である場合を含む)
      dest = MPI_PROC_NULL;
      sendcount = 0;
    }
    
    ierr = MPI_Sendrecv(a, sendcount, global_dist_array[proc], dest, 0, b, recvcount, local_dist_array[proc], source, 0, cart_comm, &status);

    if (ierr != 0) {
      printf("Problem with Sendrecv (Scatter).\nExiting\n");
      MPI_Abort(MPI_COMM_WORLD,3);
    }
  } // for (int proc = 0; proc < numprocs_cart; ++proc)
  
  MPI_Barrier(MPI_COMM_WORLD);
  return ierr;
}

// 送受信(gather)ルーチン
int Distributed_matrix::gather()
{
  char c[ROWS*COLS];
  
  if (myrank_cart == 0) {
    for (int ii=0; ii<ROWS*COLS; ++ii) {
      c[ii] = (char)(-1);
    }
  }

  int rank_send = 0;  // プロセスrank_sendに集約
  int sendcount, recvcount;
  for (int proc = 0; proc < numprocs_cart; ++proc) {    
    int rank_recv = proc;

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
    ierr = MPI_Sendrecv(b, sendcount, local_dist_array[proc], source, 0, c, recvcount, global_dist_array[proc], dest, 0, cart_comm, &status);

    if (ierr != 0) {
      printf("Problem with Sendrecv (Scatter).\nExiting\n");
      MPI_Abort(MPI_COMM_WORLD,3);
    }
  } // for (int proc = 0; proc < numprocs_cart; ++proc)
  
  MPI_Barrier(MPI_COMM_WORLD);
  
  /* each proc prints it's local "b" out, in order */
  if (myrank_cart == 0) {
    printf("Global matrix: \n");
    for (int ii=0; ii<ROWS; ++ii) {
      for (int jj=0; jj<COLS; ++jj) {
	printf("%3d ",(int)(unsigned char)c[ii*COLS+jj]);
      }
      printf("\n");
    }
  }
  
  MPI_Barrier(MPI_COMM_WORLD);
  return ierr;
}

int Distributed_matrix::create_global_matrix()
{
  a = new char[ROWS*COLS];
  
  if (myrank_cart == 0) {
    for (int ii=0; ii<ROWS*COLS; ++ii) {
      a[ii] = (char)ii;
    }
  }
}

// ローカル行列の作成
int Distributed_matrix::create_local_matrix()
{
  int tmp = ROWS / blockrows;
  int num_block_rows = (tmp + nprows-1 - myrow) / nprows;
  
  tmp = COLS / blockcols;
  int num_block_cols = (tmp + npcols-1 - mycol) / npcols;
  int amari_proc_row = ((ROWS / blockrows) % nprows + 2) % nprows; // 最後のブロックを持つプロセスの次のプロセス
  int rest_rows = ROWS % blockrows;
  int my_rest_block_rows = 0;
  if (myrow == amari_proc_row)
    my_rest_block_rows = rest_rows;
  int amari_proc_col = ((COLS / blockcols) % npcols + 2) % npcols; // 最後のブロックを持つプロセスの次のプロセス
  int rest_cols = COLS % blockcols;
  int my_rest_block_cols = 0;
  if (mycol == amari_proc_col)
    my_rest_block_cols = rest_cols;
  int local_matrix_rows = num_block_rows * blockrows + my_rest_block_rows;
  int local_matrix_cols = num_block_cols * blockcols + my_rest_block_cols;
  //cout << "myrank=" << myrank_cart << "  local_matrix_cols=" << local_matrix_cols << "  local_matrix_rows=" << local_matrix_rows << endl;
  
  b = new char[local_matrix_rows * local_matrix_cols];
  for (int ii=0; ii<local_matrix_rows * local_matrix_cols; ii++) b[ii] = -1; //rank;
}



void Distributed_matrix::print_matrix()
{
  /* each proc prints it's local "b" out, in order */
  for (int proc=0; proc<numprocs_cart; ++proc) {
    if (proc == myrank_cart) {
      printf("Rank = %d  myrow=%d mycol=%d\n", myrank_cart, myrow, mycol);
      if (myrank_cart == 0) {
	printf("Global matrix: \n");
	for (int ii=0; ii<ROWS; ++ii) {
	  for (int jj=0; jj<COLS; ++jj) {
	    printf("%3d ",(int)(unsigned char)a[ii*COLS+jj]);
	  }
	  printf("\n");
	}
      }
      printf("Local Matrix:\n");
      for (int ii=0; ii<local_matrix_rows; ++ii) {
	for (int jj=0; jj<local_matrix_cols; ++jj) {
	  printf("%3d ",(int)(unsigned char)b[ii*local_matrix_cols+jj]);
	}
	printf("\n");
      }
      printf("\n");
    }
    MPI_Barrier(MPI_COMM_WORLD);
  }
}


int main(int argc, char **argv)
{
  MPI_Init(&argc, &argv);
  
  Distributed_matrix Dist;
  Dist.create_global_matrix();
  Dist.create_local_matrix();
  
  Dist.create_struct();

  Dist.scatter();
  Dist.gather();

  Dist.print_matrix();
  
  
  MPI_Finalize();
  return 0;
}

