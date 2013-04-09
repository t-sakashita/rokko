#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

using namespace std;

#define NDIMS  2
#define COLS  7
#define ROWS  7

int main(int argc, char **argv) {

    MPI_Init(&argc, &argv);
    int p, rank, my_rank, pool_size;
    int array_of_psizes[NDIMS];

    /* Prepare for calling MPI_Type_create_darray */
    //MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    array_of_psizes[0] = 0;   array_of_psizes[1] = 0;
    MPI_Comm_size(MPI_COMM_WORLD, &pool_size);

    MPI_Dims_create(pool_size, NDIMS, array_of_psizes);
    int periodic[2] = {0, 0};
    int nprows = 2, npcols = 2;   // 例として
    array_of_psizes[0] = nprows;
    array_of_psizes[1] = npcols;
    MPI_Comm cart_comm;
    int coords[2];
    MPI_Cart_create(MPI_COMM_WORLD, NDIMS, array_of_psizes, periodic, false, &cart_comm);
    int my_rank_cart, pool_size_cart;

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &p);

    MPI_Comm_rank(cart_comm, &my_rank_cart);
    MPI_Comm_size(cart_comm, &pool_size_cart);
    MPI_Cart_coords(cart_comm, my_rank_cart, NDIMS, coords);


    char a[ROWS*COLS];
    const int BLOCKROWS = ROWS / nprows;  /* number of rows in _block_ */
    const int BLOCKCOLS = COLS / npcols; /* number of cols in _block_ */

    if (rank == 0) {
        for (int ii=0; ii<ROWS*COLS; ++ii) {
            a[ii] = (char)ii;
        }
    }

    int blockrows = 2;  //BLOCKROWS/2;  // ブロック長(行方向)任意の長さ
    int blockcols = 3;  //BLOCKCOLS/2;  // ブロック長(列方向)任意の長さ

    MPI_Datatype * global_dist_array = new MPI_Datatype[pool_size_cart];
    MPI_Datatype * local_dist_array = new MPI_Datatype[pool_size_cart];
    for (int no = 0; no < pool_size_cart; ++no) {
      MPI_Cart_coords(cart_comm, no, NDIMS, coords);
      //int myrow = 0, mycol = 0;
      int myrow = coords[0], mycol = coords[1];
      cout << endl << endl << "!!my_rank_cart=" << my_rank_cart << "  myrow=" << myrow << "  mycol=" << mycol << endl;

    int tmp = ROWS / blockrows;
    int num_block_rows = (tmp + nprows-1 - myrow) / nprows;
    //cout << "tmp=" << tmp << "  num_block_rows=" << num_block_rows << endl;
    //getchar();

    tmp = COLS / blockcols;
    int num_block_cols = (tmp + npcols-1 - mycol) / npcols;
    cout << "num_block_rows=" << num_block_rows << endl;
    cout << "num_block_cols=" << num_block_cols << endl;

    int amari_proc_row = ((ROWS / blockrows) % nprows + 2) % nprows; // 最後のブロックを持つプロセスの次のプロセス  % rank番号は0から始まることに注意．割り切れている場合は，my_rest_block_rows=0となるので問題ない．
    int rest_rows = ROWS % blockrows;
    int my_rest_block_rows = 0;
    if (myrow == amari_proc_row)
      my_rest_block_rows = rest_rows;

    int amari_proc_col = ((COLS / blockcols) % npcols + 2) % npcols; // 最後のブロックを持つプロセスの次のプロセス
    int rest_cols = COLS % blockcols;
    int my_rest_block_cols = 0;
    if (mycol == amari_proc_col)
      my_rest_block_cols = rest_cols;

    cout << "my_rest_block_cols=" << my_rest_block_cols << "  my_rest_block_rows=" << my_rest_block_rows << endl;
    int num_block_cols_r = (COLS + npcols * blockcols - 1) /  (npcols * blockcols);
    cout << "num_block_cols_r=" << num_block_cols_r << endl;

    int local_matrix_rows = num_block_rows * blockrows + my_rest_block_rows;
    int local_matrix_cols = num_block_cols * blockcols + my_rest_block_cols;
    cout << "local_matrix_rows=" << local_matrix_rows << " local_matrix_cols=" << local_matrix_cols << endl;

    int count_max = local_matrix_rows * num_block_cols_r;
    int* array_of_blocklengths = new int[count_max];
    MPI_Aint* array_of_displacements = new MPI_Aint[count_max];
    MPI_Datatype* array_of_types = new MPI_Datatype[count_max];

    // 送信データの構造体
    int count = 0;

    for (int i=0; i<num_block_rows; ++i) {
      for (int k=0; k<blockrows; ++k) {
	for (int j=0; j<num_block_cols; ++j) {
	  array_of_blocklengths[count] = blockcols;
	  array_of_displacements[count] = ((i*nprows + myrow)*blockrows+k) * COLS + (j * npcols + mycol) * blockcols;
	  array_of_types[count] = MPI_CHAR;
	  count++;
	}
	if (my_rest_block_cols != 0) {
	  array_of_blocklengths[count] = my_rest_block_cols;
	  array_of_displacements[count] = ((i*nprows + myrow)*blockrows+k) * COLS + (num_block_cols * npcols + mycol) * blockcols;
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
	array_of_displacements[count] = ((num_block_rows * nprows + myrow)*blockrows+k) * COLS + (j * npcols + mycol) * blockcols;
	array_of_types[count] = MPI_CHAR;
	count++;
      }
      if (my_rest_block_cols != 0) {
	array_of_blocklengths[count] = my_rest_block_cols;
	array_of_displacements[count] = ((num_block_rows * nprows + myrow)*blockrows+k) * COLS + (num_block_cols * npcols + mycol) * blockcols;
	cout << "amari: count=" << count << "  disp="  << array_of_displacements[count] << endl;
	array_of_types[count] = MPI_CHAR;
	count++;
      }
    }

    MPI_Barrier(MPI_COMM_WORLD);
    cout << "rank" << rank << "  imano_count=" << count << "   num_block_cols_r=" << num_block_cols_r << "  num_block_rows=" << num_block_rows << "  num_block_cols=" << num_block_cols << endl;
    MPI_Barrier(MPI_COMM_WORLD);

    // 作成した送信データTypeの表示
    if (rank == 0) {
      for (int i=0; i<count; ++i) {
	printf("send Type count %d:  length=%3d  disp=%3d\n", i, array_of_blocklengths[i], array_of_displacements[i]);
      }
    }

    MPI_Type_struct(count, array_of_blocklengths, array_of_displacements, array_of_types, &global_dist_array[no]);
    MPI_Type_commit(&global_dist_array[no]);

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
    if (rank == 0) {
      for (int i=0; i<count; ++i) {
	printf("recv Type count %d:  length=%3d  disp=%3d\n", i, array_of_blocklengths[i], array_of_displacements[i]);
      }
    }

    MPI_Type_struct(count, array_of_blocklengths, array_of_displacements, array_of_types, &local_dist_array[no]);
    MPI_Type_commit(&local_dist_array[no]);

    } // end of "for (int no = 0; no < pool_size_cart; ++no)"

    MPI_Barrier(MPI_COMM_WORLD);

    // ローカル行列の作成
    MPI_Cart_coords(cart_comm, my_rank_cart, NDIMS, coords);
    int myrow = coords[0], mycol = coords[1];
    int tmp = ROWS / blockrows;
    int num_block_rows = (tmp + nprows-1 - myrow) / nprows;

    //cout << "myrank=" << my_rank_cart << " num_block_rows0=" << num_block_rows << endl;
    tmp = COLS / blockcols;
    int num_block_cols = (tmp + npcols-1 - mycol) / npcols;
    //cout << "myrank=" << my_rank_cart << " num_block_cols0=" << num_block_cols << endl;
    //cout << "myrank=" << my_rank_cart << " num_block_rows++=" << num_block_rows << endl;
    //cout << "myrank=" << my_rank_cart << " num_block_cols++=" << num_block_cols << endl;
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
    //cout << "myrank=" << my_rank_cart << "  my_rest_block_cols=" << my_rest_block_cols << "  my_rest_block_rows=" << my_rest_block_rows << endl;
    int local_matrix_rows = num_block_rows * blockrows + my_rest_block_rows;
    int local_matrix_cols = num_block_cols * blockcols + my_rest_block_cols;
    //cout << "myrank=" << my_rank_cart << "  local_matrix_cols=" << local_matrix_cols << "  local_matrix_rows=" << local_matrix_rows << endl;

    char* b = new char[local_matrix_rows * local_matrix_cols];
    for (int ii=0; ii<local_matrix_rows * local_matrix_cols; ii++) b[ii] = -1; //rank;

    MPI_Status  status;
    int ierr;

    //　送受信(Scatter)
    for (int no = 0; no < pool_size_cart; ++no) {
      MPI_Cart_coords(cart_comm, no, NDIMS, coords);
      //int myrow = 0, mycol = 0;
      int myrow = coords[0], mycol = coords[1];
      //cout << endl << endl << "##my_rank_cart=" << my_rank_cart << "  myrow=" << myrow << "  mycol=" << mycol << endl;

    int rank_send = 0;
    int rank_recv = no;

    int sendcount = 1; //BLOCKROWS/2;
    int recvcount = 1; //BLOCKROWS/2;
    int dest, source;
    if (rank == rank_recv) {
      source = rank_send;
    } 
    else {  // 自プロセスが受信者ではない場合(送信者である場合を含む)
      source = MPI_PROC_NULL;
      recvcount = 0;
    }
    if (rank == rank_send)  dest = rank_recv; 
    else {  // 自プロセスが送信者ではない場合(受信者である場合を含む)
      dest = MPI_PROC_NULL;
      sendcount = 0;
    }

    ierr = MPI_Sendrecv(a, sendcount, global_dist_array[no], dest, 0, b, recvcount, local_dist_array[no], source, 0, MPI_COMM_WORLD, &status);
    //ierr = MPI_Sendrecv(a, sendcount, blocktype_sender, dest, 0, b, recvcount, blocktype_receiver, source, 0, MPI_COMM_WORLD, &status);
    //ierr = MPI_Sendrecv(a, sendcount, blocktype, dest, 0, b, recvcount, blocktype_local, source, 0, MPI_COMM_WORLD, &status);
    //ierr = MPI_Sendrecv(a, sendcount, MPI_CHAR, dest, 0, b, recvcount, MPI_CHAR, source, 0, MPI_COMM_WORLD, &status);
    if (ierr != 0) {
        printf("Problem with Sendrecv (Scatter).\nExiting\n");
        MPI_Abort(MPI_COMM_WORLD,3);
    }

    //cout << "end sendrecv" << endl;
    } // 二回目:for (int no = 0; no < pool_size_cart; ++no)

    MPI_Barrier(MPI_COMM_WORLD);
    cout.flush();
    MPI_Barrier(MPI_COMM_WORLD);

    /* each proc prints it's local "b" out, in order */
    for (int proc=0; proc<p; proc++) {
        if (proc == rank) {
	  printf("Rank = %d  myrow=%d mycol=%d\n", rank, myrow, mycol);
            if (rank == 0) {
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

    // 以下，送受信(gather)
    char c[ROWS*COLS];

    if (rank == 0) {
        for (int ii=0; ii<ROWS*COLS; ++ii) {
	  c[ii] = (char)(-1);
        }
    }

    for (int no = 0; no < pool_size_cart; ++no) {
      MPI_Cart_coords(cart_comm, no, NDIMS, coords);
      //int myrow = 0, mycol = 0;
      int myrow = coords[0], mycol = coords[1];
      //cout << endl << endl << "##my_rank_cart=" << my_rank_cart << "  myrow=" << myrow << "  mycol=" << mycol << endl;
      
      int rank_send = 0;
      int rank_recv = no;
      
      int sendcount = 1; //BLOCKROWS/2;
      int recvcount = 1; //BLOCKROWS/2;
      int dest, source;
      if (rank == rank_recv) {
	source = rank_send;
      } 
      else {  // 自プロセスが受信者ではない場合(送信者である場合を含む)
	source = MPI_PROC_NULL;
	recvcount = 0;
      }
      if (rank == rank_send)  dest = rank_recv; 
      else {  // 自プロセスが送信者ではない場合(受信者である場合を含む)
	dest = MPI_PROC_NULL;
	sendcount = 0;
      }
      ierr = MPI_Sendrecv(b, recvcount, local_dist_array[no], source, 0, c, sendcount, global_dist_array[no], dest, 0, MPI_COMM_WORLD, &status);

      if (ierr != 0) {
        printf("Problem with Sendrecv (Scatter).\nExiting\n");
        MPI_Abort(MPI_COMM_WORLD,3);
      }
    } // 二回目:for (int no = 0; no < pool_size_cart; ++no)

    MPI_Barrier(MPI_COMM_WORLD);

    /* each proc prints it's local "b" out, in order */
    if (rank == 0) {
      printf("Global matrix: \n");
      for (int ii=0; ii<ROWS; ++ii) {
	for (int jj=0; jj<COLS; ++jj) {
	  printf("%3d ",(int)(unsigned char)c[ii*COLS+jj]);
	}
	printf("\n");
      }
    }

    MPI_Barrier(MPI_COMM_WORLD);

    MPI_Finalize();
    return 0;
}

