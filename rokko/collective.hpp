#ifndef ROKKO_COLLECTIVE_H
#define ROKKO_COLLECTIVE_H

#include <mpi.h>

#include <Eigen/Dense>

namespace rokko {

// ローカル行列のサイズ，ブロック数の計算
template<typename MATRIX_MAJOR>
void calculate_local_matrix_size(const rokko::distributed_matrix<MATRIX_MAJOR>& mat, int proc_row, int proc_col, int& local_num_block_rows, int& local_num_block_cols, int& local_matrix_rows, int& local_matrix_cols, int& local_rest_block_rows, int& local_rest_block_cols)
{
  int m_global = mat.m_global;  int n_global = mat.n_global;  int mb = mat.mb;  int nb = mat.nb;
  int m_local = mat.m_local;  int n_local = mat.n_local;
  //int n_local = mat.m_local;  int m_local = mat.n_local;

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


template<typename MATRIX_MAJOR>
void create_struct_local_eigenK(const rokko::distributed_matrix<MATRIX_MAJOR>& mat, MPI_Datatype& local_array_type)
{
  int m_global = mat.m_global;  int n_global = mat.n_global;  int mb = mat.mb;  int nb = mat.nb;
  //int m_local = mat.n_local;  int n_local = mat.m_local;
  int m_local = mat.m_local;  int n_local = mat.n_local;
  int myrank = mat.myrank;
  int myrow = mat.myrow;  int mycol = mat.mycol; int nprow = mat.nprow;  int npcol = mat.npcol;
  int lld = mat.lld;

  int count_max = m_local;
  int*          array_of_blocklengths = new int[count_max];
  MPI_Aint*     array_of_displacements = new MPI_Aint[count_max];
  MPI_Datatype* array_of_types = new MPI_Datatype[count_max];

  int count = 0;
  for (int i=0; i<m_local; ++i) {
    array_of_blocklengths[count] = n_local;  //m_local;
    array_of_displacements[count] = i * lld * sizeof(double);
    array_of_types[count] = MPI_DOUBLE;
    //if (myrank == 0)
    //  std::cout << "verify: count=" << count << "  length=" << array_of_blocklengths[count] << "  disp="  << (int)array_of_displacements[count] << std::endl;
    ++count;
  }

  // print out struct of local matrix
  if (myrank == 0) {
    for (int i=0; i<count; ++i) {
      //cout << "proc=" << proc << "  count=" << count << " lengthhhhh=" << array_of_blocklengths[i] << std::endl;
      printf("eigen_K local Type proc=%d count=%d:  length=%3d  disp=%3d\n", myrank, i, array_of_blocklengths[i], (int)array_of_displacements[i]/8);
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


template<typename MATRIX_MAJOR>
void create_struct_global_eigenK(const rokko::distributed_matrix<MATRIX_MAJOR>& mat, MPI_Datatype& local_array_type, int proc)
{
  int m_global = mat.m_global;  int n_global = mat.n_global;  int mb = mat.mb;  int nb = mat.nb;
  //int m_local = mat.n_local;  int n_local = mat.m_local;
  int m_local = mat.m_local;  int n_local = mat.n_local;
  //int myrow = mat.myrow;  int mycol = mat.mycol;
  int nprow = mat.nprow;  int npcol = mat.npcol;
  //int lld = mat.lld;

  int myrank = mat.myrank;
  int myrow = mat.g.calculate_grid_row(proc);
  int mycol = mat.g.calculate_grid_col(proc);

  int count_max = m_local * n_local;
  int*          array_of_blocklengths = new int[count_max];
  MPI_Aint*     array_of_displacements = new MPI_Aint[count_max];
  MPI_Datatype* array_of_types = new MPI_Datatype[count_max];

  int count = 0;

  // C言語と同様にi(行), j(列)の順番にfor文を回す．
  for (int i=0; i<m_local; ++i) {
    for (int j=0; j<n_local; ++j) {
      array_of_blocklengths[count] = 1; //mb;
      // C言語と同様にrow-major
      array_of_displacements[count] = ( (myrow + nprow * i) * m_global + (mycol + npcol * j) ) * sizeof(double);
      //array_of_displacements[count] = j * m_local * sizeof(double);
      array_of_types[count] = MPI_DOUBLE;
      //if (myrank == 0)
      //  std::cout << "verify: count=" << count << "  length=" << array_of_blocklengths[count] << "  disp="  << (int)array_of_displacements[count] << std::endl;
      ++count;
    }
  }
  /*
  for (int j=0; j<n_local; ++j) {
    array_of_blocklengths[count] = m_local;
    array_of_displacements[count] = j * lld * sizeof(double);
    array_of_types[count] = MPI_DOUBLE;
    //if (myrank == 0)
    //  std::cout << "verify: count=" << count << "  length=" << array_of_blocklengths[count] << "  disp="  << (int)array_of_displacements[count] << std::endl;
    ++count;
  }
  */

  // print out struct of local matrix
  if (myrank == 0) {
    for (int i=0; i<count; ++i) {
      //cout << "proc=" << proc << "  count=" << count << " lengthhhhh=" << array_of_blocklengths[i] << std::endl;
      printf("eigen_K global Type proc=%d count=%d:  length=%3d  disp=%3d\n", myrank, i, array_of_blocklengths[i], (int)array_of_displacements[i]/8);
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

template<typename MATRIX_MAJOR>
void create_struct_local_general(const rokko::distributed_matrix<MATRIX_MAJOR>& mat, MPI_Datatype& local_array_type)
{
  int m_global = mat.m_global;  int n_global = mat.n_global;  int mb = mat.mb;  int nb = mat.nb;
  //int m_local = mat.n_local;  int n_local = mat.m_local;
  int m_local = mat.m_local;  int n_local = mat.n_local;
  int myrank = mat.myrank;
  int myrow = mat.myrow;  int mycol = mat.mycol; int nprow = mat.nprow;  int npcol = mat.npcol;
  int lld = mat.lld;

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

  // print out struct of local matrix
  if (myrank == 0) {
    for (int i=0; i<count; ++i) {
      //cout << "proc=" << proc << "  count=" << count << " lengthhhhh=" << array_of_blocklengths[i] << std::endl;
      printf("local Type proc=%d count=%d:  length=%3d  disp=%3d\n", myrank, i, array_of_blocklengths[i], (int)array_of_displacements[i]/8);
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

// struct of global matrix
template<typename MATRIX_MAJOR>
void create_struct_global_general(const rokko::distributed_matrix<MATRIX_MAJOR>& mat, MPI_Datatype& global_array_type, int proc)
{
  int m_global = mat.m_global;  int n_global = mat.n_global;  int mb = mat.mb;  int nb = mat.nb;
  int m_local = mat.m_local;  int n_local = mat.n_local;
  int nprow = mat.nprow;  int npcol = mat.npcol;

  int myrank = mat.myrank;
  int myrow, mycol;
  myrow = mat.g.calculate_grid_row(proc);
  mycol = mat.g.calculate_grid_col(proc);

  int num_block_rows, num_block_cols, local_matrix_rows, local_matrix_cols, rest_num_block_rows, rest_num_block_cols;
  calculate_local_matrix_size(mat, myrow, mycol, num_block_rows, num_block_cols, local_matrix_rows, local_matrix_cols, rest_num_block_rows, rest_num_block_cols);

  const int type_block_rows = (m_local + mb - 1) / mb;   // 切り上げ
  int count_max = type_block_rows * n_local;

  //cout << "proc=" << myrank << "count_max=" << count_max << std::endl;
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
	  if (myrank == 0)
	    std::cout << "verify: count=" << count << "  length=" << array_of_blocklengths[count] << "  disp="  << (int)array_of_displacements[count] << std::endl;
	  ++count;
      }
      if (rest_num_block_rows != 0) {
	array_of_blocklengths[count] = rest_num_block_rows;
	array_of_displacements[count] = ( ((i*npcol + mycol)*nb+k) * m_global + (num_block_rows * nprow + myrow) * mb ) * sizeof(double);
	if (myrank == 0)
	    std::cout << "amari: count=" << count << "  length=" << array_of_blocklengths[count] << "  disp="  << array_of_displacements[count] << std::endl;
	array_of_types[count] = MPI_DOUBLE;
	++count;
      }
    }
  }
  //cout << "before amari: count=" << count << std::endl;
  // 列のあまり
  for (int k=0; k<rest_num_block_cols; ++k) {
    std::cout << "ddddddddddddddddddddddddddddddd" << std::endl;
    for (int j=0; j<num_block_rows; ++j) {
      array_of_blocklengths[count] = mb;
      array_of_displacements[count] = ( ((num_block_cols * npcol + mycol)*nb+k) * m_global + (j * nprow + myrow) * mb ) * sizeof(double);
      array_of_types[count] = MPI_DOUBLE;
      ++count;
    }
    if (rest_num_block_rows != 0) {
      std::cout << "rest: count=" << count << "  disp="  << array_of_displacements[count-1] << std::endl;
      array_of_blocklengths[count] = rest_num_block_rows;
      array_of_displacements[count] = ( ((num_block_cols * npcol + mycol)*nb+k) * m_global + (num_block_rows * nprow + myrow) * mb ) * sizeof(double);
      if (myrank == 0)
	std::cout << "rest_amari: count=" << count << "  disp="  << array_of_displacements[count] << std::endl;
      array_of_types[count] = MPI_DOUBLE;
      ++count;
    }
  }

  //if (myrank == 0)
  //cout << "myrank" << myrank << "  imano_count=" << count << "   num_block_rows_r=" << num_block_rows_r << "  num_block_cols=" << num_block_cols << "  num_block_rows=" << num_block_rows << std::endl;

  // print out struct of global matrix
  if (myrank == 0) {
    for (int i=0; i<count; ++i) {
	//cout << "proc=" << proc << "  count=" << count << " lengthhhhh=" << array_of_blocklengths[i] << std::endl;
      printf("global Type proc=%d count=%d:  length=%3d  disp=%3d\n", proc, i, array_of_blocklengths[i], (int)array_of_displacements[i]/8);
    }
    std::cout << "num_block_rows=" << num_block_rows << std::endl;
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

template<typename MATRIX_MAJOR>
void copy_g2l_root(Eigen::MatrixXd& mat_global, const rokko::distributed_matrix<MATRIX_MAJOR>& mat)
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
  // 列のあまり
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

template<typename MATRIX_MAJOR>
void copy_l2g_root(const rokko::distributed_matrix<MATRIX_MAJOR>& mat, Eigen::MatrixXd& mat_global)
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
	  std::cout << "ROW=" << (num_block_rows*nprow+myrow)*mb+c << "  COL=" << (i*npcol+mycol)*nb+k << "rest_num_block_rows=" << rest_num_block_rows << std::endl;
	  mat_global( (num_block_rows*nprow+myrow)*mb+c, (i*npcol+mycol)*nb+k ) = mat.array[(i*nb+k) * lld + num_block_rows * mb + c];
	}
	++count;
      }
    }
  }
  // 列のあまり
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

void create_struct_local(const rokko::distributed_matrix<rokko::matrix_col_major>& mat, MPI_Datatype& local_array_type)
{
  create_struct_local_general(mat, local_array_type);
  //if ((mat.mb == 1) && (mat.nb == 1)) {
  //  create_struct_local_eigenK(mat, local_array_type);
  //}
  //else {
  //  create_struct_local_general(mat, local_array_type);
  //}
}


void create_struct_global(const rokko::distributed_matrix<rokko::matrix_col_major>& mat, MPI_Datatype& global_array_type, int proc)
{
  create_struct_global_general(mat, global_array_type, proc);
  //if ((mat.mb == 1) && (mat.nb == 1)) {
  //  create_struct_global_eigenK(mat, global_array_type, proc);
  //}
  //else {
  //  create_struct_global_general(mat, global_array_type, proc);
  //}
}

void create_struct_local(const rokko::distributed_matrix<rokko::matrix_row_major>& mat, MPI_Datatype& local_array_type)
{
  if ((mat.mb == 1) && (mat.nb == 1)) {
    create_struct_local_eigenK(mat, local_array_type);
  }
  else {
    create_struct_local_general(mat, local_array_type);
  }
}

void create_struct_global(const rokko::distributed_matrix<rokko::matrix_row_major>& mat, MPI_Datatype& global_array_type, int proc)
{
  if ((mat.mb == 1) && (mat.nb == 1)) {
    create_struct_global_eigenK(mat, global_array_type, proc);
  }
  else {
    create_struct_global_general(mat, global_array_type, proc);
  }
}



template<typename MATRIX_MAJOR>
int gather(rokko::distributed_matrix<MATRIX_MAJOR>& mat, Eigen::MatrixXd& mat_global, int root)
{
  double* global_array;

  if (mat.myrank == root) {
    mat_global.resize(mat.m_global, mat.n_global);
    global_array = &mat_global(0,0);  // 本当に、内部では連続な配列になっているか？
    //global_array = mat_global.data();  // 本当に、内部では連続な配列になっているか
    //for (int ii=0; ii<mat.m_global * mat.n_global; ++ii) {
    //  global_array[ii] = 100*mat.myrank + ii;
    //}
  }

  MPI_Status  status;
  int ierr;

  int local_matrix_rows, local_matrix_cols;
  int local_num_block_rows, local_num_block_cols;

  int local_rest_block_rows, local_rest_block_cols;
  int count_max;

  MPI_Comm cart_comm = MPI_COMM_WORLD;

  int m_global = mat.m_global;  int n_global = mat.n_global;  int mb = mat.mb;  int nb = mat.nb;
  int m_local = mat.m_local;  int n_local = mat.n_local;
  int myrow = mat.myrow;  int mycol = mat.mycol; int nprow = mat.nprow;  int npcol = mat.npcol;
  int myrank = mat.myrank; int nprocs = mat.nprocs;

  double* local_array = mat.get_array_pointer();

  MPI_Datatype local_array_type, global_array_type;
  create_struct_local(mat, local_array_type);

  int rank_recv = root;  // プロセスrank_recvに集約
  int sendcount = 1;
  int recvcount = 1;

  for (int proc = 0; proc < nprocs; ++proc) {
    if ((myrank == proc) && (myrank != root)) {
      ierr = MPI_Send(local_array, sendcount, local_array_type, root, 0, cart_comm);
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
      global_array_type = NULL;
    }

    if ((proc == root) && (myrank == root)) {
      create_struct_global(mat, global_array_type, root);
      ierr = MPI_Sendrecv(local_array, sendcount, local_array_type, root, 0, global_array, recvcount, global_array_type, root, 0, cart_comm, &status);
      if (ierr != 0) {
      	printf("Error with Sendrecv (Gather). ierr=%d\nExiting\n", ierr);
      	MPI_Abort(MPI_COMM_WORLD,78);
      	exit(78);
      }
      MPI_Type_free(&global_array_type);
      global_array_type = NULL;
      //copy_l2g_root(mat, mat_global);
    }

  } // for (int proc = 0; proc < nprocs; ++proc)

  MPI_Type_free(&local_array_type);
  local_array_type = NULL;

  return ierr;
}


template<typename MATRIX_MAJOR>
int scatter(rokko::distributed_matrix<MATRIX_MAJOR>& mat, Eigen::MatrixXd& mat_global, int root)
{
  double* global_array;

  if (mat.myrank == root) {
    // mat_globalは初期化，値の代入後に呼び出されているから，resizeはしない
    //mat_global.resize(mat.m_global, mat.n_global);
    global_array = &mat_global(0,0);  // 本当に、内部では連続な配列になっているか？
    //global_array = mat_global.data();  // 本当に、内部では連続な配列になっているか
    for (int ii=0; ii<mat.m_global * mat.n_global; ++ii) {
      global_array[ii] = ii;
    }
  }

  MPI_Status  status;
  int ierr;

  int local_matrix_rows, local_matrix_cols;
  int local_num_block_rows, local_num_block_cols;

  int local_rest_block_rows, local_rest_block_cols;
  int count_max;

  MPI_Comm cart_comm = MPI_COMM_WORLD;

  int m_global = mat.m_global;  int n_global = mat.n_global;  int mb = mat.mb;  int nb = mat.nb;
  int m_local = mat.m_local;  int n_local = mat.n_local;
  int myrow = mat.myrow;  int mycol = mat.mycol; int nprow = mat.nprow;  int npcol = mat.npcol;
  int myrank = mat.myrank;
  int nprocs = mat.nprocs;

  double* local_array = mat.get_array_pointer();

  MPI_Datatype local_array_type, global_array_type;
  create_struct_local(mat, local_array_type);

  int sendcount = 1, recvcount = 1;

  for (int proc = 0; proc < nprocs; ++proc) {
    if ((proc != root) &&  (myrank == root)) {
      create_struct_global(mat, global_array_type, proc);
      ierr = MPI_Send(global_array, sendcount, global_array_type, proc, 0, cart_comm);
      if (ierr != 0) {
	printf("Error with Recv (Scatter). ierr=%d\nExiting\n", ierr);
	MPI_Abort(MPI_COMM_WORLD,68);
	exit(78);
      }
      MPI_Type_free(&global_array_type);
      global_array_type = NULL;
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
      ierr = MPI_Sendrecv(global_array, sendcount, global_array_type, root, 0, local_array, recvcount, local_array_type, root, 0, cart_comm, &status);
      if (ierr != 0) {
      printf("Error with Sendrecv (Scatter). ierr=%d\nExiting\n", ierr);
      	MPI_Abort(MPI_COMM_WORLD,70);
	exit(78);
      }
      MPI_Type_free(&global_array_type);
      global_array_type = NULL;
      //copy_g2l_root(mat_global, mat);
    }

  } // for (int proc = 0; proc < nprocs; ++proc)

  MPI_Type_free(&local_array_type);
  local_array_type = NULL;

  return ierr;
}


} // namespace rokko

#endif // ROKKO_COLLECTIVE_H


