#include <mpi.h>
#include <rokko/rokko.h>
#include <stdio.h>
#include <stdlib.h>

using matrix_major = rokko::matrix_col_major;

template<typename T, typename MATRIX_MAJOR>
void function_matrix(Eigen::VectorXd const& eigval_tmp, rokko::distributed_matrix<T, MATRIX_MAJOR> const& eigvec, rokko::distributed_matrix<T, MATRIX_MAJOR>& result, rokko::distributed_matrix<T, MATRIX_MAJOR>& tmp) {
  for (int local_j=0; local_j<eigvec.get_n_local(); ++local_j) {
    int global_j = eigvec.translate_l2g_col(local_j);
    double coeff = eigval_tmp(global_j);
    for (int local_i=0; local_i<eigvec.get_m_local(); ++local_i) {
      double value = eigvec.get_local(local_i, local_j);
      tmp.set_local(local_i, local_j, coeff * value); 
    }
  }
  product(1, tmp, false, eigvec, true, 0, result);
}

template<typename T, typename MATRIX_MAJOR>
void diagonalize_fixedB(rokko::parallel_dense_ev& solver, rokko::distributed_matrix<T, MATRIX_MAJOR>& A, rokko::distributed_matrix<T, MATRIX_MAJOR>& B, Eigen::VectorXd& eigval, rokko::distributed_matrix<T, MATRIX_MAJOR>& eigvec, T tol = 0) {
  rokko::distributed_matrix<double, matrix_major> tmp(A.get_mapping()), Binvroot(A.get_mapping()), mat(A.get_mapping());
  rokko::parameters params;
  params.set("routine", "");
  solver.diagonalize(B, eigval, eigvec, params);
  // computation of B^{-1/2}
  for(int i=0; i<eigval.size(); ++i)
    eigval(i) = (eigval(i) > tol) ? sqrt(1/eigval(i)) : 0;
  function_matrix(eigval, eigvec, Binvroot, tmp);
  
  // computation of B^{-1/2} A B^{-1/2}
  product(1, Binvroot, false, A, false, 0, tmp);
  product(1, tmp, false, Binvroot, false, 0, mat);
  // diagonalization of B^{-1/2} A B^{-1/2}
  solver.diagonalize(mat, eigval, tmp, params);

  // computation of {eigvec of Ax=lambda Bx} = B^{-1/2} {eigvec of B^{-1/2} A B^{-1/2}}
  product(1, Binvroot, false, tmp, false, 0, eigvec);
}

template<typename T, int MATRIX_MAJOR>
void set_A_B(Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic,MATRIX_MAJOR>& locA, Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic,MATRIX_MAJOR>& locB) {
  if ((locA.rows() != 4) || (locA.cols() != 4) || (locB.rows() != 4) || (locB.cols() != 4)) {
    std::cerr << "error: size must be 4!" << std::endl;
    throw;
  }
  locA << 0.24, 0.39, 0.42, -0.16,
          0.39, -0.11, 0.79, 0.63,
          0.42, 0.79, -0.25, 0.48,
         -0.16, 0.63, 0.48, -0.03;

  locB << 4.16, -3.12, 0.56, -0.10,
         -3.12, 5.03, -0.83, 1.09,
          0.56, -0.83, 0.76, 0.34,
         -0.10, 1.09, 0.34, 1.18;
}

