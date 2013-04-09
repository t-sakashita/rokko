#include <mpi.h>
#include <rokko/scalapack.hpp>

int main(int argc, char* argv[])
{
  MPI_Init(&argc, &argv); 
  MPI_Comm comm = MPI_COMM_WORLD;
  rokko::Grid g(comm);
  
  int dim = 10;

  rokko::distributed_matrix mat(dim, dim, g);
  rokko::generate_Frank_matrix_local(mat);
  print_matrix(mat);

  eigen::VectorXd eigvals(dim);
  rokko::distributed_matrix eigvecs(dim, dim, g);
  rokko::scalapack::diagonalize(mat, eigvals, eigvecs);
    
  eigen::MatrixXd eigvecs_global;
  const int root = 0;
  gather(eigvecs, eigvecs_global, root);
  rokko::print_matrix(eigvecs);
  if (myrank == root) {
    std::cout << eigvecs_global << std::endl;
  }

  // 固有値の絶対値の降順に固有値(と対応する固有ベクトル)をソート
  // ソート後の固有値の添字をベクトルqに求める．
  double absmax;
  int qq;

  VectorXi q(dim);
  
  VectorXd eigval_sorted(dim);
  cout << "aadim=" << dim << endl;
  MatrixXd eigvec_sorted(dim,dim);
  
  if(myrank == root) {
    // 固有値・固有ベクトルを絶対値の降順にソート
    for (int i=0; i<eigval.size(); ++i) q[i] = i;
    for (int m=0; m<eigval.size(); ++m) {
      absmax = abs(eigval[q[m]]);
      for (int i=m+1; i<eigvec.rows(); ++i) {
	if (absmax < abs(eigval[q[i]])) {
	  absmax = eigval[q[i]];
	  qq = q[m];
	  q[m] = q[i];
	  q[i] = qq;
	}
      }
      eigval_sorted(m) = eigval(q[m]);
      eigvec_sorted.col(m) = eigvec.col(q[m]);    
    }

    cout << "Computed eigenvalues= " << eigval_sorted.transpose() << endl;
    
    cout << "Eigenvector:" << endl << eigvec_sorted << endl << endl;
    cout << "Check the orthogonality of eigenvectors:" << endl
	 << eigvec_sorted * eigvec_sorted.transpose() << endl;   // Is it equal to indentity matrix?
    //cout << "residual := A x - lambda x = " << endl
    //	 << A_global_matrix * eigvec_sorted.col(0)  -  eigval_sorted(0) * eigvec_sorted.col(0) << endl;
    //cout << "Are all the following values equal to some eigenvalue = " << endl
    //	 << (A_global_matrix * eigvec_sorted.col(0)).array() / eigvec_sorted.col(0).array() << endl;
    //cout << "A_global_matrix=" << endl << A_global_matrix << endl; 

  }

  rokko::scalapack::Finzalize();
  MPI_Finalize();
}
