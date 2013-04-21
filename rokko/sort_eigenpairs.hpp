#ifndef ROKKO_SORT_EIGENPAIRS_H
#define ROKKO_SORT_EIGENPAIRS_H

namespace rokko {

int sort_eigenpairs(const Eigen::VectorXd& eigval, const Eigen::MatrixXd& eigvec_global, Eigen::VectorXd& eigval_sorted, Eigen::MatrixXd& eigvec_sorted)
{
  int dim = eigval.size();
  int* q = new int[dim];

  // 固有値を（絶対値のではなく）昇順に並べる
  if (q==NULL) {
    cerr << "error: q" << endl;
    return 1;
  }

  double emax;
  for (int i=0; i<dim; ++i) q[i] = i;
  for (int k=0; k<dim; ++k) {
    emax = eigval(q[k]);
    for (int i=k+1; i<dim; ++i) {
      if (emax < eigval(q[i])) {       // 昇順になっていないとき、交換
	emax = eigval(q[i]);
	int qq = q[k];
	q[k] = q[i];
	q[i] = qq;
      }
    }
    eigval_sorted(k) = eigval(q[k]);
    eigvec_sorted.col(k) = eigvec_global.col(q[k]);
  }
  delete[] q;
  return 0;
}

} // namespace rokko

#endif // ROKKO_SORT_EIGENPAIRS_H
