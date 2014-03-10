#ifndef ROKKO_UTILITY_SORT_EIGENPAIRS_H
#define ROKKO_UTILITY_SORT_EIGENPAIRS_H

namespace rokko {

template<typename MATRIX_MAJOR>  
void sort_eigenpairs(const localized_vector& eigval, const localized_matrix<MATRIX_MAJOR>& eigvec,
  localized_vector& eigval_sorted, localized_matrix<MATRIX_MAJOR>& eigvec_sorted) {
  int dim = eigval.size();
  int* q = new int[dim];

  // 固有値を（絶対値ではなく）昇順に並べる
  if (q==NULL) {
    std::cerr << "error: q" << std::endl;
    return;
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
    eigvec_sorted.col(k) = eigvec.col(q[k]);
  }
  delete[] q;
}

} // namespace rokko

#endif // ROKKO_UTILITY_SORT_EIGENPAIRS_H
