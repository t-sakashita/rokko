#ifndef ROKKO_UTILITY_SORT_EIGENPAIRS_H
#define ROKKO_UTILITY_SORT_EIGENPAIRS_H

#include <vector>
#include <rokko/localized_vector.hpp>
#include <rokko/localized_matrix.hpp>

namespace rokko {

// sort eivanvalue and eigenvectors 

template<typename MATRIX_MAJOR>  
void sort_eigenpairs(const localized_vector& eigval, const localized_matrix<MATRIX_MAJOR>& eigvec,
                     localized_vector& eigval_sorted, localized_matrix<MATRIX_MAJOR>& eigvec_sorted,
                     bool ascending = true) 
{
  int dim = eigval.size();
  std::vector<int> q(dim);
  double tmp;
  for (int i=0; i<dim; ++i) q[i] = i;
  for (int k=0; k<dim; ++k) {
    tmp = eigval(q[k]);
    for (int i=k+1; i<dim; ++i) {
      if(ascending){
        if (tmp > eigval(q[i])) { 
          tmp = eigval(q[i]);
          int qq = q[k];
          q[k] = q[i];
          q[i] = qq;
        }
      }else{
        if (tmp < eigval(q[i])) { 
          tmp = eigval(q[i]);
          int qq = q[k];
          q[k] = q[i];
          q[i] = qq;
        }
      }
    }
    eigval_sorted(k) = eigval(q[k]);
    if(eigvec_sorted.is_col_major()){
      eigvec_sorted.col(k) = eigvec.col(q[k]);
    }else{
      eigvec_sorted.row(k) = eigvec.row(q[k]);
    }
  }
}

} // namespace rokko

#endif // ROKKO_UTILITY_SORT_EIGENPAIRS_H
