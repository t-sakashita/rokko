// calculatong a product of hamiltonian and vector in quantum spin model
// Copyright (C) 2013 by Tatsuya Sakashita

#include <vector>

namespace rokko {

namespace spin_hamiltonian {

void multiply(int L, std::vector<std::pair<int, int> >& lattice, const double* v, double* w) {
  int N = 1 << L;
  for (int k=0; k<N; ++k) {
    for (int l=0; l<lattice.size(); ++l) {
      int i = lattice[l].first;
      int j = lattice[l].second;
      //cout << "k=" << k << " i=" << i << " j=" << j << endl;
      int m1 = 1 << i;
      int m2 = 1 << j;
      int m3 = m1 + m2;
      if (((k & m3) == m1)) {  // when (bit i == 1, bit j == 0) or (bit i == 0, bit j == 1)
        w[k] += 0.5 * v[k^m3] - 0.25 * v[k];
        //cout << "if" << endl;
      }
      else if ((k & m3) == m2) {
        w[k] += 0.5 * v[k^m3] - 0.25 * v[k];
      }
      else {
        w[k] += 0.25 * v[k];
        //      cout << "else" << endl;
      }
    }
  }
}

void multiply(int L, std::vector<std::pair<int, int> >& lattice, const std::vector<double>& v, std::vector<double>& w) {
  multiply(L, lattice, &v[0], &w[0]);
}

void fill_diagonal(int L, std::vector<std::pair<int, int> >& lattice, double* w) {
  int N = 1 << L;
  for (int k=0; k<N; ++k) {
    w[k] = 0;
    for (int l=0; l<lattice.size(); ++l) {
      int i = lattice[l].first;
      int j = lattice[l].second;
      //cout << "k=" << k << " i=" << i << " j=" << j << endl;
      int m1 = 1 << i;
      int m2 = 1 << j;
      int m3 = m1 + m2;
      if (((k & m3) == m1)) {  // when (bit i == 1, bit j == 0) or (bit i == 0, bit j == 1)
        w[k] -= 0.25;
        //cout << "if" << endl;
      }
      else if ((k & m3) == m2) {
        w[k] -= 0.25;
      }
      else {
        w[k] += 0.25;
        //      cout << "else" << endl;
      }
    }
  }
}

void fill_diagonal(int L, std::vector<std::pair<int, int> >& lattice, std::vector<double>& w) {
  fill_diagonal(L, lattice, &w[0]);
}

}

} 


