// calculatong a product of hamiltonian and vector in quantum spin model
// Copyright (C) 2013 by Tatsuya Sakashita

#include <vector>


namespace spin_hamiltonian {

void multiply(int L, std::vector<std::pair<int, int> >& lattice, const std::vector<double>& v, std::vector<double>& w) {
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

} 


