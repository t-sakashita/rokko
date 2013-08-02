// calculatong a product of hamiltonian and vector in quantum spin model
// Copyright (C) 2013 by Tatsuya Sakashita

#include <iostream>
#include <vector>
#include <rokko/utility/spin_hamiltonian.hpp>

int main()
{
  int L = 5; //4;
  int N = 1 << L;
  std::vector<std::pair<int, int> > lattice;
  for (int i=0; i<L-1; ++i) {
    lattice.push_back(std::make_pair(i, i+1));
  }

  for (int i=0; i<N; ++i) {
    std::vector<double> v, w;
    v.assign(N, 0);
    v[i] = 1;
    w.assign(N, 0);
    rokko::spin_hamiltonian::multiply(L, lattice, v, w);
    for (int j=0; j<N; ++j) {
      std::cout << w[j] << " ";
    }
    std::cout << std::endl;
  }

  std::vector<double> v(N);
  rokko::spin_hamiltonian::fill_diagonal(L, lattice, v);
  for (int j=0; j<N; ++j) {
    std::cout << v[j] << " ";
  }
  std::cout << std::endl;

}




