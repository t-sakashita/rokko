#include <mpi.h>
#include <rokko/grid.hpp>

namespace rokko {

class distributed_matrix
{
public:
  distributed_matrix(int m_global, int n_global, const grid& g)
    : m_global(m_global), n_global(n_global), g(g), myrank(g.myrank), nprocs(g.nprocs), myrow(g.myrow), mycol(g.mycol), nprow(g.nprow), npcol(g.npcol), ictxt(g.ictxt)
  {
    // ローカル行列の形状を指定
    mb = m_global / g.nprow;
    if (mb == 0) mb = 1;
    //mb = 10;        
    nb = n_global / g.npcol;
    if (nb == 0) nb = 1;
    //nb = 10;
    // mbとnbを最小値にそろえる．（注意：pdsyevではmb=nbでなければならない．）
    mb = min(mb, nb);
    nb = mb;

    const int ZERO=0, ONE=1;
    m_local = numroc_( m_global, mb, myrow, ZERO, nprow );
    n_local = numroc_( n_global, nb, mycol, ZERO, npcol );
    int lld = m_local;
    if (lld == 0) lld = 1;

    for (int proc=0; proc<nprocs; ++proc) {
      if (proc == g.myrank) {
	cout << "proc=" << proc << endl;
	cout << "  mb=" << mb << "  nb=" << nb << endl;
	cout << "  mA=" << m_local << "  nprow=" << g.nprow << endl;
	cout << "  nA=" << n_local << "  npcol=" << g.npcol << endl;
	cout << "  lld=" << lld << " m_local=" << m_local << " n_local=" << n_local << endl;
      }
      MPI_Barrier(MPI_COMM_WORLD);
    }

    descinit_(desc, m_global, n_global, mb, nb, ZERO, ZERO, ictxt, lld, info);
    if (info) {
      cerr << "error " << info << " at descinit function of descA " << "mA=" << m_local << "  nA=" << n_local << "  lld=" << lld << "." << endl;
      exit(1);
    }

    array = new double[m_local * n_local];
    if (array == NULL) {
      cerr << "failed to allocate array." << endl;
      exit(1);
    }
  }

  /*
  ~distributed_matrix()
  {
    cout << "Destructor ~Distributed_Matrix()" << endl;
    delete[] array;
    array = NULL;
  }
  */

  int translate_l2g_row(const int& local_i) const
  {
    return (g.myrow * mb) + local_i + (local_i / mb) * mb * (g.nprow - 1);
  }

  int translate_l2g_col(const int& local_j) const
  {
    return (g.mycol * nb) + local_j + (local_j / nb) * nb * (g.npcol - 1);
  }

  int translate_g2l_row(const int& global_i) const
  {
    const int local_offset_block = global_i / mb;
    return (local_offset_block - g.myrow) / g.nprow * mb + global_i % nb;
  }

  int translate_g2l_col(const int& global_j) const
  {
    const int local_offset_block = global_j / nb;
    return (local_offset_block - g.mycol) / g.npcol * nb + global_j % nb;
  }

  bool is_gindex_myrow(const int& global_i) const
  {
    int local_offset_block = global_i / mb;
    return (local_offset_block % g.nprow) == g.myrow;
  }

  bool is_gindex_mycol(const int& global_j) const
  {
    int local_offset_block = global_j / nb;
    return (local_offset_block % g.npcol) == g.mycol;
  }

  void print_matrix() const
  {
    /* each proc prints it's local_array out, in order */
    for (int proc=0; proc<nprocs; ++proc) {
      if (proc == myrank) {
	printf("Rank = %d  myrow=%d mycol=%d\n", myrank, myrow, mycol);
	printf("Local Matrix:\n");
	for (int ii=0; ii<m_local; ++ii) {
	  for (int jj=0; jj<n_local; ++jj) {
	    printf("%3.2f ",array[ii * n_local + jj]);
	  }
	  printf("\n");
	}
	printf("\n");
      }
      MPI_Barrier(MPI_COMM_WORLD);
    }
}

  int m_global, n_global;
  double* array;
  int mb, nb;
  int m_local, n_local;
  // variables of class Grid
  int myrank, nprocs;
  int myrow, mycol;
  int ictxt, nprow, npcol;
  Grid g;

private:
  int info;
};

} // namespace rokko

