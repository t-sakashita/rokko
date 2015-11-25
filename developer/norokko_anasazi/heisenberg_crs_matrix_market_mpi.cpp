#include "AnasaziConfigDefs.hpp"
#include "AnasaziBasicEigenproblem.hpp"
//#include "AnasaziSimpleLOBPCGSolMgr.hpp"
#include "AnasaziBlockKrylovSchurSolMgr.hpp"
#include "AnasaziBasicOutputManager.hpp"
#include "AnasaziEpetraAdapter.hpp"
#include "Epetra_CrsMatrix.h"
#include "Teuchos_CommandLineProcessor.hpp"

#ifdef HAVE_MPI
#include "Epetra_MpiComm.h"
#include <mpi.h>
#else
#include "Epetra_SerialComm.h"
#endif
#include "Epetra_Map.h"

using namespace Anasazi;

int main(int argc, char *argv[]) {

  using std::endl;

#ifdef HAVE_MPI
  // Initialize MPI
  MPI_Init(&argc,&argv);
#endif

  // Create an Epetra communicator
#ifdef HAVE_MPI
  Epetra_MpiComm Comm(MPI_COMM_WORLD);
#else
  Epetra_SerialComm Comm;
#endif

  // Create an Anasazi output manager
  BasicOutputManager<double> printer;
  //printer.stream(Errors) << Anasazi_Version() << endl << endl;
  // Get the sorting string from the command line
  std::string which("LM");
  Teuchos::CommandLineProcessor cmdp(false,true);
  cmdp.setOption("sort",&which,"Targetted eigenvalues (SM or LM).");
  if (cmdp.parse(argc,argv) != Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL) {
#ifdef HAVE_MPI
    MPI_Finalize();
#endif
    return -1;
  }

  int L = 8;
  cmdp.setOption("L", &L ,"Lattice size.");
  int N = 1 << L;
  std::vector<std::pair<int, int> > lattice;
  for (int i=0; i<L; ++i) {
    lattice.push_back(std::make_pair(i, (i+1)%L));
  }

  // Construct a Map that puts approximately the same number of
  // equations on each processor.
  Epetra_Map Map(N, 0, Comm);

  // Get update list and number of local equations from newly created Map.
  int NumMyElements = Map.NumMyElements();
  std::vector<int> MyGlobalElements(NumMyElements);
  Map.MyGlobalElements(&MyGlobalElements[0]);

  // Create an Epetra_Matrix
  Teuchos::RCP<Epetra_CrsMatrix> A = Teuchos::rcp( new Epetra_CrsMatrix(Copy, Map, N) );  // fix me: NumEntriesPerRow

  // Compute coefficients for hamiltonian matrix of quantum Heisenberg model
  std::vector<double> values;
  std::vector<int> cols;

  for (int l=0; l<lattice.size(); ++l) {
    for (int row=0; row<NumMyElements; ++row) {
      cols.clear();
      values.clear();
      int k = MyGlobalElements[row];
      int i = lattice[l].first;
      int j = lattice[l].second;
      int m1 = 1 << i;
      int m2 = 1 << j;
      int m3 = m1 + m2;
      if (((k & m3) == m1) || ((k & m3) == m2)) {  // when (bit i == 1, bit j == 0) or (bit i == 0, bit j == 1)
        cols.push_back(k^m3);
        values.push_back(0.5);
        cols.push_back(k);
        values.push_back(-0.25);
      } else {
        cols.push_back(k);
        values.push_back(0.25);
      }
      int info = A->InsertGlobalValues(k, cols.size(), &values[0], &cols[0]);
      //cout << "info=" << info << endl;
      assert( info==0 );
    }
  }

  // Finish up
  int info = A->FillComplete();
  assert( info==0 );
  A->SetTracebackMode(1); // Shutdown Epetra Warning tracebacks
  int num_cols;
  double* values_;
  int* cols_;
  int local_row = 0;
  std::cout << "%%MatrixMarket matrix coordinate real general" << std::endl;
  std::cout << N << " " << N << " " << A->NumGlobalNonzeros() << std::endl;
  for (int global_row=0; global_row<N; ++global_row) {
    if (global_row == MyGlobalElements[local_row]) {
      A->ExtractMyRowView(local_row, num_cols, values_, cols_);
      for (int i=0; i<num_cols; ++i) {
	std::cout << global_row + 1 << " " << cols_[i] + 1 << " " << values_[i] << std::endl;
      }
      ++local_row;
    }
    MPI_Barrier(MPI_COMM_WORLD);
  }

#ifdef HAVE_MPI
  MPI_Finalize();
#endif
  return 0;
}

