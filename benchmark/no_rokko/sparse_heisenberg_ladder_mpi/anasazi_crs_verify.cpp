#include "AnasaziConfigDefs.hpp"
#include "AnasaziBasicEigenproblem.hpp"
#include "AnasaziLOBPCGSolMgr.hpp"
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

#include <boost/lexical_cast.hpp>
#include <rokko/utility/lattice.hpp>
#include <rokko/utility/machine_info.hpp>

using namespace Anasazi;

int main(int argc, char *argv[]) {
  double init_tick, gen_tick, diag_tick, end_tick;

#ifdef HAVE_MPI
  // Initialize MPI
  int provided;
  MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &provided);
#endif
  init_tick = MPI_Wtime();
  // Create an Epetra communicator
#ifdef HAVE_MPI
  Epetra_MpiComm Comm(MPI_COMM_WORLD);
#else
  Epetra_SerialComm Comm;
#endif

  // Create an Anasazi output manager
  BasicOutputManager<double> printer;
  printer.stream(Errors) << Anasazi_Version() << std::endl << std::endl;

  int len_ladder = 5;
  if (argc >= 2) len_ladder = boost::lexical_cast<int>(argv[1]);
  int L = 2 * len_ladder;
  std::vector<std::pair<int, int>> lattice;
  rokko::create_ladder_lattice_1dim(len_ladder, lattice);
  rokko::output_lattice(printer.stream(Anasazi::Errors), lattice);
  int N = 1 << L;

  gen_tick = MPI_Wtime();
  // Construct a Map that puts approximately the same number of
  // equations on each processor.
  Epetra_Map Map(N, 0, Comm);

  // Get update list and number of local equations from newly created Map.
  int NumMyElements = Map.NumMyElements();
  std::vector<int> MyGlobalElements(NumMyElements);
  Map.MyGlobalElements(MyGlobalElements.data());

  // Create an Epetra_Matrix
  int NumEntriesPerRow = 2 * L;
  Teuchos::RCP<Epetra_CrsMatrix> A = Teuchos::rcp( new Epetra_CrsMatrix(Copy, Map, NumEntriesPerRow) );

  // Compute coefficients for hamiltonian matrix of quantum Heisenberg model
  std::vector<double> values;
  std::vector<int> cols;

  for (int local_row=0; local_row<NumMyElements; ++local_row) {
    int row = MyGlobalElements[local_row];
    cols.clear();
    values.clear();
    double diag = 0.;
    for (int l=0; l<lattice.size(); ++l) {
      int i = lattice[l].first;
      int j = lattice[l].second;
      int m1 = 1 << i;
      int m2 = 1 << j;
      int m3 = m1 + m2;
      if (((row & m3) == m1) || ((row & m3) == m2)) {  // when (bit i == 1, bit j == 0) or (bit i == 0, bit j == 1)
        cols.push_back(row^m3);
        values.push_back(0.5);
        diag += -0.25;
      } else {
        diag += 0.25;
      }
    }
    if (diag != 0.) {
      cols.push_back(row);
      values.push_back(diag);
    }
    int info = A->InsertGlobalValues(row, cols.size(), values.data(), cols.data());
    //cout << "info=" << info << endl;
    assert( info==0 );
  }

  // Finish up
  int info = A->FillComplete();
  assert( info==0 );
  A->SetTracebackMode(1); // Shutdown Epetra Warning tracebacks

  //************************************
  // Call the LOBPCG solver manager
  //***********************************
  //  Variables used for the LOBPCG Method
  diag_tick = MPI_Wtime();
  std::string which("LM");
  constexpr int    nev       = 1;
  constexpr int    blockSize = nev;
  //constexpr int    maxIters  = 500;
  //constexpr double tol       = 1.0e-12;

  using MV = Epetra_MultiVector;
  using OP = Epetra_Operator;
  using MVT = Anasazi::MultiVecTraits<double, Epetra_MultiVector>;

  // Create an Epetra_MultiVector for an initial vector to start the solver.
  // Note:  This needs to have the same number of columns as the blocksize.
  Teuchos::RCP<Epetra_MultiVector> ivec = Teuchos::rcp( new Epetra_MultiVector(Map, blockSize) );
  ivec->Random();

  // Create the eigenproblem.
  Teuchos::RCP<BasicEigenproblem<double, MV, OP>> MyProblem =
    Teuchos::rcp( new BasicEigenproblem<double, MV, OP>(A, ivec) );

  // Inform the eigenproblem that the operator A is symmetric
  MyProblem->setHermitian(true);

  // Set the number of eigenvalues requested
  MyProblem->setNEV( nev );

  // Inform the eigenproblem that you are finishing passing it information
  bool boolret = MyProblem->setProblem();
  if (boolret != true) {
    printer.print(Errors,"Anasazi::BasicEigenproblem::setProblem() returned an error.\n");
#ifdef HAVE_MPI
    MPI_Finalize();
#endif
    return -1;
  }

  // Create parameter list to pass into the solver manager
  Teuchos::ParameterList MyPL;
  //std::cout << MyPL.get<std::string>("Which") << std::endl;
  MyPL.set( "Verbosity", Anasazi::Errors | Anasazi::Warnings | Anasazi::IterationDetails | Anasazi::FinalSummary | Anasazi::Debug | Anasazi::OrthoDetails );
  //MyPL.set( "Which", which );
  //MyPL.set( "Block Size", blockSize );
  //MyPL.set( "Maximum Iterations", maxIters );
  //MyPL.set( "Convergence Tolerance", tol );

  // Create the solver manager
  LOBPCGSolMgr<double, MV, OP> MySolverMan(MyProblem, MyPL);

  // Solve the problem
  ReturnType returnCode = MySolverMan.solve();
  end_tick = MPI_Wtime();

  // Get the eigenvalues and eigenvectors from the eigenproblem
  Eigensolution<double,MV> sol = MyProblem->getSolution();
  std::vector<Value<double>> evals = sol.Evals;
  Teuchos::RCP<MV> evecs = sol.Evecs;

  // Compute residuals.
  std::vector<double> normR(sol.numVecs);
  if (sol.numVecs > 0) {
    Teuchos::SerialDenseMatrix<int,double> T(sol.numVecs, sol.numVecs);
    Epetra_MultiVector tempAevec( Map, sol.numVecs );
    T.putScalar(0.0); 
    for (int i=0; i<sol.numVecs; i++) {
      T(i,i) = evals[i].realpart;
    }
    A->Apply( *evecs, tempAevec );
    MVT::MvTimesMatAddMv( -1.0, *evecs, T, 1.0, tempAevec );
    MVT::MvNorm( tempAevec, normR );
  }

  // Print the results
  std::ostringstream os;
  os.setf(std::ios_base::right, std::ios_base::adjustfield);
  os << "Solver manager returned " << (returnCode == Converged ? "converged." : "unconverged.") << std::endl;
  os << std::endl;
  os << "------------------------------------------------------"<< std::endl;
  os << std::setw(16) << "Eigenvalue"
     << std::setw(18) << "Direct Residual"
     << std::endl;
  os<<"------------------------------------------------------" << std::endl;
  for (int i=0; i<sol.numVecs; i++) {
    os<< std::setw(16)<<evals[i].realpart
      << std::setw(18)<<normR[i]/evals[i].realpart
      << std::endl;
  }
  os<<"------------------------------------------------------" << std::endl;
  os << "gen_time = " << diag_tick - gen_tick << std::endl
     << "diag_time = " << end_tick - diag_tick << std::endl;
  rokko::machine_info(os);
  os << "the number of eigenvalues to compute=" << MyProblem->getNEV() << std::endl;
  printer.print(Errors, os.str());

#ifdef HAVE_MPI
  MPI_Finalize();
#endif
  return 0;
}
