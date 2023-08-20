#include "AnasaziConfigDefs.hpp"
#include "AnasaziBasicEigenproblem.hpp"
#include "AnasaziLOBPCGSolMgr.hpp"
//#include "AnasaziBlockKrylovSchurSolMgr.hpp"
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

#include <rokko/utility/heisenberg_hamiltonian_mpi.hpp>
#include <rokko/utility/lattice.hpp>
#include <rokko/utility/math.hpp>
#include <rokko/utility/machine_info.hpp>

class HeisenbergOp : public Epetra_Operator {
 public:
  //! @name Constructor/Destructor
  //@{

  //! Basic constructor.  Accepts reference-counted pointer to an Epetra_Operator.
  HeisenbergOp(const MPI_Comm& comm, int L, const std::vector<std::pair<int, int>>& lattice) : comm_(comm), L_(L), lattice_(lattice), ep_comm(comm), ep_map(1 << L_, 0, ep_comm) {
    int nproc;
    MPI_Comm_size(MPI_COMM_WORLD, &nproc);
    const int p = rokko::find_power_of_two(nproc);
    int local_N = 1 << (L-p);
    buffer_.assign(local_N, 0);
  }

  //! Destructor
  ~HeisenbergOp() {}
  //@}

  virtual int SetUseTranspose(bool UseTranspose) { return 0; };
  //@}
  
  //! @name Mathematical functions
  //@{ 

    //! Returns the result of a Epetra_Operator applied to a Epetra_MultiVector X in Y.
    /*! 
    \param In
	   X - A Epetra_MultiVector of dimension NumVectors to multiply with matrix.
    \param Out
	   Y -A Epetra_MultiVector of dimension NumVectors containing result.

    \return Integer error code, set to 0 if successful.
  */
  virtual int Apply(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const {
    const int numvectors = X.NumVectors();
    //std::cout << "numvectors=" << numvectors << std::endl;
    //std::cout << "X.GlobalLength()=" << X.GlobalLength() << std::endl;
    //std::cout << "X.MyLength()=" << X.MyLength() << std::endl;

    Y.PutScalar(0);
    for (int i=0; i<numvectors; ++i) {
      const double* x = X[i];
      double* y = Y[i];
      rokko::heisenberg_hamiltonian::multiply(comm_, L_, lattice_, x, y, buffer_.data());
    }
    //std::cout << "X=" << X << std::endl;
    return 0;
  }

    //! Returns the result of a Epetra_Operator inverse applied to an Epetra_MultiVector X in Y.
    /*! 
    \param In
	   X - A Epetra_MultiVector of dimension NumVectors to solve for.
    \param Out
	   Y -A Epetra_MultiVector of dimension NumVectors containing result.

    \return Integer error code, set to 0 if successful.

    \warning In order to work with AztecOO, any implementation of this method must 
              support the case where X and Y are the same object.
  */
  virtual int ApplyInverse(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const { return 0; }

    //! Returns the infinity norm of the global matrix.
    /* Returns the quantity \f$ \| A \|_\infty\f$ such that
       \f[\| A \|_\infty = \max_{1\lei\lem} \sum_{j=1}^n |a_{ij}| \f].

       \warning This method must not be called unless HasNormInf() returns true.
    */ 
  virtual double NormInf() const { return 0; }
  //@}
  
  //! @name Attribute access functions
  //@{ 

    //! Returns a character string describing the operator
  virtual const char * Label() const { return "Heisenberg Hamiltonian"; }

    //! Returns the current UseTranspose setting.
  virtual bool UseTranspose() const { return false; }

    //! Returns true if the \e this object can provide an approximate Inf-norm, false otherwise.
  virtual bool HasNormInf() const  { return false; }

    //! Returns a pointer to the Epetra_Comm communicator associated with this operator.
  virtual const Epetra_Comm & Comm() const { return ep_comm; }

    //! Returns the Epetra_Map object associated with the domain of this operator.
  virtual const Epetra_Map & OperatorDomainMap() const { return ep_map; }

    //! Returns the Epetra_Map object associated with the range of this operator.
  virtual const Epetra_Map & OperatorRangeMap() const { return ep_map; }
  //@}

  //! @name Operator application method
  //@{

  /*! \brief This method takes the Anasazi::MultiVec \c X and
      applies the operator to it resulting in the Anasazi::MultiVec \c Y.
  */
  //@}

 private:
  MPI_Comm comm_;
  mutable std::vector<double> buffer_;
  int L_;
  std::vector<std::pair<int, int>> lattice_;
  Epetra_MpiComm ep_comm;
  Epetra_Map ep_map;
};

int main(int argc, char *argv[]) {

#ifdef HAVE_MPI
  // Initialize MPI
  int provided;
  MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &provided);
#endif
  MPI_Barrier(MPI_COMM_WORLD);
  const auto init_tick = MPI_Wtime();
  // Create an Epetra communicator
  //
#ifdef HAVE_MPI
  Epetra_MpiComm Comm(MPI_COMM_WORLD);
#else
  Epetra_SerialComm Comm;
#endif
  MPI_Barrier(MPI_COMM_WORLD);
  const auto initend_tick = MPI_Wtime();

  // Create an Anasazi output manager
  Anasazi::BasicOutputManager<double> printer;
  printer.stream(Anasazi::Errors) << Anasazi::Anasazi_Version() << std::endl << std::endl;

  const int len_ladder = (argc >= 2) ? std::stoi(argv[1]) : 5;
  const auto L = 2 * len_ladder;
  const auto lattice = rokko::create_ladder_lattice_1dim(len_ladder);
  rokko::output_lattice(printer.stream(Anasazi::Errors), lattice);
  const auto N = 1 << L;

  MPI_Barrier(MPI_COMM_WORLD);
  const auto gen_tick = MPI_Wtime();
  // Construct a Map that puts approximately the same number of
  // equations on each processor.
  Epetra_Map Map(N, 0, Comm);

  // Get update list and number of local equations from newly created Map.
  int NumMyElements = Map.NumMyElements();

  std::vector<int> MyGlobalElements(NumMyElements);
  Map.MyGlobalElements(MyGlobalElements.data());

  // Create an Epetra_Matrix
  Teuchos::RCP<HeisenbergOp> A = Teuchos::rcp( new HeisenbergOp(MPI_COMM_WORLD, L, lattice) );

  //************************************
  // Call the LOBPCG solver manager
  //***********************************
  //  Variables used for the LOBPCG Method
  MPI_Barrier(MPI_COMM_WORLD);
  const auto diag_tick = MPI_Wtime();
  //std::string which("LM");
  constexpr int    nev       = 1;
  constexpr int    blockSize = nev;
  constexpr int    maxIters  = 200;
  constexpr double tol       = 1.0e-8;

  using MV = Epetra_MultiVector;
  using OP = Epetra_Operator;
  using MVT = Anasazi::MultiVecTraits<double, Epetra_MultiVector>;

  // Create an Epetra_MultiVector for an initial vector to start the solver.
  // Note:  This needs to have the same number of columns as the blocksize.
  Teuchos::RCP<Epetra_MultiVector> ivec = Teuchos::rcp( new Epetra_MultiVector(Map, blockSize) );
  ivec->Random();

  // Create the eigenproblem.
  Teuchos::RCP<Anasazi::BasicEigenproblem<double, MV, OP>> MyProblem =
    Teuchos::rcp( new Anasazi::BasicEigenproblem<double, MV, OP>(A, ivec) );

  // Inform the eigenproblem that the operator A is symmetric
  MyProblem->setHermitian(true);

  // Set the number of eigenvalues requested
  MyProblem->setNEV( nev );

  // Inform the eigenproblem that you are finishing passing it information
  bool boolret = MyProblem->setProblem();
  if (boolret != true) {
    printer.print(Anasazi::Errors,"Anasazi::BasicEigenproblem::setProblem() returned an error.\n");
#ifdef HAVE_MPI
    MPI_Finalize();
#endif
    return -1;
  }

  // Create parameter list to pass into the solver manager
  Teuchos::ParameterList MyPL;
  //MyPL.set( "Which", which );
  //MyPL.set( "Block Size", blockSize );
  MyPL.set( "Maximum Iterations", maxIters );
  MyPL.set( "Convergence Tolerance", tol );

  // Create the solver manager
  Anasazi::LOBPCGSolMgr<double, MV, OP> MySolverMan(MyProblem, MyPL);
  //Anasazi::BlockKrylovSchurSolMgr<double, MV, OP> MySolverMan(MyProblem, MyPL);

  // Solve the problem
  Anasazi::ReturnType returnCode = MySolverMan.solve();
  MPI_Barrier(MPI_COMM_WORLD);
  const auto end_tick = MPI_Wtime();

  // Get the eigenvalues and eigenvectors from the eigenproblem
  Anasazi::Eigensolution<double,MV> sol = MyProblem->getSolution();
  std::vector<Anasazi::Value<double>> evals = sol.Evals;
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
  os << "Solver manager returned " << (returnCode == Anasazi::Converged ? "converged." : "unconverged.") << std::endl;
  os << std::endl;
  os << "------------------------------------------------------" << std::endl;
  os << std::setw(16) << "Eigenvalue"
     << std::setw(18) << "Direct Residual"
     << std::endl;
  os << "------------------------------------------------------" << std::endl;
  for (int i=0; i<sol.numVecs; i++) {
    os << std::setw(16) << evals[i].realpart
       << std::setw(18) << normR[i]/evals[i].realpart
       << std::endl;
  }
  os << "------------------------------------------------------" << std::endl;
  os << "init_time = " << initend_tick - init_tick << std::endl
     << "gen_time = " << diag_tick - gen_tick << std::endl
     << "diag_time = " << end_tick - diag_tick << std::endl;
  rokko::machine_info(os);
  printer.print(Anasazi::Errors,os.str());

#ifdef HAVE_MPI
  MPI_Finalize();
#endif
  return 0;
}
