#include "AnasaziConfigDefs.hpp"
#include "AnasaziBasicEigenproblem.hpp"
#include "AnasaziSimpleLOBPCGSolMgr.hpp"
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

#include <rokko/utility/laplacian_matrix.hpp>
#include <vector>

class LaplacianOp : public Epetra_Operator {
 public:
  //! @name Constructor/Destructor
  //@{

  //! Basic constructor.  Accepts reference-counted pointer to an Epetra_Operator.
  LaplacianOp(int dim, const MPI_Comm& comm = MPI_COMM_WORLD) : dim_(dim), comm_(comm), ep_comm(comm), ep_map(dim, 0, ep_comm) {
    //comm_ = MPI_COMM_WORLD;
    MPI_Comm_size(comm_, &nprocs);
    MPI_Comm_rank(comm_, &myrank);

    int tmp = dim_ / nprocs;
    int rem = dim % nprocs;
    num_local_rows_ = (dim + nprocs - myrank - 1) / nprocs;
    start_row_ = tmp * myrank + std::min(rem, myrank);
    end_row_ = start_row_ + num_local_rows_ - 1;

    if (start_row_ == 0)  is_first_proc = true;
    else is_first_proc = false;
    
    if (end_row_ == (dim-1))  is_last_proc = true;
    else is_last_proc = false;

    end_k_ = num_local_rows_ - 1;
    
    std::cout << "myrank=" << myrank << " start_row=" << start_row_ << " end_row=" << end_row_ << std::endl;
    std::cout << "myrank=" << myrank << " num_local_rows_=" << num_local_rows_ << std::endl;
  }

  //! Destructor
  ~LaplacianOp() {};
  //@}

  virtual int SetUseTranspose(bool UseTranspose) { return 0; };
  //@}

  void multiply(const double* x, double* y) const {
    MPI_Status *status;
    if (!is_first_proc) {
      //std::cout << "recv myrank=" << myrank << std::endl;
      MPI_Send(&x[0], 1, MPI_DOUBLE, myrank-1, 0, comm_);
      MPI_Recv(&buf1, 1, MPI_DOUBLE, myrank-1, 0, comm_, status);
      y[0] = - buf1 + 2 * x[0] - x[1];
    }
    else { // for the first process 0
      y[0] = x[0] - x[1];
    }
    if (!is_last_proc) {
      //std::cout << "send myrank=" << myrank << std::endl;
      MPI_Send(&x[end_k_], 1, MPI_DOUBLE, myrank+1, 0, comm_);
      MPI_Recv(&buf2, 1, MPI_DOUBLE, myrank+1, 0, comm_, status);
      y[end_k_] = - x[end_k_ - 1] + 2 * x[end_k_] - buf2;
    }
    else { // for the last process
      y[end_k_] = 2 * x[end_k_] - x[end_k_ - 1];
    }
    
    for (int k=1; k<end_k_; ++k) {  //  from 1 to end-1
      y[k] = - x[k-1] + 2 * x[k] - x[k+1];
    }
  }
  
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

    Y.PutScalar(0.);
    for (int i=0; i<numvectors; ++i) {
      const double* x = X[i];
      double* y = Y[i];
      multiply(x, y);
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
  virtual const char * Label() const { return "Laplacian matrix"; }

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
  Epetra_MpiComm ep_comm;
  Epetra_Map ep_map;
  int nprocs, myrank;
  int dim_, local_offset_, num_local_rows_;
  int start_row_, end_row_;
  int start_k_, end_k_;
  mutable double buf1, buf2;
  bool is_first_proc, is_last_proc;
};

int main(int argc, char *argv[]) {

  using std::endl;

#ifdef HAVE_MPI
  // Initialize MPI
  //
  MPI_Init(&argc,&argv);
#endif

  // Create an Epetra communicator
  //
#ifdef HAVE_MPI
  Epetra_MpiComm Comm(MPI_COMM_WORLD);
#else
  Epetra_SerialComm Comm;
#endif

  // Create an Anasazi output manager
  Anasazi::BasicOutputManager<double> printer;
  printer.stream(Anasazi::Errors) << Anasazi::Anasazi_Version() << endl << endl;

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

  int dim = 20;

  // Construct a Map that puts approximately the same number of
  // equations on each processor.
  Epetra_Map Map(dim, 0, Comm);

  // Get update list and number of local equations from newly created Map.
  int NumMyElements = Map.NumMyElements();

  std::vector<int> MyGlobalElements(NumMyElements);
  Map.MyGlobalElements(&MyGlobalElements[0]);

  // Create an Epetra_Matrix
  Teuchos::RCP<LaplacianOp> A = Teuchos::rcp( new LaplacianOp(dim, MPI_COMM_WORLD) );

  //************************************
  // Call the LOBPCG solver manager
  //***********************************
  //  Variables used for the LOBPCG Method
  const int    nev       = 10;
  const int    blockSize = 5;
  const int    maxIters  = 500;
  const double tol       = 1.0e-8;

  typedef Epetra_MultiVector MV;
  typedef Epetra_Operator OP;
  typedef Anasazi::MultiVecTraits<double, Epetra_MultiVector> MVT;

  // Create an Epetra_MultiVector for an initial vector to start the solver.
  // Note:  This needs to have the same number of columns as the blocksize.
  Teuchos::RCP<Epetra_MultiVector> ivec = Teuchos::rcp( new Epetra_MultiVector(Map, blockSize) );
  ivec->Random();

  // Create the eigenproblem.
  Teuchos::RCP<Anasazi::BasicEigenproblem<double, MV, OP> > MyProblem =
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
  MyPL.set( "Which", which );
  MyPL.set( "Block Size", blockSize );
  MyPL.set( "Maximum Iterations", maxIters );
  MyPL.set( "Convergence Tolerance", tol );

  // Create the solver manager
  Anasazi::SimpleLOBPCGSolMgr<double, MV, OP> MySolverMan(MyProblem, MyPL);
  //Anasazi::BlockKrylovSchurSolMgr<double, MV, OP> MySolverMan(MyProblem, MyPL);

  // Solve the problem
  Anasazi::ReturnType returnCode = MySolverMan.solve();

  // Get the eigenvalues and eigenvectors from the eigenproblem
  Anasazi::Eigensolution<double,MV> sol = MyProblem->getSolution();
  std::vector<Anasazi::Value<double> > evals = sol.Evals;
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
  os<<"Solver manager returned " << (returnCode == Anasazi::Converged ? "converged." : "unconverged.") << endl;
  os<<endl;
  os<<"------------------------------------------------------"<<endl;
  os<<std::setw(16)<<"Eigenvalue"
    <<std::setw(18)<<"Direct Residual"
    <<endl;
  os<<"------------------------------------------------------"<<endl;
  for (int i=0; i<sol.numVecs; i++) {
    os<<std::setw(16)<<evals[i].realpart
      <<std::setw(18)<<normR[i]/evals[i].realpart
      <<endl;
  }
  os<<"------------------------------------------------------"<<endl;

  os << "theoretical eigenvalues:" << std::endl;
  for(int k=0; k<dim; ++k) {
    os << rokko::laplacian_matrix::eigenvalue(dim, k) << " ";
  }
  os << std::endl;
  printer.print(Anasazi::Errors,os.str());

  os << "Apply_test" << std::endl;
  Epetra_MultiVector X( Map, 2 ), Y( Map, 2 );
  X.PutScalar(0.);
  X.ReplaceGlobalValue(4, 1, 1.);
  X.Print(std::cout);
  A->Apply( X, Y );
  Y.Print(std::cout);
  printer.print(Anasazi::Debug, os.str());

#ifdef HAVE_MPI
  MPI_Finalize();
#endif
  return 0;
}
