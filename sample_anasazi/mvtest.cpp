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

#include <vector>


//=======================================================================
void myprint(Epetra_MultiVector& vec)  {
  int MyPID = vec.Map().Comm().MyPID();
  int NumProc = vec.Map().Comm().NumProc();

  for (int iproc=0; iproc < NumProc; iproc++) {
    if (MyPID==iproc) {
      int NumVectors1 = vec.NumVectors();
      int NumMyElements1 = vec.Map(). NumMyElements();
      int MaxElementSize1 = vec.Map().MaxElementSize();
      int * FirstPointInElementList1 = NULL;
      if (MaxElementSize1!=1) FirstPointInElementList1 = vec.Map().FirstPointInElementList();
      double ** A_Pointers = vec.Pointers();

      if (MyPID==0) {
  std::cout.width(8);
  std::cout <<  "     MyPID"; std::cout << "    ";
  std::cout.width(12);
  if (MaxElementSize1==1)
    std::cout <<  "GID  ";
  else
    std::cout <<  "     GID/Point";
  for (int j = 0; j < NumVectors1 ; j++)
    {
      std::cout.width(20);
      std::cout <<  "Value  ";
    }
  std::cout << std::endl;
      }
      for (int i=0; i < NumMyElements1; i++) {
  for (int ii=0; ii< vec.Map().ElementSize(i); ii++) {
       int iii;
    std::cout.width(10);
    std::cout <<  MyPID; std::cout << "    ";
    std::cout.width(10);
    if (MaxElementSize1==1) {
      if(vec.Map().GlobalIndicesInt())
      {
#ifndef EPETRA_NO_32BIT_GLOBAL_INDICES
      int * MyGlobalElements1 = vec.Map().MyGlobalElements();
      std::cout << MyGlobalElements1[i] << "    ";
#else
      throw vec.ReportError("Epetra_MultiVector::Print: ERROR, GlobalIndicesInt but no API for it.",-1);
#endif
      }
      else if(vec.Map().GlobalIndicesLongLong())
      {
#ifndef EPETRA_NO_64BIT_GLOBAL_INDICES
      long long * MyGlobalElements1 = vec.Map().MyGlobalElements64();
      std::cout << MyGlobalElements1[i] << "    ";
#else
      throw vec.ReportError("Epetra_MultiVector::Print: ERROR, GlobalIndicesLongLong but no API for it.",-1);
#endif
      }
      else
        throw vec.ReportError("Epetra_MultiVector::Print ERROR, Don't know map global index type.",-1);

       iii = i;
       }
    else {
      if(vec.Map().GlobalIndicesInt())
      {
#ifndef EPETRA_NO_32BIT_GLOBAL_INDICES
      int * MyGlobalElements1 = vec.Map().MyGlobalElements();
      std::cout <<  MyGlobalElements1[i]<< "/" << ii << "    ";
#else
      throw vec.ReportError("Epetra_MultiVector::Print: ERROR, GlobalIndicesInt but no API for it.",-1);
#endif
      }
      else if(vec.Map().GlobalIndicesLongLong())
      {
#ifndef EPETRA_NO_64BIT_GLOBAL_INDICES
      long long * MyGlobalElements1 = vec.Map().MyGlobalElements64();
      std::cout <<  MyGlobalElements1[i]<< "/" << ii << "    ";
#else
      throw vec.ReportError("Epetra_MultiVector::Print: ERROR, GlobalIndicesLongLong but no API for it.",-1);
#endif
      }
      else
	throw vec.ReportError("Epetra_MultiVector::Print ERROR, Don't know map global index type.",-1);

         iii = FirstPointInElementList1[i]+ii;
       }
    for (int j = 0; j < NumVectors1 ; j++)
      {
        std::cout.width(20);
        std::cout <<  A_Pointers[j][iii];
      }
    std::cout << std::endl;
  }
      }
      std::cout << std::flush;
    }

    // Do a few global ops to give I/O a chance to complete
    vec.Map().Comm().Barrier();
    vec.Map().Comm().Barrier();
    vec.Map().Comm().Barrier();
  }
  return;
}

int main(int argc, char *argv[]) {

  using std::endl;

#ifdef HAVE_MPI
  // Initialize MPI
  //
  MPI_Init(&argc,&argv);
  int my_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

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

  int N = 8;
  Epetra_Map Map(N, 0, Comm);
  int NumMyElements = Map.NumMyElements();
  int blockSize = 2;
  std::vector<int> MyGlobalElements(NumMyElements);
  Map.MyGlobalElements(&MyGlobalElements[0]);

  Teuchos::RCP<Epetra_MultiVector> ivec = Teuchos::rcp( new Epetra_MultiVector(Map, blockSize) );
  ivec->Random();

  ivec->PutScalar(my_rank);
  ivec->Print(std::cout);
  myprint(*ivec);

#ifdef HAVE_MPI
  MPI_Finalize();
#endif
  return 0;
}
