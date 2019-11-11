*     This is a modified version of the example on the Netlib,
*     which can be found at
*     http://www.netlib.org/scalapack/examples/sample_pdsyev_call.f

*     The license of the original one is the modified BSD license.
*     http://www.netlib.org/scalapack/LICENSE
*
      PROGRAM SAMPLE_PDSYEVD_CALL
*
*  -- ScaLAPACK routine (version 1.2) --
*     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
*     and University of California, Berkeley.
*     May 10, 1996
*
*     This routine contains a sample call to PDSYEV.
*     When compiled and run, it produces output which can be
*     pasted directly into matlab.
*
      IMPLICIT NONE
*     .. Parameters ..
      INTEGER            N
      character(len=10) :: tmp_str
      integer arg_len, status
*     ..
*     .. Local Scalars ..
      INTEGER            IERR
      INTEGER            LWORK, LIWORK, LDA
      INTEGER            CONTEXT, I, INFO, MYCOL, MYROW,
     $                   NPCOL, NPROCS, NPROW
      INTEGER            MB, NB, BB, M_LOCAL, N_LOCAL
      DOUBLE PRECISION   INIT_TICK, GEN_TICK, DIAG_TICK, END_TICK
*     ..
*     .. Local Arrays ..
      INTEGER            DESC( 50 )
      DOUBLE PRECISION, pointer :: A(:,:),W(:),Z(:,:)
      DOUBLE PRECISION   TMP_WORK(1)
      INTEGER            TMP_IWORK(1)
      DOUBLE PRECISION, allocatable :: WORK(:)
      INTEGER, allocatable :: IWORK(:)
*     ..
*     .. External Subroutines ..
      EXTERNAL           BLACS_EXIT, BLACS_GET, BLACS_GRIDEXIT,
     $                   BLACS_GRIDINFO, BLACS_GRIDINIT, BLACS_PINFO,
     $                   BLACS_SETUP, DESCINIT,
     $                   PDSYEV
      INTEGER            NUMROC
*     ..
*     .. Executable Statements ..
*
*
*     Set up the problem
*
      include 'mpif.h'

      call mpi_init_thread( MPI_THREAD_MULTIPLE, i, ierr )
      if (command_argument_count().eq.1) then
         call get_command_argument(1, tmp_str, arg_len, status)
         read(tmp_str, *) n
      else
         write(*,'(A)') "Error: argument dimension is needed"
      stop
      endif
      call mpi_barrier(mpi_comm_world, ierr)
      INIT_TICK = MPI_WTIME()
      call MPI_Comm_Size(mpi_comm_world,nprocs,ierr)
      NPROW = INT(SQRT(NPROCS + 0.5))
      DO WHILE ( (NPROW .ne. 1) .AND. (MOD(NPROCS,NPROW).ne.0) )
         NPROW = NPROW - 1
      ENDDO
      NPCOL = NPROCS / NPROW
*     
*     Initialize a single BLACS context
*
      CALL BLACS_GET( -1, 0, CONTEXT )
      CALL BLACS_GRIDINIT( CONTEXT, 'R', NPROW, NPCOL )
      CALL BLACS_GRIDINFO( CONTEXT, NPROW, NPCOL, MYROW, MYCOL )
*     
*     These are basic array descriptors
*
      call mpi_barrier(mpi_comm_world, ierr)
      GEN_TICK = MPI_WTIME()
      MB = N / NPROW
      IF( MB.EQ.0 ) THEN
         MB = 1
      ENDIF
      NB = N / NPCOL
      IF( NB.EQ.0 ) THEN
         NB = 1
      ENDIF
      BB = MIN(MB, NB)
      M_LOCAL = NUMROC(N, BB, MYROW, 0, NPROW)
      N_LOCAL = NUMROC(N, BB, MYCOL, 0, NPCOL)
!      print*, "M_LOCAL=", M_LOCAL, " N_LOCAL=", N_LOCAL
      LDA = MAX(M_LOCAL, 1)
      allocate( A(M_LOCAL, N_LOCAL), Z(M_LOCAL, N_LOCAL), W(N) )
      CALL DESCINIT( DESC, N, N, BB, BB, 0, 0, CONTEXT, LDA, INFO )
      CALL MINIJ_MATRIX( N, A, DESC, INFO )
*      CALL PDLAPRNT( N, N, A, 1, 1, DESC, 0, 0, 'A', 6, PRNWORK )
*     
*     Ask PDSYEV to compute the entire eigendecomposition
*
      call mpi_barrier(mpi_comm_world, ierr)
      DIAG_TICK = MPI_WTIME()
      CALL PDSYEVD( 'V', 'U', N, A, 1, 1, DESC, W, Z, 1, 1,
     $             DESC, TMP_WORK, -1, TMP_IWORK, -1, INFO )
      LWORK = INT(TMP_WORK(1))
      allocate ( WORK(LWORK) )
      LIWORK = TMP_IWORK(1)
      allocate ( IWORK(LIWORK) )
      CALL PDSYEVD( 'V', 'U', N, A, 1, 1, DESC, W, Z, 1, 1,
     $             DESC, WORK, LWORK, IWORK, LIWORK, INFO )
      call mpi_barrier(mpi_comm_world, ierr)
      END_TICK = MPI_WTIME()
*     
*     Print out the eigenvalues and eigenvectors
*
      IF( MYROW.EQ.0 .AND. MYCOL.EQ.0 ) THEN
         write(*,'(" Eigenvalues:",5E12.5e2)') W(1:N)
         write(*,*) "init_time = ", GEN_TICK - INIT_TICK
         write(*,*) "gen_time = ", DIAG_TICK - GEN_TICK
         write(*,*) "diag_time = ", END_TICK - DIAG_TICK
      END IF

      CALL BLACS_GRIDEXIT( CONTEXT )
      call mpi_finalize(ierr)
      STOP
      END
*

*
      SUBROUTINE MINIJ_MATRIX( N, A, DESCA, INFO )
*  This is originally PDLAMODHILB in ScaLAPACK
*  -- ScaLAPACK routine (version 1.2) --
*     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
*     and University of California, Berkeley.
*     May 10, 1996
*     ..
*     .. Scalar Arguments ..
      INTEGER            N
*     ..
*     .. Array Arguments ..
      INTEGER            DESCA( * )
      DOUBLE PRECISION   A( * )
*     ..
*     .. Local Scalars ..
      INTEGER            I, J
*     ..
*     .. External Subroutines ..
      EXTERNAL           PDELSET
*
*     Create minij matrix
      DO 20 J = 1, N
         DO 10 I = 1, N
            CALL PDELSET( A, I, J, DESCA, dble(MIN(I,J)))
   10    CONTINUE
   20 CONTINUE

*     
*
      RETURN
*
*     End of MINIJ_MATRIX
*
      END
