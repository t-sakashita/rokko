*
      PROGRAM SAMPLE_PDSYEV_CALL
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
*     .. Parameters ..
      INTEGER            LWORK, MAXN
      DOUBLE PRECISION   ZERO
      PARAMETER          ( LWORK = 264, MAXN = 100, ZERO = 0.0D+0 )
      INTEGER            LDA
      DOUBLE PRECISION   MONE
      INTEGER            MAXPROCS
      PARAMETER          ( LDA = MAXN, MONE = -1.0D+0, MAXPROCS = 512 )
*     ..
*     .. Local Scalars ..
      INTEGER            CONTEXT, I, IAM, INFO, MYCOL, MYROW, N, NB,
     $                   NPCOL, NPROCS, NPROW
*     ..
*     .. Local Arrays ..
      INTEGER            DESCA( 50 ), DESCZ( 50 )
      DOUBLE PRECISION   A( LDA, LDA ), W( MAXN ),
     $                   WORK( LWORK ), Z( LDA, LDA )
*     ..
*     .. External Subroutines ..
      EXTERNAL           BLACS_EXIT, BLACS_GET, BLACS_GRIDEXIT,
     $                   BLACS_GRIDINFO, BLACS_GRIDINIT, BLACS_PINFO,
     $                   BLACS_SETUP, DESCINIT, PDLAMODHILB, PDLAPRNT,
     $                   PDSYEV
*     ..
*     .. Executable Statements ..
*
*
*     Set up the problem
*
      N = 4
      NB = 1
      NPROW = 2
      NPCOL = 2
*
*     Initialize the BLACS
*
      CALL BLACS_PINFO( IAM, NPROCS )
      IF( ( NPROCS.LT.1 ) ) THEN
         CALL BLACS_SETUP( IAM, NPROW*NPCOL )
      END IF
*
*     Initialize a single BLACS context
*
      CALL BLACS_GET( -1, 0, CONTEXT )
      CALL BLACS_GRIDINIT( CONTEXT, 'R', NPROW, NPCOL )
      CALL BLACS_GRIDINFO( CONTEXT, NPROW, NPCOL, MYROW, MYCOL )
*
*     These are basic array descriptors
*
      CALL DESCINIT( DESCA, N, N, NB, NB, 0, 0, CONTEXT, LDA, INFO )
      CALL DESCINIT( DESCZ, N, N, NB, NB, 0, 0, CONTEXT, LDA, INFO )
*
*     Build a matrix that you can create with
*     a one line matlab command:  hilb(n) + diag([1:-1/n:1/n])
*
      CALL PDLAMODHILB( N, A, 1, 1, DESCA, INFO )
*
*     Ask PDSYEV to compute the entire eigendecomposition
*
      CALL PDSYEV( 'V', 'U', N, A, 1, 1, DESCA, W, Z, 1, 1,
     $             DESCZ, WORK, LWORK, INFO )
*
*     Print out the eigenvalues and eigenvectors
*
      IF( MYROW.EQ.0 .AND. MYCOL.EQ.0 ) THEN
         write(*,'(" Eigenvalues:",5f10.5)') W(1:N)
      END IF
      CALL PDLAPRNT( N, N, Z, 1, 1, DESCZ, 0, 0, 'Z', 6, WORK )
*
      CALL BLACS_GRIDEXIT( CONTEXT )
      CALL BLACS_EXIT( 0 )
*
      STOP
      END
*
      SUBROUTINE PDLAMODHILB( N, A, IA, JA, DESCA, INFO )
*
*  -- ScaLAPACK routine (version 1.2) --
*     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
*     and University of California, Berkeley.
*     May 10, 1996
*
*
*
*
*     .. Parameters ..
      INTEGER            BLOCK_CYCLIC_2D, DLEN_, DT_, CTXT_, M_, N_,
     $                   MB_, NB_, RSRC_, CSRC_, LLD_
      PARAMETER          ( BLOCK_CYCLIC_2D = 1, DLEN_ = 9, DT_ = 1,
     $                   CTXT_ = 2, M_ = 3, N_ = 4, MB_ = 5, NB_ = 6,
     $                   RSRC_ = 7, CSRC_ = 8, LLD_ = 9 )
      DOUBLE PRECISION   ONE
      PARAMETER          ( ONE = 1.0D+0 )
*     ..
*     .. Scalar Arguments ..
      INTEGER            IA, INFO, JA, N
*     ..
*     .. Array Arguments ..
      INTEGER            DESCA( * )
      DOUBLE PRECISION   A( * )
*     ..
*     .. Local Scalars ..
      INTEGER            I, J, MYCOL, MYROW, NPCOL, NPROW
*     ..
*     .. External Subroutines ..
      EXTERNAL           BLACS_GRIDINFO, PDELSET
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          DBLE
*     ..
*     .. Executable Statements ..
*
*
*     The matlab code for a real matrix is:
*         hilb(n) + diag( [ 1:-1/n:1/n ] )
*     The matlab code for a complex matrix is:
*         hilb(N) + toeplitz( [ 1 (1:(N-1))*i ] )
*
*       This is just to keep ftnchek happy
      IF( BLOCK_CYCLIC_2D*CSRC_*CTXT_*DLEN_*DT_*LLD_*MB_*M_*NB_*N_*
     $    RSRC_.LT.0 )RETURN
*
      INFO = 0
*
      CALL BLACS_GRIDINFO( DESCA( CTXT_ ), NPROW, NPCOL, MYROW, MYCOL )
*
*
      IF( IA.NE.1 ) THEN
         INFO = -3
      ELSE IF( JA.NE.1 ) THEN
         INFO = -4
      END IF
*
      DO 20 J = 1, N
         DO 10 I = 1, N
            CALL PDELSET( A, I, J, DESCA, dble(N - MAX(I,J) + 1))
   10    CONTINUE
   20 CONTINUE
*
*
      RETURN
*
*     End of PDLAMODHLIB
*
      END
