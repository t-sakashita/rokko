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

      n = 4

      do j = 1, n
         do i = 1, n
            a(i,j)= n - max(i,j) + 1
         enddo
      enddo

      call dsyev('v','u', n, a, lda, w, work, n*n, info)

      write(*,'(" Eigenvalues:",4f10.5)') w(1:n)

      write(*,'(" Eigenvectors:")') 
      write(*,'(f10.5)') a(1:n,1:n)

      end
