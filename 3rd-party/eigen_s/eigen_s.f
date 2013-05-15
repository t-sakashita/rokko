!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  n   : dimension of matrix
!  a   : matrix
!  lda : size of array "a"
!  w   : array for eigen value
!  z   : array for eigen vector
!  ldz : size of array "z"
!  d   : array to store a main diagonal element
!  e   : array to store a sub diagonal element
!  m0  : coefficient of blocking
!  ifl : switch flag ( 0:eigen-value and eigen-vector )
!                    ( 1:eigen-value only )
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine eigen_s(n, a, lda, w, z, ldz, m0, ifl)
      implicit none

      integer, intent(in)    :: n, lda, ldz, m0, ifl
      real(8), intent(inout) :: a(lda,*), w(*), z(ldz,*)
      real(8), pointer       :: d(:), e(:), work(:)
      real(8)                :: hs0, hs1
      integer                :: m_forward, m_backward, my_rank,
     &                          info, ierr
      include 'mpif.h'

      allocate(d(1:n), e(1:n), work(1:n))

      hs0 = mpi_wtime()
      call mpi_comm_rank(mpi_comm_world, my_rank, ierr)
!----
      m_forward = m0
      call eigen_trd(n, a(1,1), lda, d(1), e(1), m_forward)
!----
      work(1:n) = e(1:n)
      w(1:n)  = d(1:n)
!----
      call eigen_dc(n, w(1), work(2), z(1,1), ldz, info)
!----
      if(ifl==0) then
         m_backward = 128
         call eigen_tbk(n, a(1,1), lda, z(1,1), ldz, e(1), m_backward)
      endif
!----
      hs1 = mpi_wtime()
c #ifdef TIMER
c       if(my_rank==0)then
c          print*," "
c          print*,"Total Execution Time of eigen_s = ",hs1-hs0,"(sec)"
c       endif
c       call flush(6)
c #endif
!----
      deallocate(d, e, work)
!----
      end subroutine

