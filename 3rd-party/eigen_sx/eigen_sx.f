!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  n   : dimension of matrix
!  a   : matrix
!  nma : size of array "a"
!  w   : array for eigen value
!  z   : array for eigen vector
!  d   : array to store a main diagonal element
!  e   : array to store a sub diagonal element
!  nme : size of array "e"
!  m0  : coefficient of blocking
!  ifl : switch flag ( 0:eigen-value and eigen-vector )
!                    ( 1:eigen-value only )
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine eigen_sx(n,a,nma,w,z,d,e,nme,m0,ifl)
      implicit none

      real(8) :: a(*),z(*),d(*),e(*),w(*)
      real(8),allocatable :: work(:)
      real(8) :: hs0, hs1
      integer :: n, nma, nme, m0, ifl
      integer :: m_forward, m_backward, istat
!
      include 'mpif.h'
      include 'trd.h'
!
      hs0=MPI_Wtime()

      m_forward = m0
      call eigen_prd(n, a(1), nma, d(1), e(1), nme, m_forward)

      allocate(work(2*nme), stat=istat)
      if(istat.ne.0) then
         print*,"Memory exhausted in eigen_sx routine"
         call flush(6)
         stop
      endif

      work(0*nme+1:0*nme+n-1) = e(0*nme+2:0*nme+n)
      work(0*nme+n) = 0
      work(1*nme+1:1*nme+n-2) = e(1*nme+3:1*nme+n)
      work(1*nme+n-1) = 0
      work(1*nme+n) = 0
      w(1:n)=d(1:n)

      call eigen_dcx(n, w(1), work(1), nme, z(1), nma)
      deallocate(work)

      if(ifl==0) then
         m_backward = 128
         call eigen_pbk(n, a(1), nma, z(1), nma, e(1+nme), m_backward)
      endif

c #ifdef TIMER
c       hs1=MPI_Wtime()
c       if(myrank==1) then
c          print*," "
c          print*,"Total Execution Time of eigen_sx =",hs1-hs0,"(sec)"
c       endif
c       call flush(6)
c #endif
*
      return
      end

