!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! n      : dimension of matrix
! a      : matrix
! nma    : size of array "a"
! z      : array for eigen vector
! nmz    : size of array "z"
! e      : array to store a sub diagonal element
! m0     : coefficient of blocking
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine eigen_pbk(n, a, nma, z, nmz, e, m0)
!----
      use communication_sx, only : eigen_init
     &                           , eigen_free
!----
      implicit double precision(a-h,o-z),integer(i-n)
      real(8) :: a(*)
      real(8) :: z(*)
      real(8) :: e(*)

      real(8) , pointer :: d(:)
      real(8) , pointer :: v(:)
      real(8) , pointer :: tt(:), ss(:)

      real(8) :: d1, d2

      include 'mpif.h'
!$    include 'omp_lib.h'
      include 'trd.h'
      include 'cstab.h'
!----
      iblk = 2

      nma  = nma
      nmz  = nmz
      m    = min(nsm,m0)
      if ( m < 1 ) m = 1
!----
      call eigen_init(2)
!----
      na   = (n-1)/size_of_row+1
      na   = na  +mod(na  -1,2)
      call cstab_get_optdim(nma,9,16*4,16*6,nm)
!----
      allocate(
     &          d(1:n),
     &          v(1:max(nm*m,n)+n_columns),
     &          ss(1:na*m+nsm*nsm+n_columns),
     &          tt(1:na*m+nsm*nsm+n_columns),
     &          stat=i_stat)
      if(i_stat/=0)then
         if(myrank==1)print*,"memory allocation error."
         call mpi_abort(mpi_comm_world,1,ierr)
      endif
!----
      call cstab_adjust_base(v(1), z(1), i_v)
      call cstab_adjust_base(ss(1), z(1), i_s)
      call cstab_adjust_base(tt(1), z(1), i_t)
          kx = (l1_window/4)
!----&           +(l1_window)
!----&           +(l1_lsize/8)
     &           +(l1_lsize)
     &           +(l2_lsize/4)
      i_v = i_v + kx*1
      i_s = i_s + kx*2
      i_t = i_t + kx*3
      call cstab_round_offset(i_v)
      call cstab_round_offset(i_s)
      call cstab_round_offset(i_t)
!----
      call mpi_barrier(mpi_comm_world,ierr)
      d1 = mpi_wtime()
!$omp parallel
      call eigen_pbk_main1( n,
     $      a(1), nma,
     $      z(1), nmz,
     $      e(1), d(1), v(1+i_v), nm, m, iblk,
     $      ss(1+i_s), tt(1+i_t) )
!$omp end parallel
      call mpi_barrier(mpi_comm_world,ierr)
      d2 = mpi_wtime()
!----
      deallocate( d )
      deallocate( v )
      deallocate( ss )
      deallocate( tt )
!----
      call eigen_free(3)

#ifdef TIMER
      if(myrank==1) then
         print*,"Exectime of \"eigen_pbk\" routine  =",d2-d1,"(sec)"
      endif
#endif
!----
      return
      end subroutine
