!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! n     : dimension of matrix
! a     : matrix
! nma   : size of array "a" 
! d_out : array to store a main diagonal element 
! e_out : array to store a main diagonal element 
! nme   : size of array "e_out" 
! m0    : coefficient of blocking 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine eigen_prd(n, a, nma,
     &                 d_out, e_out, nme, m0)
!----
      use communication_sx, only : eigen_init
     &                           , eigen_free
     &                           , datacast_dbl
     &                           , datacast_dbl2
!----
      implicit none

      include 'mpif.h'

      integer                :: n, nma, nme, m0
      integer                :: nm, ne, m
      real(8)                :: a(*)
      real(8)                :: d_out(*), e_out(*)
      real(8)                :: d1, d2

      include 'trd.h'
      include 'cstab.h'

      integer                :: nx

      integer, parameter     :: nm_max_l1 = 16*4
      integer, parameter     :: nm_max_l2 = 16*4*2

      integer :: lda, ldz, nz

      nm = nma
      ne = nme
      m  = m0

      call eigen_init(2)

      lda = nma
      ldz = (n-1)/size_of_col+1
      call cstab_get_optdim(ldz,9,nm_max_l1,nm_max_l2,nx)
      ldz = nx
      nz = (n-1)/size_of_row+1

      d1=mpi_wtime()
      call eigen_prd_main1(a, lda, d_out, e_out, ne, n, m)
      d2=mpi_wtime()

      call eigen_free(1)
#ifdef TIMER
!$omp master
      if(myrank==1)then
         print*,"Exectime of \"eigen_prd\" routine  =",(d2-d1),"(sec)"
      endif
      call flush(6)
!$omp end master
#endif

      return
      end subroutine ! eigen_prd

      subroutine eigen_prd_main1(a, nm, d_out, e_out, ne, n, m0)
      implicit none

      include 'mpif.h'

      integer                :: n, nm, ne, nv, m0
      real(8)                :: a(1:nm,*)
      real(8)                :: d_out(1:n)
      real(8)                :: e_out(1:ne,*)

      real(8), pointer       :: u_t(:), v_t(:), d_t(:)
      real(8), pointer       :: w(:)
      real(8), pointer       :: uv_col(:)
      real(8), pointer       :: uv_row(:)

      include 'trd.h'
      include 'cstab.h'
      include 'param.h'

      integer                :: nx

      integer                :: ierr, kx

      integer                :: offset1, offset2, offset3
      integer                :: offset4, offset5, offset6
      integer                :: offset7

      integer, parameter     :: nm_max_l1 = 16*4
      integer, parameter     :: nm_max_l2 = 16*6


      nx = (n-1)/size_of_col+1 +2
      call cstab_get_optdim(nx, 6, nm_max_l1, nm_max_l2, nv)
!---- if(myrank==1)then
!----    print*,"nv=",nv
!---- endif

      allocate(
     &         w(1:nm*m0+n_columns),
     &         v_t(1:max(2*nv,2*(n+3),2*8*m0)+n_columns),
     &         u_t(1:max(2*nv,2*(n+3),2*8*m0)+n_columns),
     &         d_t(1:nv+n_columns),
     &         uv_col(1:nv*2*m0+2*n_columns),
     &         uv_row(1:nv*2*m0+2*n_columns),
     &         stat=ierr)
      if ( ierr /= 0 ) then
         if ( myrank == 1 ) print*,"memory allocation error."
         call mpi_abort(mpi_comm_world, 1, ierr)
      end if

      w = zero
      v_t = zero
      u_t = zero
      d_t = zero
      uv_col = zero
      uv_row = zero

      kx = nv*m0+n_columns
      call cstab_adjust_base(uv_col(1), a(1,1),offset1)
      call cstab_adjust_base(uv_col(1+kx), a(1,1),offset3)
      call cstab_adjust_base(uv_row(1), a(1,1),offset2)
      call cstab_adjust_base(uv_row(1+kx), a(1,1),offset4)
      call cstab_adjust_base(u_t(1),a(1,1),offset5)
      call cstab_adjust_base(v_t(1),a(1,1),offset6)
      call cstab_adjust_base(w(1),a(1,1),offset7)
      kx = (l1_window/8)
!----&           +(l1_window)
!----&           +(l1_lsize/8)
     &           +(l1_lsize)
     &           +(l2_lsize/8)
      offset1 = offset1+kx*1
      offset2 = offset2+kx*2
      offset3 = offset3+kx*3
      offset4 = offset4+kx*4
      offset5 = offset5+kx*5
      offset6 = offset6+kx*6
      offset7 = offset7+kx*7
      call cstab_round_offset(offset1)
      call cstab_round_offset(offset2)
      call cstab_round_offset(offset3)
      call cstab_round_offset(offset4)
      call cstab_round_offset(offset5)
      call cstab_round_offset(offset6)
      call cstab_round_offset(offset7)

!----
      kx = nv*m0+n_columns

!$omp parallel
      call eigen_prd_main2(a, nm, d_out, e_out, ne,
     &              n, m0,
     &              w(1+offset7),
     &              uv_col(1   +offset1),    ! u_x(1:nv,m)
     &              uv_row(1   +offset2),    ! u_y(1:nv,m)
     &              uv_col(1+kx+offset3),    ! v_x(1:nv,m)
     &              uv_row(1+kx+offset4),    ! v_y(1:nv,m)
     &              u_t(1+offset5),
     &              v_t(1+offset6),
     &              d_t(1),
     &              nv)
!$omp end parallel

      deallocate(w)
      deallocate(v_t)
      deallocate(u_t)
      deallocate(d_t)
      deallocate(uv_col)
      deallocate(uv_row)

      return
      end subroutine ! eigen_prd_main1

      subroutine eigen_prd_main2(a, nm, d_out, e_out, ne,
     &              n, m0, w, u_x, u_y, v_x, v_y, u_t, v_t, d_t, nv)
      implicit none

      integer, intent(in)    :: n, nm, ne, nv, m0
      real(8)                :: a(1:nm,*)
      real(8)                :: d_out(1:n), e_out(1:ne,*)
      real(8)                :: w(1:nm,*)
      real(8)                :: u_x(1:nv,*), u_y(1:nv,*)
      real(8)                :: v_x(1:nv,*), v_y(1:nv,*)
      real(8)                :: u_t(1:2*nv), v_t(1:2*nv)
      real(8)                :: d_t(1:nv)

      include 'trd.h'

      real(8)                :: c(mband,mband)
      integer                :: i, m

      integer                :: k_1, k_2, k_3, k_4

      integer                :: m1, mm
      integer                :: i_block, i_base

      include 'mpif.h'
!$    include 'omp_lib.h'

      real(8)                :: d1,d2,dd(100)

!----
! initialization
!----
      dd = 0.0d+00
      m  = m0

!$omp master

#ifdef DETAIL
      if ( myrank == 1 ) then
         print*,"num.of.process=",nprocs,"(",size_of_col,size_of_row,")"
!$       print*,"num.of.threads=",omp_get_num_threads()
      endif
      call flush(6)
#endif

      call eigen_prd_init(a(1,1), nm, n,
     &        d_out(1), e_out(1,1), ne,
     &        u_t(1), v_t(1), nv)

!$omp end master

      if ( m > n ) then
         m = n
      end if
!----
! main iteration to reduce the input matrix into a penta-diagonal
!----
      mm = ((n-(mband+mod(n,mband)))-1)/m+1+1

!$omp barrier

      do i_block=mm,2,-1

!$omp barrier

         i_base = (i_block-2)*m+(mband+mod(n,mband))
         m1     = min(m,n-i_base)

!$omp master

         call eigen_prd_d_preserve(a(1,1), w(1,1), nm,
     &           d_t(1), u_x(1,1), u_y(1,1), v_x(1,1), v_y(1,1), nv,
     &           m1, i_base, i_block)

!$omp end master

         k_2 = m1
         k_3 = max(1,3*(2-i_block))

!$omp barrier

         do k_1=k_2,k_3,-mband; k_4=k_1-mband+1

!$omp barrier

            i = i_base+k_1
!----
! u=...
!----
!$omp master

            d1=mpi_wtime()
            call eigen_prd_u(
     &                   w(1,k_4), nm,
     &                   u_x(1,k_4), u_y(1,k_4), nv,
     &                   u_t(1), v_t(1), i,
     &                   c(1,1), e_out(1,1), ne)
            d2=mpi_wtime()
            dd(1)=dd(1)+(d2-d1)

!$omp end master

!$omp barrier

!$omp master
            d1=mpi_wtime()
!$omp end master
!----
! v:=Au
!----
            call eigen_prd_Au(
     &                    a(1,1), nm,
     &                    u_x(1,k_4), u_y(1,k_4), v_x(1,k_4), nv,
     &                    u_t(1), v_t(1), d_t(1),
     &                    i, i_base, m1)

!$omp barrier

!$omp master
            d2=mpi_wtime()
            dd(2)=dd(2)+(d2-d1)
!----
! v=v-(uv+vu)u
!----
            d1=mpi_wtime()
            call eigen_prd_2update_v(
     &                        u_x(1,1), v_x(1,1), nv,
     &                        u_t(1), v_t(1),
     &                        i, i_base, m1)
            d2=mpi_wtime()
            dd(5)=dd(5)+(d2-d1)

!----
! v':= v-((u,v)/2|u|^2)u
!----
            d1=mpi_wtime()
            call eigen_prd_v(
     &                   u_x(1,k_4), v_x(1,k_4), v_y(1,k_4), nv,
     &                   u_t(1), v_t(1), c(1,1), i)
            d2=mpi_wtime()
            dd(4)=dd(4)+(d2-d1)
!----
! w':= w-uv-vu
!----
            d1=mpi_wtime()
            call eigen_prd_2update_local(
     &              w(1,1), nm,
     &              u_x(1,1), u_y(1,1), v_x(1,1), v_y(1,1), nv,
     &              i_base, i)
            d2=mpi_wtime()
            dd(6)=dd(6)+(d2-d1)

!$omp end master

!$omp barrier

         end do! k_1

!$omp barrier

!$omp master

         call eigen_prd_d_restore(a(1,1), w(1,1), nm,
     &           d_t(1), m1, i_base)
!$omp end master

!$omp barrier

!$omp master
         d1=mpi_wtime()
!$omp end master
!----
! a:=a-v^tu-ut^v
!----
         call eigen_prd_2update(
     &           a(1,1), nm,
     &           u_x(1,1), u_y(1,1), v_x(1,1), v_y(1,1), nv,
     &           m1, i_base)

!$omp barrier

!$omp master
         d2=mpi_wtime()
         dd(3)=dd(3)+(d2-d1)

!$omp end master

!$omp barrier

      end do! i_block

!$omp barrier

!$omp master

      call eigen_prd_final(a(1,1), nm, n,
     &        d_out(1), e_out(1,1), ne, u_t(1))

#ifdef DETAIL
      if(myrank==1)then
         print*," "
         print*,"detail of exectime in \"eigen_prd\""
         print*,"   calc (u,beta)    ",dd(1),"(sec)"
         print*,"   mat-vec (au)     ",dd(2),"(sec)"
         print*,"   2update (a-uv-vu)",dd(3),"(sec)"
         print*,"   calc v           ",dd(4),"(sec)"
         print*,"   v=v-(uv+vu)u     ",dd(5),"(sec)"
         print*,"   uv post reduction",dd(6),"(sec)"
      endif
      call flush(6)
#endif

!$omp end master

      return
      end subroutine ! eigen_prd_main2

