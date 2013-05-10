      subroutine eigen_pbk_main2(n, z, nmz,
     &           d, v, nm, m, iblk,
     &           i, ss, tt)
!----
      use communication_sx, only : reduce_dbl
     &                           , get_loop_start
     &                           , get_loop_end
!----
      implicit double precision(a-h,o-z),integer(i-n)
!----
      integer :: n, nmz, nm, m, iblk, i
      real(8) :: z(nmz,*),d(*),v(nm,*)
      real(8) :: ss(*), tt(*)
!----
      include 'mpif.h'
!$    include 'omp_lib.h'
      include 'trd.h'
!----
      i_2 = 1
      i_3 = n
!----
      j_2 = get_loop_start(1,         size_of_col,my_col)
      j_3 = get_loop_end  (i+m-1-iblk,size_of_col,my_col)
!----
!$omp master
! sm:= 0
! ss:= 0
      ss(1:(i_3-i_2+1)*m+nsm*nsm) = 0.0d+0
!$omp end master
!----
!$omp barrier
! sm:= sm+lower(v^tv)
! ss(1:m,1:n):= ss+v(j_2:j_3,1:m)^t*z(j_2:j_3,1:n)
      call eigen_pbk_step1(z, nmz,
     &           v, nm, m, ss(1+nsm*nsm), ss(1),
     &           i_2,i_3,j_2,j_3)
!$omp barrier

!$omp master
!----
      call reduce_dbl(ss, tt, (i_3-i_2+1)*m+nsm*nsm,
     &                1, mpi_comm_col)
!----
! sm:= d*sm
      do m_0=1,m-1
!dir$ ivdep
         do i_0=m_0+1,m
            ss(i_0+(m_0-1)*nsm) = ss(i_0+(m_0-1)*nsm) * d(i+i_0-1)
         end do! i_0
      end do! m_0
!----
! ss:= d*ss
      i_4 = mod(i_3-i_2+1,4)+i_2
      do i_1=i_2,i_4-1                       ! 0
!dir$ ivdep
         do m_0=1,m
            ss((i_1+0-i_2)*m+m_0+nsm*nsm) =
     &      ss((i_1+0-i_2)*m+m_0+nsm*nsm) * d(i+m_0-1)
         end do! m_0
      end do! i_1
      do i_1=i_4,i_3,4                  ! 3
!dir$ ivdep
         do m_0=1,m
            ss((i_1+0-i_2)*m+m_0+nsm*nsm) =
     &      ss((i_1+0-i_2)*m+m_0+nsm*nsm) * d(i+m_0-1)
            ss((i_1+1-i_2)*m+m_0+nsm*nsm) =
     &      ss((i_1+1-i_2)*m+m_0+nsm*nsm) * d(i+m_0-1)
            ss((i_1+2-i_2)*m+m_0+nsm*nsm) =
     &      ss((i_1+2-i_2)*m+m_0+nsm*nsm) * d(i+m_0-1)
            ss((i_1+3-i_2)*m+m_0+nsm*nsm) =
     &      ss((i_1+3-i_2)*m+m_0+nsm*nsm) * d(i+m_0-1)
         end do! m_0
      end do! i_1
!$omp end master
!----
!$omp barrier
! v:= v*(i-sm)^{-1}
! z(j_2:j_3,1:n):= z + v(j_2:j_3,1:m)*ss(1:m,1:n)
      call eigen_pbk_step2( z, nmz,
     &            v, nm, m,
     &            ss(1+nsm*nsm), ss(1),
     &            i_2,i_3,j_2,j_3 )
!$omp barrier
!----
      return
      end subroutine ! eigen_pbk_main2
