      subroutine eigen_tbk_main2(n, z, nmz,
     $           d, v, nm, m, i, ss, tt, dcom, dx, dy)
!----
      use communication_s, only : reduce_dbl
     $                          , get_loop_start
     $                          , get_loop_end
!----
      implicit double precision(a-h,o-z),integer(i-n)
*
      integer :: n, nmz, nm, i
      real(8) :: z(nmz,*),d(*),v(nm,*)
      real(8) :: ss(*), tt(*), dcom, dx, dy, ds, de
*
      include 'mpif.h'
!$    include 'omp_lib.h'
      include 'trd.h'
*
*-
      i_2 = 1
      i_3 = n
*
      j_2 = get_loop_start(1,    size_of_col,my_col)
      j_3 = get_loop_end  (i+m-2,size_of_col,my_col)
*
!$omp master
! sm:= 0
! ss:= 0
          ss(1:(i_3-i_2+1)*m+nsm*nsm) = 0.0d+0
!$omp end master
*-
!$omp barrier
! sm:= sm+lower(v^tv)
! ss(1:m,1:n):= ss+v(j_2:j_3,1:m)^t*z(j_2:j_3,1:n)
      ds=mpi_wtime()
      call eigen_tbk_step1(z, nmz,
     $           v, nm, m, ss(1+nsm*nsm), ss(1),
     $           i_2,i_3,j_2,j_3)
!$omp barrier
      de=mpi_wtime()

!$omp master
      dx = dx + (de-ds)
*-
c #ifdef TIMER
c       call mpi_barrier(mpi_comm_col,ierr)
c       ds=mpi_wtime()
c #endif
      call reduce_dbl(ss, tt, (i_3-i_2+1)*m+nsm*nsm,
     $                1, mpi_comm_col)
c #ifdef TIMER
c       de=mpi_wtime()
c       dcom=dcom+(de-ds)
c #endif
*-
! sm:= d*sm
      do m_0=1,m-1
!dir$ ivdep
         do i_0=m_0+1,m
            ss(i_0+(m_0-1)*nsm) = ss(i_0+(m_0-1)*nsm) * d(i+i_0-1)
         end do! i_0
      end do! m_0
*-
! ss:= d*ss
      i_4 = mod(i_3-i_2+1,4)+i_2
      do i_1=i_2,i_4-1                       ! 0
!dir$ ivdep
         do m_0=1,m
            ss((i_1+0-i_2)*m+m_0+nsm*nsm) =
     $      ss((i_1+0-i_2)*m+m_0+nsm*nsm) * d(i+m_0-1)
         end do! m_0
      end do! i_1
      do i_1=i_4,i_3,4                  ! 3
!dir$ ivdep
         do m_0=1,m
            ss((i_1+0-i_2)*m+m_0+nsm*nsm) =
     $      ss((i_1+0-i_2)*m+m_0+nsm*nsm) * d(i+m_0-1)
            ss((i_1+1-i_2)*m+m_0+nsm*nsm) =
     $      ss((i_1+1-i_2)*m+m_0+nsm*nsm) * d(i+m_0-1)
            ss((i_1+2-i_2)*m+m_0+nsm*nsm) =
     $      ss((i_1+2-i_2)*m+m_0+nsm*nsm) * d(i+m_0-1)
            ss((i_1+3-i_2)*m+m_0+nsm*nsm) =
     $      ss((i_1+3-i_2)*m+m_0+nsm*nsm) * d(i+m_0-1)
         end do! m_0
      end do! i_1
!$omp end master
*-
!$omp barrier
! v:= v*(i-sm)^{-1}
! z(j_2:j_3,1:n):= z + v(j_2:j_3,1:m)*ss(1:m,1:n)
      ds=mpi_wtime()
      call eigen_tbk_step2( z,
     $            nmz, v, nm, m,
     $            ss(1+nsm*nsm), ss(1),
     $            i_2,i_3,j_2,j_3 )

!$omp barrier
      de=mpi_wtime()

!$omp master
      dy = dy + (de-ds)
!$omp end master
*-

      return
      end subroutine  eigen_tbk_main2

