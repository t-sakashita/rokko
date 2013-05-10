      subroutine eigen_tbk_step1(z, nmz,
     $           v, nm, m, ss, sm, i_2,i_3,j_2,j_3)
      implicit double precision(a-h,o-z),integer(i-n)
*
      integer :: nmz, nm, m
      real(8) :: z(nmz,*),v(nm,*)
      real(8) :: ss(m,*)
*
      include 'mpif.h'
!$    include 'omp_lib.h'
      include 'trd.h'
*
      real(8) :: sm(nsm,nsm)
      integer :: thread_rank, thread_size
*
      thread_size = 1
      thread_rank = 0
!$    thread_size = omp_get_num_threads()
!$    thread_rank = omp_get_thread_num()
*
*-
      do j_0=j_2,j_3,nsx; j_4=min(j_0+nsx-1,j_3)

         do m_0=1+thread_rank,m-1,thread_size

            i_0 = m_0 + 1
            i_4 = mod(m-m_0-1+1,4)+i_0
            if ( i_0 == i_4 - 1 ) then
               s0 = sm(i_0+0,m_0)
!dir$ ivdep
!dir$ vector aligned
               do j_1=j_0,j_4
                  t0 = v(j_1,m_0)
                  v0 = v(j_1,i_0+0)
                  s0 = s0+v0*t0
               end do! k
               sm(i_0+0,m_0) = s0
            end if
            if ( i_0 == i_4 - 2 ) then
               s0 = sm(i_0+0,m_0)
               s1 = sm(i_0+1,m_0)
!dir$ ivdep
!dir$ vector aligned
               do j_1=j_0,j_4
                  t0 = v(j_1,m_0)
                  v0 = v(j_1,i_0+0)
                  v1 = v(j_1,i_0+1)
                  s0 = s0+v0*t0
                  s1 = s1+v1*t0
               end do! k
               sm(i_0+0,m_0) = s0
               sm(i_0+1,m_0) = s1
            end if
            if ( i_0 == i_4 - 3 ) then
               s0 = sm(i_0+0,m_0)
               s1 = sm(i_0+1,m_0)
               s2 = sm(i_0+2,m_0)
!dir$ ivdep
!dir$ vector aligned
               do j_1=j_0,j_4
                  t0 = v(j_1,m_0)
                  v0 = v(j_1,i_0+0)
                  v1 = v(j_1,i_0+1)
                  v2 = v(j_1,i_0+2)
                  s0 = s0+v0*t0
                  s1 = s1+v1*t0
                  s2 = s2+v2*t0
               end do! k
               sm(i_0+0,m_0) = s0
               sm(i_0+1,m_0) = s1
               sm(i_0+2,m_0) = s2
            end if
            do i_0=i_4,m,4              ! 3
               s0 = sm(i_0+0,m_0)
               s1 = sm(i_0+1,m_0)
               s2 = sm(i_0+2,m_0)
               s3 = sm(i_0+3,m_0)
!dir$ ivdep
!dir$ vector aligned
               do j_1=j_0,j_4
                  t0 = v(j_1,m_0)
                  v0 = v(j_1,i_0+0)
                  s0 = s0+v0*t0
                  v1 = v(j_1,i_0+1)
                  s1 = s1+v1*t0
                  v2 = v(j_1,i_0+2)
                  s2 = s2+v2*t0
                  v3 = v(j_1,i_0+3)
                  s3 = s3+v3*t0
               end do! k
               sm(i_0+0,m_0) = s0
               sm(i_0+1,m_0) = s1
               sm(i_0+2,m_0) = s2
               sm(i_0+3,m_0) = s3
            end do! i_0

         end do! m_0

      end do! j_0
*-
      j_5 = j_3 - j_2 + 1

      do ii_2 = i_2, i_3, (1024*thread_size)

         ii_3 = min(ii_2+(1024*thread_size)-1,i_3)
         ii_4 = ii_2-i_2+1

         i_0 = (ii_3-ii_2) / thread_size + 1
         i_4 = i_0 * thread_rank
         i_5 = min(i_0, ii_3-(ii_2+i_4)+1)

         if ( j_5 > 0 .and. i_5 > 0 ) then
            call dgemm('t','n',
     &                 m, i_5, j_5,
     &                 1.0d+00, v (j_2    ,1       ), nm,
     &                          z (j_2    ,ii_2+i_4), nmz,
     &                 1.0d+00, ss(1      ,ii_4+i_4), m)
         endif
      end do

*-
      return
      end subroutine  eigen_tbk_step1
*-
      subroutine eigen_tbk_step2( z,
     $           nmz,
     $           v, nm, m, ss, sm, i_2,i_3, j_2,j_3 )
      implicit double precision(a-h,o-z),integer(i-n)
*
      integer :: nmz, nm, m
      real(8) :: z(nmz,*)
      real(8) :: v(nm,*)
      real(8) :: ss(m,*)
*
      include 'mpif.h'
!$    include 'omp_lib.h'
      include 'trd.h'

*
      real(8) :: sm(nsm,*)
      integer :: thread_rank, thread_size

      thread_size = 1
      thread_rank = 0
!$    thread_size = omp_get_num_threads()
!$    thread_rank = omp_get_thread_num()
*
*-
      j_5 = (j_3-j_2) / thread_size + 1
      j_5 = ((j_5-1)/2+1)*2
      j_4 = j_5 * thread_rank
      j_5 = min(j_5, j_3-(j_2+j_4)+1)
      j_6 = j_2+j_4
      j_7 = j_6+j_5-1
         
      do j_0=j_6,j_7,nsx; j_8=min(j_0+nsx-1,j_7)

         do m_0=m,1,-2

            i_0 = m
            i_4 = mod(m-m_0, 3)
            if ( i_4 == 1 ) then
!dir$ ivdep
!dir$ vector aligned
               do j_1=j_0,j_8
                  v0 = v(j_1,i_0-0)
                  u0 =
     $                 + sm(i_0-0,m_0  ) * v0
                  u1 =
     $                 + sm(i_0-0,m_0-1) * v0
                  v(j_1,m_0  ) = v(j_1,m_0  ) + u0
                  v(j_1,m_0-1) = v(j_1,m_0-1) + u1
               end do! j_1
            end if
            if ( i_4 == 2 ) then
!dir$ ivdep
!dir$ vector aligned
               do j_1=j_0,j_8
                  v0 = v(j_1,i_0-0)
                  v1 = v(j_1,i_0-1)
                  u0 =
     $                + sm(i_0-0,m_0  ) * v0
     $                + sm(i_0-1,m_0  ) * v1
                  u1 =
     $                + sm(i_0-0,m_0-1) * v0
     $                + sm(i_0-1,m_0-1) * v1
                  v(j_1,m_0  ) = v(j_1,m_0  ) + u0
                  v(j_1,m_0-1) = v(j_1,m_0-1) + u1
               end do! j_1
            end if

            do i_0=m-i_4,m_0+1,-3
!dir$ ivdep
!dir$ vector aligned
               do j_1=j_0,j_8
                  v0 = v(j_1,i_0-0)
                  v1 = v(j_1,i_0-1)
                  v2 = v(j_1,i_0-2)
                  u0 =
     $                + sm(i_0-0,m_0  ) * v0
     $                + sm(i_0-1,m_0  ) * v1
     $                + sm(i_0-2,m_0  ) * v2
                  u1 =
     $                + sm(i_0-0,m_0-1) * v0
     $                + sm(i_0-1,m_0-1) * v1
     $                + sm(i_0-2,m_0-1) * v2
                  v(j_1,m_0  ) = v(j_1,m_0  ) + u0
                  v(j_1,m_0-1) = v(j_1,m_0-1) + u1
               end do! j_1
            end do! i_0

            i_0=m_0
!dir$ ivdep
!dir$ vector aligned
            do j_1=j_0,j_8
               v0 = v(j_1,i_0-0)
               u1 =
     $             + sm(i_0-0,m_0-1) * v0
               v(j_1,m_0-1) = v(j_1,m_0-1) + u1
            end do! j_1

         end do! m_0

      end do! j_0

!$omp barrier
*-
      j_5 = j_3 - j_2 + 1

      i_5 = (i_3-i_2) / thread_size + 1
      i_4 = i_5 * thread_rank
      i_5 = min(i_5, i_3-(i_2+i_4)+1)

      if ( j_5 > 0 .and. i_5 > 0 ) then
         call dgemm('n','n',
     &               j_5, i_5, m,
     &               1.0d+00, v (j_2    ,1      ), nm,
     &                        ss(1      ,1  +i_4), m,
     &               1.0d+00, z (j_2    ,i_2+i_4), nmz)
      end if
*-


      return
      end subroutine  eigen_tbk_step2

