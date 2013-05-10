      subroutine eigen_pbk_main1(n, a, nma, z, nmz, e, d,
     &             v, nm, m, iblk, ss, tt)
!----
      use communication_sx, only : reduce_dbl
     &                           , datacast_dbl
     &                           , bcast_dbl
     &                           , get_loop_start
     &                           , get_loop_end
     &                           , get_loop_node
     &                           , translate_g2l
     &                           , translate_l2g
!----
      implicit double precision(a-h,o-z),integer(i-n)

      real(8) :: a(1:nma,*)
      real(8) :: z(1:nmz,*)
      real(8) :: e(1:n)
      real(8) :: d(1:n)

      real(8) :: v(1:nm,*)

      real(8) :: ss(*), tt(*)

      include 'mpif.h'
!$    include 'omp_lib.h'
      include 'trd.h'
!----
      integer :: nodes(0:nsm-1)

      real(8), pointer :: wk(:)

!----
!$omp master
      call mpi_barrier(mpi_comm_world,ierr)
      d1=mpi_wtime()

      lwk=((m-1)/size_of_row+1)*((n-1)/size_of_col+1)
      allocate(wk(lwk))
!$omp end master
!----
      i_start = get_loop_start(1+iblk, size_of_row,my_row)
      i_end   = get_loop_end  (n,      size_of_row,my_row)
!----
!$omp master
      d(1:n)=0.0
      do i0=i_start, i_end
!----    i = translate_l2g(i0, size_of_row,my_row)
         i = (i0-1)*size_of_row+my_row
         l = i-iblk
!----    owner_node_col = get_loop_node(l, size_of_col,my_col)
         owner_node_col = mod(l-1,size_of_col)+1
         if ( owner_node_col == my_col ) then
!----       j0 = translate_g2l(l, size_of_col,my_col)
            j0  = (l-1)/size_of_col+1
            d(i) = a(j0,i0)
         end if
      end do! i0
      call reduce_dbl( d, v, n, 1, mpi_comm_world )
      do i0=1+iblk,n
         if ( e(i0)*d(i0) == 0.0d+00 ) then
            s0 = 0.0
         else
            s0 = (1.0d0/d(i0))/e(i0)
         end if
         d(i0) = s0
      end do! i0
!$omp end master
!----
      i_2 = get_loop_start(1, size_of_row,my_row)
      i_3 = get_loop_end  (n, size_of_row,my_row)
!----
      nx = min(mod(n-(1+iblk)+1,m)+(1+iblk)-1,n)
!----
      do i=(1+iblk),nx
!----
         if ( e(i) == 0.0 ) cycle
!----
         j_2 = get_loop_start(1,   size_of_col,my_col)
         j_3 = get_loop_end  (i-iblk, size_of_col,my_col)
         i_4=mod(i_3-i_2+1,4)+i_2
!----
!$omp master
         nodes(0) = get_loop_node(i, size_of_row, my_row)
         if ( nodes(0) == my_row ) then
            i_1 = translate_g2l(i, size_of_row, my_row)
            do j_1=j_2,j_3
               v(j_1,1) = a(j_1,i_1)
            end do! j_1
         end if
         call bcast_dbl(v(j_2,1), j_3-j_2+1, nodes(0), mpi_comm_row)
!$omp end master
!----
!$omp barrier
!----
!$omp master
         if ( i_4 == i_2 + 1 ) then
            i_1 = i_2
            s0=0.0d+00
            do j_1=j_2,j_3
               s0=s0+v(j_1,1)*z(j_1,i_1+0)
            end do! j_1
            ss(i_1+0)=s0
         end if
         if ( i_4 == i_2 + 2 ) then
            i_1 = i_2
            s0=0.0d+00
            s1=0.0d+00
            do j_1=j_2,j_3
               s0=s0+v(j_1,1)*z(j_1,i_1+0)
               s1=s1+v(j_1,1)*z(j_1,i_1+1)
            end do! j_1
            ss(i_1+0)=s0
            ss(i_1+1)=s1
         end if
         if ( i_4 == i_2 + 3 ) then
            i_1 = i_2
            s0=0.0d+00
            s1=0.0d+00
            s2=0.0d+00
            do j_1=j_2,j_3
               s0=s0+v(j_1,1)*z(j_1,i_1+0)
               s1=s1+v(j_1,1)*z(j_1,i_1+1)
               s2=s2+v(j_1,1)*z(j_1,i_1+2)
            end do! j_1
            ss(i_1+0)=s0
            ss(i_1+1)=s1
            ss(i_1+2)=s2
         end if
!$omp end master
!$omp do
         do i_1=i_4,i_3,4
            s0=0.0d+00
            s1=0.0d+00
            s2=0.0d+00
            s3=0.0d+00
            do j_1=j_2,j_3
               s0=s0+v(j_1,1)*z(j_1,i_1+0)
               s1=s1+v(j_1,1)*z(j_1,i_1+1)
               s2=s2+v(j_1,1)*z(j_1,i_1+2)
               s3=s3+v(j_1,1)*z(j_1,i_1+3)
            end do! j_1
            ss(i_1+0)=s0
            ss(i_1+1)=s1
            ss(i_1+2)=s2
            ss(i_1+3)=s3
         end do! i_1
!$omp end do
!----
!$omp barrier
!----
!$omp master
         call reduce_dbl(ss(i_2),tt, i_3-i_2+1, 1, mpi_comm_col)

         s0 = d(i)
         do i_1=i_2,i_3
            ss(i_1) = ss(i_1) * s0
         end do! i_1
!$omp end master

!$omp barrier

         j_2 = get_loop_start(1,      size_of_col,my_col)
         j_3 = get_loop_end  (i-iblk, size_of_col,my_col)
         i_4=mod(i_3-i_2+1,4)+i_2

!$omp master
         if ( i_4 == i_2 + 1 ) then
            i_1 = i_2
            s0 = ss(i_1+0)
            do j_1=j_2,j_3
               z(j_1,i_1+0) = z(j_1,i_1+0) + s0 * v(j_1,1)
            end do! j_1
         end if
         if ( i_4 == i_2 + 2 ) then
            i_1 = i_2
            s0 = ss(i_1+0)
            s1 = ss(i_1+1)
            do j_1=j_2,j_3
               z(j_1,i_1+0) = z(j_1,i_1+0) + s0 * v(j_1,1)
               z(j_1,i_1+1) = z(j_1,i_1+1) + s1 * v(j_1,1)
            end do! j_1
         end if
         if ( i_4 == i_2 + 3 ) then
            i_1 = i_2
            s0 = ss(i_1+0)
            s1 = ss(i_1+1)
            s2 = ss(i_1+2)
            do j_1=j_2,j_3
               z(j_1,i_1+0) = z(j_1,i_1+0) + s0 * v(j_1,1)
               z(j_1,i_1+1) = z(j_1,i_1+1) + s1 * v(j_1,1)
               z(j_1,i_1+2) = z(j_1,i_1+2) + s2 * v(j_1,1)
            end do! j_1
         end if
!$omp end master
!$omp do
         do i_1=i_4,i_3,4
            s0 = ss(i_1+0)
            s1 = ss(i_1+1)
            s2 = ss(i_1+2)
            s3 = ss(i_1+3)
            do j_1=j_2,j_3
               z(j_1,i_1+0) = z(j_1,i_1+0) + s0 * v(j_1,1)
               z(j_1,i_1+1) = z(j_1,i_1+1) + s1 * v(j_1,1)
               z(j_1,i_1+2) = z(j_1,i_1+2) + s2 * v(j_1,1)
               z(j_1,i_1+3) = z(j_1,i_1+3) + s3 * v(j_1,1)
            end do! j_1
         end do! i_1
!$omp enddo
!----
!$omp barrier
!----
      end do !i
!----
      do i=nx+1, n, m
!----
!$omp master
         ds=mpi_wtime()

         if ( m > size_of_row ) then

            do j=0,m-1
               nodes(j) = get_loop_node(i+j, size_of_row, my_row)
            enddo

            do iy=1,size_of_row

               k0=0
               do j=0,m-1
                  if ( nodes(j) == iy ) then
                     i_1 = translate_g2l(i+j, size_of_row, iy)
                     j_2 = get_loop_start(1,    size_of_col,my_col)
                     j_3 = get_loop_end  (i+m-1-iblk,size_of_col,my_col)
                     if ( my_row == iy ) then
                        do j_1=j_2,j_3
                           wk(k0+j_1) = a(j_1,i_1)
                        end do! k
                     endif
                     k0=k0+(j_3-j_2+1)
                  endif
               enddo

               call bcast_dbl(wk, k0, iy, mpi_comm_row)

               k0=0
               do j=0,m-1
                  if ( nodes(j) == iy ) then
                     i_1 = translate_g2l(i+j, size_of_row, iy)
                     j_2 = get_loop_start(1,    size_of_col,my_col)
                     j_3 = get_loop_end  (i+m-1-iblk,size_of_col,my_col)
                     do j_1=j_2,j_3
                        v(j_1,j+1) = wk(k0+j_1)
                     end do! k
                     k0=k0+(j_3-j_2+1)
                     j_2 = get_loop_start(i+j,  size_of_col,my_col)
                     j_3 = get_loop_end  (i+m-1-iblk,size_of_col,my_col)
                     do j_1=j_2,j_3
                        v(j_1,j+1) = 0.0d+00
                     end do
                  endif
               enddo

            enddo

         else

            do j=0,m-1
               nodes(j) = get_loop_node(i+j, size_of_row, my_row)
               if ( nodes(j) == my_row ) then
                  i_1 = translate_g2l(i+j, size_of_row, my_row)
                  j_2 = get_loop_start(1,    size_of_col,my_col)
                  j_3 = get_loop_end  (i+m-1-iblk,size_of_col,my_col)
                  do j_1=j_2,j_3
                     v(j_1,j+1) = a(j_1,i_1)
                  end do! k
                  j_2 = get_loop_start(i+j,  size_of_col,my_col)
                  j_3 = get_loop_end  (i+m-1-iblk,size_of_col,my_col)
!cdir novector
                  do j_1=j_2,j_3
                     v(j_1,j+1) = 0.0d+00
                  end do
               end if
            end do
            j_2 = get_loop_start(1,     size_of_col,my_col)
            j_3 = get_loop_end  (i+m-1-iblk, size_of_col,my_col)
            do j=0,m-1
               call bcast_dbl(v(1,j+1), j_3-j_2+1,
     &                        nodes(j), mpi_comm_row)
            enddo

         endif
!$omp end master
!----
!$omp barrier
!----
         call eigen_pbk_main2(i_3, z, nmz,
     $           d, v, nm, m, iblk, i, ss, tt)
!----
!$omp barrier
!----
      end do !i
!----
!$omp master
      deallocate(wk)
      call mpi_barrier(mpi_comm_world,ierr)
!$omp end master
!----
      return
      end subroutine
