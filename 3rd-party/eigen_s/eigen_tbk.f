
      subroutine eigen_tbk(n, a, nma0, z, nmz0, e, m0)
!----
      use communication_s, only : eigen_init
     &                          , eigen_free
!----
      implicit double precision(a-h,o-z),integer(i-n)
      real(8) :: a(*)
      real(8) :: z(*)
      real(8) :: e(*)
      real(8) :: hs0, hs1

      real(8) , pointer :: d(:)
      real(8) , pointer :: v(:)
      real(8) , pointer :: tt(:), ss(:)

      include 'mpif.h'
!$    include 'omp_lib.h'
      include 'trd.h'
      include 'CSTAB.h'
*
      nma  = nma0
      nmz  = nmz0
      m    = min(nsm,m0)
      if ( m < 1 ) m = 1
*
      call eigen_init(2)
*
      na   = (n-1)/size_of_row+1
      na   = na  +mod(na  -1,2)
      call cstab_get_optdim(nma,9,16*4,16*6,nm)
*
      allocate(
     &         d(1:n),
     &         v(1:max(nm*m,n)+n_columns),
     &         ss(1:na*m+nsm*nsm+n_columns),
     &         tt(1:na*m+nsm*nsm+n_columns),
     &         stat=i_stat)
      if(i_stat/=0)then
         if(myrank==1)print*,"memory allocation error."
         call mpi_abort(mpi_comm_world,1,ierr)
      endif
*
      call cstab_adjust_base(v(1), z(1), i_v)
      call cstab_adjust_base(ss(1), z(1), i_s)
      call cstab_adjust_base(tt(1), z(1), i_t)
      kx = (l1_window/4)
     &    +(l1_lsize)
     &    +(l2_lsize/4)
      i_v = i_v + kx*1
      i_s = i_s + kx*2
      i_t = i_t + kx*3
      call cstab_round_offset(i_v)
      call cstab_round_offset(i_s)
      call cstab_round_offset(i_t)

       hs0=mpi_wtime()
*
!$omp parallel
      call eigen_tbk_main1( n,
     $     a(1), nma,
     $     z(1), nmz,
     $     e(1), d(1), v(1+i_v), nm, m, ss(1+i_s), tt(1+i_t) )
!$omp end parallel
*
       hs1=mpi_wtime()

      deallocate( d )
      deallocate( v )
      deallocate( ss )
      deallocate( tt )
*
      call eigen_free(3)
*
#ifdef TIMER
      if ( myrank == 1 ) then
         print *,"Exectime of \"eigen_tbk\" routine =",hs1-hs0,"(sec)"
      endif
#endif
*
      return
      end subroutine
*
*
*
      subroutine eigen_tbk_main1(n, a,nma, z,nmz, e,d,
     &             v,nm, m, ss,tt )
!----
      use communication_s, only : bcast_dbl
     &                          , reduce_dbl
     &                          , get_loop_start
     &                          , get_loop_end
     &                          , translate_l2g
     &                          , translate_g2l
     &                          , get_owner_node
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
*
      integer :: nodes(0:nsm-1)

      real(8), pointer :: wk(:)

*
!$omp master
#ifdef TIMER
      call mpi_barrier(mpi_comm_world,ierr)
      d1=mpi_wtime()
      dcom=0.0d0+00
#endif
      dx = 0.0; dy = 0.0

      lwk=((m-1)/size_of_row+1)*((n-1)/size_of_col+1)
      allocate(wk(lwk))
!$omp end master
*
      i_2 = get_loop_start(2, size_of_row,my_row)
      i_3 = get_loop_end  (n, size_of_row,my_row)
*
      nx = min(mod(n-2+1,m)+2-1,n)
*
!$omp master
      d(1:n)=0.0
      do i_1=i_2,i_3
!----    i = translate_l2g(i_1, size_of_row,my_row)
         i = (i_1-1)*size_of_row+my_row
         l = i-1
!----    owner_node_col = get_owner_node(l, size_of_col,my_col)
         owner_node_col = mod(l-1,size_of_col)+1
         if ( owner_node_col == my_col ) then
!----       j_1 = translate_g2l(l, size_of_col,my_col)
            j_1  = (l-1)/size_of_col+1
            d(i) = a(j_1,i_1)
         end if
      end do! i_1
      call reduce_dbl( d, v, n, 1, mpi_comm_world )
      do i=2,n
         if ( e(i) == 0.0d+00 ) then
            s0 = 0.0
         else
            s0 = (1.0d0/d(i))/e(i)
         end if
         d(i) = s0
      end do! i
!$omp end master
*
      i_2 = get_loop_start(1, size_of_row,my_row)
      i_3 = get_loop_end  (n, size_of_row,my_row)

      do i=2,nx
*
!$omp barrier
*
         if ( e(i) == 0.0 ) cycle
*
         j_2 = get_loop_start(1,   size_of_col,my_col)
         j_3 = get_loop_end  (i-1, size_of_col,my_col)
         i_4=mod(i_3-i_2+1,4)+i_2
*
!$omp master
#ifdef TIMER
         call mpi_barrier(mpi_comm_row,ierr)
         ds=mpi_wtime()
#endif
         nodes(0) = get_owner_node(i, size_of_row, my_row)
         if ( nodes(0) == my_row ) then
            i_1 = translate_g2l(i, size_of_row, my_row)
            do j_1=j_2,j_3
               v(j_1,1) = a(j_1,i_1)
            end do! j_1
         end if
         call bcast_dbl(v(j_2,1), j_3-j_2+1, nodes(0), mpi_comm_row)
         de=mpi_wtime()
         dcom=dcom+(de-ds)
!$omp end master
*
!$omp barrier
*
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
         do i_1=i_4,i_3,4 !             3
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

!$omp barrier

!$omp master
#ifdef TIMER
         call mpi_barrier(mpi_comm_col,ierr)
         ds=mpi_wtime()
#endif
         call reduce_dbl(ss(i_2),tt, i_3-i_2+1, 1, mpi_comm_col)

         de=mpi_wtime()
         dcom=dcom+(de-ds)

         s0 = d(i)
         do i_1=i_2,i_3
            ss(i_1) = ss(i_1) * s0
         end do! i_1
!$omp end master

!$omp barrier

         j_2 = get_loop_start(1,   size_of_col,my_col)
         j_3 = get_loop_end  (i-1, size_of_col,my_col)
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
         do i_1=i_4,i_3,4 !             3
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
*
!$omp barrier
*
      end do
*
      do i=nx+1, n, m

*
!$omp barrier
*
!$omp master
#ifdef TIMER
         call mpi_barrier(mpi_comm_row,ierr)
         ds=mpi_wtime()
#endif
         if ( m > size_of_row ) then

            do j=0,m-1
               nodes(j) = get_owner_node(i+j, size_of_row, my_row)
            enddo

            do iy=1,size_of_row

               k0=0
               do j=0,m-1
                  if ( nodes(j) == iy ) then
                     i_1 = translate_g2l(i+j, size_of_row, iy)
                     j_2 = get_loop_start(1,    size_of_col,my_col)
                     j_3 = get_loop_end  (i+m-2,size_of_col,my_col)
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
                     j_3 = get_loop_end  (i+m-2,size_of_col,my_col)
                     do j_1=j_2,j_3
                        v(j_1,j+1) = wk(k0+j_1)
                     end do! k
                     k0=k0+(j_3-j_2+1)
                     j_2 = get_loop_start(i+j,  size_of_col,my_col)
                     j_3 = get_loop_end  (i+m-2,size_of_col,my_col)
                     do j_1=j_2,j_3
                        v(j_1,j+1) = 0.0d+00
                     end do
                  endif
               enddo

            enddo

         else

            do j=0,m-1
               nodes(j) = get_owner_node(i+j, size_of_row, my_row)
               if ( nodes(j) == my_row ) then
                  i_1 = translate_g2l(i+j, size_of_row, my_row)
                  j_2 = get_loop_start(1,    size_of_col,my_col)
                  j_3 = get_loop_end  (i+m-2,size_of_col,my_col)
                  do j_1=j_2,j_3
                     v(j_1,j+1) = a(j_1,i_1)
                  end do! k
                  j_2 = get_loop_start(i+j,  size_of_col,my_col)
                  j_3 = get_loop_end  (i+m-2,size_of_col,my_col)
!cdir novector
                  do j_1=j_2,j_3
                     v(j_1,j+1) = 0.0d+00
                  end do
               end if
            end do
            j_2 = get_loop_start(1,     size_of_col,my_col)
            j_3 = get_loop_end  (i+m-2, size_of_col,my_col)
            do j=0,m-1
               call bcast_dbl(v(1,j+1), j_3-j_2+1,
     &                        nodes(j), mpi_comm_row)
            enddo

         endif

         de=mpi_wtime()
         dcom=dcom+(de-ds)
!$omp end master
*
!$omp barrier
*
         call eigen_tbk_main2(i_3, z, nmz,
     $        d, v, nm, m, i, ss, tt, dcom, dx, dy)
*
!$omp barrier
*
      end do
*
!$omp master
      deallocate(wk)
      call mpi_barrier(mpi_comm_world,ierr)
#ifdef DETAIL
      d2=mpi_wtime()
      if ( 1 == myrank ) then
         print*," "
         print*,"detail of exectime in eigen_tbk "
         print*,"   time of eigen_tbk_main1=",(d2-d1),"(sec)"
         print*,"   communication in eigen_tbk_main1 =",dcom,"(sec)"
!----    print*,"   ",(2d0*n*n*n)/(d2-d1)*1d-9,"gflops"
!----    print*,"   ",(1d0*n*n*n)/(dx)*1d-9,"gflops"
!----    print*,"   ",(1d0*n*n*n)/(dy)*1d-9,"gflops"
      end if
#endif
!$omp end master
*
*
      return
      end subroutine

