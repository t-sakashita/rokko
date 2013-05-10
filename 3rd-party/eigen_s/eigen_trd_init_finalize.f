      subroutine  eigen_trd_init(a, nm, n,
     &            d_out, e_out,
     &            u_t, v_t, nv)
!----
      use communication_s, only : get_loop_start
     &                          , get_loop_end
     &                          , translate_l2g
!----
      implicit none

      integer, intent(in)    ::  nm, n, nv
      real(8), intent(inout) ::  a(1:nm, *)
      real(8), intent(out)   ::  d_out(*)
      real(8), intent(out)   ::  e_out(*)
      real(8), intent(in)    ::  u_t(*), v_t(*)

      integer                ::  i_1, i_2, i_3, i_4
      integer                ::  j_1, j_2, j_3, j_4
      integer                ::  k_1, k_2, k_3
      integer                ::  j

      integer                ::  thread_size

!$    include 'omp_lib.h'
      include 'trd.h'
      include 'trd_au.h'
      include 'CSTAB.h'
      include 'param.h'


      d_out(1:n) = zero
      e_out(1:n) = zero

      if ( diag_0 > 0 ) then
         i_2 = diag_0
         i_3 = get_loop_end  (n, size_of_row, my_row)
         if ( i_2 <= i_3 ) then
            i_4 = size_of_col/n_common
            j_2 = diag_1
            j_3 = get_loop_end  (n, size_of_col, my_col)
            j_4 = size_of_row/n_common
            k_2 = 0
            k_3 = (i_3-i_2)/i_4
            do  k_1 = k_2, k_3
               i_1 = i_2 + k_1 * i_4
               j_1 = j_2 + k_1 * j_4
               j   = (i_1-1)*size_of_row+my_row
               d_out(j) = a(j_1, i_1)
            end do ! k_1
         end if
      end if

      i_2 = get_loop_start(1, size_of_row, my_row)
      i_3 = get_loop_end  (n, size_of_row, my_row)
      do  i_1 = i_2, i_3
         j   = translate_l2g(i_1, size_of_row, my_row)
         j_2 = get_loop_start    (j+1, size_of_col, my_col)
         if ( j <= n ) then
            a(j_2:nm, i_1) = zero
         else
            a(1:nm, i_1) = zero
         end if
      end do ! i_1

      thread_size = 1
!$    thread_size = omp_get_num_threads()

      if ( thread_size > 1 ) then

         allocate(u0_z(nv*thread_size+n_columns),
     &            v0_z(nv*thread_size+n_columns))
         call cstab_adjust_base(u0_z(1), u_t(1), offset1)
         call cstab_adjust_base(v0_z(1), v_t(1), offset2)

      end if


      return
      end subroutine  eigen_trd_init


      subroutine  eigen_trd_finalize(a, nm, n, d_out, e_out, u_t)
!----
      use communication_s, only : bcast_dbl
     &                          , reduce_dbl
     &                          , get_loop_end
     &                          , get_owner_node
     &                          , translate_g2l
!----
      implicit none

      integer, intent(in)    ::  nm, n
      real(8), intent(inout) ::  a(1:nm, *)
      real(8), intent(out)   ::  d_out(*)
      real(8), intent(out)   ::  e_out(*)
      real(8), intent(out)   ::  u_t(*)

      integer                ::  owner_node_col
      integer                ::  owner_node_row
      integer                ::  i_1, i_2, i_3, i_4
      integer                ::  j_1, j_2, j_3, j_4
      integer                ::  k_1, k_2, k_3
      integer                ::  i, j, k, l
      real(8)                ::  t

      integer                ::  thread_size

!$    include 'omp_lib.h'
      include 'mpif.h'
      include 'trd.h'
      include 'trd_au.h'
      include 'CSTAB.h'


      if ( n >= 2 ) then
         i = 2; l = i-1
         owner_node_row = get_owner_node (i, size_of_row, my_row)
         owner_node_col = get_owner_node (l, size_of_col, my_col)
         i_1 = translate_g2l(i, size_of_row, owner_node_row)
         j_1 = translate_g2l(l, size_of_col, owner_node_col)
         if ( owner_node_row == my_row
     &         .and. owner_node_col == my_col ) then
            e_out(i)    = -a(j_1, i_1)
            a(j_1, i_1) =2*a(j_1, i_1)
         end if
         j = size_of_col * (owner_node_row-1) + (owner_node_col-1) + 1
         call bcast_dbl( e_out(2), 1, j, mpi_comm_world )
      end if

      if ( diag_0 > 0 ) then
         i_2 = diag_0
         i_3 = get_loop_end  (n, size_of_row, my_row)
         if ( i_2 <= i_3 ) then
            i_4 = size_of_col/n_common
            j_2 = diag_1
            j_3 = get_loop_end  (n, size_of_col, my_col)
            j_4 = size_of_row/n_common
            k_2 = 0
            k_3 = (i_3-i_2)/i_4
            do  k_1 = k_2, k_3
               i_1 = i_2 + k_1 * i_4
               j_1 = j_2 + k_1 * j_4
               j   = (i_1-1)*size_of_row+my_row
               t           = d_out(j)
               d_out(j)    = a(j_1, i_1)
               a(j_1, i_1) = t
            end do ! k_1
         end if
      end if

      call reduce_dbl(d_out(1),   u_t(1), n, 1, mpi_comm_world)

      thread_size = 1
!$    thread_size = omp_get_num_threads()

      if ( thread_size > 1 ) then

         deallocate(u0_z,v0_z)

      end if

      return
      end subroutine eigen_trd_finalize

