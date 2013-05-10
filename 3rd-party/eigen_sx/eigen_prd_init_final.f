      subroutine eigen_prd_init(a, nm, n,
     &              d_out, e_out, ne,
     &              u_t, v_t, nv)
!----
      use communication_sx, only : get_loop_start
     &                           , get_loop_end
     &                           , translate_g2l
     &                           , translate_l2g
!----
!$    use omp_lib
      implicit none

      integer, intent(in)    :: nm, n, ne, nv
      real(8), intent(inout) :: a(1:nm, *)
      real(8), intent(out)   :: d_out(*)
      real(8), intent(out)   :: e_out(1:ne, *)
      real(8), intent(in)    :: u_t(*), v_t(*)

      integer                ::  i_1, i_2, i_3, i_4
      integer                ::  j_1, j_2
      integer                ::  j

      integer                :: kx
      integer                :: thread_size, thread_rank

      real(8), pointer       :: u0_z(:),v0_z(:)
      real(8), pointer       :: u1_z(:),v1_z(:)
      integer                :: offset1, offset2
      integer                :: offset3, offset4
      common/tred_au_common/    u0_z, v0_z, u1_z, v1_z,
     &                          offset1, offset2, offset3, offset4

      include 'trd.h'
      include 'param.h'
      include 'cstab.h'


      d_out(1:n)    = zero
      e_out(1:n,1)  = zero
      e_out(1:n,2)  = zero

      if ( diag_0 > 0 ) then
         i_2 = diag_0
         i_3 = get_loop_end  (n, size_of_row,my_row)
         i_4 = size_of_col/n_common
         do i_1=i_2,i_3,i_4
!----       j   = translate_l2g(i_1, size_of_row,my_row)
!----       j_1 = translate_g2l(j, size_of_col,my_col)
            j   = (i_1-1)*size_of_row+my_row
            j_1 = (j-1)/size_of_col+1
            d_out(j) = a(j_1,i_1)
         end do! i_1
      end if

      i_2 = get_loop_start(1, size_of_row,my_row)
      i_3 = get_loop_end  (n, size_of_row,my_row)
      do i_1=i_2,i_3
         j   = translate_l2g(i_1, size_of_row,my_row)
         j_2 = get_loop_start(j+1, size_of_col,my_col)
         do j_1=j_2,nm
            a(j_1,i_1) = zero
         end do! j_1
         if ( j > n ) then
            do j_1=j_2,nm
               a(1:nm,i_1) = zero
            end do! j_1
         end if
      end do! i_1

      thread_rank = 0
      thread_size = 1
!$    thread_rank = omp_get_thread_num()
!$    thread_size = omp_get_num_threads()

      allocate(u0_z(nv*thread_size+n_columns),
     &         v0_z(nv*thread_size+n_columns))
      allocate(u1_z(nv*thread_size+n_columns),
     &         v1_z(nv*thread_size+n_columns))
      call cstab_adjust_base(u0_z(1), u_t(1), offset1)
      call cstab_adjust_base(v0_z(1), v_t(1), offset2)
      call cstab_adjust_base(u1_z(1), u_t(1), offset3)
      call cstab_adjust_base(v1_z(1), v_t(1), offset4)
      kx =  (l1_window/8)
!----&           +(l1_window)
!----&           +(l1_lsize/8)
     &           +(l1_lsize)
     &           +(l2_lsize/8)
      offset1 = offset1 + kx * 1
      offset2 = offset2 + kx * 2
      offset3 = offset3 + kx * 3
      offset4 = offset4 + kx * 4
      call cstab_round_offset(offset1)
      call cstab_round_offset(offset2)
      call cstab_round_offset(offset3)
      call cstab_round_offset(offset4)

      return
      end subroutine eigen_prd_init


      subroutine eigen_prd_final(a, nm, n, d_out, e_out, ne, u_t)
!----
      use communication_sx, only : reduce_dbl
     &                           , get_loop_node
     &                           , get_loop_end
     &                           , translate_g2l
     &                           , translate_l2g
!----
      implicit none

      integer, intent(in)    :: nm, ne, n
      real(8), intent(inout) :: a(1:nm, *)
      real(8), intent(out)   :: d_out(*)
      real(8), intent(out)   :: e_out(1:ne, *)
      real(8), intent(out)   :: u_t(*)

      integer                ::  owner_node_col
      integer                ::  owner_node_row
      integer                ::  i_1, i_2, i_3, i_4
      integer                ::  j_1
      integer                ::  i, j, l
      real(8)                ::  t

      real(8), pointer       :: u0_z(:),v0_z(:)
      real(8), pointer       :: u1_z(:),v1_z(:)
      integer                :: offset1, offset2
      integer                :: offset3, offset4
      common/tred_au_common/    u0_z, v0_z, u1_z, v1_z,
     &                          offset1, offset2, offset3, offset4

      include 'mpif.h'
      include 'trd.h'
      include 'param.h'


      if ( n >= 2 ) then
         do i = 2+mod(n,2), 2, -1
            l = i-1
            owner_node_row = get_loop_node (i, size_of_row,my_row)
            owner_node_col = get_loop_node (l, size_of_col,my_col)
            if ( owner_node_row == my_row
     &         .and. owner_node_col == my_col ) then
               i_1 = translate_g2l(i, size_of_row,my_row)
               j_1 = translate_g2l(l, size_of_col,my_col)
               e_out(i,1)  =  a(j_1,i_1)
               a(j_1,i_1)  =  zero
            end if
            l = i-2; if ( l < 1 ) cycle
            owner_node_row = get_loop_node (i, size_of_row,my_row)
            owner_node_col = get_loop_node (l, size_of_col,my_col)
            if ( owner_node_row == my_row
     &         .and. owner_node_col == my_col ) then
               i_1 = translate_g2l(i, size_of_row,my_row)
               j_1 = translate_g2l(l, size_of_col,my_col)
               e_out(i,2)  =  a(j_1,i_1)
               a(j_1,i_1)  =  zero
            end if
         end do
      end if

      if ( diag_0 > 0 ) then
         i_2 = diag_0
         i_3 = get_loop_end  (n, size_of_row,my_row)
         i_4 = size_of_col/n_common
         do i_1=i_2,i_3,i_4
            j   = translate_l2g(i_1, size_of_row,my_row)
            j_1 = translate_g2l(j, size_of_col,my_col)
            t          = d_out(j)
            d_out(j)   = a(j_1,i_1)
            a(j_1,i_1) = t
         end do! i_1
      end if

      call reduce_dbl(d_out(1),   u_t(1), n, 1, mpi_comm_world)
      call reduce_dbl(e_out(1,1), u_t(1), n, 1, mpi_comm_world)
      call reduce_dbl(e_out(1,2), u_t(1), n, 1, mpi_comm_world)

      deallocate(u0_z, v0_z)
      deallocate(u1_z, v1_z)


      return
      end subroutine eigen_prd_final

