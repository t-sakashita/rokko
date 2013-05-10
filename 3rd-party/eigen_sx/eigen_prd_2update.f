c======================================================================c
      subroutine eigen_prd_2update(
     &      ar, nma,
     &      ur, uyr, vr, vyr, nmv,
     &      m_size, i_base)
c======================================================================c
!----
      use communication_sx, only : get_loop_start
     &                           , get_loop_end
     &                           , translate_l2g
!----
!$    use omp_lib
      implicit none
!----
      real(8), intent(inout) :: ar (nma, *)
      real(8), intent(in)    :: ur (nmv, *)
      real(8), intent(in)    :: uyr(nmv, *)
      real(8), intent(in)    :: vr (nmv, *)
      real(8), intent(in)    :: vyr(nmv, *)
      integer, intent(in)    :: nma, nmv, m_size, i_base
!----
      include 'param.h'
!----
      include 'trd.h'

      integer ::  k0, k1, k2, m, n

      integer ::  i_1,i_2,i_3,i_4
      integer ::  j_1,j_2,j_3,j_4

      integer ::  ii_2,ii_3
      integer ::  jj_2,jj_3

      integer ::  blk_size1, blk_size2
      integer ::  ii_step
!----
      integer ::  thread_rank, thread_size
!----
      integer, parameter :: blas_chunk = 7*8
      integer, parameter :: blas_nvect = 256*3
!----
      intrinsic :: min, max
      external  :: dgemm
!----
c======================================================================c
!----
      if ( i_base <= 0 ) return
!----
c======================================================================c
!----
      thread_rank = 0
      thread_size = 1
!$    thread_rank = omp_get_thread_num()
!$    thread_size = omp_get_num_threads()
!----
      m  = m_size                               ! preserve the argument
                                                ! variable 
      n  = get_loop_end(i_base, size_of_row, my_row) ! local matrix size
      k0 = translate_l2g(n, size_of_row, my_row)  ! translation to the 
                                                ! global index

      jj_2 = 1                                  ! beggining of loop
      jj_3 = get_loop_end  (k0, size_of_col, my_col)   ! end of loop

      ii_step = 0
      do j_1 = jj_2, jj_3, blas_nvect
         j_4 = min(j_1+blas_nvect-1, jj_3)      ! [j_1:j_4] available 
                                                ! on this iteration 

         k1 = translate_l2g(j_1, size_of_col, my_col) ! translation to 
                                                      ! the global index
         ii_2 = get_loop_start(k1, size_of_row, my_row) ! beggining of
                                                        ! loop
         ii_2 = max(1, ii_2)                          ! ** should be 
                                                      ! .ge. 1
         ii_3 = n                                     ! end of loop

         do i_1 = ii_2, ii_3, blas_chunk
            i_4 = min(i_1+blas_chunk-1, ii_3)   ! [i_1:i_4] available 
                                                ! on this iteration

            k2  = translate_l2g(i_4, size_of_row, my_row) ! translation 
                                                   ! to the global index
            j_3 = get_loop_end  (k2, size_of_col, my_col) ! end of loop
            j_3 = min(j_4, j_3)

            i_2 = i_1; i_3 = i_4
            j_2 = j_1

            blk_size1 = j_3-j_2+1
            blk_size2 = i_3-i_2+1

            if ( blk_size1 > 0 .and. blk_size2 > 0 ) then

               if ( mod(ii_step, thread_size) == thread_rank ) then

                  call dgemm('n','t',
     &                 blk_size1, blk_size2, m,
     &                 m_one, ur (j_1, 1), nmv,
     &                        vyr(i_1, 1), nmv,
     &                 one,   ar (j_1, i_1), nma)
                  call dgemm('n','t',
     &                 blk_size1, blk_size2, m,
     &                 m_one, vr (j_1, 1), nmv,
     &                        uyr(i_1, 1), nmv,
     &                 one,   ar (j_1, i_1), nma)

               end if

               ii_step = ii_step+1

            end if

         end do ! i_1

      end do ! j_1
!----
      return
      end subroutine eigen_prd_2update
c======================================================================c

