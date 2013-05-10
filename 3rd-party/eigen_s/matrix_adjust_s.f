       subroutine matrix_adjust_s(n, a, b, nm)
!----
      use communication_s, only : get_loop_start
     &                          , get_loop_end
     &                          , translate_l2g
     &                          , eigen_init
     &                          , eigen_free
!----
       implicit double precision(a-h,o-z),integer(i-n)
       real(8)                :: a(n,n), b(1:nm,*)

       include 'mpif.h'
       include 'trd.h'

       call eigen_init(2)

       j_2=get_loop_start(1, size_of_col, my_col)
       j_3=get_loop_end  (n, size_of_col, my_col)
       i_2=get_loop_start(1, size_of_row, my_row)
       i_3=get_loop_end  (n, size_of_row, my_row)

       do i_1=i_2,i_3
          i = translate_l2g(i_1, size_of_row, my_row)
          do j_1=j_2,j_3
             j = translate_l2g(j_1, size_of_col, my_col)
             b(j_1,i_1)=a(j,i)
          enddo
       enddo

       call eigen_free(0)

       return
       end

