      subroutine eigen_s_init_wrapper(ndim, size_of_col_local,
     &size_of_row_local )

      use communication_s, only : eigen_init

      implicit none

      include 'mpif.h'
      include 'trd.h'

      integer ::  ndim
      integer ::  size_of_col_local
      integer ::  size_of_row_local

      call eigen_init(ndim)
      size_of_col_local = size_of_col
      size_of_row_local = size_of_row

      return
      end subroutine            ! eigen_s_init_wrapper


      subroutine eigen_s_free_wrapper(flag)
      use communication_s, only : eigen_free

      implicit none

      include 'mpif.h'
      include 'trd.h'

      integer :: flag

      call eigen_free(flag)

      return
      end subroutine            ! eigen_s_free_wrapper

