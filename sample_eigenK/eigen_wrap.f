      subroutine eigen_init_wrapper(ndim, size_of_col_local,
     &size_of_row_local )

      use communication_s, only : eigen_init

      implicit none

      include 'mpif.h'
!      include 'trd.h'

      integer                :: my_col, size_of_col, mpi_comm_col,
     &                          my_row, size_of_row, mpi_comm_row,
     &                          p0_(1024),q0_(1024), n_common,
     &                          diag_0, diag_1
      common /CYCL2D/           my_col, size_of_col, mpi_comm_col,
     &                          my_row, size_of_row, mpi_comm_row,
     &                          p0_      ,q0_      , n_common,
     &                          diag_0, diag_1

      integer ::  ndim
      integer ::  size_of_col_local
      integer ::  size_of_row_local

      call eigen_init(ndim)
      size_of_col_local = size_of_col
      size_of_row_local = size_of_row

      return
      end subroutine            ! eigen_init_wrapper    


      subroutine eigen_free_wrapper(flag)
      use communication_s, only : eigen_free

      implicit none

      include 'mpif.h'
!      include 'trd.h'

      integer :: flag

      call eigen_free(flag)

      return
      end subroutine            ! eigen_free_wrapper
