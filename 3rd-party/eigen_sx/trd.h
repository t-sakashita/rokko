       integer                :: myrank,nprocs
       common /USEMPI/           myrank,nprocs

       integer                :: my_col, size_of_col, mpi_comm_col,
     &                           my_row, size_of_row, mpi_comm_row,
     &                           p0_(1024),q0_(1024), n_common,
     &                           diag_0, diag_1
       common /CYCL2D/           my_col, size_of_col, mpi_comm_col,
     &                           my_row, size_of_row, mpi_comm_row,
     &                           p0_      ,q0_      , n_common,
     &                           diag_0, diag_1

!      external               :: get_loop_range
!      integer, external      :: get_loop_start, get_loop_end
!      integer, external      :: get_loop_node
!      integer, external      :: translate_g2l
!      integer, external      :: translate_l2g

       integer, parameter :: nsx = 480
       integer, parameter :: nsm = 256 ! 64 +32

       integer, parameter :: mband = 2  ! band width
