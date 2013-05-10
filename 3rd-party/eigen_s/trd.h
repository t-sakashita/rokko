!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  nprocs          : the number of all processes
!  myrank          : number of the process
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      integer                :: nprocs,myrank
      common /USEMPI/           nprocs,myrank

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  my_col          : number of the column
!  size_of_col     : size of column direction
!  mpi_comm_col    : comunicator for column
!  my_row          : number of the row
!  size_of_row     : size of row direction
!  mpi_comm_row    : comunicator for row
!  p0_             : data for diagonal elements calculation
!  q0_             : data for diagonal elements calculation
!  n_common        : data for diagonal elements calculation
!  diag_0          : data for diagonal elements calculation
!  diag_1          : data for diagonal elements calculation
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      integer                :: my_col, size_of_col, mpi_comm_col,
     &                          my_row, size_of_row, mpi_comm_row,
     &                          p0_(1024),q0_(1024), n_common,
     &                          diag_0, diag_1
      common /CYCL2D/           my_col, size_of_col, mpi_comm_col,
     &                          my_row, size_of_row, mpi_comm_row,
     &                          p0_      ,q0_      , n_common,
     &                          diag_0, diag_1

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  get_loop_start  : calculate a local loop start position from the 
!                    global loop start position.
!  get_loop_end    : calculate a local loop end position from the 
!                    global loop end position.
!  get_owner_node  : calculate number of the specified rank.
!  translate_g2l   : calculate a local array index from a global array
!                    index.
!  translate_l2g   : calculate a global array index from a local array
!                    index.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     integer, external      :: get_loop_start, get_loop_end
!     integer, external      :: get_owner_node
!     integer, external      :: translate_g2l
!     integer, external      :: translate_l2g

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  nsx             : vector length of middle data.
!  nsm             : maximum number of the block.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      integer, parameter :: nsx = 480
      integer, parameter :: nsm = 256 

