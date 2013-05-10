! power 4
!     integer, parameter     :: l1_size   = 32*1024
!     integer, parameter     :: l1_way    = 2
!     integer, parameter     :: l2_size   = 512*1024
!     integer, parameter     :: l2_way    = 8
 
! pentium 4(northwood or prior)
!     integer, parameter     :: l1_size   = 8*1024
!     integer, parameter     :: l1_way    = 4
!     integer, parameter     :: l2_size   = 512*1024
!     integer, parameter     :: l2_way    = 8
 
! pentium 4(prescott)
!     integer, parameter     :: l1_size   = 16*1024
!     integer, parameter     :: l1_way    = 8
!     integer, parameter     :: l2_size   = 1024*1024
!     integer, parameter     :: l2_way    = 8
 
! cerelon d(prescott)
!     integer, parameter     :: l1_size   = 16*1024
!     integer, parameter     :: l1_way    = 8
!     integer, parameter     :: l2_size   = 256*1024
!     integer, parameter     :: l2_way    = 4

! itanium 2
!     integer, parameter     :: l1_size   = 16*1024
!     integer, parameter     :: l1_way    = 8
!     integer, parameter     :: l2_size   = 256*1024
!     integer, parameter     :: l2_way    = 8

! core2
      integer, parameter     :: l1_size   = 32*1024
      integer, parameter     :: l1_way    = 8
      integer, parameter     :: l2_size   = 2*1024*1024
      integer, parameter     :: l2_way    = 16

      integer, parameter     :: l1_lsize  = (l1_size/l1_way)/8
!      integer, parameter     :: l1_window = 8
      integer, parameter     :: l1_window = 16
      integer, parameter     :: l2_lsize  = (l2_size/l2_way)/8
      integer, parameter     :: l2_window = 64

      integer, parameter     :: n_columns = l2_lsize
      integer, parameter     :: prefetch_size = 512/8

! fr ia32-linux
      integer, parameter     :: page_size  = 4096
      integer, parameter     :: page_lsize = page_size/8

