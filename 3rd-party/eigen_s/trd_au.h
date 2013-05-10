!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  u0_z           : temporary array for thread
!  v0_z           : temporary array for thread
!  u1_z           : temporary array for thread
!  v1_z           : temporary array for thread
!  offset1        : store the data alignment.
!  offset2        : store the data alignment.
!  offset3        : store the data alignment.
!  offset4        : store the data alignment.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      real(8), pointer       ::  u0_z(:),v0_z(:)
      real(8), pointer       ::  u1_z(:),v1_z(:)
      integer                ::  offset1, offset2
      integer                ::  offset3, offset4
      common  /tred_au_common/   u0_z, v0_z, u1_z, v1_z,
     &                           offset1, offset2, offset3, offset4

