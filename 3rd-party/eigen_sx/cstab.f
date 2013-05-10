      subroutine cstab_get_optdim(n_min, n_unroll,
     &                            delta_l1, delta_l2, n_opt)
      implicit none

      integer                :: n_min, n_unroll, delta_l1, delta_l2
      integer                :: n_opt, n_delta

      include 'cstab.h'

      integer                :: i0, k0

      n_opt=n_min
      do

         n_opt   = (n_opt-1)/l1_window+1
         n_opt   = (n_opt/2)*2+1
         n_opt   = n_opt*l1_window

         n_delta = 0

         do i0=1,int((n_unroll*1.2-1.0)/l1_way+1)

            k0=mod(i0*n_opt+l1_lsize/2,l1_lsize)-l1_lsize/2
            if(abs(k0)<=delta_l1/2)then
               n_delta=(delta_l1/2-k0-1)/i0+1
               goto 10000
            end if

         end do

         do i0=1,int((n_unroll*1.2-1.0)/l2_way+1)

            k0=mod(i0*n_opt+l2_lsize/2,l2_lsize)-l2_lsize/2
            if(abs(k0)<=delta_l2/2)then
               n_delta=(delta_l2/2-k0-1)/i0+1
               goto 10000
            end if

         end do

10000    continue

         if(n_delta==0)exit
         n_opt = n_opt + n_delta

      end do

      return
      end subroutine

      subroutine cstab_adjust_base(a,b,offset)
      implicit none

      real(8)                :: a(*), b(*)
      integer                :: offset

      include 'cstab.h'


      call get_delta(a(1),b(1),offset)
      offset=(offset/8)

      if(offset>0)then
         offset=mod(l2_lsize-mod(+offset,l2_lsize),l2_lsize)
      else
         offset=mod(l2_lsize+mod(-offset,l2_lsize),l2_lsize)
      endif

      return
      end subroutine
!----
      subroutine cstab_round_offset(offset)
      implicit none

      integer                :: offset

      include 'cstab.h'

      offset=mod(offset,l2_lsize)

      return
      end subroutine
