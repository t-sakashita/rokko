!
! utility for measuring elapse time
!

subroutine clock(t)
  implicit none
  real*8 t
  integer*4 ti, t_rate, t_max, diff
  call system_clock(ti, t_rate, t_max)
  t = ti / dble(t_rate)
end subroutine
