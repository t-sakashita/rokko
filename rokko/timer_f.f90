!/*****************************************************************************
!*
!* Rokko: Integrated Interface for libraries of eigenvalue decomposition
!*
!* Copyright (C) 2012-2013 by Tatsuya Sakashita <t-sakashita@issp.u-tokyo.ac.jp>,
!*                            Synge Todo <wistaria@comp-phys.org>,
!*                            Tsuyoshi Okubo <t-okubo@issp.u-tokyo.ac.jp>
!*
!* Distributed under the Boost Software License, Version 1.0. (See accompanying
!* file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
!*
!*****************************************************************************/

! still writing, not yet available

module rokko_timer
  use iso_c_binding
  implicit none

!timer routines
  subroutine set_timer(timer_)
    implicit none
    type(timer),intent(out)::timer_

    call initialize_timer(timer_%ptr_timer)
  end subroutine set_timer

  subroutine del_timer(timer_)
    implicit none
    type(timer),intent(out)::timer_

    call delete_timer(timer_%ptr_timer)
  end subroutine del_timer

  subroutine start_timer(timer_, id)
    type(timer),intent(inout)::timer_
    integer,intent(in):: id
    call timer_start(timer_%ptr_timer, id)
  end subroutine start_timer

  subroutine end_timer(timer_, id)
    type(timer),intent(inout)::timer_
    integer,intent(in):: id
    call timer_stop(timer_%ptr_timer, id)
  end subroutine end_timer

  subroutine registrate_timer(timer_, id, label)
    implicit none
    type(timer),intent(inout)::timer_
    integer,intent(in)::id
    character(*),intent(in)::label
    call timer_registrate(timer_%ptr_timer, id, label)
  end subroutine registrate_timer

  real(kind(1d0)) function get_count_timer(timer_, id)
    type(timer),intent(inout)::timer_
    integer,intent(in):: id
    real(kind(1d0))::timer_get_count
    external timer_get_count

    get_count_timer = timer_get_count(timer_%ptr_timer, id)
  end function get_count_timer

  real(kind(1d0)) function get_average_timer(timer_, id)
    type(timer),intent(inout)::timer_
    integer,intent(in):: id
    real(kind(1d0))::timer_get_average
    external timer_get_average

    get_average_timer = timer_get_average(timer_%ptr_timer, id)
  end function get_average_timer

end module rokko_timer
