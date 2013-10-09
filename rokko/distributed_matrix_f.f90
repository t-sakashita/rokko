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

module rokko_distributed_matrix
  use iso_c_binding
  implicit none

  type::rokko_distributed_matrix
     type(c_ptr) ptr_distributed_matrix
  end type rokko_distributed_matrix

  interface
     subroutine generate_distributed_matrix_function_col_major_fortran(mat, f) bind(c)
       use iso_c_binding
       type(c_ptr), value :: mat
       type(c_funptr), value :: f
     end subroutine generate_distributed_matrix_function_col_major_fortran

     subroutine generate_distributed_matrix_function_row_major_fortran(mat, f) bind(c)
       use iso_c_binding
       type(c_ptr), value :: mat
       type(c_funptr), value :: f
     end subroutine generate_distributed_matrix_function_row_major_fortran

     subroutine generate_distributed_matrix_array_col_major_fortran(mat, array, rows, cols, ld) &
          bind(c)
       use iso_c_binding
       type(c_ptr), value :: mat
       type(c_ptr), value, intent(in) :: array
       integer(c_int), value, intent(in) :: rows, cols, ld
     end subroutine generate_distributed_matrix_array_col_major_fortran

     subroutine generate_distributed_matrix_array_row_major_fortran(mat, array, rows, cols, ld) &
          bind(c)
       use iso_c_binding
       type(c_ptr), value :: mat
       type(c_ptr), value, intent(in) :: array
       integer(c_int), value, intent(in) :: rows, cols, ld
     end subroutine generate_distributed_matrix_array_row_major_fortran

     subroutine set_distributed_matrix_local_col_major(matrix, i, j, val) bind(c)
       use iso_c_binding
       type(c_ptr), intent(in), value :: matrix
       integer(c_int), intent(in), value :: i, j
       real(c_double), intent(in), value :: val
     end subroutine set_distributed_matrix_local_col_major

     subroutine set_distributed_matrix_local_row_major(matrix, i, j, val) bind(c)
       use iso_c_binding
       type(c_ptr), intent(in), value :: matrix
       integer(c_int), intent(in), value :: i, j
       real(c_double), intent(in), value :: val
     end subroutine set_distributed_matrix_local_row_major

     real(c_double) function get_distributed_matrix_local_col_major(matrix, i, j) bind(c)
       use iso_c_binding
       type(c_ptr), intent(in), value :: matrix
       integer(c_int), intent(in), value :: i, j
     end function get_distributed_matrix_local_col_major

     real(c_double) function get_distributed_matrix_local_row_major(matrix, i, j) bind(c)
       use iso_c_binding
       type(c_ptr), intent(in), value :: matrix
       integer(c_int), intent(in), value :: i, j
     end function get_distributed_matrix_local_row_major

  end interface

contains

  subroutine set_distributed_matrix(matrix, dim1, dim2, g, solver_, matrix_major)
    implicit none
    type(rokko_distributed_matrix), intent(inout) :: matrix
    integer, intent(in) :: dim1, dim2
    type(grid), intent(in) :: g
    type(solver), intent(in) :: solver_
    character(*), intent(in), optional :: matrix_major

    if (present(matrix_major)) then
       if (matrix_major == "matrix_col_major") then
          call initialize_distributed_matrix_col_major(matrix%ptr_distributed_matrix, dim1, dim2, &
               g%ptr_grid, solver_%ptr_solver)
       else if (matrix_major == "matrix_row_major") then
          call initialize_distributed_matrix_row_major(matrix%ptr_distributed_matrix, dim1, dim2, &
               g%ptr_grid, solver_%ptr_solver)
       else
          write(0,*) "Incorrect matrix_major. matrix_col_major or matrix_row_major is accepted"
       endif
    else
       call initialize_distributed_matrix_col_major(matrix%ptr_distributed_matrix, dim1, dim2, &
            g%ptr_grid, solver_%ptr_solver)
    endif
  end subroutine set_distributed_matrix

  subroutine del_distributed_matrix(matrix, matrix_major)
    implicit none
    type(distributed_matrix)::matrix
    character(*),intent(in),optional::matrix_major

    if (present(matrix_major)) then
      if (matrix_major == "matrix_col_major") then

        call delete_distributed_matrix_col_major(matrix&
           &%ptr_distributed_matrix)
      else if (matrix_major == "matrix_row_major") then
        call delete_distributed_matrix_row_major(matrix&
           &%ptr_distributed_matrix)
      else
        write(0,*) "Incorrect matrix_major. matrix_col_major or &
           &matrix_row_major is accepted"
      endif
    else
      call delete_distributed_matrix_col_major(matrix&
         &%ptr_distributed_matrix)
    endif
  end subroutine del_distributed_matrix

  subroutine set_localized_vector(vector, dim)
    implicit none
    type(localized_vector),intent(out)::vector
    integer,intent(in)::dim

    call initialize_localized_vector(vector%ptr_localized_vector, dim)
  end subroutine set_localized_vector

  subroutine del_localized_vector(vector)
    implicit none
    type(localized_vector),intent(inout)::vector

    call delete_localized_vector(vector%ptr_localized_vector)
  end subroutine del_localized_vector

  real(kind(1d0)) function get_element_localized_vector(w_, i)
    implicit none
    type(localized_vector),intent(in)::w_
    integer,intent(in)::i
    real(kind(1d0))::localized_vector_get_element
    external localized_vector_get_element

    get_element_localized_vector = localized_vector_get_element(w_&
       &%ptr_localized_vector, i-1)
  end function get_element_localized_vector

  subroutine copy_localized_vector(w,w_)
    implicit none
    real(kind(1d0)),intent(out)::w(:)
    type(localized_vector),intent(in)::w_
    integer::dim,i

    dim = size(w)

    do i=1,dim
      w(i) = get_element_localized_vector(w_, i)
    enddo
  end subroutine copy_localized_vector


  subroutine diagonalize(solver_in, matrix, w, Z, tm_in, matrix_major)
    implicit none
    type(solver), intent(inout) :: solver_in
    type(distributed_matrix), intent(inout) :: matrix, Z
    real(kind(1d0)) :: w(:)
    type(timer), intent(inout), optional :: tm_in
    character(*), intent(in), optional :: matrix_major
    type(timer) :: tm

    type(localized_vector) :: w_
    integer :: dim

    dim = size(w)
    call set_localized_vector(w_, dim)

    if (present(tm_in)) then
       if (present(matrix_major)) then
          if (matrix_major == "matrix_col_major") then
             call solver_diagonalize_matrix_col_major(solver_in%ptr_solver, &
                  matrix%ptr_distributed_matrix, w_%ptr_localized_vector, &
                  Z%ptr_distributed_matrix, tm_in%ptr_timer)
             call copy_localized_vector(w, w_)
          else if (matrix_major == "matrix_row_major") then
             call solver_diagonalize_matrix_row_major(solver_in%ptr_solver, &
                  matrix%ptr_distributed_matrix, w_%ptr_localized_vector, &
                  Z%ptr_distributed_matrix, tm_in%ptr_timer)
             call copy_localized_vector(w, w_)
          else
             write(0,*) "Incorrect matrix_major. matrix_col_major or m&
                  &atrix_row_major is accepted"
          endif
       else
          call solver_diagonalize_matrix_col_major(solver_in%ptr_solver, &
               matrix%ptr_distributed_matrix, w_%ptr_localized_vector, &
               Z%ptr_distributed_matrix, tm_in%ptr_timer)
          call copy_localized_vector(w, w_)
       endif
    else
       if (present(matrix_major)) then
          if (matrix_major == "matrix_col_major") then
             call solver_diagonalize_matrix_col_major(solver_in%ptr_solver, &
                  matrix%ptr_distributed_matrix, w_%ptr_localized_vector, &
                  Z%ptr_distributed_matrix, tm%ptr_timer)
             call copy_localized_vector(w, w_)
          else if (matrix_major == "matrix_row_major") then
             call solver_diagonalize_matrix_row_major(solver_in%ptr_solver, &
                  matrix%ptr_distributed_matrix, w_%ptr_localized_vector, &
                  Z%ptr_distributed_matrix, tm%ptr_timer)
             call copy_localized_vector(w, w_)
          else
             write(0,*) "Incorrect matrix_major. matrix_col_major or m&
                  &atrix_row_major is accepted"
          endif
       else
          call solver_diagonalize_matrix_col_major(solver_in%ptr_solver, &
               matrix%ptr_distributed_matrix, w_%ptr_localized_vector, &
               Z%ptr_distributed_matrix, tm%ptr_timer)
          call copy_localized_vector(w, w_)
       endif
    endif
    call del_localized_vector(w_)
  end subroutine diagonalize

! frank matrix generator
  subroutine generate_frank_matrix_distributed(matrix,matrix_major)
    implicit none
    type(distributed_matrix),intent(inout)::matrix
    character(*),intent(in),optional::matrix_major

    if (present(matrix_major)) then
      if (matrix_major == "matrix_col_major") then

        call frank_generate_distributed_matrix_col_major(matrix&
           &%ptr_distributed_matrix)
      else if (matrix_major == "matrix_row_major") then
        call frank_generate_distributed_matrix_row_major(matrix&
           &%ptr_distributed_matrix)
      else
        write(0,*) "Incorrect matrix_major. matrix_col_major or m&
           &atrix_row_major is accepted"
      endif
    else
      call frank_generate_distributed_matrix_col_major(matrix&
         &%ptr_distributed_matrix)
    endif
  end subroutine generate_frank_matrix_distributed

  subroutine generate_frank_matrix_localized(matrix,matrix_major)
    implicit none
    type(localized_matrix),intent(inout)::matrix
    character(*),intent(in),optional::matrix_major

    if (present(matrix_major)) then
      if (matrix_major == "matrix_col_major") then
        call frank_generate_localized_matrix_col_major(matrix&
           &%ptr_localized_matrix)
      else if (matrix_major == "matrix_row_major") then
        call frank_generate_localized_matrix_row_major(matrix&
           &%ptr_localized_matrix)
      else
        write(0,*) "Incorrect matrix_major. matrix_col_major or m&
           &atrix_row_major is accepted"
      endif
    else
      call frank_generate_localized_matrix_col_major(matrix&
         &%ptr_localized_matrix)
    endif
  end subroutine generate_frank_matrix_localized

  ! matrix generator according to a given function pointer
  subroutine generate_distributed_matrix_function(matrix, f, matrix_major)
    implicit none
    type(distributed_matrix), intent(inout) :: matrix
    interface
       real(c_double) function f(i, j) bind(c)
         use iso_c_binding
         integer(c_int), intent(in) :: i, j
       end function f
    end interface
    character(*),intent(in),optional :: matrix_major

    if (present(matrix_major)) then
       if (matrix_major == "matrix_col_major") then
          call generate_distributed_matrix_function_col_major_fortran(matrix%ptr_distributed_matrix, &
               c_funloc(f))
       else if (matrix_major == "matrix_row_major") then
          call generate_distributed_matrix_function_row_major_fortran(matrix%ptr_distributed_matrix, &
               c_funloc(f))
       else
          write(0,*) "Incorrect matrix_major. matrix_col_major or matrix_row_major&
               & is accepted"
       endif
    else
       call generate_distributed_matrix_function_col_major_fortran(matrix%ptr_distributed_matrix, &
            c_funloc(f))
    endif
  end subroutine generate_distributed_matrix_function

  ! matrix generator according to a given fortran array
  ! it is inefficient way. do not use in pracitcal calculation.
  subroutine generate_distributed_matrix_array(matrix, array, rows, cols, ld, matrix_major)
    use iso_c_binding
    implicit none
    type(distributed_matrix), intent(inout) :: matrix
    integer(c_int), intent(in) :: rows, cols, ld
    real(c_double), intent(in), target :: array(ld, cols)
    character(*), intent(in), optional :: matrix_major

    if (present(matrix_major)) then
       if (matrix_major == "matrix_col_major") then
          call generate_distributed_matrix_array_col_major_fortran(matrix%ptr_distributed_matrix, &
               c_loc(array(1,1)), rows, cols, ld)
       else if (matrix_major == "matrix_row_major") then
          call generate_distributed_matrix_array_row_major_fortran(matrix%ptr_distributed_matrix, &
               c_loc(array(1,1)), rows, cols, ld)
       else
          write(0,*) "Incorrect matrix_major. matrix_col_major or matrix_row_major&
               &is accepted"
       endif
    else
       call generate_distributed_matrix_array_col_major_fortran(matrix%ptr_distributed_matrix, &
            c_loc(array(1,1)), rows, cols, ld)
    endif
  end subroutine generate_distributed_matrix_array

! set a value matrix element to distributed_matrix
  subroutine set_distributed_matrix_local(matrix, i, j, val, matrix_major)
    implicit none
    type(distributed_matrix), intent(inout) :: matrix
    integer(c_int), intent(in), value :: i, j
    real(kind(1d0)), intent(in), value :: val
    character(*), intent(in), optional :: matrix_major

    if (present(matrix_major)) then
       if (matrix_major == "matrix_col_major") then
          call set_distributed_matrix_local_col_major(matrix&
               &%ptr_distributed_matrix, i-1, j-1, val)
       else if (matrix_major == "matrix_row_major") then
          call set_distributed_matrix_local_row_major(matrix&
               &%ptr_distributed_matrix, i-1, j-1, val)
       else
          write(0,*) "Incorrect matrix_major. matrix_col_major or m&
               &atrix_row_major is accepted"
       endif
    else
       call set_distributed_matrix_local_col_major(matrix&
            &%ptr_distributed_matrix, i-1, j-1, val)
    endif
  end subroutine set_distributed_matrix_local

! get a value matrix element to distributed_matrix
  real(c_double) function get_distributed_matrix_local(matrix, i, j, matrix_major)
    implicit none
    type(distributed_matrix), intent(inout) :: matrix
    integer(c_int), intent(in), value :: i, j
    character(*), intent(in), optional :: matrix_major

    if (present(matrix_major)) then
       if (matrix_major == "matrix_col_major") then
          get_distributed_matrix_local = get_distributed_matrix_local_col_major(matrix&
               &%ptr_distributed_matrix, i-1, j-1)
       else if (matrix_major == "matrix_row_major") then
          get_distributed_matrix_local = get_distributed_matrix_local_row_major(matrix&
               &%ptr_distributed_matrix, i-1, j-1)
       else
          write(0,*) "Incorrect matrix_major. matrix_col_major or m&
               &atrix_row_major is accepted"
       endif
    else
       get_distributed_matrix_local = get_distributed_matrix_local_col_major(matrix&
            &%ptr_distributed_matrix, i-1, j-1)
    endif
  end function get_distributed_matrix_local


! set a value matrix element to distributed_matrix
  subroutine print_distributed_matrix(matrix, matrix_major)
    implicit none
    type(distributed_matrix), intent(inout) :: matrix
    character(*), intent(in), optional :: matrix_major

    if (present(matrix_major)) then
       if (matrix_major == "matrix_col_major") then
          call print_distributed_matrix_col_major(matrix&
               &%ptr_distributed_matrix)
       else if (matrix_major == "matrix_row_major") then
          call print_distributed_matrix_row_major(matrix&
               &%ptr_distributed_matrix)
       else
          write(0,*) "Incorrect matrix_major. matrix_col_major or m&
               &atrix_row_major is accepted"
       endif
    else
       call print_distributed_matrix_col_major(matrix&
            &%ptr_distributed_matrix)
    endif
  end subroutine print_distributed_matrix

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

end module rokko
