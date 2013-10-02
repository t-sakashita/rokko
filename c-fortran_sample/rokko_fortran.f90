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

module rokko
  use ISO_C_BINDING
  implicit none
  
  type::distributed_matrix
     Type(c_ptr) ptr_distributed_matrix
  end type distributed_matrix

  type::localized_matrix
     Type(c_ptr) ptr_localized_matrix
  end type localized_matrix

  type::grid
     Type(c_ptr) ptr_grid
  end type grid

  type::solver
     Type(c_ptr) ptr_solver
  end type solver

  type::localized_vector
     Type(c_ptr) ptr_localized_vector
  end type localized_vector

  type::timer
     Type(c_ptr) ptr_timer
  end type timer

  interface generate_frank_matrix
     module procedure generate_frank_matrix_distributed, generate_frank_matrix_localized
  end interface
contains
  subroutine set_solver(solver_, solver_name_)
    implicit none
    Type(solver),intent(inout)::solver_
    character(*),intent(in):: solver_name_
    integer::argc
    character(len=1)::argv
    
    !    integer::argc,lenght_solver_name
    !    character::argv
    !    Type(c_ptr)::initialize_solver    
    !    external initialize_solver
    character(len=1),parameter::NULL =char(0)

    write(*,*) trim(solver_name_), len(solver_name_)
    !solver_%ptr_solver = initialize_solver(trim(solver_name_)//NULL)
    call initialize_solver(solver_%ptr_solver, trim(solver_name_)&
         &//NULL, argc, argv)
    
  end subroutine set_solver
  subroutine del_solver(solver_)
    implicit none
    Type(solver),intent(inout)::solver_
    
    call delete_solver(solver_%ptr_solver)
  end subroutine del_solver
  
  
  subroutine set_grid(grid_, comm_, grid_major_type)
    implicit none
    Type(grid),intent(out)::grid_
    integer,intent(in)::comm_
    Type(c_ptr)::test
    character(len=100),intent(in),optional::grid_major_type
!    Type(c_ptr)::initialize_grid_col_major,initialize_grid_row_major
!    external initialize_grid_col_major,initialize_grid_row_major


    write(*,*) "MPI_in = ",comm_
    if (present(grid_major_type)) then
      if (grid_major_type == "grid_col_major") then
        
!        grid_%ptr_grid = initialize_grid_col_major(comm_)
        call initialize_grid_col_major(grid_%ptr_grid, comm_)
        
      else if (grid_major_type == "grid_row_major") then

!        grid_%ptr_grid = initialize_grid_row_major(comm_)
        call initialize_grid_row_major(grid_%ptr_grid, comm_)
        
      else
        write(0,*) "Incorrect grid_major_type. grid_col_major or grid&
           &_row_major is accepted"
      endif
    else
!      grid_%ptr_grid = initialize_grid_col_major(comm_)
        call initialize_grid_col_major(test, comm_)
        grid_%ptr_grid = test
    end if
  end subroutine set_grid

  subroutine del_grid(grid_)
    implicit none
    Type(grid),intent(out)::grid_
    
    call delete_grid(grid_%ptr_grid)
  end subroutine del_grid

  integer function get_myrank_grid(grid_)
    implicit none
    Type(grid),intent(in)::grid_
    integer::grid_get_myrank
    external grid_get_myrank

    get_myrank_grid = grid_get_myrank(grid_%ptr_grid)
  end function get_myrank_grid

  integer function get_nprocs_grid(grid_)
    implicit none
    Type(grid),intent(in)::grid_
    integer::grid_get_nprocs
    external grid_get_nprocs

    get_nprocs_grid = grid_get_nprocs(grid_%ptr_grid)
  end function get_nprocs_grid
    
  subroutine set_distributed_matrix(matrix, dim1, dim2, g,&
     & solver_, matrix_major_type)
    implicit none
    Type(distributed_matrix),intent(inout)::matrix
    integer,intent(in)::dim1,dim2
    Type(grid),intent(in)::g
    Type(solver),intent(in)::solver_
    character(*),intent(in),optional::matrix_major_type

    if (present(matrix_major_type)) then
      if (matrix_major_type == "matrix_col_major") then

        call initialize_distributed_matrix_col_major(matrix&
           &%ptr_distributed_matrix, dim1, dim2, g&
           &%ptr_grid, solver_%ptr_solver)
      else if (matrix_major_type == "matrix_row_major") then
        call initialize_distributed_matrix_row_major(matrix&
           &%ptr_distributed_matrix, dim1, dim2, g&
           &%ptr_grid, solver_%ptr_solver)
      else
        write(0,*) "Incorrect matrix_major_type. matrix_col_major or m&
           &atrix_row_major is accepted"
      endif
    else
      call initialize_distributed_matrix_col_major(matrix&
         &%ptr_distributed_matrix, dim1, dim2, g&
         &%ptr_grid, solver_%ptr_solver)
    endif
  end subroutine set_distributed_matrix

  subroutine del_distributed_matrix(matrix, matrix_major_type)
    implicit none
    Type(distributed_matrix)::matrix
    character(*),intent(in),optional::matrix_major_type
    
    if (present(matrix_major_type)) then
      if (matrix_major_type == "matrix_col_major") then

        call delete_distributed_matrix_col_major(matrix&
           &%ptr_distributed_matrix)
      else if (matrix_major_type == "matrix_row_major") then
        call delete_distributed_matrix_row_major(matrix&
           &%ptr_distributed_matrix)
      else
        write(0,*) "Incorrect matrix_major_type. matrix_col_major or m&
           &atrix_row_major is accepted"
      endif
    else
      call delete_distributed_matrix_col_major(matrix&
         &%ptr_distributed_matrix)
    endif
  end subroutine del_distributed_matrix

  


  subroutine set_localized_vector(vector, dim)
    implicit none
    Type(localized_vector),intent(out)::vector
    integer,intent(in)::dim
    
    call initialize_localized_vector(vector%ptr_localized_vector, dim)
  end subroutine set_localized_vector

  subroutine del_localized_vector(vector)
    implicit none
    Type(localized_vector),intent(inout)::vector

    call delete_localized_vector(vector%ptr_localized_vector)
  end subroutine del_localized_vector

  real(kind(1d0)) function get_element_localized_vector(w_, i)
    implicit none
    Type(localized_vector),intent(in)::w_
    integer,intent(in)::i
    real(kind(1d0))::localized_vector_get_element
    external localized_vector_get_element
    
    get_element_localized_vector = localized_vector_get_element(w_&
       &%ptr_localized_vector, i-1)
  end function get_element_localized_vector
    
  subroutine copy_localized_vector(w,w_)
    implicit none
    real(kind(1d0)),intent(out)::w(:)
    Type(localized_vector),intent(in)::w_
    integer::dim,i
    
    dim = size(w)

    do i=1,dim
      w(i) = get_element_localized_vector(w_, i)
    enddo
  end subroutine copy_localized_vector


  subroutine Diagonalize(solver_, matrix, w, Z, timer_, matrix_major_type)
    implicit none
    Type(solver),intent(inout)::solver_
    Type(distributed_matrix),intent(inout)::matrix,Z
    real(kind(1d0))::w(:)
    Type(timer),intent(inout)::timer_
    character(*),intent(in),optional::matrix_major_type

    Type(localized_vector)::w_
    integer::dim
        
    dim = size(w)
    call set_localized_vector(w_, dim)

    if (present(matrix_major_type)) then
      if (matrix_major_type == "matrix_col_major") then

        call solver_diagonalize_matrix_col_major(solver_%ptr_solver, matrix&
           &%ptr_distributed_matrix, w_&
           &%ptr_localized_vector, Z%ptr_distributed_matrix, timer_%ptr_timer)

        call copy_localized_vector(w,w_)       

      else if (matrix_major_type == "matrix_row_major") then
        call solver_diagonalize_matrix_row_major(solver_%ptr_solver, matrix&
           &%ptr_distributed_matrix, w_&
           &%ptr_localized_vector, Z%ptr_distributed_matrix, timer_%ptr_timer)

        call copy_localized_vector(w,w_)

      else
        write(0,*) "Incorrect matrix_major_type. matrix_col_major or m&
           &atrix_row_major is accepted"
      endif
    else
      call solver_diagonalize_matrix_col_major(solver_%ptr_solver, matrix&
         &%ptr_distributed_matrix, w_&
         &%ptr_localized_vector, Z%ptr_distributed_matrix, timer_%ptr_timer)

      call copy_localized_vector(w,w_)
    endif
    call del_localized_vector(w_)
  end subroutine Diagonalize    

! frank matrix generator
  subroutine generate_frank_matrix_distributed(matrix,matrix_major_type)
    implicit none
    Type(distributed_matrix),intent(inout)::matrix
    character(*),intent(in),optional::matrix_major_type

    if (present(matrix_major_type)) then
      if (matrix_major_type == "matrix_col_major") then
        
        call frank_generate_distributed_matrix_col_major(matrix&
           &%ptr_distributed_matrix)
      else if (matrix_major_type == "matrix_row_major") then
        call frank_generate_distributed_matrix_row_major(matrix&
           &%ptr_distributed_matrix)
      else
        write(0,*) "Incorrect matrix_major_type. matrix_col_major or m&
           &atrix_row_major is accepted"
      endif
    else
      call frank_generate_distributed_matrix_col_major(matrix&
         &%ptr_distributed_matrix)
    endif
  end subroutine generate_frank_matrix_distributed
    
  subroutine generate_frank_matrix_localized(matrix,matrix_major_type)
    implicit none
    Type(localized_matrix),intent(inout)::matrix
    character(*),intent(in),optional::matrix_major_type

    if (present(matrix_major_type)) then
      if (matrix_major_type == "matrix_col_major") then
        
        call frank_generate_localized_matrix_col_major(matrix&
           &%ptr_localized_matrix)
      else if (matrix_major_type == "matrix_row_major") then
        call frank_generate_localized_matrix_row_major(matrix&
           &%ptr_localized_matrix)
      else
        write(0,*) "Incorrect matrix_major_type. matrix_col_major or m&
           &atrix_row_major is accepted"
      endif
    else
      call frank_generate_localized_matrix_col_major(matrix&
         &%ptr_localized_matrix)
    endif
  end subroutine generate_frank_matrix_localized


!timer routines
  subroutine set_timer(timer_)
    implicit none
    Type(timer),intent(out)::timer_

    call initialize_timer(timer_%ptr_timer)
  end subroutine set_timer

  subroutine del_timer(timer_)
    implicit none
    Type(timer),intent(out)::timer_

    call delete_timer(timer_%ptr_timer)
  end subroutine del_timer

  subroutine start_timer(timer_, id)
    Type(timer),intent(inout)::timer_
    integer,intent(in):: id
    call timer_start(timer_%ptr_timer, id)
  end subroutine start_timer

  subroutine end_timer(timer_, id)
    Type(timer),intent(inout)::timer_
    integer,intent(in):: id
    call timer_stop(timer_%ptr_timer, id)
  end subroutine end_timer

  subroutine registrate_timer(timer_, id, label)
    implicit none
    Type(timer),intent(inout)::timer_
    integer,intent(in)::id
    character(*),intent(in)::label
    call timer_registrate(timer_%ptr_timer, id, label)
    
  end subroutine registrate_timer
    
  real(kind(1d0)) function get_count_timer(timer_, id)
    Type(timer),intent(inout)::timer_
    integer,intent(in):: id
    real(kind(1d0))::timer_get_count
    external timer_get_count

    get_count_timer = timer_get_count(timer_%ptr_timer, id)

  end function get_count_timer

  real(kind(1d0)) function get_average_timer(timer_, id)
    Type(timer),intent(inout)::timer_
    integer,intent(in):: id
    real(kind(1d0))::timer_get_average
    external timer_get_average

    get_average_timer = timer_get_average(timer_%ptr_timer, id)
  end function get_average_timer
   
end module rokko
