subroutine func(i, j)
  use iso_c_binding
  integer(c_int) :: i
  integer(c_int) :: j
  i = j
end subroutine

program main
  call func(1, 2)
end program

