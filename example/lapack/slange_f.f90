program main
  implicit none
  integer, parameter :: n = 5
  real :: a(n, n), norm
  integer :: i, j
  double precision :: slange
  do j = 1, n
     do i = 1, n
        a(i, j) = n - 0.253E0 * max(i-1, j-1)
     end do
  end do
  print *, "Matrix A: ", n, n
  do j = 1, n
     print *, a(1:n, j)
  end do

  norm = slange('F', n, n, a, n)
  print *, "norm =", norm
end program main
