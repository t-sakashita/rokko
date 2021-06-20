program main
  implicit none
  integer, parameter :: n = 5
  real :: a(n, n), norm, work(n)
  integer :: i, j
  real :: slange
  do j = 1, n
     do i = 1, n
        a(i, j) = min(i, j)
     end do
  end do
  print *, "Matrix A: ", n, n
  do j = 1, n
     print *, a(1:n, j)
  end do

  norm = slange('M', n, n, a, n, work)
  print *, "element of largest absolute value =", norm
  norm = slange('1', n, n, a, n, work)
  print *, "one norm =", norm
  norm = slange('I', n, n, a, n, work)
  print *, "infinity norm =", norm
  norm = slange('F', n, n, a, n, work)
  print *, "Frobenius norm =", norm
end program main
