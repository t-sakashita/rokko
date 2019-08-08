module my_elpa_utility
  implicit none

contains

  subroutine prepare_matrix_minij(na, a, nblk, np_rows, np_cols, my_prow, my_pcol)
    use iso_c_binding
    implicit none
    
    integer, intent(in)           :: na, nblk, np_rows, np_cols, my_prow, my_pcol
    real(kind=c_double)    :: a(:,:)
    integer                       :: i, j, rowLocal, colLocal
    
    do i = 1, na
       do j = 1, na
          if (map_global_array_index_to_local_index(i, j, rowLocal, colLocal, nblk, np_rows, np_cols, my_prow, my_pcol)) then
             a(rowLocal,colLocal) = min(i,j)
          endif
       enddo
    enddo
    
  end subroutine prepare_matrix_minij
  
  !Processor col for global col number
  pure function pcol(global_col, nblk, np_cols) result(local_col)
    use iso_c_binding, only : c_int
    implicit none
    integer(kind=c_int), intent(in) :: global_col, nblk, np_cols
    integer(kind=c_int)             :: local_col
    local_col = MOD((global_col-1)/nblk,np_cols)
  end function pcol
  
  !Processor row for global row number
  pure function prow(global_row, nblk, np_rows) result(local_row)
    use iso_c_binding, only : c_int
    implicit none
    integer(kind=c_int), intent(in) :: global_row, nblk, np_rows
    integer(kind=c_int)             :: local_row
    local_row = MOD((global_row-1)/nblk,np_rows)
  end function prow
  
  function map_global_array_index_to_local_index(iGLobal, jGlobal, iLocal, jLocal , nblk, np_rows, np_cols, my_prow, my_pcol) &
       result(possible)
    use iso_c_binding, only : c_int
    implicit none
    
    integer(kind=c_int)              :: pi, pj, li, lj, xi, xj
    integer(kind=c_int), intent(in)  :: iGlobal, jGlobal, nblk, np_rows, np_cols, my_prow, my_pcol
    integer(kind=c_int), intent(out) :: iLocal, jLocal
    logical                       :: possible
    
    possible = .true.
    iLocal = 0
    jLocal = 0
    
    pi = prow(iGlobal, nblk, np_rows)
    
    if (my_prow .ne. pi) then
       possible = .false.
       return
    endif
    
    pj = pcol(jGlobal, nblk, np_cols)
    
    if (my_pcol .ne. pj) then
       possible = .false.
       return
    endif
    li = (iGlobal-1)/(np_rows*nblk) ! block number for rows
    lj = (jGlobal-1)/(np_cols*nblk) ! block number for columns
    
    xi = mod( (iGlobal-1),nblk)+1   ! offset in block li
    xj = mod( (jGlobal-1),nblk)+1   ! offset in block lj
    
    iLocal = li * nblk + xi
    jLocal = lj * nblk + xj
  end function map_global_array_index_to_local_index

end module my_elpa_utility
