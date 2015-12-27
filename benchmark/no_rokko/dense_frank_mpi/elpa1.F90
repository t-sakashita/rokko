!    This file is modified by T. Sakashita.
!    It is originally part of ELPA.
!
!    The ELPA library was originally created by the ELPA consortium,
!    consisting of the following organizations:
!
!    - Rechenzentrum Garching der Max-Planck-Gesellschaft (RZG),
!    - Bergische Universität Wuppertal, Lehrstuhl für angewandte
!      Informatik,
!    - Technische Universität München, Lehrstuhl für Informatik mit
!      Schwerpunkt Wissenschaftliches Rechnen ,
!    - Fritz-Haber-Institut, Berlin, Abt. Theorie,
!    - Max-Plack-Institut für Mathematik in den Naturwissenschaften,
!      Leipzig, Abt. Komplexe Strukutren in Biologie und Kognition,
!      and
!    - IBM Deutschland GmbH
!
!
!    More information can be found here:
!    http://elpa.rzg.mpg.de/
!
!    ELPA is free software: you can redistribute it and/or modify
!    it under the terms of the version 3 of the license of the
!    GNU Lesser General Public License as published by the Free
!    Software Foundation.
!
!    ELPA is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU Lesser General Public License for more details.
!
!    You should have received a copy of the GNU Lesser General Public License
!    along with ELPA.  If not, see <http://www.gnu.org/licenses/>
!
!    ELPA reflects a substantial effort on the part of the original
!    ELPA consortium, and we ask you to respect the spirit of the
!    license that we chose: i.e., please contribute any changes you
!    may have back to the original ELPA library distribution, and keep
!    any derivatives of ELPA under the same license that we chose for
!    the original distribution, the GNU Lesser General Public License.
!
!
!>
!> By calling executable [arg1] [arg2] [arg3]
!> one can define the size (arg1), the number of
!> Eigenvectors to compute (arg2), and the blocking (arg3).
!> If these values are not set default values (4000, 1500, 16)
!> are choosen.
!>
!> The real ELPA 2 kernel is set as the default kernel.
!> However, this can be overriden by setting
!> the environment variable "REAL_ELPA_KERNEL" to an
!> appropiate value.
!>

subroutine read_input_parameters(na, nev, nblk)
  implicit none

  integer, intent(out) :: na, nev, nblk
  ! Command line arguments
  character(len=128)   :: arg
  
  ! default parameters
  na = 4000
  nev = 1500
  nblk = 16

  if (COMMAND_ARGUMENT_COUNT() >= 1) then
     call GET_COMMAND_ARGUMENT(1, arg)
     read(arg, *) na
     if (COMMAND_ARGUMENT_COUNT() >= 2) then
        call GET_COMMAND_ARGUMENT(2, arg)
        read(arg, *) nev
     else
        nev = na
     endif
  endif

  if (COMMAND_ARGUMENT_COUNT() >= 3) then
     call GET_COMMAND_ARGUMENT(3, arg)     
     read(arg, *) nblk
  endif

end subroutine read_input_parameters

subroutine set_up_blacsgrid(mpi_comm_world, my_blacs_ctxt, np_rows, &
                                np_cols, nprow, npcol, my_prow, my_pcol)

  implicit none
  integer, intent(in)     :: mpi_comm_world
  integer, intent(inout)  :: my_blacs_ctxt, np_rows, &
       np_cols, nprow, npcol, my_prow, my_pcol
  
  my_blacs_ctxt = mpi_comm_world
  call BLACS_Gridinit(my_blacs_ctxt, 'C', np_rows, np_cols)
  call BLACS_Gridinfo(my_blacs_ctxt, nprow, npcol, my_prow, my_pcol)

end subroutine set_up_blacsgrid
    
subroutine set_up_blacs_descriptor(na, nblk, my_prow, my_pcol, &
     np_rows, np_cols, na_rows,  &
     na_cols, sc_desc, my_blacs_ctxt, info)

  use elpa_utilities, only : error_unit
  use mpi
  implicit none
  
  integer, intent(inout)  :: na, nblk, my_prow, my_pcol, np_rows,   &
       np_cols, na_rows, na_cols, sc_desc(1:9), &
       my_blacs_ctxt, info
  
  integer, external       :: numroc
  integer                 :: mpierr
  
  ! determine the neccessary size of the distributed matrices,
  ! we use the scalapack tools routine NUMROC
  
  na_rows = numroc(na, nblk, my_prow, 0, np_rows)
  na_cols = numroc(na, nblk, my_pcol, 0, np_cols)
  
  ! set up the scalapack descriptor for the checks below
  ! For ELPA the following restrictions hold:
  ! - block sizes in both directions must be identical (args 4 a. 5)
  ! - first row and column of the distributed matrix must be on
  !   row/col 0/0 (arg 6 and 7)
  
  call descinit(sc_desc, na, na, nblk, nblk, 0, 0, my_blacs_ctxt, na_rows, info)

  if (info .ne. 0) then
     write(error_unit,*) 'Error in BLACS descinit! info=',info
     write(error_unit,*) 'Most likely this happend since you want to use'
     write(error_unit,*) 'more MPI tasks than are possible for your'
     write(error_unit,*) 'problem size (matrix size and blocksize)!'
     write(error_unit,*) 'The blacsgrid can not be set up properly'
     write(error_unit,*) 'Try reducing the number of MPI tasks...'
     call MPI_ABORT(mpi_comm_world, 1, mpierr)
  endif
  
end subroutine set_up_blacs_descriptor


subroutine generate_matrix( n, a, desca, info )
  !     ..
  !     .. Scalar Arguments ..
  integer            n
  !     ..
  !     .. Array Arguments ..
  integer            desca( * )
  double precision   a( * )
  !     ..
  !     .. local Scalars ..
  integer            i, j
  !     ..
  !     .. External Subroutines ..
  external           pdelset

  ! Create Frank matrix
  do j = 1, n
     do i = 1, n
        call pdelset( a, i, j, desca, dble(n - max(i,j) + 1))
     enddo
  enddo
 
end subroutine generate_matrix

      
program test_real1

!-------------------------------------------------------------------------------
! Standard eigenvalue problem - REAL version
!
! This program demonstrates the use of the ELPA module
! together with standard scalapack routines
!
! Copyright of the original code rests with the authors inside the ELPA
! consortium. The copyright of any additional modifications shall rest
! with their original authors, but shall adhere to the licensing terms
! distributed along with the original code in the file "COPYING".
!
!-------------------------------------------------------------------------------

   use ELPA1
 !  use ELPA2
   use elpa_utilities, only : error_unit
   use mpi
   implicit none

   !-------------------------------------------------------------------------------
   ! Please set system size parameters below!
   ! na:   System size
   ! nev:  Number of eigenvectors to be calculated
   ! nblk: Blocking factor in block cyclic distribution
   !-------------------------------------------------------------------------------
   integer :: nblk, na, nev

   !-------------------------------------------------------------------------------
   !  Local Variables
   integer :: omp_get_max_threads, required_mpi_thread_level, provided_mpi_thread_level
   integer :: np_rows, np_cols, na_rows, na_cols
   integer :: myid, nprocs, my_prow, my_pcol, mpi_comm_rows, mpi_comm_cols
   integer :: i, mpierr, my_blacs_ctxt, sc_desc(9), info, nprow, npcol
   logical :: success

   double precision, allocatable :: a(:,:), z(:,:), ev(:)
   double precision :: init_tick, gen_tick, diag_tick, end_tick

   !  MPI Initialization
   call mpi_init_thread(MPI_THREAD_MULTIPLE, provided_mpi_thread_level, mpierr)
   call mpi_comm_rank(mpi_comm_world,myid,mpierr)
   call mpi_comm_size(mpi_comm_world,nprocs,mpierr)

   call read_input_parameters(na, nev, nblk)

#ifdef WITH_OPENMP
   if (myid .eq. 0) then
      print *,"Threaded version of test program"
      print *,"Using ",omp_get_max_threads()," threads"
      print *," "
   endif
#endif

   if (myid .eq. 0) then
      print *," "
      print *,"This ELPA1 is build with"
#ifdef WITH_REAL_AVX_BLOCK2_KERNEL
      print *,"AVX optimized kernel (2 blocking) for real matrices"
#endif
#ifdef WITH_REAL_AVX_BLOCK4_KERNEL
      print *,"AVX optimized kernel (4 blocking) for real matrices"
#endif
#ifdef WITH_REAL_AVX_BLOCK6_KERNEL
      print *,"AVX optimized kernel (6 blocking) for real matrices"
#endif

#ifdef WITH_REAL_GENERIC_KERNEL
     print *,"GENERIC kernel for real matrices"
#endif
#ifdef WITH_REAL_GENERIC_SIMPLE_KERNEL
     print *,"GENERIC SIMPLE kernel for real matrices"
#endif
#ifdef WITH_REAL_SSE_KERNEL
     print *,"SSE ASSEMBLER kernel for real matrices"
#endif
#ifdef WITH_REAL_BGP_KERNEL
     print *,"BGP kernel for real matrices"
#endif
#ifdef WITH_REAL_BGQ_KERNEL
     print *,"BGQ kernel for real matrices"
#endif
   endif

   !-------------------------------------------------------------------------------
   ! Selection of number of processor rows/columns
   ! We try to set up the grid square-like, i.e. start the search for possible
   ! divisors of nprocs with a number next to the square root of nprocs
   ! and decrement it until a divisor is found.

   do np_cols = nint(sqrt(real(nprocs))),2,-1
      if(mod(nprocs,np_cols) == 0 ) exit
   enddo
   ! at the end of the above loop, nprocs is always divisible by np_cols

   np_rows = nprocs/np_cols

   if(myid==0) then
      print *
      print '(a)','Standard eigenvalue problem - REAL version'
      print *
      print '(3(a,i0))','Matrix size=',na,', Number of eigenvectors=',nev,', Block size=',nblk
      print '(3(a,i0))','Number of processor rows=',np_rows,', cols=',np_cols,', total=',nprocs
      print *
   endif

   !-------------------------------------------------------------------------------
   ! Set up BLACS context and MPI communicators
   !
   ! The BLACS context is only necessary for using Scalapack.
   !
   ! For ELPA, the MPI communicators along rows/cols are sufficient,
   ! and the grid setup may be done in an arbitrary way as long as it is
   ! consistent (i.e. 0<=my_prow<np_rows, 0<=my_pcol<np_cols and every
   ! process has a unique (my_prow,my_pcol) pair).
   init_tick = mpi_wtime()
   call set_up_blacsgrid(mpi_comm_world, my_blacs_ctxt, np_rows, np_cols, &
                         nprow, npcol, my_prow, my_pcol)

   ! All ELPA routines need MPI communicators for communicating within
   ! rows or columns of processes, these are set in get_elpa_row_col_comms.
   mpierr = get_elpa_row_col_comms(mpi_comm_world, my_prow, my_pcol, &
                                   mpi_comm_rows, mpi_comm_cols)

   gen_tick = mpi_wtime()
   call set_up_blacs_descriptor(na ,nblk, my_prow, my_pcol, np_rows, np_cols, &
                                na_rows, na_cols, sc_desc, my_blacs_ctxt, info)

   !-------------------------------------------------------------------------------
   ! Allocate matrices and set up a test matrix for the eigenvalue problem
   allocate(a(na_rows,na_cols))
   allocate(z(na_rows,na_cols))
   allocate(ev(na))

   CALL generate_matrix( na, A, SC_DESC, INFO )

   !-------------------------------------------------------------------------------
   ! Calculate eigenvalues/eigenvectors
   diag_tick = mpi_wtime()
   success = solve_evp_real(na, nev, a, na_rows, ev, z, na_rows, nblk, &
                            na_cols, mpi_comm_rows, mpi_comm_cols)
   end_tick = mpi_wtime()

   if (.not.(success)) then
      write(error_unit,*) "solve_evp_real_2stage produced an error! Aborting..."
      call MPI_ABORT(mpi_comm_world, 1, mpierr)
   endif

   if (myid == 0) then
      print *, i,ev(1:na)
      print *, "init_time = ", gen_tick - init_tick
      print *, "gen_time = ", diag_tick - gen_tick
      print *, "diag_time = ", end_tick - diag_tick
   endif

   deallocate(a)
   deallocate(z)
   deallocate(ev)

   call blacs_gridexit(my_blacs_ctxt)
   call mpi_finalize(mpierr)

end program test_real1

