module input_parameters_mod
  implicit none

contains

  subroutine read_input_parameters(na, nev, nblk)
    implicit none

    integer, intent(out) :: na, nev, nblk
    ! Command line arguments
    character(len=128)   :: arg

    ! default parameters
    na = 1000
    nev = na
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

end module input_parameters_mod
