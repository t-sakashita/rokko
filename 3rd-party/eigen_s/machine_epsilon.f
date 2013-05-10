*
*
* NOTICE : BECAUSE OF THE ACCURACY PROBLEMS,
*          DO NOT USE ANY OPTIMIZATION
*          OPTION ON THE COMPILATION PHASE.
*          IN MOST CASES, YOU MAY USE -O0.
*
*
       real(8) function machine_epsilon()
       implicit NONE

       real(8)                :: epsilon, tmp

       include 'param.h'


       epsilon = ONE
       do
          tmp = epsilon
          epsilon = epsilon * HALF
          if ( ONE + epsilon == ONE ) exit
       end do
       machine_epsilon = tmp


       return
       end function machine_epsilon
