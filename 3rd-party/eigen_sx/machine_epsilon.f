*
*
* NOTICE : BECAUSE OF THE ACCURACY PROBLEMS, DO NOT USE ANY OPTIMIZATION
*          OPTION ON THE COMPILATION PHASE
*
*
       real(8) function machine_epsilon()

       implicit NONE
       real(8)                :: epsilon, tmp

       epsilon = 1.0D+00
       do
          tmp = epsilon
          epsilon = epsilon*5.0D-01
          if ( 1.0D+00+epsilon == 1.0D+00 ) exit
       enddo
       machine_epsilon = tmp

       return
       end function ! machine_epsilon
