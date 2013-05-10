       subroutine matrix_set(n, a)
       implicit double precision(a-h,o-z),integer(i-n)
       real(8)   :: a(n,n)

       do i=1,n
          do j=1,n
             a(j,i)=(n+1-Max(n+1-i,n+1-j))*1.0D+00
ccc          CALL RANDOM_NUMBER(t)
ccc          a(j,i)=1.D0-2*t
ccc          if(i==j)then
ccc             a(j,i)=-7.2D0
ccc          else
ccc             a(j,i)=-3.0D0/(i-j)**2
ccc          endif
          enddo
       enddo

       return
       end

