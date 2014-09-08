c>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
c        ACKNOWLEDGE THE USE OF THIS PACKAGE
c         WHEN YOU PUBLISH YOUR RESULTS !!!
c>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
c
c************ Sample main program #11 *****************
c  1d Heisenberg antiferromagnet with 8 spins
c    Eigenvalues and an eigenvector by diag
c*************************************************
      parameter (n=14,idim=3432,ibond=n)
      implicit real*8 (a-h,o-z)
      dimension E(4)
      dimension elemnt(idim,idim)
      dimension list1(idim),list2(2,0:2**15)
      dimension bondwt(ibond),ipair(2*ibond),zrtio(ibond)
      dimension npair(2),sxx(1),szz(1)
      dimension wk(idim,8),iwk(idim),v(idim)
c
      data bondwt/ibond*-1.0d0/
      data zrtio/ibond*1.0d0/
      data ipair/1,2, 2,3, 3,4, 4,5, 5,6, 6,7, 7,8, 8,9,
     &     9,10, 10,11, 11,12, 12,13, 13,14, 14,1/
c
      call clock(t1)
      call sz(n,idim,0.0d0,list1,list2)
c- You may alternatively use szdy or sztn for faster processing -
c      call szdy(n,idim,0.0d0,list1,list2)
c  or
c      call sztn(n,idim,0.0d0,list1,list2)
c------------------------------------------------------------
      call elm3(n,idim,ipair,bondwt,zrtio,ibond,
     &          elemnt,idim,list1,list2)
c
      call clock(t2)
      eps=1.d-13
      nvec=1
      ne=4
      call diag(elemnt,idim,idim,E,v,ne,nvec,eps,wk,iwk)
      call clock(t3)
      print 100,E
  100 format(/' [Eigenvalues]  '/2x,4f14.8)
c*** Do not forget to call elm3 again before calling check3
      call elm3(n,idim,ipair,bondwt,zrtio,ibond,
     &          elemnt,idim,list1,list2)
      call check3(elemnt,idim,idim,v,wk,Hexpec)
      call clock(t4)
c
      npair(1)=1
      npair(2)=2
      call xcorr(n,idim,npair,1,v,sxx,list1,list2)
      print *,' sxx:',sxx
      call zcorr(n,idim,npair,1,v,szz,list1)
      print *,' szz:',szz
      call clock(t5)
      write(6,'(A,f8.3,A)') 'initialize      ', (t2-t1), ' sec'
      write(6,'(A,f8.3,A)') 'diagonalization ', (t3-t2), ' sec'
      write(6,'(A,f8.3,A)') 'check           ', (t4-t3), ' sec'
      write(6,'(A,f8.3,A)') 'correlation     ', (t5-t4), ' sec'
      end
