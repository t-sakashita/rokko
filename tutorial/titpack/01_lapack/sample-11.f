c>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
c        ACKNOWLEDGE THE USE OF THIS PACKAGE
c         WHEN YOU PUBLISH YOUR RESULTS !!!
c>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
c
c************ Sample main program #11 *****************
c  1d Heisenberg antiferromagnet with 8 spins
c    Eigenvalues and an eigenvector by diag
c*************************************************
      parameter (n=8,idim=70,ibond=n)
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
      data ipair/1,2, 2,3, 3,4, 4,5, 5,6, 6,7, 7,8, 8,1/
c
      call sz(n,idim,0.0d0,list1,list2)
c- You may alternatively use szdy or sztn for faster processing -
c      call szdy(n,idim,0.0d0,list1,list2)
c  or
c      call sztn(n,idim,0.0d0,list1,list2)
c------------------------------------------------------------
      call elm3(n,idim,ipair,bondwt,zrtio,ibond,
     &          elemnt,idim,list1,list2)
c
      nvec=1
      ne=4
      call diag(elemnt,idim,idim,E,v,ne,nvec,wk,iwk)
      print 100,E
  100 format(/' [Eigenvalues]  '/2x,4f14.8)
c*** Do not forget to call elm3 again before calling check3
      call elm3(n,idim,ipair,bondwt,zrtio,ibond,
     &          elemnt,idim,list1,list2)
      call check3(elemnt,idim,idim,v,wk,Hexpec)
c
      npair(1)=1
      npair(2)=2
      call xcorr(n,idim,npair,1,v,sxx,list1,list2)
      print *,' sxx:',sxx
      call zcorr(n,idim,npair,1,v,szz,list1)
      print *,' szz:',szz
      end
