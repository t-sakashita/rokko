c>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
c        ACKNOWLEDGE THE USE OF THIS PACKAGE
c         WHEN YOU PUBLISH YOUR RESULTS !!!
c>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
c
c************* Sample main program #6 *****************
c  1d Heisenberg antiferromagnet with 16 spins
c   Eigenvalues and an eigenvector/ lnc2, lncv2
c    Precision check and correlation functions
c***********************************************
c
      parameter (n=16,idim=12870,ibond=n,ic=ibond+1)
      implicit real*8 (a-h,o-z)
      dimension E(4)
      dimension elemnt(idim,ic),loc(idim,ic)
      dimension list1(idim),list2(2,0:2**15)
      dimension bondwt(ibond),ipair(2*ibond),zrtio(ibond)
      dimension npair(2)
      dimension wk(idim,2)
      dimension x(idim),sxx(1),szz(1)
c
      data bondwt/ibond*-1.0d0/
      data zrtio/ibond*1.0d0/
      data ipair/1,2, 2,3,  3,4,  4,5,  5,6,  6,7,  7,8, 8,9,
     &     9,10, 10,11, 11,12, 12,13, 13,14, 14,15, 15,16, 16,1/
      nvec=1
      iv=idim/5
c
      call sz(n,idim,0.0d0,list1,list2)
c- You may alternatively use szdy or sztn for faster processing -
c      call szdy(n,idim,0.0d0,list1,list2)
c  or
c      call sztn(n,idim,0.0d0,list1,list2)
c------------------------------------------------------------
      call elm2(n,idim,ipair,bondwt,zrtio,ibond,
     &          elemnt,loc,idim,ic,list1,list2)
      call lnc2(elemnt,loc,idim,idim,ic,nvec,iv,E,itr,wk)
      print 100,e,itr
  100 format(/' [Eigenvalues]  '/2x,4f14.8
     &       /' [Iteration number]'/i8)
      call lncv2(elemnt,loc,idim,idim,ic,nvec,iv,x,itr,wk)
c- You may alternatively use inv2 / Note: dimension wk(idim,4)-
c    call inv2(elemnt,loc,idim,idim,ic,iv,E(1),x,wk)
c---------------------------------------------------------------
      print *,'[Eigenvector components (selected)]'
      print 120,(x(j),j=13,idim,idim/20)
  120 format(4d18.9)
c
      call check2(elemnt,loc,idim,idim,ic,x,wk,Hexpec)
      npair(1)=1
      npair(2)=2
      call xcorr(n,idim,npair,1,x,sxx,list1,list2)
      call zcorr(n,idim,npair,1,x,szz,list1)
      print 130,sxx,szz
  130  format(/' [Nearest neighbor correlation functions]'/
     &       '    sxx :',d18.10,',    szz :',d18.10)
      end
