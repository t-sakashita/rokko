c>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
c        ACKNOWLEDGE THE USE OF THIS PACKAGE
c         WHEN YOU PUBLISH YOUR RESULTS !!!
c>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
c
c************ Sample main program #3 **************
c  1d Heisenberg antiferromagnet with 16 spins
c    Eigenvectors of excited states by lncv1
c*********************************************
      parameter (n=16,idim=12870,ibond=n)
      parameter (nvec=3,nbond=1)
      implicit real*8 (a-h,o-z)
      dimension E(4)
      dimension list1(idim),list2(2,0:2**15)
      dimension bondwt(ibond),ipair(2*ibond),zrtio(ibond)
      dimension wk(idim,2)
      dimension x(idim,nvec)
c
      data bondwt/ibond*-1.0d0/
      data zrtio/ibond*1.0d0/
      data ipair/1,2, 2,3,  3,4,  4,5,  5,6,  6,7,  7,8, 8,9,
     &     9,10, 10,11, 11,12, 12,13, 13,14, 14,15, 15,16, 16,1/
      iv=idim/3
      call sz(n,idim,0.0d0,list1,list2)
c- You may alternatively use szdy or sztn for faster processing -
c      call szdy(n,idim,0.0d0,list1,list2)
c  or
c      call sztn(n,idim,0.0d0,list1,list2)
c------------------------------------------------------------
c
c*** Eigenvalues
      call lnc1(n,idim,ipair,bondwt,zrtio,ibond,
     &                 nvec,iv,E,itr,wk,idim,list1,list2)
      print 100,e,itr
 100  format(/' [Eigenvalues]  '/2x,4f14.8
     &       /' [Iteration number]'/i8)
c
c*** Eigenvectors
      call lncv1(n,idim,ipair,bondwt,zrtio,ibond,
     &            nvec,iv,x,itr,wk,idim,list1,list2)
      print *,'[Eigenvector components (selected)]'
      print 120,(x(j,nvec),j=13,idim,idim/20)
 120  format(4d18.9)
c
c*** Precision check and correlation functions
      call check1(n,idim,ipair,bondwt,zrtio,ibond,
     &            x(1,nvec),wk,Hexpec,list1,list2)
      end
