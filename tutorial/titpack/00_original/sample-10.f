c>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
c        ACKNOWLEDGE THE USE OF THIS PACKAGE
c         WHEN YOU PUBLISH YOUR RESULTS !!!
c>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
c
c************* Sample main program #10 *****************
c  1d Heisenberg antiferromagnet with 15 spins
c Degeneracy check by various initial vectors
c***********************************************
c
      parameter (n=15,idim=6435,ibond=n,ic=ibond+1)
      implicit real*8 (a-h,o-z)
      dimension E(4)
      dimension elemnt(idim,ic),loc(idim,ic)
      dimension list1(idim),list2(2,0:2**15)
      dimension bondwt(ibond),ipair(2*ibond),zrtio(ibond)
      dimension wk(idim,4)
      dimension x(idim)
      dimension v(idim,3),norm(3)
c
      data bondwt/ibond*-1.0d0/
      data zrtio/ibond*1.0d0/
      data ipair/1,2, 2,3,  3,4,  4,5,  5,6,  6,7,  7,8, 8,9,
     &     9,10, 10,11, 11,12, 12,13, 13,14, 14,15, 15,1/
      nvec=1
      iv=idim/5
c
      call sz(n,idim,0.5d0,list1,list2)
c- You may alternatively use szdy or sztn for faster processing -
c      call szdy(n,idim,0.5d0,list1,list2)
c  or
c      call sztn(n,idim,0.5d0,list1,list2)
c------------------------------------------------------------
      call elm2(n,idim,ipair,bondwt,zrtio,ibond,
     &          elemnt,loc,idim,ic,list1,list2)
      k=0
c *** Two different initial conditions
      do 10 iv=91,idim,idim/3
        k=k+1
        call lnc2(elemnt,loc,idim,idim,ic,nvec,iv,E,itr,wk)
        print 100,k,e
  100    format(/' #',i2,' [Eigenvalues]  '/2x,4f14.8)
        call inv2(elemnt,loc,idim,idim,ic,iv,E(1),x,wk)
        print *,'[Eigenvector components (selected)]'
        print 120,(x(j),j=13,idim,idim/20)
  120    format(4d18.9)
        call check2(elemnt,loc,idim,idim,ic,x,wk,Hexpec)
        do 20 j=1,idim
  20     v(j,k)=x(j)
  10   continue
c *** Degeneracy check
      call orthg(idim,idim,v,norm,idgn,3)
      print 110,idgn
  110  format(/' [Degeneracy]'/i6)
      end
