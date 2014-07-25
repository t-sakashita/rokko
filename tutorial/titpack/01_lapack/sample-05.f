c>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
c        ACKNOWLEDGE THE USE OF THIS PACKAGE
c         WHEN YOU PUBLISH YOUR RESULTS !!!
c>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
c
c************ Sample main program #5 **************
c  1d Heisenberg antiferromagnet with 16 spins
c Degeneracy check by various initial vectors
c*********************************************
      parameter (n=16,idim=12870,ibond=n)
      parameter (nbond=1)
      implicit real*8 (a-h,o-z)
      dimension E(4)
      dimension list1(idim),list2(2,0:2**15)
      dimension bondwt(ibond),ipair(2*ibond),zrtio(ibond)
      dimension wk(idim,2)
      dimension x(idim),v(idim,2),norm(2)
c
      data bondwt/ibond*-1.0d0/
      data zrtio/ibond*1.0d0/
      data ipair/1,2, 2,3,  3,4,  4,5,  5,6,  6,7,  7,8, 8,9,
     &     9,10, 10,11, 11,12, 12,13, 13,14, 14,15, 15,16, 16,1/
c
      nvec=1
      call sz(n,idim,0.0d0,list1,list2)
c- You may alternatively use szdy or sztn for faster processing -
c      call szdy(n,idim,0.0d0,list1,list2)
c  or
c      call sztn(n,idim,0.0d0,list1,list2)
c------------------------------------------------------------
c
      k=0
c*** Two different initial conditions
      do 10 iv=21,idim,idim/2
        k=k+1
        call lnc1(n,idim,ipair,bondwt,zrtio,ibond,
     &                   nvec,iv,E,itr,wk,idim,list1,list2)
        print 100,k,e
 100    format(/' #',i2,' [Eigenvalues]  '/2x,4f14.8)
        call lncv1(n,idim,ipair,bondwt,zrtio,ibond,
     &              nvec,iv,x,itr,wk,idim,list1,list2)
        do 20 j=1,idim
 20     v(j,k)=x(j)
 10   continue
c
c*** Degeneracy check
      call orthg(idim,idim,v,norm,idgn,2)
      print 110,idgn,norm
 110  format(/' [Degeneracy] :',i4, '   Norm :',2i4)
      end
