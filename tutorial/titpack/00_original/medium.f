c***********************************************************
c*                                                         *
c*                   TITPACK Ver. 2                        *
c*                                                         *
c*                   February, 1991                        *
c*                                                         *
c*          Copyright (C) Hidetoshi Nishimori              *
c*                                                         *
c***********************************************************
c
c===========  SUBROUTINES / MID-SIZE MATRICES ==============
c
c*** variables marked @ should be given from the main program
c*** variables marked # are evaluated and returned
c*** The following variables are common to all routines
c    elemnt    @ nonzero elements
c    loc       @ location of nonzero elements
c    idim      @ matrix dimension
c    ideclr    @ matrix dimension declared in the main program
c    ic        @ max number of nonzero elements in a row
c
c
c********* matrix elements for general bond weights ***********
c
      subroutine elm2(n,idim,ipair,bondwt,zrtio,ibond,
     &                elemnt,loc,ideclr,ic,list1,list2)
c
c    n           @ lattice size
c    idim        @ matrix dimension
c    ipair       @ pairs of sites connected by bonds
c    bondwt      @ exchange interaction of each bond Jxy
c    zrtio       @ ratio of Jz to Jxy
c    ibond       @ number of bonds
c    list1,list2 @ @spin configurations generated in 'sz'
c
      implicit real*8 (a-h,o-z)
      dimension elemnt(ideclr,ic),loc(ideclr,ic)
      dimension list1(idim),list2(2,0:2**15)
      dimension bondwt(ibond),ipair(2*ibond),zrtio(ibond)
c
      ihfbit=2**((n+1)/2)
      irght=2**((n+1)/2)-1
      ilft=ieor(2**n-1,irght)
c
      call datack(ipair,ibond,n)
c*** initialization
      do 10 i=1,ic
      do 10 j=1,idim
 10   elemnt(j,i)=0.0d0
c
c*** diagonal elements
      do 20 k=1,ibond
          isite1=ipair(k*2-1)-1
          isite2=ipair(k*2  )-1
          is1=2**isite1
          is2=2**isite2
          is=is1+is2
          wght=-bondwt(k)*zrtio(k)*0.5d0
          do 20 i=1,idim
            ibit=iand(list1(i),is)
            if(ibit.eq.0.or.ibit.eq.is)then
               elemnt(i,ic)=elemnt(i,ic)+wght
              else
               elemnt(i,ic)=elemnt(i,ic)-wght
            end if
            loc(i,ic)=i
 20   continue
c
c*** off-diagonal elements
      do 30 k=1,ibond
          isite1=ipair(k*2-1)-1
          isite2=ipair(k*2  )-1
          is1=2**isite1
          is2=2**isite2
          is=is1+is2
          wght=-bondwt(k)
          do 30 i=1,idim
            ibit=iand(list1(i),is)
            if(ibit.eq.0.or.ibit.eq.is)then
               elemnt(i,k)=0.0d0
               loc(i,k)=i
             else
               elemnt(i,k)=wght
               iexchg=ieor(list1(i),is)
               ia=iand(iexchg,irght)
               ib=iand(iexchg,ilft)/ihfbit
               loc(i,k)=list2(1,ia)+list2(2,ib)
            end if
 30   continue
      return
      end
c**** eigenvalues by the Lanczos method / matrix element given *****
c           --- dummy routine for simple working area
c
      subroutine lnc2(elemnt,loc,idim,ideclr,ic,nvec,iv,E,itr,wk)
c
c    nvec      @ number of eigenvectors to calculate in lncv2
c    iv        @ location of the nonzero element of the initial vector
c    E         # eigenvalues
c    itr       # number of iterations required for convergence
c    wk          working area
c
      implicit real*8(a-h,o-z)
      dimension elemnt(ideclr,ic),loc(ideclr,ic)
      dimension E(4)
      dimension wk(ideclr,2)
c
      if(iv.le.0.or.iv.gt.idim)then
          print *,' #(E08)# Incorrect iv given to lnc2'
          stop
      end if
      if(nvec.lt.0.or.nvec.gt.4)then
          print *,' #(W12)# Wrong value given to nvec in lnc2'
          print *,'         Only the eigenvalues are calculated'
          nvec=0
      end if
      call lnc2z(elemnt,loc,idim,ideclr,ic,nvec,iv,E,itr,
     &           wk(1,1),wk(1,2))
      end
c
c**** eigenvalues by the Lanczos method / matrix element given
c
      subroutine lnc2z(elemnt,loc,idim,ideclr,ic,nvec,iv,E,itr,v1,v0)
c
      implicit real*8(a-h,o-z)
      dimension elemnt(ideclr,ic),loc(ideclr,ic)
      dimension E(4)
      dimension v0(idim),v1(idim)
      dimension wk(150,5),iwk(150)
      common /vecdat/alpha(150),beta(150),coef(150,5)
c
c*** initialization
      do 10 i=1,idim
          v0(i)=0.0d0
          v1(i)=0.0d0
 10   continue
      v1(iv)=1.0d0
c
c*** alpha(1) and beta(1)
      prdct=0.d0
      do 20 k=1,ic
      do 20 j=1,idim
         temp=elemnt(j,k)*v1(loc(j,k))
         v0(j)=v0(j)+temp
         prdct=prdct+v1(j)*temp
 20   continue
      alpha1=prdct
      alpha(1)=alpha1
      beta1=0.d0
      do 50 i=1,idim
 50   beta1=beta1+(v0(i)-alpha1*v1(i))**2
      beta1=sqrt(beta1)
      beta(1)=beta1
c
c*** iteration
      do 100 i=2,150
          do 110 j=1,idim
            temp1=v1(j)
            temp2=(v0(j)-alpha1*v1(j))/beta1
            v0(j)=-beta1*temp1
            v1(j)=temp2
 110      continue
c
          prdct=0.d0
          do 120 k=1,ic
          do 120 j=1,idim
             temp=elemnt(j,k)*v1(loc(j,k))
             v0(j)=v0(j)+temp
             prdct=prdct+v1(j)*temp
 120      continue
          alpha1=prdct
          alpha(i)=alpha1
          beta1=0.d0
          do 130 j=1,idim
 130      beta1=beta1+(v0(j)-alpha1*v1(j))**2
          beta1=sqrt(beta1)
          beta(i)=beta1
          if(beta(i).lt.0.5d-30)then
            print *,' #(E09)# Tridiagonalization unsuccessful in lnc2'
            print *,'         Beta(i) is too small at i=  ',i
            stop
          end if
c
c*** convergence check
          if(i.gt.20.and.mod(i,5).eq.0)then
             call bisec(alpha,beta,i,E,4,eps)
             if(abs((ebefor-E(2))/E(2)).lt.1.0d-13)then
                if(nvec.gt.0)call vec12(E,i,nvec,wk(1,1),wk(1,2),
     &                            wk(1,3),wk(1,4),wk(1,5),iwk)
                itr=i
                return
             end if
             ebefor=E(2)
          end if
          if(i.eq.20)then
             eps=1.d-10
             call bisec(alpha,beta,20,E,4,eps)
             ebefor=E(2)
          end if
 100  continue
c
      print *,' #(W13)# lnc2 did not converge within 150 steps'
      itr=150
      return
      end
c
c*** eigenvector by the Lanczos method / matrix elements given ***
c        --- dummy routine for simple working area
c
      subroutine lncv2(elemnt,loc,idim,ideclr,ic,nvec,iv,x,itr,wk)
c
c    nvec      @ number of eigenvectors to calculate in lncvec
c    iv        @ location of the nonzero element of the initial vector
c    x         # eigenvectors
c    itr       # number of iterations required for convergence
c    wk          working area
c
      implicit real*8(a-h,o-z)
      dimension elemnt(ideclr,ic),loc(ideclr,ic)
      dimension x(ideclr,nvec)
      dimension wk(ideclr,2)
c
      if(nvec.le.0.or.nvec.gt.4)then
          print *,'#(W14)# nvec given to lncv2 out of range'
          return
      end if
      call lncv2z(elemnt,loc,idim,ideclr,ic,nvec,iv,x,itr,
     &            wk(1,1),wk(1,2))
      end
c
c*** eigenvector by the Lanczos method / matrix elements given
c
      subroutine lncv2z(elemnt,loc,idim,ideclr,ic,nvec,iv,x,itr,
     &                  v1,v0)
c
      implicit real*8(a-h,o-z)
      dimension elemnt(ideclr,ic),loc(ideclr,ic)
      dimension x(ideclr,nvec)
      dimension v0(idim),v1(idim)
      common /vecdat/alpha(150),beta(150),coef(150,5)
c
c*** initialization
      do 10 i=1,idim
          v0(i)=0.0d0
          v1(i)=0.0d0
 10   continue
      do 15 k=1,nvec
      do 15 i=1,idim
 15   x(i,k)=0.d0
c
      v1(iv)=1.0d0
      do 20 k=1,nvec
 20   x(iv,k)=coef(1,k)
c
c*** alpha(1) and beta(1)
      prdct=0.d0
      do 30 k=1,ic
        do 30 j=1,idim
          temp=elemnt(j,k)*v1(loc(j,k))
          v0(j)=v0(j)+temp
          prdct=prdct+v1(j)*temp
 30   continue
      alpha1=alpha(1)
      beta1=beta(1)
      do 40 k=1,nvec
      do 40 j=1,idim
 40   x(j,k)=x(j,k)+coef(2,k)*(v0(j)-alpha1*v1(j))/beta1
c
c*** iteration
      do 100 i=2,itr-1
          do 110 j=1,idim
            temp1=v1(j)
            temp2=(v0(j)-alpha1*v1(j))/beta1
            v0(j)=-beta1*temp1
            v1(j)=temp2
 110      continue
c
          prdct=0.d0
          do 120 k=1,ic
            do 120 j=1,idim
              temp=elemnt(j,k)*v1(loc(j,k))
              v0(j)=v0(j)+temp
              prdct=prdct+v1(j)*temp
 120      continue
          alpha1=alpha(i)
          beta1=beta(i)
          do 130 k=1,nvec
          do 130 j=1,idim
 130      x(j,k)=x(j,k)+coef(i+1,k)*(v0(j)-alpha1*v1(j))/beta1
 100  continue
c
c*** normalization
      do 200 k=1,nvec
          dnorm=0.d0
          do 210 j=1,idim
 210      dnorm=dnorm+x(j,k)**2
          dnorm=sqrt(dnorm)
          do 220 j=1,idim
 220      x(j,k)=x(j,k)/dnorm
 200  continue
      return
      end
c
c*************** check of the eigenvector and eigenvalue ************
c
      subroutine check2(elemnt,loc,idim,ideclr,ic,x,v,Hexpec)
c
c    x         @ eigenvector to be checked
c    v         # H*x
c    Hexpec    # <x*H*x>
c
      implicit real*8(a-h,o-z)
      dimension elemnt(ideclr,ic),loc(ideclr,ic)
      dimension x(idim),v(idim)
c
      dnorm=0.d0
      do 5 j=1,idim
 5    dnorm=dnorm+x(j)**2
      if(dnorm.lt.1.d-30)then
         print *,' #(W15)# Null vector given to check2'
         return
       end if
      do 10 i=1,idim
 10   v(i)=0.0d0
      do 20 k=1,ic
      do 20 j=1,idim
      v(j)=v(j)+elemnt(j,k)*x(loc(j,k))
 20   continue
c
      prd=0.0d0
      do 30 i=1,idim
 30   prd=prd+v(i)*x(i)
      Hexpec=prd
      print *
      print 200
 200  format(' -------------------------- Information from check2')
      print 210,prd
 210  format(' <x*H*x> =',1pd16.8)
      print 220
 220  format(' H*x(j)/x(j) (j=min(idim/3,13),idim,max(1,idim/20))')
      print 230,(v(i)/x(i),i=min(idim,13),idim,max(1,idim/20))
 230  format(4d18.9)
      print 240
 240  format(' --------------------------------------------------')
      return
      end
c
c****************** inverse iteration ************************
c      --- dummy routine for simple working area
c
      subroutine inv2(elemnt,loc,idim,ideclr,ic,iv,Eig,x,wk)
c
c    iv        @ location of nonzero element in the initial vector
c    Eig       @ approximate eigenvalue
c    x         # eigenvector
c    wk         working area
c
      implicit real*8(a-h,o-z)
      dimension elemnt(ideclr,ic),loc(ideclr,ic)
      dimension x(idim),wk(ideclr,4)
c
      call inv2z(elemnt,loc,idim,ideclr,ic,iv,Eig,x,
     &           wk(1,1),wk(1,2),wk(1,3),wk(1,4))
      end
c
c****************** inverse iteration
c
      subroutine inv2z(elemnt,loc,idim,ideclr,ic,iv,Eig,x,b,p,r,y)
c
c    iv        @ location of nonzero element in the initial vector
c    Eig       @ approximate eigenvalue
c    x         # eigenvector
c    b           working area for the rhs of (H-E(approx))*x=b
c    p,r,y       working area used in the routine cg
c
      implicit real*8(a-h,o-z)
      dimension elemnt(ideclr,ic),loc(ideclr,ic)
      dimension x(idim),b(idim),p(idim),r(idim),y(idim)
c
      do 10 i=1,idim
 10   b(i)=0.0d0
      b(iv)=1.0d0
      do 20 itr=1,20
        call cg2(elemnt,loc,idim,ideclr,ic,Eig,x,b,p,r,y,iterat)
        if(iterat.ge.idim.or.iterat.ge.500)then
          xnorm=0.0d0
          do 22 i=1,idim
 22       xnorm=xnorm+x(i)**2
          xnorm=sqrt(xnorm)
          do 24 i=1,idim
 24       x(i)=x(i)/xnorm
          print *,' #(W16)# itr in cg2 exceeds idim or 500'
          print *,'         Approximate eigenvector returned'
          print *,'         Iteration number in inv2 is',itr
          return
        end if
        xnorm=0.0d0
        do 30 i=1,idim
 30     xnorm=xnorm+x(i)**2
        xnorm=sqrt(xnorm)
        do 40 i=1,idim
 40     x(i)=x(i)/xnorm
        xb=0.0d0
        do 50 i=1,idim
 50     xb=xb+x(i)*b(i)
        if(abs(abs(xb)-1.0d0).lt.1.0d-12)then
c         print 100,itr
 100      format('       number of iterations in inv2 :',i5)
          return
        end if
        do 60 i=1,idim
 60     b(i)=x(i)
 20   continue
      print *,' #(W17)# inv2 iteration did not converge'
      return
      end
c
c************** solution of linear equations -- cg method *************
c
      subroutine cg2(elemnt,loc,idim,ideclr,ic,Eig,x,b,p,r,y,itr)
c
c    Eig       @ approximate eigenvalue
c    x         # solution
c    b         @ right hand side of the equation : (H-Eig)*x=b
c    p,r,y       working area
c    itr       # number of iterations to be returned
c
      implicit real*8(a-h,o-z)
      dimension elemnt(ideclr,ic),loc(ideclr,ic)
      dimension b(idim),x(idim),p(idim),r(idim),y(idim)
c
c*** initialization
      bnorm=0.0d0
      do 10 i=1,idim
          bnorm=bnorm+b(i)**2
          r(i)=b(i)
          p(i)=b(i)
          x(i)=0.0d0
 10   continue
c
c*** iteration
      do 20 itr=1,min(500,idim)
        do 30 i=1,idim
 30     y(i)=-Eig*p(i)
        do 40 k=1,ic
        do 40 i=1,idim
 40     y(i)=y(i)+elemnt(i,k)*p(loc(i,k))
        rp=0.0d0
        yp=0.0d0
        do 50 i=1,idim
          rp=rp+r(i)*p(i)
          yp=yp+y(i)*p(i)
 50     continue
        alpha=rp/yp
        rnorm=0.0d0
        do 60 i=1,idim
          x(i)=x(i)+alpha*p(i)
          rnorm=rnorm+r(i)**2
 60     continue
        rnorm2=0.0d0
        do 70 i=1,idim
          r(i)=r(i)-alpha*y(i)
          rnorm2=rnorm2+r(i)**2
 70     continue
        beta=rnorm2/rnorm
        do 90 i=1,idim
 90     p(i)=r(i)+beta*p(i)
        if(mod(itr,5).ne.0)go to 20
        if(sqrt(rnorm2).lt.1.0d-09*sqrt(bnorm))then
c         write(6,150)itr
 150      format('       number of iterations in cg2    :',i5)
          return
        end if
 20   continue
c     print *,' #(Wxx)# cg2 did not converge'
c     print *,'         Approximate solution returned'
      return
      end
