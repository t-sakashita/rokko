c***********************************************************
c*                                                         *
c*                   TITPACK Ver. 2                        *
c*                                                         *
c*                   February, 1991                        *
c*                                                         *
c*         Copyright (C) Hidetoshi Nishimori               *
c*                                                         *
c***********************************************************
c
c============  SUBROUTINES / LARGE MATRICES =================
c
c*** variables marked @ should be given from the main program
c*** variables marked # are evaluated and returned
c*** The following variables are common to all routines
c    n         @ lattice size
c    idim      @ matrix dimension
c    ipair     @ pairs of sites connected by bonds
c    bondwt    @ exchange interaction of each bond Jxy
c    zrtio     @ ratio of Jz to Jxy
c    ibond     @ number of bonds
c
c************ eigenvalues by the Lanczos method **************
c         --- dummy routine for simple working area
c
      subroutine lnc1(n,idim,ipair,bondwt,zrtio,ibond,
     &                nvec,iv,E,itr,wk,ideclr,list1,list2)
c
c    nvec      @ number of eigenvectors to calculate in lncvec
c    iv        @ location of the nonzero element of the initial vector
c    E         # eigenvalues
c    itr     # number of iterations required for convergence
c    wk       working areas
c    list1,list2 @ spin configurations generated in 'sz'
c
      implicit real*8(a-h,o-z)
      dimension E(4)
      dimension ipair(ibond*2),bondwt(ibond),zrtio(ibond)
      dimension list1(idim),list2(2,0:2**15)
      dimension wk(ideclr,2)
c
      if(iv.le.0.or.iv.gt.idim)then
          print *,' #(E06)# Incorrect iv given to lnc1'
          return
      end if
      if(nvec.lt.0.or.nvec.gt.4)then
          print *,' #(W06)# Wrong value given to nvec in lnc1'
          print *,'         Only the eigenvalues are calculated'
          nvec=0
      end if
c
      call lnc1z(n,idim,ipair,bondwt,zrtio,ibond,
     &     nvec,iv,E,itr,wk(1,1),wk(1,2),list1,list2)
      end
c
c************ eigenvalues by the Lanczos method
c
      subroutine lnc1z(n,idim,ipair,bondwt,zrtio,ibond,
     &                nvec,iv,E,itr,v1,v0,list1,list2)
c
      implicit real*8(a-h,o-z)
      dimension E(4)
      dimension ipair(ibond*2),bondwt(ibond),zrtio(ibond)
      dimension list1(idim),list2(2,0:2**15)
      dimension v0(idim),v1(idim)
      dimension wk(150,5),iwk(150)
      common /vecdat/alpha(150),beta(150),coef(150,5)
c
c*** initialization
      do 10 i=1,idim
          v0(i)=0.0d0
          v1(i)=0.0d0
 10   continue
c
      v1(iv)=1.0d0
c
      call datack(ipair,ibond,n)
c
c*** alpha(1) and beta(1)
      call mltply(n,idim,ipair,bondwt,zrtio,ibond,
     &            v1,v0,prdct,list1,list2)
      alpha1=prdct
      alpha(1)=alpha1
      beta1=0.d0
      do 50 i=1,idim
 50   beta1=beta1+(v0(i)-alpha1*v1(i))**2
      beta1=sqrt(beta1)
      beta(1)=beta1
c
C*** iteration
      do 100 i=2,150
          do 110 j=1,idim
            temp1=v1(j)
            temp2=(v0(j)-alpha1*v1(j))/beta1
            v0(j)=-beta1*temp1
            v1(j)=temp2
 110      continue
c
          call mltply(n,idim,ipair,bondwt,zrtio,ibond,
     &                v1,v0,prdct,list1,list2)
          alpha1=prdct
          alpha(i)=alpha1
          beta1=0.d0
          do 120 j=1,idim
 120      beta1=beta1+(v0(j)-alpha1*v1(j))**2
          beta1=sqrt(beta1)
          beta(i)=beta1
          if(beta(i).lt.0.5d-30)then
            print *,' #(E07)# Tridiagonalization unsuccessful in lnc1'
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
      print *,' #(W07)# lnc1 did not converge within 150 steps'
      itr=150
      return
      end
c
c************ eigenvector by the Lanczos method *************
c
      subroutine lncv1(n,idim,ipair,bondwt,zrtio,ibond,
     &                 nvec,iv,x,itr,wk,ideclr,list1,list2)
c
c    nvec      @ number of eigenvectors to be calculated
c    iv        @ location of the nonzero element of the initial vector
c    x         # eigenvector
c    itr     @ number of interations for convergence
c    wk         working area
c    ideclr    @ declared array size in the main program
c    list1,list2 @ spin configurations generated in 'sz'
c
      implicit real*8(a-h,o-z)
      dimension x(ideclr,nvec)
      dimension ipair(ibond*2),bondwt(ibond),zrtio(ibond)
      dimension list1(idim),list2(2,0:2**15)
      dimension wk(ideclr,2)
c
      if(nvec.le.0.or.nvec.gt.4)then
          print *,'#(W08)# nvec given to lncv1 out of range'
          return
      end if
      call lncv1z(n,idim,ipair,bondwt,zrtio,ibond,
     &       nvec,iv,x,ideclr,itr,wk(1,1),wk(1,2),list1,list2)
      end
c
c************ eigenvector by the Lanczos method
c
      subroutine lncv1z(n,idim,ipair,bondwt,zrtio,ibond,
     &                 nvec,iv,x,ideclr,itr,v1,v0,list1,list2)
c
      implicit real*8(a-h,o-z)
      dimension x(ideclr,nvec)
      dimension ipair(ibond*2),bondwt(ibond),zrtio(ibond)
      dimension list1(idim),list2(2,0:2**15)
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
      call mltply(n,idim,ipair,bondwt,zrtio,ibond,
     &            v1,v0,prdct,list1,list2)

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
          call mltply(n,idim,ipair,bondwt,zrtio,ibond,
     &                v1,v0,prdct,list1,list2)
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
c*************** matrix multiplication *****************
c
      subroutine mltply(n,idim,ipair,bondwt,zrtio,ibond,
     &                  v1,v0,prdct,list1,list2)
c
c    v1          @  input vector
c    v0          @# output vector : H*v1+v0(input)
c    prdct       #  <v1*H*v1>
c    list1,list2 @  spin configurations generated in 'sz'
c
      implicit real*8(a-h,o-z)
      dimension ipair(ibond*2),bondwt(ibond),zrtio(ibond)
      dimension list1(idim),list2(2,0:2**15)
      dimension v0(idim),v1(idim)
c
      ihfbit=2**((n+1)/2)
      irght=2**((n+1)/2)-1
      ilft=ieor(2**n-1,irght)
c
      prdct=0.d0
      do 10 k=1,ibond
          isite1=ipair(k*2-1)-1
          isite2=ipair(k*2  )-1
          is1=2**isite1
          is2=2**isite2
          is=is1+is2
          wght=bondwt(k)*zrtio(k)*0.5d0
          do 10 j=1,idim
            ibit=iand(list1(j),is)
            if(ibit.eq.0.or.ibit.eq.is)then
               v0(j)=v0(j)-wght*v1(j)
               factor=1.d0
               offdg=0.d0
            else
               v0(j)=v0(j)+wght*v1(j)
               iexchg=ieor(list1(j),is)
               ia=iand(iexchg,irght)
               ib=iand(iexchg,ilft)/ihfbit
               temp=v1(list2(1,ia)+list2(2,ib))*bondwt(k)
               v0(j)=v0(j)-temp
               factor=-1.d0
               offdg=-temp*v1(j)
            end if
            prdct=prdct-factor*wght*v1(j)**2+offdg
 10   continue
      end
c
c*************** check of the eigenvector and eigenvalue ************
c
      subroutine check1(n,idim,ipair,bondwt,zrtio,ibond,
     &                  x,v,Hexpec,list1,list2)
c
c    x           @ eigenvector to be checked
c    v           # H*x
c    Hexpec      # <x*H*x>
c    list1,list2 @ spin configurations generated in 'sz'
c
      implicit real*8(a-h,o-z)
      dimension x(idim),v(idim)
      dimension list1(idim),list2(2,0:2**15)
      dimension ipair(ibond*2),bondwt(ibond),zrtio(ibond)
c
      ihfbit=2**((n+1)/2)
      irght=2**((n+1)/2)-1
      ilft=ieor(2**n-1,irght)
c
      dnorm=0.d0
      do 5 i=1,idim
 5    dnorm=dnorm+x(i)**2
      if(dnorm.lt.1.d-30)then
         print *,' #(W09)# Null vector given to check1'
         return
       end if
      do 10 i=1,idim
 10   v(i)=0.0d0
      do 20 k=1,ibond
          isite1=ipair(k*2-1)-1
          isite2=ipair(k*2  )-1
          is1=2**isite1
          is2=2**isite2
          is=is1+is2
          wght=bondwt(k)
          diag=wght*0.5d0*zrtio(k)
          do 20 i=1,idim
            ibit=iand(list1(i),is)
             if(ibit.eq.0.or.ibit.eq.is)then
                v(i)=v(i)-diag*x(i)
               else
                v(i)=v(i)+diag*x(i)
                iexchg=ieor(list1(i),is)
                ia=iand(iexchg,irght)
                ib=iand(iexchg,ilft)/ihfbit
                v(i)=v(i)-
     &               x(list2(1,ia)+list2(2,ib))*wght
             end if
 20   continue
c
      prd=0.0d0
      do 30 i=1,idim
 30   prd=prd+v(i)*x(i)
      Hexpec=prd
      print *
      print 200
 200  format(' ---------------------------- Information from check1')
      print 210,prd
 210  format(' <x*H*x> =',1pd16.8)
      print 220
 220  format(' H*x(j)/x(j) (j=min(idim/3,13),idim,max(1,idim/20))')
      print 230,(v(i)/x(i),i=min(idim/3,13),idim,max(1,idim/20))
 230  format(4d18.9)
      print 240
 240  format(' ---------------------------------------------------')
      return
      end
c
c****************** inverse iteration ************************
c      --- dummy routine for simple working area
c
      subroutine inv1(n,idim,ipair,bondwt,zrtio,ibond,
     &                  Eig,iv,x,wk,ideclr,list1,list2)
c
c    Eig         @ eigenvalue
c    iv          @ location of the nonzero element of the initial vector
c    x           # eigenvector
c    wk            working area
c    ideclr      @ declared array size in the main program
c    list1,list2 @ spin configurations generated in 'sz'
c
      implicit real*8(a-h,o-z)
      dimension ipair(ibond*2),bondwt(ibond),zrtio(ibond)
      dimension list1(idim),list2(2,0:2**15)
      dimension x(idim)
      dimension wk(ideclr,4)
c
      call inv1z(n,idim,ipair,bondwt,zrtio,ibond,
     & Eig,iv,x,wk(1,1),wk(1,2),wk(1,3),wk(1,4),list1,list2)
      end
c
c****************** inverse iteration
c
      subroutine inv1z(n,idim,ipair,bondwt,zrtio,ibond,
     &                  Eig,iv,x,b,p,r,y,list1,list2)
c
c    b            working area for the rhs of (H-E(approx))*x=b
c    p,r,y        working area used in the routine cg1
c


      implicit real*8(a-h,o-z)
      dimension ipair(ibond*2),bondwt(ibond),zrtio(ibond)
      dimension list1(idim),list2(2,0:2**15)
      dimension b(idim),x(idim),r(idim),y(idim),p(idim)
c
      do 10 i=1,idim
 10   b(i)=0.0d0
      b(iv)=1.0d0
      do 20 itr=1,20
        call cg1(n,idim,ipair,bondwt,zrtio,ibond,
     &          Eig,x,b,p,r,y,iterat,list1,list2)
        if(iterat.gt.idim)then
          xnorm=0.0d0
          do 22 i=1,idim
 22       xnorm=xnorm+x(i)**2
          xnorm=sqrt(xnorm)
          do 24 i=1,idim
 24       x(i)=x(i)/xnorm
          print *,' #(W10)# Iterat in cg1 exceeds idim or 500'
          print *,'         Approximate eigenvector returned'
          print *,'         Itration number in inv1 is',itr
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
c          print 100,itr
 100      format('       number of iterations in inv1 :',i5)
          return
        end if
        do 60 i=1,idim
 60     b(i)=x(i)
 20   continue
      print *,' #(W11)# inv1 did not converge'
      return
      end
c
c************** solution of linear equations -- cg method ************
c
      subroutine cg1(n,idim,ipair,bondwt,zrtio,ibond,
     &              Eig,x,b,p,r,y,itr,list1,list2)
c
c    Eig         @ eigenvalue
c    x           # eigenvector
c    b             working area for the rhs of (H-E(approx))*x=b
c    p,r,y         working area used in the routine cg
c    itr         # number of iterations required for convergence
c    list1,list2 @ spin configurations generated in 'sz'
c
      implicit real*8(a-h,o-z)
      dimension ipair(ibond*2),bondwt(ibond),zrtio(ibond)
      dimension list1(idim),list2(2,0:2**15)
      dimension b(idim),x(idim),r(idim),y(idim),p(idim)
c
      ihfbit=2**((n+1)/2)
      irght=2**((n+1)/2)-1
      ilft=ieor(2**n-1,irght)
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
 30     y(i)=0.0d0
        do 40 k=1,ibond
          isite1=ipair(k*2-1)-1
          isite2=ipair(k*2  )-1
          is1=2**isite1
          is2=2**isite2
          is=is1+is2
          eperbd=Eig/float(ibond)
          wght=bondwt(k)
          diag1= wght*0.5d0*zrtio(k)+eperbd
          diag2=-wght*0.5d0*zrtio(k)+eperbd
          do 40 i=1,idim
            ibit=iand(list1(i),is)
            if(ibit.eq.0.or.ibit.eq.is)then
                y(i)=y(i)-diag1*p(i)
            else
                iexchg=ieor(list1(i),is)
                ia=iand(iexchg,irght)
                ib=iand(iexchg,ilft)/ihfbit
                y(i)=y(i)-diag2*p(i)-p(list2(1,ia)+list2(2,ib))*wght
            end if
 40     continue
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
        if(sqrt(rnorm2).lt.1.0d-9*sqrt(bnorm))then
c          print 150,itr
 150      format('       number of iterations in cg1     :',i5)
          return
        end if
 20   continue
c     print *,' #(Wxx)# cg1 did not converge'
      return
      end
