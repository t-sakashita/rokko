c***********************************************************
c*                                                         *
c*                    TITPACK Ver. 2                       *
c*                                                         *
c*                    February, 1991                       *
c*                                                         *
c*          Copyright (C) Hidetoshi Nishimori              *
c*                                                         *
c***********************************************************
c
c============== SUBROUTINES / SMALL MATRICES ===============
c
c*** variables marked @ should be given from the main program
c*** variables marked # are evaluated and returned
c
c******************* matrix elements *********************
c
      subroutine elm3(n,idim,ipair,bondwt,zrtio,ibond,
     &                elemnt,ideclr,list1,list2)
c
c    n           @ lattice size
c    idim        @ matrix dimension
c    ipair       @ pairs of sites connected by bonds
c    bondwt      @ exchange interaction of each bond Jxy
c    zrtio       @ ratio of Jz to Jxy
c    ibond       @ number of bonds
c    ideclr      @ declared array size in the main program
c    elemnt      # matrix elements
c    list1,list2 @ spin configurations generated in 'sz'
c
      implicit real*8 (a-h,o-z)
      dimension elemnt(ideclr,idim)
      dimension ipair(ibond*2),bondwt(ibond),zrtio(ibond)
      dimension list1(idim),list2(2,0:2**15)
c
      ihfbit=2**((n+1)/2)
      irght=2**((n+1)/2)-1
      ilft=ieor(2**n-1,irght)
c
      call datack(ipair,ibond,n)
c*** initialization
      do 10 i=1,idim
      do 10 j=1,idim
 10   elemnt(j,i)=0.0d0
c
c*** elements
      do 30 k=1,ibond
          isite1=ipair(k*2-1)-1
          isite2=ipair(k*2  )-1
          is1=2**isite1
          is2=2**isite2
          is=is1+is2
          wght=bondwt(k)
          diag=wght*0.5d0*zrtio(k)
*voption vec
          do 30 i=1,idim
            ibit=iand(list1(i),is)
            if(ibit.eq.0.or.ibit.eq.is)then
               elemnt(i,i)=elemnt(i,i)-diag
              else
               elemnt(i,i)=elemnt(i,i)+diag
               iexchg=ieor(list1(i),is)
               ia=iand(iexchg,irght)
               ib=iand(iexchg,ilft)/ihfbit
               newcfg=list2(1,ia)+list2(2,ib)
               elemnt(i,newcfg)=-wght
            end if
 30   continue
      return
      end
c
c************ eigenvalues of a small matrix *************
c
      subroutine diag(elemnt,ideclr,idim,E,v,ne,nvec,eps,wk,iwk)
c
c    elemnt      @ matrix elements
c    ideclr      @ declared array size in the main program
c    idim        @ matrix dimension
c    E           # eigenvalues
c    v           # eigenvector
c    ne          @ number of eigenvalues to calculate
c    nvec        @ number of eigenvectors to calculate
c    eps         @ limit of error
c    wk,iwk        working areas
c
      implicit real*8(a-h,o-z)
      dimension E(ne),v(ideclr,nvec),elemnt(ideclr,idim)
      dimension wk(ideclr,8),iwk(ideclr)
c
      if(nvec.lt.0.or.nvec.gt.ne)then
           print *,' #(E10)# nvec given to diag out of range'
           stop
      end if
c
      call hshldr(elemnt,ideclr,idim,wk(1,1),wk(1,2),wk(1,3),wk(1,4),
     &            wk(1,5),wk(1,6))
      call bisec(wk(1,1),wk(1,2),idim,E,ne,eps)
      if(nvec.eq.0)return
      call vec3(E,elemnt,ideclr,idim,ne,nvec,wk(1,4),wk(1,5),wk(1,6),
     &          wk(1,7),wk(1,8),iwk,wk(1,1),wk(1,2),wk(1,3),v)
      end
c
c************ Householder tridiagonalization *************
c
      subroutine hshldr(elemnt,ideclr,idim,alpha,beta,c,w,p,q)
c
c    elemnt     @ matrix
c    ideclr     @ declared array size in the main program
c    idim       @ matrix dimension
c    alpha      # diagonal element
c    beta       # subdiagonal element
c    c          # normalization factor for eigenvector calculation
c    w,p,q        working areas
c
      implicit real*8(a-h,o-z)
      dimension elemnt(ideclr,idim),alpha(idim),beta(idim),c(idim)
      dimension w(idim),p(idim),q(idim)
c
      do 10 k=1,idim-2
        s=0.d0
        do 20 i=k+1,idim
 20     s=s+elemnt(i,k)**2
        s=sqrt(s)
        if(elemnt(k+1,k).lt.0.d0)s=-s
c
        alpha(k)=elemnt(k,k)
        beta(k)=-s
        c(k)=0.0d0
        if(s**2.lt.1.d-26)goto 10
c
        c(k)=1.d0/(s**2+elemnt(k+1,k)*s)
        w(k+1)=elemnt(k+1,k)+s
        do 30 i=k+2,idim
 30     w(i)=elemnt(i,k)
        elemnt(k+1,k)=w(k+1)
c
        do 40 i=k+1,idim
          t=0.d0
          do 50 j=k+1,i
 50       t=t+elemnt(i,j)*w(j)
          do 55 j=i+1,idim
 55       t=t+elemnt(j,i)*w(j)
          p(i)=c(k)*t
 40     continue
c
        t=0.d0
        do 60 i=k+1,idim
 60     t=t+p(i)*w(i)
        s=0.5d0*c(k)*t
        do 70 i=k+1,idim
 70     q(i)=p(i)-s*w(i)
c
        do 80 j=k+1,idim
        do 80 i=j,idim
 80     elemnt(i,j)=elemnt(i,j)-w(i)*q(j)-q(i)*w(j)
c
 10   continue
c
      alpha(idim-1)=elemnt(idim-1,idim-1)
      alpha(idim)=elemnt(idim,idim)
      beta(idim-1)=elemnt(idim,idim-1)
      return
      end
c
c**** eigenvector of a tridiagonal matrix by inverse iteration ****
c
      subroutine vec3(E,elemnt,ideclr,idim,ne,nvec,di,bl,bu,bv,cm,lex,
     &                alpha,beta,c,v)
c
c    E        @ eigenvalues
c    elemnt   @ matrix elements
c    ideclr   @ declared array size in the main program
c    idim     @ matrix dimension
c    ne       @ number of eigenvalues
c    nvec     @ number of eigenvectors
c    di - lex   working areas
c    alpha    @ diagonal element in the tridiagonal form
c    beta     @ subdiagonal element in the tridiagonal form
c    c        @ normalization factor given from hshldr
c    v        # eigenvectors
c
      implicit real*8 (a-h,o-z)
      dimension E(ne),elemnt(ideclr,idim)
      dimension di(idim),bl(idim),bu(idim),bv(idim),cm(idim),lex(idim)
      dimension alpha(idim),beta(idim),c(idim),v(ideclr,nvec)
c
      if(nvec.gt.ne)then
          print *,' #(E11)# nvec given to vec3 out of range'
          stop
      end if
      do 10 k=1,nvec
c
        do 100 j=1,idim
          di(j)=E(k)-alpha(j)
          bl(j)=-beta(j)
          bu(j)=-beta(j)
 100    continue
c
c*** LU decomposition
        do 110 j=1,idim-1
         if(abs(di(j)).gt.abs(bl(j)))then
c--- non pivoting
           lex(j)=0
           if(abs(di(j)).lt.1.d-13)di(j)=1.d-13
           cm(j+1)=bl(j)/di(j)
           di(j+1)=di(j+1)-cm(j+1)*bu(j)
           bv(j)=0.d0
         else
c--- pivoting
           lex(j)=1
           cm(j+1)=di(j)/bl(j)
           di(j)=bl(j)
           s=bu(j)
           bu(j)=di(j+1)
           bv(j)=bu(j+1)
           di(j+1)=s-cm(j+1)*bu(j)
           bu(j+1)= -cm(j+1)*bv(j)
         end if
 110    continue
        if(abs(di(idim)).lt.1.d-13)di(idim)=1.d-13
c
c--- initial vector
        do 120 j=1,idim
 120    v(j,k)=1.d0/(float(j)*5.d0)
c
c*** degeneracy check up
        if(k.eq.1)then
           km=k
        else if(abs(E(k)-E(km)).gt.1.d-13)then
           km=k
        else
           do 130 i=km,k-1
             prd=0.d0
             do 140 j=1,idim
 140         prd=prd+v(j,i)*v(j,k)
             do 150 j=1,idim
 150         v(j,k)=v(j,k)-prd*v(j,i)
 130       continue
        end if
c
c*** inverse iteration
        do 160 l=1,k-km+3
          if((l.ne.1).or.(k.ne.km))then
c--- forward substitution
            do 170 j=1,idim-1
              if(lex(j).eq.0)then
                v(j+1,k)=v(j+1,k)-cm(j+1)*v(j,k)
              else
                s=v(j,k)
                v(j,k)=v(j+1,k)
                v(j+1,k)=s-cm(j+1)*v(j,k)
              end if
 170        continue
          end if
c--- backward substitution
          do 180 j=idim,1,-1
            s=v(j,k)
            if(j.le.idim-1)s=s-bu(j)*v(j+1,k)
            if(j.le.idim-2)s=s-bv(j)*v(j+2,k)
            v(j,k)=s/di(j)
 180      continue
c
c*** normalization
          dnorm=0.d0
          do 190 j=1,idim
 190      dnorm=dnorm+v(j,k)**2
          if(dnorm.gt.1.d-13)dnorm=1./sqrt(dnorm)
          do 200 j=1,idim
 200      v(j,k)=v(j,k)*dnorm
 160    continue
c
 10   continue
c
c*** back transformation to the original representation
      do 210 k=1,nvec
c
        do 220 i=idim-2,1,-1
          prd=0.d0
          do 230 j=i+1,idim
 230      prd=prd+elemnt(j,i)*v(j,k)
          s=prd*c(i)
          do 240 j=i+1,idim
 240      v(j,k)=v(j,k)-s*elemnt(j,i)
 220    continue
 210  continue
c
c*** orthogonalization for degenerate case
      km=1
      do 250 k=2,nvec
        if(abs(E(k)-E(km)).ge.1.0d-13)then
           km=k
        else
           do 260 i=km,k-1
             prd=0.d0
             do 270 j=1,idim
 270         prd=prd+v(j,i)*v(j,k)
             do 280 j=1,idim
 280         v(j,k)=v(j,k)-prd*v(j,i)
 260       continue
c
           dnorm=0.0d0
           do 290 j=1,idim
 290       dnorm=dnorm+v(j,k)**2
           s=1.d0/sqrt(dnorm)
           do 300 j=1,idim
 300       v(j,k)=v(j,k)*s
c
        end if
 250  continue
      end
c
c*************** check of the eigenvector and eigenvalue ************
c
      subroutine check3(elemnt,idim,ideclr,x,v,Hexpec)
c
c    elemnt    @ nonzero elements
c    idim      @ matrix dimension
c    ideclr    @ declared array size in the main program
c    x         @ eigenvector to be checked
c    v         # H*x
c    Hexpec    # <x*H*x>
c
      implicit real*8(a-h,o-z)
      dimension elemnt(ideclr,idim)
      dimension x(idim),v(idim)
c
      dnorm=0.d0
      do 5 j=1,idim
 5    dnorm=dnorm+x(j)**2
      if(dnorm.lt.1.d-30)then
         print *,' #(W18)# Null vector given to check3'
         return
       end if
      do 10 i=1,idim
 10   v(i)=0.0d0
      do 20 j=1,idim
      do 20 k=1,idim
      v(j)=v(j)+elemnt(j,k)*x(k)
 20   continue
c
      prd=0.0d0
      do 30 i=1,idim
 30   prd=prd+v(i)*x(i)
      Hexpec=prd
      print *
      print 200
 200  format(' ---------------------------- Information from check3')
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
