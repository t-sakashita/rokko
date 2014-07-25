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
      subroutine diag(elemnt, ideclr, idim, E, v, ne, nvec, wk, iwk)
c
c    elemnt      @ matrix elements
c    ideclr      @ declared array size in the main program
c    idim        @ matrix dimension
c    E           # eigenvalues
c    v           # eigenvector
c    ne          @ number of eigenvalues to calculate
c    nvec        @ number of eigenvectors to calculate
c    wk,iwk        working areas
c
      implicit none
      real*8 elemnt(ideclr,idim), E(ne), v(ideclr,nvec), wk(ideclr,8)
      integer*4 ideclr, idim, ne, nvec, iwk(ideclr)
      integer*4 i, j, info
c
      call dsyev('V', 'U', idim, elemnt, ideclr, wk(1, 1), wk(1, 2),
     &     7 * ideclr, info)
      do j = 1, ne
         E(j) = wk(j, 1)
      end do
      do j = 1, nvec
         do i = 1, idim
            v(i, j) = elemnt(i, j)
         end do
      end do
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
