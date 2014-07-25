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
c================== COMMON SUBROUTINES ====================
c
c*** variables marked @ should be given from the main program
c*** variables marked # are evaluated and returned
c
c******** configurations with the specified sz ***********
c
      subroutine sz(n,idim,szval,list1,list2)
c
c    n          @  lattice size
c    idim       @  dimension of the matrix
c    szval      @  total sz
c    list1(i)   #  i-th spin configuration
c    list2      #  inverse list of list1 expressed by the
c                  2-dim search method of M. Ogata and H.Q. Lin.
c
      real*8 szval
      dimension list1(idim),list2(2,0:2**15)
c
      if(szval.lt.-1.0d-13.or.szval.gt.n/2.d0-1.d0+1.d-13)then
        print *,' #(E01)# Variable szval given to sz out of range'
        stop
      end if
      if(idim.lt.3)then
        print *,' #(E02)# Incorrect idim or n given to sz'
        stop
      end if
c
c* initialization
      ihfbit=2**((n+1)/2)
      irght=2**((n+1)/2)-1
      ilft=ieor(2**n-1,irght)
      iupspn=n/2+mod(n,2)+int(szval+0.001d0)
      icnt=0
      ja=0
      jb=0
      ibpatn=0
c
c* main loop
      do 10 i=1,2**n
        isz=0
        do 20 j=0,n-1
 20     isz=isz+mod(i/2**j,2)
        if(isz.ne.iupspn)go to 10
        icnt=icnt+1
        if(icnt.gt.idim)then
          print *,' #(E02)# Incorrect idim or n given to sz'
          stop
        end if
        ia=iand(i,irght)
        ib=iand(i,ilft)/ihfbit
        if(ib.eq.ibpatn)then
          ja=ja+1
        else
          ibpatn=ib
          ja=1
          jb=icnt-1
        end if
        list1(icnt)=i
        list2(1,ia)=ja
        list2(2,ib)=jb
 10   continue
      if(icnt.eq.idim)return
      print *,' #(E02)# Incorrect idim or n given to sz'
      stop
      end
c
c************* data check of pairs of sites ************
c
      subroutine datack(ipair,ibond,n)
      dimension ipair(ibond*2)
c
      do 10 k=1,ibond
         isite1=ipair(k*2-1)
         isite2=ipair(k*2  )
         if(isite1.le.0.or.isite2.le.0.or.
     &      isite1.gt.n.or.isite2.gt.n)then
            print *,' #(E03)# Incorrect data in ipair'
            print *,'         Location :  ',k*2-1,k*2
            stop
         end if
 10    continue
       end
c
c********* eigenvalues by the bisection method **********
c
      subroutine bisec(alpha,beta,ndim,E,ne,eps)
c
c    alpha  @ diagonal element
c    beta   @ subdiagonal element
c    ndim   @ matrix dimension
c    E      # eigenvalues
c    ne     @ number of eigenvalues to calculate
c    eps    @ limit of error

      implicit real*8 (a-h,o-z)
      dimension alpha(ndim),beta(ndim),E(ne),b2(2000)
c
      if(ndim.gt.2000)then
          print *,' #(E04)# ndim given to bisec exceeds 2000'
          stop
      end if
      if(ne.gt.ndim.or.ne.le.0)then
          print *,' #(E05)# ne given to bisec out of range'
          stop
      end if
c
c*** initial bound
      range=abs(alpha(1))+abs(beta(1))
      do 10 k=2,ndim-1
 10   range=max(range,abs(beta(k-1))+abs(alpha(k))+abs(beta(k)))
      range=max(range,abs(beta(ndim-1))+abs(alpha(ndim)))
      range=-range
c
      b2(1)=0.d0
      do 20 i=2,ndim
 20   b2(i)=beta(i-1)**2
c
      epsabs=abs(range)*eps
      do 30 i=1,ne
 30   E(i)=-range
      b=range
c
c*** bisection method
      do 100 k=1,ne
        a=E(k)
        do 110 j=1,100
          c=(a+b)/2.d0
          if(abs(a-b).lt.epsabs)goto 100
          numneg=0
          g=1.d0
          ipass=0
          do 120 i=1,ndim
            if(ipass.eq.0)then
              g=c-alpha(i)-b2(i)/g
              else if(ipass.eq.1)then
                ipass=2
              else
                g=c-alpha(i)
                ipass=0
            end if
c
            if(ipass.eq.0)then
              if(g.le.0.d0)numneg=numneg+1
              if(abs(g).le.abs(b2(i)*epsabs*eps))ipass=1
            end if
 120      continue
          numneg=ndim-numneg
          if(numneg.lt.k)then
            b=c
          else
            a=c
            do 130 i=k,min(numneg,ne)
 130        E(i)=c
          end if
 110    continue
 100   continue
       end
c
c*** eigenvector of a tridiagonal matrix by inverse iteration ***
c                 for the large/medium routines
c
      subroutine vec12(E,ndim,nvec,di,bl,bu,bv,cm,lex)
c
c    E(4)       @  4 lowest eigenvalues
c    ndim       @  matrix dimension
c    nvec       @  number of vectors to calculate
c    di - lex      working areas
c
      implicit real*8 (a-h,o-z)
      dimension E(4)
      dimension di(ndim),bl(ndim),bu(ndim),bv(ndim),cm(ndim),lex(ndim)
      common /vecdat/alpha(150),beta(150),v(150,5)
c
      do 10 k=1,nvec
c
        do 100 j=1,ndim
          di(j)=E(k)-alpha(j)
          bl(j)=-beta(j)
          bu(j)=-beta(j)
 100    continue
c
c*** LU decomposition
        do 110 j=1,ndim-1
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
        if(abs(di(ndim)).lt.1.d-13)di(ndim)=1.d-13
c
c--- initial vector
        do 120 j=1,ndim
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
             do 140 j=1,ndim
 140         prd=prd+v(j,i)*v(j,k)
             do 150 j=1,ndim
 150         v(j,k)=v(j,k)-prd*v(j,i)
 130       continue
        end if
c
c*** inverse iteration
        do 160 l=1,k-km+3
          if((l.ne.1).or.(k.ne.km))then
c--- forward substitution
            do 170 j=1,ndim-1
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
          do 180 j=ndim,1,-1
            s=v(j,k)
            if(j.le.ndim-1)s=s-bu(j)*v(j+1,k)
            if(j.le.ndim-2)s=s-bv(j)*v(j+2,k)
            v(j,k)=s/di(j)
 180      continue
c
c*** normalization
          dnorm=0.d0
          do 190 j=1,ndim
 190      dnorm=dnorm+v(j,k)**2
          if(dnorm.gt.1.d-13)dnorm=1./sqrt(dnorm)
          do 200 j=1,ndim
 200      v(j,k)=v(j,k)*dnorm
 160    continue
c
 10   continue
      end
c
c************* xx correlation function **************
c
      subroutine xcorr(n,idim,npair,nbond,x,sxx,list1,list2)
c
c    n           @ lattice size
c    idim        @ matrix dimension
c    npair       @ pair of sites (k,l) <Sx(k)Sx(l)>
c    nbond       @ number of bonds to be investigated
c    x           @ eigenvetor
c    sxx         # xx correlation function
c    list1,list2 @ spin configurations generated in 'sz'
c
      implicit real*8 (a-h,o-z)
      dimension list1(idim),list2(2,0:2**15)
      dimension x(idim)
      dimension npair(nbond*2),sxx(nbond)
c
      ihfbit=2**((n+1)/2)
      irght=2**((n+1)/2)-1
      ilft=ieor(2**n-1,irght)
c
      do 10 k=1,nbond
        i1=npair(k*2-1)-1
        i2=npair(k*2  )-1
        if(i1.lt.0.or.i1.ge.n.or.i2.lt.0.or.i2.ge.n.or.i1.eq.i2)then
          print *,' #(W01)# Wrong site number given to xcorr'
          return
        end if
        corr=0.d0
        is=2**i1+2**i2
        do 20 j=1,idim
          ibit=iand(list1(j),is)
          if(ibit.eq.0.or.ibit.eq.is)goto 20
          iexchg=ieor(list1(j),is)
          ia=iand(iexchg,irght)
          ib=iand(iexchg,ilft)/ihfbit
          corr=corr+x(j)*x(list2(1,ia)+list2(2,ib))
 20     continue
        sxx(k)=corr/4.d0
 10   continue
      return
      end
c
c************* zz correlation function **************
c
      subroutine zcorr(n,idim,npair,nbond,x,szz,list1)
c
c    n           @ lattice size
c    idim        @ matrix dimension
c    npair       @ pair of sites (k,l) <Sz(k)Sz(l)>
c    nbond       @ number of bonds to be investigated
c    x           @ eigenvetor
c    szz         # zz correlation function
c    list1,list2 @ spin configurations generated in 'sz'
c
      implicit real*8 (a-h,o-z)
      dimension list1(idim)
      dimension x(idim)
      dimension npair(nbond*2),szz(nbond)
c
      do 10 k=1,nbond
        i1=npair(k*2-1)-1
        i2=npair(k*2  )-1
        if(i1.lt.0.or.i1.ge.n.or.i2.lt.0.or.i2.ge.n.or.i1.eq.i2)then
          print *,' #(W02)# Wrong site number given to zcorr'
          return
        end if
        corr=0.d0
        is=2**i1+2**i2
        do 20 j=1,idim
          ibit=iand(list1(j),is)
          if(ibit.eq.0.or.ibit.eq.is)then
            factor=1.d0
          else
            factor=-1.d0
          end if
          corr=corr+factor*x(j)**2
 20     continue
        szz(k)=corr/4.d0
 10   continue
      return
      end
c
c******* Orthogonalization of the eigenvectors ************
c
      subroutine orthg(idim,ideclr,ev,norm,idgn,numvec)
c
c   idim    @  matrix dimension
c   ideclr  @  declared array size in the main program
c   ev      @# vectors to be orthogonalized / orthogonalized vectors
c   norm(j) #  norm of the j-th vector returned
c   idgn    #  degree of degenearcy
c   numvec  @  number of vectors to be checked
c
      implicit real*8(a-h,o-z)
      dimension ev(ideclr,numvec),norm(numvec)
      if(numvec.le.1)then
         print *,' #(W03)# Number of vectors is less than 2 in orthg'
         return
      end if
      do 10 i=1,numvec
        dnorm=0.0d0
        do 20 j=1,idim
 20     dnorm=dnorm+ev(j,i)**2
        if(dnorm.lt.1.0d-20)then
           print *,' #(W04)# Null vector given to orthg. Location is',i
           return
        end if
        dnorm=1.0d0/sqrt(dnorm)
        do 25 j=1,idim
 25     ev(j,i)=ev(j,i)*dnorm
 10   continue
      idgn=numvec
      norm(1)=1
c*** orthogonalization
      do 30 i=2,numvec
       norm(i)=1
       do 40 j=1,i-1
         prjct=0.0d0
         do 50 l=1,idim
 50      prjct=prjct+ev(l,i)*ev(l,j)
         do 60 l=1,idim
 60      ev(l,i)=ev(l,i)-prjct*ev(l,j)
 40    continue
       vnorm=0.0d0
       do 70 l=1,idim
 70    vnorm=vnorm+ev(l,i)**2
       if(vnorm.gt.1.0d-15)then
         vnorm=1.0d0/sqrt(vnorm)
         do 80 l=1,idim
 80      ev(l,i)=ev(l,i)*vnorm
        else
         do 90 l=1,idim
 90      ev(l,i)=0.0d0
         idgn=idgn-1
         norm(i)=0
       end if
 30   continue
c*** check orthogonality
      do 100 i=2,numvec
       do 100 j=1,i-1
       prd=0.0d0
       do 110 l=1,idim
 110   prd=prd+ev(l,i)*ev(l,j)
       if(abs(prd).lt.1.0d-10)go to 100
         print 200,i,j
 200     format(' #(W05)# Non-orthogonal vectors at',2i4)
         print 210,prd
 210     format('         Overlap : ',d14.7)
         print *,'         Unsuccessful orthogonalization'
         return
 100  continue
      return
      end
c
c******** configurations with the specified sz ***********
      subroutine szdy(n,idim,szval,list1,list2)
c
c    n          @  lattice size
c    idim       @  dimension of the matrix
c    szval      @  total sz
c    list1(i)   #  i-th spin configuration
c    list2      #  inverse list of list1 expressed by the
c                  2-dim search method of M. Ogata and H.Q. Lin.
c
c==============================================================
c     This routine is equivalent to sz but is faster than sz.
c     This routine has been developed by Daijiro Yoshioka,
c     University of Tokyo.  The copyright of szdy belongs to him.
c                                             1993/5/10
c==============================================================
      real*8 szval
      dimension list1(idim),list2(2,0:2**15),j1(33)
c
      if(szval.lt.-1.0d-13.or.szval.gt.n/2.d0-1.d0+1.d-13)then
        print *,' #(E01)# Variable szval given to sz out of range'
        stop
      end if
      if(idim.lt.3)then
        print *,' #(E02)# Incorrect idim or n given to sz'
        stop
      end if
c
c* initialization
      ihfbit=2**((n+1)/2)
      irght=2**((n+1)/2)-1
      ilft=ieor(2**n-1,irght)
      iupspn=n/2+mod(n,2)+int(szval+0.001d0)
      icnt=0
      ja=0
      jb=0
      ibpatn=0
c
c* main loop
      icnt=0
      do 1 n1=1,iupspn
 1    j1(n1)=n1
      j1(iupspn+1)=n+1
 2    continue
      icnt=icnt+1
        if(icnt.gt.idim)then
          print *,' #(E02)# Incorrect idim or n given to sz'
          stop
        end if
      i=0
      do 4 k=1,iupspn
 4    i=ibset(i,j1(k)-1)
        ia=iand(i,irght)
        ib=iand(i,ilft)/ihfbit
        if(ib.eq.ibpatn)then
          ja=ja+1
        else
          ibpatn=ib
          ja=1
          jb=icnt-1
        end if
        list1(icnt)=i
        list2(1,ia)=ja
        list2(2,ib)=jb
      do 5 k=1,iupspn
 5    if(j1(k+1) .gt. j1(k)+1) go to 6
      if(icnt.eq.idim) return
      print *,' #(E02)# Incorrect idim or n given to sz'
      stop
 6    j1(k)=j1(k)+1
      if(k .eq. 1) go to 2
      do 7 k1=1,k-1
 7    j1(k1)=k1
      go to 2
      end
c******** configurations with the specified sz ***********
c
      subroutine sztn(n,idim,szval,list1,list2)
c
c    n          @  lattice size
c    idim       @  dimension of the matrix
c    szval      @  total sz
c    list1(i)   #  i-th spin configuration
c    list2      #  inverse list of list1 expressed by the
c                  2-dim search method of M. Ogata and H.Q. Lin.
c
c==============================================================
c     This routine is equivalent to sz but is faster than sz
c     and szdy (by D. Yoshioka).
c     This routine has been developed by Tota Nakamura,
c     Univesity of Tsukuba. The copyright of sztn belongs to him.
c                                             1994/02/02
c==============================================================
      real*8 szval
      parameter (ibin18=2**18)
      dimension list1(idim),list2(2,0:2**15)
      dimension listr1(12870,0:16),listl1(2**16),icomb(0:32),j1(33)
c
      iupspn=n/2+mod(n,2)+int(szval+0.001d0)
c---check idim and szval
      if(szval.lt.-1.0d-13.or.szval.gt.n/2.d0-1.d0+1.d-13)then
       print *,' #(E01)# Variable szval given to sz out of range'
       stop
      end if
      if(idim.lt.3)then
       print *,' #(E02)# Incorrect idim or n given to sz'
       stop
      end if
      dummy=1.0d0
       do 30 j=n,n-iupspn+1,-1
 30    dummy=dummy*j/(n-j+1)
       idum=dummy+0.2d0
      if(idim.ne.idum)then
       print *,' #(E02)# Incorrect idim or n given to sz'
       stop
      end if
c* initialization
      nr=(n+1)/2
      nl=n-nr
      do 500 i=0,2**15
       list2(1,i)=0
       list2(2,i)=0
 500  continue
      icomb(0)=1
      do 10 i=1,min(iupspn,nr)
       dummy=1.0d0
       do 20 j=nr,nr-i+1,-1
 20    dummy=dummy*j/(nr-j+1)
       icomb(i)=dummy+0.2d0
 10   continue
c---Generate listr1
      listr1(1,0)=0
      list2(1,0)=1
      do 100 i=iupspn-min(iupspn,nl),min(iupspn,nr)
       if(icomb(i).ge.2)then
c---The following are taken from 'szdy' programmed by Prof. D. Yoshioka.
        idimdy=icomb(i)
        iupdy=i
        icntdy=0
        do 1 n1=1,iupdy
 1       j1(n1)=n1
        j1(iupdy+1)=nr+1
 2      continue
        icntdy=icntdy+1
        if(icntdy.gt.idimdy)then
         print *,' #(E02)# Incorrect idim or n given to sz !!'
         stop
        end if
        idy=0
        do 4 k=1,iupdy
 4      idy=ibset(idy,j1(k)-1)
        listr1(icntdy,i)=idy
        do 5 k=1,iupdy
 5      if(j1(k+1) .gt. j1(k)+1) go to 6
        if(icntdy.eq.idimdy) goto 9
        print *,' #(E02)# Incorrect idim or n given to sz'
        stop
 6      j1(k)=j1(k)+1
        if(k .eq. 1) go to 2
        do 7 k1=1,k-1
 7      j1(k1)=k1
        go to 2
 9      continue
       else
        listr1(1,i)=2**i-1
       end if
 100  continue
c---Generate list2(1,)
      do 110 i=iupspn-min(iupspn,nl),min(iupspn,nr)
      do 110 j=1,icomb(i)
 110   list2(1,listr1(j,i))=j
c---Generate listl1,list2(2,)
      icnt=0
      lold=0
      do 200 i=0,2**nl-1
       isz=0
       do 210 j=0,nl-1
 210   if(btest(i,j))isz=isz+1
       if(isz.le.min(iupspn,nl).and.isz.ge.(iupspn-min(iupspn,nr)))then
        icnt=icnt+1
        listl1(icnt)=i+ishft(isz,18)
        list2(2,i)=lold
        lold=lold+icomb(iupspn-isz)
       end if
 200  continue
c---Generate list1
      do 300 ib=1,icnt
       isb=iand(listl1(ib),ibin18-1)
       isz=ishft(listl1(ib),-18)
       j=list2(2,isb)
*voption indep(list1)
      do 300 ia=1,icomb(iupspn-isz)
       isa=listr1(ia,iupspn-isz)
       list1(list2(1,isa)+j)=ishft(isb,nr)+isa
 300  continue
      end
