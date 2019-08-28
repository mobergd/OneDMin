      subroutine pot(symb,x,y,z,v,dvdx,dvdy,dvdz,natom,maxatom)

c NEW "AUTOFIT" PESs

      implicit real*8(a-h,o-z)
      dimension x(maxatom),y(maxatom),z(maxatom)
      dimension x2(maxatom),y2(maxatom),z2(maxatom)
      dimension x3(maxatom),y3(maxatom),z3(maxatom)
      dimension dvdx(maxatom),dvdy(maxatom),dvdz(maxatom)
      dimension dvdx2(maxatom),dvdy2(maxatom),dvdz2(maxatom)
      dimension dvdx3(maxatom),dvdy3(maxatom),dvdz3(maxatom)
      parameter(autoev=27.2113961d0)
      parameter(autocmi=219474.63067d0)
      parameter(autoang=0.529177249d0)
      character*2 symb(maxatom),symb2(maxatom),symb3(maxatom)
      integer at(maxatom),at2(maxatom),at3(maxatom)

      dimension tmpprint(50)
      common/tmp/tmpprint


      v=0.d0
      v2=0.d0
      v3=0.d0
      do i=1,natom
      symb2(i)="xx"
      symb3(i)="xx"
      x2(i)=0.d0
      y2(i)=0.d0
      z2(i)=0.d0
      x3(i)=0.d0
      y3(i)=0.d0
      z3(i)=0.d0
      dvdz2(i)=0.d0
      dvdy2(i)=0.d0
      dvdx2(i)=0.d0
      dvdz3(i)=0.d0
      dvdy3(i)=0.d0
      dvdx3(i)=0.d0
      dvdz(i)=0.d0
      dvdy(i)=0.d0
      dvdx(i)=0.d0
      enddo

      do i=1,natom
      at(i)=0
      if ((symb(i).eq."H").or.
     &    (symb(i).eq."h"))  at(i)=1       ! hydrogen
      if ((symb(i).eq."C").or.
     &    (symb(i).eq."c"))  at(i)=2       ! carbon
      if ((symb(i).eq."N").or.
     &    (symb(i).eq."n"))  at(i)=2       ! carbon
      if ((symb(i).eq."O").or.
     &    (symb(i).eq."o"))  at(i)=4       ! oxygen
      if ((symb(i).eq."HO").or.
     &    (symb(i).eq."Ho"))  at(i)=5       ! OH H
      if ((symb(i).eq."He").or.
     &    (symb(i).eq."he").or.
     &    (symb(i).eq."HE")) at(i)=21      ! helium
      if ((symb(i).eq."Ar").or.
     &    (symb(i).eq."ar").or.
     &    (symb(i).eq."AR")) at(i)=23      ! argon
      if ((symb(i).eq."N2").or.
     &    (symb(i).eq."n2")) at(i)=26      ! N2  ! label your bath atoms N2, not N
      if (at(i).eq.0) then ! atom not found
           write(6,*)"Atom # ",i," (",symb(i),") not found"
           stop
      endif
      enddo

      natom2=0
      natom3=0
      do i=1,natom
        if (at(i).le.20) then ! collect target atoms
          natom2=natom2+1
          x2(natom2)=x(i)
          y2(natom2)=y(i)
          z2(natom2)=z(i)
          symb2(natom2)=symb(i)
          at2(natom2)=at(i)
        else ! collect bath atoms
          natom3=natom3+1
          x3(natom3)=x(i)
          y3(natom3)=y(i)
          z3(natom3)=z(i)
          symb3(natom3)=symb(i)
          at3(natom3)=at(i)
        endif
      enddo
c      print *,natom,natom2,natom3
c      if (natom3.ne.0) 
c     &   call bath(at3,x3,y3,z3,v3,dvdx3,dvdy3,dvdz3,natom3,maxatom)
c      if (natom2.ne.0) 
c     &   call tinkerpot(symb2,x2,y2,z2,v2,dvdx2,dvdy2,dvdz2,
c     &          natom2,maxatom)
c      if (natom3.ne.0) 
c     &   call rgexp(at,x,y,z,v,dvdx,dvdy,dvdz,natom,maxatom)
      if (natom3.ne.0) 
     &   call lsfit(at,x,y,z,v,dvdx,dvdy,dvdz,natom,maxatom)

      v=v+v2+v3

      tmpprint(1)=v    ! interaction

      natom2=0
      natom3=0
      do i=1,natom
        if (at(i).le.20) then
        natom2=natom2+1
c        print *,i,dvdx2(i),dvdy2(i),dvdz2(i)
        dvdx(i)=dvdx(i)+dvdx2(natom2)
        dvdy(i)=dvdy(i)+dvdy2(natom2)
        dvdz(i)=dvdz(i)+dvdz2(natom2)
        else
        natom3=natom3+1
        dvdx(i)=dvdx(i)+dvdx3(natom3)
        dvdy(i)=dvdy(i)+dvdy3(natom3)
        dvdz(i)=dvdz(i)+dvdz3(natom3)
        endif
      enddo

      return

      end


      subroutine bath(at,x,y,z,v,dvdx,dvdy,dvdz,natom,maxatom)

      implicit real*8(a-h,o-z)
      parameter(autoev=27.2113961d0)
      parameter(autocmi=219474.63067d0)
      dimension x(maxatom),y(maxatom),z(maxatom)
      dimension dvdx(maxatom),dvdy(maxatom),dvdz(maxatom)
      integer at(maxatom)

      v=0.
      do i=1,natom
        dvdx(i)=0.d0
        dvdy(i)=0.d0
        dvdz(i)=0.d0
      enddo

      if (natom.eq.1) return

      if (natom.gt.2) then
         print *,"Can't handle more than 2 bath atoms"
         stop
      endif

      if (natom.eq.2) then
      dx=x(1)-x(2)
      dy=y(1)-y(2)
      dz=z(1)-z(2)
      rr=dsqrt(dx*dx+dy*dy+dz*dz)  ! au

        if (at(1).eq.25.and.at(2).eq.25) then
!       H2 bath
!       From Hack's fit (eq 8 in Hack, Truhlar, JCP 110, 4315 (1999))
!                        to Kolos and Wolniewicz JCP 43, 2429 (1965)
!       Rmin = 1.40121 au, Vmin = -4.74772265 eV relative to H+H
        c1=139.7160d0        ! eV
        c2=-123.8978d0       ! eV / bohr
        c3=3.4031d0          ! 1 / bohr
        c4=-6.8725d0         ! eV / bohr**2
        c5=-23.0440d0        ! eV / bohr**3
        c6=2.032d0           ! 1 / bohr

        v=(c1+c2*rr)*dexp(-c3*rr)
     &   +(c4+c5*rr)*dexp(-c6*rr)*rr**2
c       move zero from asymptote to minimum
        v=v+4.74772265
c        print *,rr,v
        v=v/autoev

        dvdr=((c1+c2*rr)*(-c3)+c2)*dexp(-c3*rr)
     &      +((c4+c5*rr)*(-c6)+c5)*dexp(-c6*rr)*rr**2
     &       +(c4+c5*rr)*dexp(-c6*rr)*rr*2.d0
        dvdr=dvdr/autoev

        elseif (at(1).eq.26.and.at(2).eq.26) then
!       N2 bath
!       fit to MRCI+Q/CBS(AQZ,A5Z) full valence
!       agrees reasonably well with more complicated form of LeRoy (JCP 125, 164310 (2006))
!       Jasper June 9, 2010
        de=79845.d0 ! exp De in cm-1 (Ronin, Luanay, Larzillier, PRL 53, 159 (1984), as quoted by LeRoy)
        re=1.097679d0 ! exp
        c1=2.68872341 ! my fit
        c2=0.240070803
        c3=0.472261727

        yy=rr*autoang-re
        beta = c1+c2*yy+c3*yy**2
        v = de*(1.d0-dexp(-beta*yy))**2-de    ! A and cm-1
        v=v/autocmi  ! convert to au

        dvdr = 2.d0*de*(1.d0-dexp(-beta*yy))*dexp(-beta*yy)*beta
        dvdr=dvdr*autocmi/autocmi  ! convert to au

!       elseif (at(1).eq.??.and.at(2).eq.??) then
!       OTHER DIATOMIC BATHS HERE
        else
        print *,"Don't know this diatomic bath"
        stop
        endif

      dvdx(1) =  dvdr*dx/rr
      dvdx(2) = -dvdr*dx/rr
      dvdy(1) =  dvdr*dy/rr
      dvdy(2) = -dvdr*dy/rr
      dvdz(1) =  dvdr*dz/rr
      dvdz(2) = -dvdr*dz/rr

      endif

      return
      end



      subroutine rgexp(at,x,y,z,v,dvdx,dvdy,dvdz,natom,maxatom)
c Rare Gas exp6 potential subroutine
c loops over geometry and looks for Rg-X interactions
c returns the full Rg-target intermolecular potential and its derivatives

      implicit real*8(a-h,o-z)
      dimension x(maxatom),y(maxatom),z(maxatom)
      dimension dvdx(maxatom),dvdy(maxatom),dvdz(maxatom)
      parameter(autoev=27.2113961d0)
      parameter(autocmi=219474.63067d0)
      parameter(autokcal=627.509d0)
      parameter(autoang=0.529177249d0)
      integer at(maxatom)
      logical troya,cutoff

      integer nfitparams,mfitparams,ix
      parameter (mfitparams=50)
      double precision ccc(mfitparams)
      common/amfit/ccc,nfitparams

      save/amfit/

      v1=0.d0
      v=0.d0
      do i=1,natom
      dvdz(i)=0.d0
      dvdy(i)=0.d0
      dvdx(i)=0.d0
      enddo

      do 1 i=1,natom
      do 2 j=i+1,natom

      m1=min(at(i),at(j))
      m2=max(at(i),at(j))
      troya=.false.   ! do or don't use Troya's form
      cutoff=.false.   ! do or don't use cutoff

      if (m1.ge.21) then ! two rare gases, skip this pair
         go to 2
      endif
      if (m2.le.20) then ! no rare gas, skip this pair
         go to 2
      endif

      if (m2.eq.m2) then  ! ANYTHING
        ix=m1
        if (ix.gt.2) ix=ix-1
        aa = ccc((ix-1)*4+1)
        bb = ccc((ix-1)*4+2)
        cc = ccc((ix-1)*4+3)
        rrc = ccc((ix-1)*4+4)
        aa=(10.d0**aa)
        cutoff=.true.
      endif

      dx=x(i)-x(j)
      dy=y(i)-y(j)
      dz=z(i)-z(j)
      rr=dsqrt(dx*dx+dy*dy+dz*dz)
      rra=rr*autoang

! NOTE CANNOT HAVE BOTH TROYA FORM AND CUTOFF FORM

      if (troya) then   ! Troya uses different form & units
        v=aa*dexp(-rra*bb)+cc/rra**6
        v=v/autokcal
        dvdr = -aa*bb*dexp(-rra*bb)-6.d0*cc/rra**7
        dvdr=dvdr/autokcal*autoang
      elseif (cutoff) then  ! cutoff 1/R**-6 at short distances
        v=aa*dexp(-rra/bb)-(cc**6/(rra**6+rrc**6))
        v=v/autocmi
        dvdr = -aa/bb*dexp(-rra/bb)
     &      +6.d0*(cc**6)*(rra**5)/(rra**6+rrc**6)**2
        dvdr=dvdr/autocmi*autoang
      else
        v=aa*dexp(-rra/bb)-(cc/rra)**6
        v=v/autocmi
        dvdr = -aa/bb*dexp(-rra/bb)+(6.d0/rra)*(cc/rra)**6
        dvdr=dvdr/autocmi*autoang
      endif
        v1=v1+v

c      print *,m1,m2,rra,v*autocmi,v1*autocmi

c derivs = sum over all bonds (DV/DRij * DRij/DXi = DV/DRij * (Xi-Xj)/Rij)
      dvdx(i) = dvdx(i) + dvdr*dx/rr
      dvdx(j) = dvdx(j) - dvdr*dx/rr
      dvdy(i) = dvdy(i) + dvdr*dy/rr
      dvdy(j) = dvdy(j) - dvdr*dy/rr
      dvdz(i) = dvdz(i) + dvdr*dz/rr
      dvdz(j) = dvdz(j) - dvdr*dz/rr

    2 continue
    1 continue
      v=v1
c      print *,'rgexp',v1
      return
      end

c **********************************************************************
c **********************************************************************
      subroutine tinkerpot(symb,xx,yy,zz,pema,dxx,dyy,dzz,nat,mnat)

      save iprepot
      entry prepot
      if (iprepot.eq.0) then
      call prepot2
      call prepot3
      iprepot=1

      return
      endif

      end
c **********************************************************************
c **********************************************************************




********************************************

      subroutine prepot2

      implicit double precision(a-h,o-z)
      include 'paramls.inc'
      dimension iagroup(maxatom),
     &  ind(maxterm,maxpair),nfound(maxatom),
     &  iatom(maxperm,maxatom),
     &  idum(maxatom),nngroup(maxatom),
     &  idum2(maxperm,maxatom),
     &  idum3(maxperm,maxatom),nperm0(maxperm),nperm1(maxperm),
     &  basis(maxterm),ibasis(maxterm),r(maxpair),
     &  rrr(maxdata,maxpair),index(maxatom,maxatom),ix(maxperm,maxpair)
      character*2 symb(maxatom),dum
      logical lreadbasis
 
      common/foox/rrr,nncoef,natom1

      save npairs,nterms,ind,ibasis

ccc GENERATE BASIS ccc
      if (.false.) then
ccc GENERATE BASIS ccc


c     generate atom permutation lists
      do i=1,natom
      nngroup(i)=0
      enddo
      do i=1,natom
      if (iagroup(i).gt.ngroup) ngroup=iagroup(i)
      nngroup(iagroup(i))=nngroup(iagroup(i))+1
      enddo

      nn=0

      do i=1,ngroup

      n=0
      do k=1,natom
      if (iagroup(k).eq.i) then
      n=n+1
      idum(n)=k
      endif
      enddo
      
      npermute=0
      call heapp(idum,n,n,idum2,npermute)
      nperm0(i)=nn+1
      nperm1(i)=nn+npermute
      do k=1,npermute
      nn=nn+1
      m=0
      do j=1,natom
      idum3(nn,j)=0
      if (iagroup(j).eq.i) then
          m=m+1
          idum3(nn,j)=idum2(k,m)
      endif
      enddo
      enddo

      enddo

      ntmp=1
      do i=1,ngroup
      idum(i)=nperm0(i)
      print *,"Group ",i," has ",(nperm1(i)-nperm0(i)+1)," permutations"
      ntmp=ntmp*(nperm1(i)-nperm0(i)+1)
      enddo
      print *,"For a total of ",ntmp," permutations"

      npermute=0
      do while (.true.)
        npermute=npermute+1
        if (npermute.gt.maxperm) then
        print *,"npermute (",npermute,") > maxperm (",maxperm,")"
        print *,"NOTE: maxperm needs to be at least npermute + 1"
        stop
        endif

        do i=1,natom
        iatom(npermute,i)=0
        do j=1,ngroup
        iatom(npermute,i)=iatom(npermute,i)+idum3(idum(j),i)
        enddo
        enddo

        idum(ngroup)=idum(ngroup)+1
 777    continue

        do i=1,ngroup
        if (idum(i).gt.nperm1(i)) then
        if (i.eq.1) go to 778
        idum(i)=nperm0(i)
        idum(i-1)=idum(i-1)+1
        go to 777
        endif
        enddo 

      enddo
 778  continue

      print *
      print *,'Atom permutations',npermute
      do i=1,min(npermute,100)
      print *,i,":",(iatom(i,j),j=1,natom)
      enddo

      ii=0
      do i=1,natom1
      do j=natom1+1,natom
      ii=ii+1
      index(i,j)=ii
      enddo
      enddo

      write(6,*)
      write(6,*)"Pair permutations"
      write(6,'(22x,100(a3,"- ",a3,4x))')
     &   ((symb(i),symb(j),j=1+natom1,natom),i=1,natom1) 
      write(6,'(21x,100(i3," -",i3,4x))')((i,j,j=1+natom1,natom),
     &   i=1,natom1) 
      do ii=1,npermute
      iix=0
      do i=1,natom1
      do j=natom1+1,natom
      iix=iix+1
      ix(ii,iix)=index(iatom(ii,i),iatom(ii,j))
      enddo
      enddo
      if (ii.le.100) print *,ii,":",(ix(ii,iix),iix=1,npairs)
      enddo

c generate terms using individual power constraints
      ii=1
      do i=1,npairs
      ind(ii,i)=0
      enddo
      do while (.true.)
        ii=ii+1
        if (ii.gt.maxterm) then
      print *,"number of terms (",ii,") > maxterm (",maxterm,")"
        stop
        endif

        do i=1,npairs
        ind(ii,i)=ind(ii-1,i)
        enddo
        ind(ii,npairs)=ind(ii,npairs)+1
 300    continue
        indtot=0
        do i=1,npairs
        indtot=indtot+ind(ii,i)
        if (ind(ii,i).gt.ipow.or.indtot.gt.ipowt) then ! ipow(i) would allow atom-atom-type-dependent limits
        if (i.eq.1) go to 400
        ind(ii,i)=0
        ind(ii,i-1)=ind(ii,i-1)+1
        go to 300
        endif
        enddo
      enddo
 400  continue
      nterms=ii-1

      print *
      print *,"Basis # (Group):  Powers"

c symmetrize
      nbasis=0
      DO ii=1,nterms
      ifail=0
      do i=1,ii-1
        do j=1,npermute
          ifail=1
          do k=1,npairs
            if (ind(i,k).ne.ind(ii,ix(j,k))) ifail=0
          enddo
          if (ifail.eq.1) go to 1010
        enddo
      enddo
 1010 continue

      if (ifail.eq.0) then
      nbasis=nbasis+1
      ibasis(ii)=nbasis
      else
      ibasis(ii)=ibasis(i)
      endif
      write(6,'(i5,"  (",i5,"):",100i8)')
     &   ii,ibasis(ii),(ind(ii,j),j=1,npairs)
      ENDDO

      nncoef=nbasis
      print *,'nncoef = ',nncoef
 
      open(55,file="basis.dat")
      write(55,*)natom1,npairs,nncoef,nterms,
     & " ! atom pairs, coefficients, terms"
      write(55,*)" TERM GROUP :     EXPONENTS"
      do ii=1,nterms
      write(55,'(2i6," : ",1000i5)')
     &   ii,ibasis(ii),(ind(ii,j),j=1,npairs)
      enddo
      close(55)

ccc READ BASIS ccc
      else
      open(55,file="basis.dat")
      read(55,*)natom1,npairs,nncoef,nterms
      read(55,*)
      do i=1,nterms
      read(55,*)k,ibasis(k),dum,(ind(k,j),j=1,npairs)
      print *,k,ibasis(k)
      enddo
      close(55)

      endif

      return

      entry funcs1(iii,basis,ncoef)

c print *,ncoef,npairs,nterms

      do j=1,ncoef
      basis(j)=0.d0
      enddo

      do j=1,npairs
      r(j)=dexp(-rrr(iii,j)*autoang)
      enddo
c      r(3)=dexp(-rrr(iii,6)*autoang)  ! AJ HACK for CH3OH+Ar
c      r(6)=dexp(-rrr(iii,3)*autoang)
c      r(8)=dexp(-rrr(iii,9)*autoang)  ! AJ HACK for C2H5OH+Ar
c      r(9)=dexp(-rrr(iii,8)*autoang)
c      r(5)=dexp(-rrr(iii,11)*autoang)  ! AJ HACK for CH3OH+N2
c      r(11)=dexp(-rrr(iii,5)*autoang)
c      r(6)=dexp(-rrr(iii,12)*autoang)  ! AJ HACK for CH3OH+N2
c      r(12)=dexp(-rrr(iii,6)*autoang)

      do i=1,nterms
      arg=1.d0
      do j=1,npairs
      arg=arg*(r(j)**ind(i,j))
c      if (iii.eq.1) print *,"lll",i,j,arg
      enddo
      basis(ibasis(i))=basis(ibasis(i))+arg
      enddo

      return 
      end

***************************************************

      recursive subroutine heapp(ia,size,n,iia,ii)

      include 'paramls.inc'
      integer i,n,size,ii
      integer ia(maxatom)
      integer iia(maxperm,maxatom)
      integer iagroup(maxatom)

      if (size.eq.1) then
         ii=ii+1
         do i=1,n
           iia(ii,i)=ia(i)
         enddo
        return
      endif

      do i=1,size
        call heapp(ia,size-1,n,iia,ii)
        if (mod(size,2).eq.1) then
          tmp=ia(1)
          ia(1)=ia(size)
          ia(size)=tmp
        else
          tmp=ia(i)
          ia(i)=ia(size)
          ia(size)=tmp
      endif

      enddo

      end subroutine

***************************************************
      subroutine lsfit(at,x,y,z,v,dvdx,dvdy,dvdz,natom,maxatomx)
c
      implicit double precision (a-h,o-z)
c
      include 'paramls.inc'
      dimension coef(maxterm),sig(maxdata)
      dimension basis(maxterm)

      dimension vv(maxdata),rrr(maxdata,maxpair)
      dimension vv2(maxdata),rrr2(maxdata,maxpair)
      dimension rcom2(maxdata),xprint(50,20)
      dimension x(maxatom),y(maxatom),z(maxatom)

      dimension dvdx(maxatom),dvdy(maxatom),dvdz(maxatom)
      parameter(autoev=27.2113961d0)
      parameter(autocmi=219474.63067d0)
      parameter(autokcal=627.509d0)
      integer at(maxatom)


      character*2 dum

      common/foox/rrr,nncoef,natom1

      save iprepot3
      entry prepot3
      if (iprepot3.eq.0) then
      open(77,file="coef.dat")
      do k=1,nncoef
      read(77,*)i,coef(k)
      print *,i,coef(k)
      enddo
      iprepot3=1
      return
      endif

c      natom1=6
      ii=0
      do j=1,natom1  ! molecule
      do k=natom1+1,natom  ! bath
      ii=ii+1
      rrr(1,ii)=dsqrt((x(j)-x(k))**2+(y(j)-y(k))**2+(z(j)-z(k))**2)
      enddo  
      enddo  

      ncoef=nncoef
      call funcs1(1,basis,ncoef) 
      v=0.d0
      do j=1,ncoef
         v=v+coef(j)*basis(j)
      enddo
      v=v/autocmi

      resp=0.0001
      ii=0
      do j=1,natom1  ! molecule
      do k=natom1+1,natom  ! bath
      dx=x(j)-x(k)
      dy=y(j)-y(k)
      dz=z(j)-z(k)
      ii=ii+1

      rrr(1,ii)=rrr(1,ii)+resp
      call funcs1(1,basis,ncoef) 
      vp=0.d0
      do l=1,ncoef
         vp=vp+coef(l)*basis(l)
      enddo
      vp=vp/autocmi

      rrr(1,ii)=rrr(1,ii)-2.d0*resp
      call funcs1(1,basis,ncoef) 
      vm=0.d0
      do l=1,ncoef
         vm=vm+coef(l)*basis(l)
      enddo
      vm=vm/autocmi
      rrr(1,ii)=rrr(1,ii)+resp

      dtmpdrr=(vp-vm)/(2.d0*resp)
      dvdx(j) = dvdx(j) + dtmpdrr*dx/rrr(1,ii)
      dvdx(k) = dvdx(k) - dtmpdrr*dx/rrr(1,ii)
      dvdy(j) = dvdy(j) + dtmpdrr*dy/rrr(1,ii)
      dvdy(k) = dvdy(k) - dtmpdrr*dy/rrr(1,ii)
      dvdz(j) = dvdz(j) + dtmpdrr*dz/rrr(1,ii)
      dvdz(k) = dvdz(k) - dtmpdrr*dz/rrr(1,ii)
      enddo
      enddo

      return
 
      end
