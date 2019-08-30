      program onedmin

c Notes: This version is compiled with the "universal" TB+exp/6 PES for CxHy+M.
c        PES evalulations are very cheap.
c        This PES can be inaccurate for unsaturated CxHy + H2,N2,O2.
c        Other PESs can be readily added, including direct dynamics.

      implicit none

c     dimensions
      integer mnat,mnsamp
      parameter(mnat=100) ! max number of atoms for each collider
      parameter(mnsamp=50000) ! max number of samples

c     units and constants
      double precision amutoau,kb,autoang,autoev,autofs,mu,autocmi,pi
      parameter(amutoau=1822.844987d0) ! amu to au
      parameter(kb=3.166829d-6)  ! Boltzmann in hartee/K
      parameter(autoang=0.52917706d0) ! au to A
      parameter(autoev=27.2113961d0) ! au to eV
      parameter(autofs=0.024189d0) ! au to fs
      parameter(autocmi=219474.63067d0) ! au to cm-1
      parameter(mu=1.d0*amutoau)  ! mass-scaling mass
      parameter(pi=3.1415926536d0)

c
      integer ntot0,ntot1,nn,is,nsamp,io,
     &    nr,ir,imin,ifail,nn0,nn1,ranseed
      character*16 filex0,filex1
      character*2 symb0(mnat),symb1(mnat),symb(mnat),symbz(mnat),tmpsymb
      double precision xx(3,mnat),mm0(mnat),mm1(mnat)
      double precision xxz(3,mnat)
      double precision xx0(mnsamp,3,mnat)
      double precision xx1(mnsamp,3,mnat)
      double precision mtot,xtot,ytot,ztot,v,vp,vm,reps,gv,rdiff,hv
      double precision x(mnat),y(mnat),z(mnat)
      double precision dx(mnat),dy(mnat),dz(mnat)
      double precision vmin,vmax,vavg,vstd
      double precision rrmin,rrmax,rravg,rrstd,rrlast
      double precision ssmin,ssmax,ssavg,ssstd
      double precision smin,smax,rmin,rmax,xr,xv,zero
      double precision mass0,mass1
      integer ncalc1,ncalc2,zeroflag
      logical lfail
      integer n,i,j,nrcl,itmp1,itmp2
      double precision rr,r

c common block used to pass the interaction potential
      double precision tmpprint(50)
      common/tmp/tmpprint

      call prepot

c READ STANDARD INPUT
c     Standard input must contain these 4 lines
      read(5,*)ranseed,nsamp
      read(5,*)filex0  !         file name of the geometry file of the target
      read(5,*)filex1  !         file name of the geometry file of the bath
      read(5,*)smin,smax
      zeroflag=0

c hard code the rest of the input for automated kinetics
      rmin = smin * (2.d0**(1.d0/6.d0))
      rmax = smax * (2.d0**(1.d0/6.d0)) 
      zero = 0.d0  ! zero for the PES used here; assume the "universal" TB+exp/6 PES is being used
      ntot0=1
      ntot1=1

c Write output msg and the input variables to the output file
      write(6,*)'---------------------------'
      write(6,*)'OneDMin'
      write(6,*)'Argonne National Laboratory'
      write(6,*)'---------------------------'
      write(6,*)
      write(6,*)'Job Parameters:'
      write(6,*)'RANSEED = ',ranseed
      write(6,*)'TARGET GEOMS: NUMBER = ',ntot0,', FILE = ',filex0
      write(6,*)'BATH GEOMS:   NUMBER = ',ntot1,', FILE = ',filex1
      write(6,*)'SAMPLES = ',nsamp
      write(6,*)'S range = ',smin,smax,' A'
      write(6,*)'R range = ',rmin,rmax,' A'
      write(6,*)
      write(6,*)

! initialize RNG
      call srand(ranseed)

! read XYZ file for the target
      open(unit=93,file=filex0)
      mass0=0.d0
      j=0
      do
        j=j+1
        read(93,*,END=321)symb0(j),xx0(1,1,j),xx0(1,2,j),xx0(1,3,j)
        xx0(1,1,j)=xx0(1,1,j)/autoang  ! A to a0
        xx0(1,2,j)=xx0(1,2,j)/autoang
        xx0(1,3,j)=xx0(1,3,j)/autoang
        call mass(symb0(j),mm0(j))
        mass0=mass0+mm0(j)
        if (symb0(j).eq."c") symb0(j)="C"
        if (symb0(j).eq."h") symb0(j)="H"
      enddo
 321  close(93)
      nn0=j-1

! read XYZ file for the bath
      open(unit=94,file=filex1)
      mass1=0.d0
      j=0
      do
        j=j+1
        read(94,*,END=322)symb1(j),xx1(1,1,j),xx1(1,2,j),xx0(1,3,j)
        xx1(1,1,j)=xx1(1,1,j)/autoang  ! A to a0
        xx1(1,2,j)=xx1(1,2,j)/autoang
        xx1(1,3,j)=xx1(1,3,j)/autoang
        call mass(symb1(j),mm1(j))
        mass1=mass1+mm1(j)
        if (symb1(j).eq."c") symb1(j)="C"
        if (symb1(j).eq."h") symb1(j)="H"
      enddo
 322  close(94)
      nn1=j-1

! Initialize averages
      vavg=0.d0
      vstd=0.d0
      vmax=-10000.d0
      vmin=0.d0
      ssmax=0.d0
      ssmin=100.d0
      rrmax=0.d0
      rrmin=100.d0
      rravg=0.d0
      ssavg=0.d0
      rrstd=0.d0
      ssstd=0.d0

      write(6,*)'Calculating Energy at Infinite Separation...'
      write(6,*)'Energy (in Eh) will be written to zero.ene'
      write(6,*)
      write(6,*)
      write(6,*)'Calculating Sigma and Epsilon for each orientation...'
      write(6,*)

c LOOP OVER ORIENTATION SAMPLES
      do is = 1,nsamp   
        
        write(6,*)'ORIENTATION',is
        write(6,*)

 22     continue

        ncalc1=0
        ncalc2=0
        rr=(rmax+rmin)/2.d0  ! first guess

!       Sample and spin the bath
        r=dble(rand())
        r=r*dble(ntot1)
        n=int(r)+1
        do j=1,nn1
          xx(1,j)=xx1(n,1,j)
          xx(2,j)=xx1(n,2,j)
          xx(3,j)=xx1(n,3,j)
        enddo
        xtot=0.d0
        ytot=0.d0
        ztot=0.d0
        mtot=0.d0
        do i=1,nn1
            xtot = xtot + mm1(i)*xx(1,i)
            ytot = ytot + mm1(i)*xx(2,i)
            ztot = ztot + mm1(i)*xx(3,i)
            mtot = mtot + mm1(i)
        enddo
        do i=1,nn1
            xx(1,i) = xx(1,i) - xtot/mtot
            xx(2,i) = xx(2,i) - ytot/mtot
            xx(3,i) = xx(3,i) - ztot/mtot
        enddo
        call spin(xx,nn1)

! move the bath to the end of the xx array
        do i=1,nn1
           xxz(1,i+nn0)=xx(1,i)
           xxz(2,i+nn0)=xx(2,i)
           xxz(3,i+nn0)=xx(3,i)
           symbz(i+nn0)=symb1(i)
        enddo
        do i=1,nn1
           xx(1,i+nn0)=xxz(1,i+nn0)
           xx(2,i+nn0)=xxz(2,i+nn0)
           xx(3,i+nn0)=xxz(3,i+nn0)
           symb(i+nn0)=symbz(i+nn0)
        enddo

! Sample and spin the target
        r=dble(rand())
        r=r*dble(ntot0)
        n=int(r)+1
        do j=1,nn0
          xx(1,j)=xx0(n,1,j)
          xx(2,j)=xx0(n,2,j)
          xx(3,j)=xx0(n,3,j)
        enddo
        xtot=0.d0
        ytot=0.d0
        ztot=0.d0
        mtot=0.d0
        do i=1,nn0
            xtot = xtot + mm0(i)*xx(1,i)
            ytot = ytot + mm0(i)*xx(2,i)
            ztot = ztot + mm0(i)*xx(3,i)
            mtot = mtot + mm0(i)
        enddo
        do i=1,nn0
            xx(1,i) = xx(1,i) - xtot/mtot
            xx(2,i) = xx(2,i) - ytot/mtot
            xx(3,i) = xx(3,i) - ztot/mtot
        enddo
        call spin(xx,nn0)

        do i=1,nn0
           symb(i)=symb0(i)
        enddo
        
        write(6,*)'Finding distance R for bottom of well [min(V(R))]'
        write(6,'(a5,2a16)')' Step','R','V(R)'

        ifail=0
 10     continue
        lfail=.false.
          
c       NR search for minimum energy
        rdiff = 100.
        do while (rdiff > .001d0)   ! Convergence tolerance of the CoM distance in A
          ncalc1=ncalc1+1
          reps = 0.001d0            ! Numerical gradient and hessian parameter

          nn=nn0+nn1

          do i=1,nn0
            x(i) = xx(1,i)
            y(i) = xx(2,i)
            z(i) = xx(3,i)
          enddo

          do i=nn0+1,nn
           x(i) = xx(1,i)+(rr+reps)/autoang
           y(i) = xx(2,i)
           z(i) = xx(3,i)
          enddo

          if (zeroflag.eq.0) then
            do i=nn0+1,nn
              x(i) = xx(1,i)+30./autoang
            enddo
            call pot(symb,x,y,z,zero,dx,dy,dz,nn,nn)
            zero=tmpprint(1)
            open(unit=82,file='zero.ene')
            write(82,*)zero
            close(82)
            zeroflag=1
            do i=nn0+1,nn
              x(i) = xx(1,i)+(rr+reps)/autoang
            enddo
          endif

          call pot(symb,x,y,z,vp,dx,dy,dz,nn,nn)
          vp=(tmpprint(1)-zero)*autocmi           ! the interaction is passed via a common block

          do i=nn0+1,nn
            x(i) = xx(1,i)+(rr-reps)/autoang
            y(i) = xx(2,i)
            z(i) = xx(3,i)
          enddo

          call pot(symb,x,y,z,vm,dx,dy,dz,nn,nn)
          vm=(tmpprint(1)-zero)*autocmi

          do i=nn0+1,nn
            x(i) = xx(1,i)+rr/autoang
            y(i) = xx(2,i)
            z(i) = xx(3,i)
          enddo
        
          call pot(symb,x,y,z,v,dx,dy,dz,nn,nn)
          v=(tmpprint(1)-zero)*autocmi

          gv=(vp-vm)/(2.d0*reps)       ! numerical gradient w.r.t. the CoM distance
          hv=(vp+vm-2.d0*v)/(reps**2)  ! numerical hessian w.r.t. the CoM distance

          rrlast=rr
          rr=rr-gv/hv               
          rdiff=abs(gv/hv)             ! optimize to minimum

          write(6,'(i5,2f16.8)')ncalc1,rr,v
c          write(6,*)'r',ncalc1,rr,rrlast,v,vp,vm

c         try to catch bad optimizations so we can start again with a different guess
          if (rr.gt.rmax) then 
            print *,"failed rmax"
               ifail=ifail+1
               lfail=.true.
          elseif (rr.lt.rmin) then 
            print *,"failed rmin"
               ifail=ifail+1
               lfail=.true.
          elseif (rr.ne.rr) then 
            print *,"failed NAN"
               ifail=ifail+1
               lfail=.true.
          elseif (v.lt.-5000) then 
            print *,"failed low v"
               ifail=ifail+1
               lfail=.true.
          elseif (v.gt.5000) then 
            print *,"failed high v"
               ifail=ifail+1
               lfail=.true.
          elseif (ncalc1.gt.200) then 
            print *,"failed steps"
               ifail=ifail+1
               lfail=.true.
          endif
          if (lfail) then
            lfail=.false.
            if (ifail.gt.20) then
              go to 999
            elseif(ifail.gt.0.and.ifail.lt.10) then
              rr=rmax-dble(ifail-1)*(rmax-rmin)/dble(20)
              go to 10
            elseif(ifail.gt.9.and.ifail.lt.20) then
              rr=rmax-dble(ifail)*(rmax-rmin)/dble(20)
              go to 10
            endif
          endif

        enddo ! END MINIMUM SEARCH LOOP
          
        write(6,*)'Found R where V(R) is minimized!!!'
        write(6,*)

        rrmax=max(rrmax,rr)
        rrmin=min(rrmin,rr)
        rravg=rravg+rr
        rrstd=rrstd+rr**2
        vmax=max(vmax,v)
        vmin=min(vmin,v)
        vavg=vavg+v
        vstd=vstd+v**2

        xr=rr
        xv=v

! Write the minimum-energy complexes to a file
        open(unit=80,file='min_geoms.out',position='append')
        write(80,*)nn
        write(80,*)is,rr,v
        do i=1,nn
          tmpsymb=symb(i)
          if (symb(i).eq."N2") tmpsymb="N"
          write(80,*)tmpsymb,x(i)*autoang,y(i)*autoang,z(i)*autoang
        enddo
        close(80) 

        write(6,*)'Finding distance R where V(R) = 0:'
        write(6,'(a5,2a16)')' Step','R','V(R)'

! NR search for inner turning point
        rr=rr/(2.d0**(1.d0/6.d0))

        ifail=0
 11     continue
        lfail=.false.

        rdiff = 100.
        do while (rdiff > .001d0) 
          ncalc2=ncalc2+1
          reps = 0.0001d0

          do i=1,nn0
            x(i) = xx(1,i)
            y(i) = xx(2,i)
            z(i) = xx(3,i)
          enddo
          do i=nn0+1,nn
            x(i) = xx(1,i)+rr/autoang
            y(i) = xx(2,i)
            z(i) = xx(3,i)
          enddo
          
          call pot(symb,x,y,z,v,dx,dy,dz,nn,nn)
          v=(tmpprint(1)-zero)*autocmi

          do i=nn0+1,nn
            x(i) = xx(1,i)+(rr+reps)/autoang
            y(i) = xx(2,i)
            z(i) = xx(3,i)
          enddo

          call pot(symb,x,y,z,vp,dx,dy,dz,nn,nn)
          vp=(tmpprint(1)-zero)*autocmi

          do i=nn0+1,nn
            x(i) = xx(1,i)+(rr-reps)/autoang
            y(i) = xx(2,i)
            z(i) = xx(3,i)
          enddo

          call pot(symb,x,y,z,vm,dx,dy,dz,nn,nn)
          vm=(tmpprint(1)-zero)*autocmi

          gv=(vp-vm)/(2.d0*reps)       ! grad
          hv=(vp+vm-2.d0*v)/(reps**2)  ! hess

          rrlast=rr
          rr=rr-v/gv                   ! optimize to inner turning point
          rdiff=abs(v/gv)

          write(6,'(i5,2f16.8)')ncalc2,rr,v
c          write(6,*)'s',ncalc2,rr,rrlast,v,vp,vm

          if (rr.gt.rmax.or.rr.ne.rr) then 
             ifail=ifail+1
             lfail=.true.
          elseif (rr.lt.smin) then 
             ifail=ifail+1
             lfail=.true.
          elseif (v.lt.-5000) then 
             ifail=ifail+1
             lfail=.true.
          elseif (v.gt.5000) then 
             ifail=ifail+1
             lfail=.true.
          elseif (ncalc2.gt.500) then 
             ifail=ifail+1
             lfail=.true.
          endif
          if (lfail) then
            lfail=.false.
            if (ifail.gt.10) then
              go to 999
              stop
            elseif(ifail.gt.0.and.ifail.lt.5) then
              rr=rmax-dble(ifail-1)*(smax-smin)/dble(10)
              go to 11
            elseif(ifail.gt.4.and.ifail.lt.10) then
              rr=rmax-dble(ifail)*(smax-smin)/dble(10)
              go to 11
            endif
          endif

        enddo ! END INNTER TURNING POINT LOOP
        
        write(6,*)'R where V(R)=0 found!!!'
        write(6,*)

        ssmax=max(ssmax,rr)
        ssmin=min(ssmin,rr)
        ssavg=ssavg+rr
        ssstd=ssstd+rr**2

        write(6,'(a17,f15.8,a11)')'Writing Epsilon =', -1.0*xv, 
     & 'to lj.out.'
        write(6,'(a15,2x,f15.8,a11)')'Writing Sigma =', rr, 
     & 'to lj.out.'
        write(6,*)
        write(6,*)
c        write(6,'(3i10,20f15.5)')is,ncalc1,ncalc2,xr,
c     & xr/(2.d0**(1.d0/6.d0)),xv,rr,v
       
!       Write the LJ parameters to lj.out
        open(unit=81,file='lj.out',position='append')
        if (is.eq.1) then
        write(81,'(a11,a15,a19)') 'Orientation', 'Sigma (Ang)',
     &  'Epsilon (cm-1)'
        endif
        write(81,'(i11,f15.5,f19.5)') is,rr,-1.0*xv
        close(81)

      enddo ! END ORIENTATION LOOP

      vstd=dsqrt(dble(nsamp)*vstd-vavg**2)/dble(nsamp)
      vavg=vavg/dble(nsamp)
      rrstd=dsqrt(dble(nsamp)*rrstd-rravg**2)/dble(nsamp)
      rravg=rravg/dble(nsamp)
      ssstd=dsqrt(dble(nsamp)*ssstd-ssavg**2)/dble(nsamp)
      ssavg=ssavg/dble(nsamp)

      open(unit=81,file='lj.out',position='append')
      write(81,*)
      write(81,*)
      write(81,*)'Simple Statistical Analysis of Calculated Parameters:'
      write(81,*)

      write(81,'(21x,4a15)')'r','s(est)','v','s'
      write(81,'(a18,3x,20f15.5)')'minimum    ',rrmin,
     & rrmin/(2.d0**(1.d0/6.d0)),
     & vmin,ssmin
      write(81,'(a18,3x,20f15.5)')'maximum    ',rrmax,
     & rrmax/(2.d0**(1.d0/6.d0)),
     & vmax,ssmax
      write(81,'(a18,3x,20f15.5)')'average    ',rravg,
     & rravg/(2.d0**(1.d0/6.d0)),
     & vavg,ssavg
      write(81,'(a18,3x,20f15.5)')'standard dev. ',rrstd,
     & rrstd/(2.d0**(1.d0/6.d0)),
     & vstd,ssstd

      write(81,*)
      write(81,*)"    SIGMA  = ",ssavg," +/- ",ssstd," A"
      write(81,*)"   EPSILON = ",vavg," +/- ",vstd," cm-1"
      write(81,*)
      close(81)

      write(6,*)'Finished sampling desired number of orientatoons'
      write(6,*)'Exiting OneDMin...' 

 188  format(1a,2f15.5)

      stop

 999  continue

      open(unit=81,file='lj.out',position='append')
      write(81,*)"      Sigma       FAILED "
      write(81,*)"      Epsilon]    FAILED "
      close(81)

      end

