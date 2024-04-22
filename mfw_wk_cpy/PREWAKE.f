c    <Copyright (c) 1998 by Mahendra Bhagwat>
c
      IMPLICIT none
c
      INTEGER nb,zif,pif,z,nz,nzt,p,np,r,nr,s,ns,w,nw,pw,nk,li,gs,k,
     &        nctrl,usaero,nfw
c
      DOUBLE PRECISION pi,ft,dp,dz,mu,muc,nt,ip,rad,flph,
     &       crd,twt,om,rada,anb,
     &       rcout,sce,vn,rcb,dcy,rpar,rrn,zrn,cnv,cvt,lin,ct0,asr0,
     &       cd0,cd1,cd2,b1c0,b1s0,t0_0,t1c0,t1s0,b0_0,
     &       cttol,fltol,lock,sig,den,bmass,ibeta,sbeta,nubeta,
     &       kbeta,betap,vinfx,vinfy,vinfz,sonic,sn,kinv,
     &       ub1,ub2,ua1,ua2,uf1,uf2,ug1,ug2,pbar,qbar
      DOUBLE PRECISION bct,alpha0
c
      CHARACTER*1 sco,trm,cvr,storetrm,method,initial
      CHARACTER*18 temp1
      CHARACTER*22 temp2
      CHARACTER*1 temp3
      INTEGER INP
      INTEGER nstaper, rotgeo
      DOUBLE PRECISION taperst,taper,altitude,isaplus
      LOGICAL outplot
c
c----------------------------------------------------------------------
c >>  DATA INPUT FILES <<
c
c >>  User defined inputs: <<
      namelist/user/nw,pw,nfw,nk,ft,bct,dp,dz,sco,sce,vn,rcb,dcy,rpar,
     &        rrn,zrn,cnv,trm,cvr,cvt,lin,li,gs,method,initial,
     &        outplot

c >>  Rotor flight and operating conditions: <<
      namelist/flight/mu,muc,pbar,qbar,ct0,t0_0,t1c0,t1s0,
     &                b0_0,b1c0,b1s0,cttol,fltol,altitude,isaplus

c >>  Rotor geometry: <<
      namelist/geometry/nr,nb,ns,asr0,rad,flph,
     &     rcout,crd,twt,om,taperst,taper,rotgeo

c >>  Rotor mass and inertial properties: <<
      namelist/rotprop/bmass,ibeta,sbeta,nubeta,kbeta,betap

c >>  Unsteady aerodynamics coefficients: <<
      namelist/usa/usaero,ub1,ub2,ua1,ua2,uf1,uf2,ug1,ug2,alpha0
c----------------------------------------------------------------------
c
      open(33,file='user.input',status='old')
      read(33,user)
      close(33)
 
      open(44,file='flight.input',status='old')
      read(44,flight)
      close(44)

      open(55,file='geometry.input',status='old')
      read(55,geometry)
      close(55)

      open(66,file='rotprop.input',status='old')
      read(66,rotprop)
      close(66)

      open(77,file='usa.input',status='old')
      read(77,usa)
      close(77)
c ----------------------------------------------------------------------
c
      pi=4.D+0*atan(1.D+0)
c ----------------------------------------------------------------------
c
      dp=dp*pi/180.D+0
      dz=dz*pi/180.D+0
      zif=1
      pif=1
      if (dz .gt. dp) zif=nint(dz/dp)
      if (dp .gt. dz) pif=nint(dp/dz)

      anb=2.D+0*pi/DBLE(nb)
      np=nint(2.D+0*pi/dp)*pif

      nt=ft+bct
c-mb      if (abs(mu).lt.0.005D+0 .AND. abs(muc).lt.0.005D+0) nt=ft+2.D+0
c-mb      at some stage the BC was used only in hover
      ip=nt*2.D+0*pi
      ip=0.D+0

      nzt=nint(nt*2.D+0*pi/dz)*zif+1

c-mb   adim - azimuthal locations
      open(11,file='adimpar.inc',status='unknown')
      read(11,20)temp1
      read(11,10)temp2,INP,temp3
      close(11)
      if(INP.ne.np) then
      open(11,file='adimpar.inc',status='unknown')
      write(11,20) '      INTEGER adim'
      write(11,10)'      parameter (adim=',np,')'
      close(11)
      endif
c-mb   zdim - points along the filament 
      open(11,file='zdimpar.inc',status='unknown')
      read(11,20)temp1
      read(11,10)temp2,INP,temp3
      close(11)
      if(INP.ne.nzt) then
      open(11,file='zdimpar.inc',status='unknown')
      write(11,20)'      INTEGER zdim'
      write(11,10)'      parameter (zdim=',nzt,')'
      close(11)
      endif
c-mb   rdim - number of rotors
      open(11,file='rdimpar.inc',status='unknown')
      read(11,20)temp1
      read(11,10)temp2,INP,temp3
      close(11)
      if(INP.ne.nr) then
      open(11,file='rdimpar.inc',status='unknown')
      write(11,20)'      INTEGER rdim'
      write(11,10)'      parameter (rdim=',nr,')'
      close(11)
      endif
c-mb   kdim - number of cores!
      open(11,file='kdimpar.inc',status='unknown')
      read(11,20)temp1
      read(11,10)temp2,INP,temp3
      close(11)
      if(INP.ne.nk) then
      open(11,file='kdimpar.inc',status='unknown')
      write(11,20)'      INTEGER kdim'
      write(11,10)'      parameter (kdim=',nk,')'
      close(11)
      endif
c-mb   wdim - number of trailers
      open(11,file='wdimpar.inc',status='unknown')
      read(11,20)temp1
      read(11,10)temp2,INP,temp3
      close(11)
      if(INP.ne.nw) then
      open(11,file='wdimpar.inc',status='unknown')
      write(11,20)'      INTEGER wdim'
      write(11,10)'      parameter (wdim=',nw,')'
      close(11)
      endif
c-mb   sdim - number of blade segments
      open(11,file='sdimpar.inc',status='unknown')
      read(11,20)temp1
      read(11,10)temp2,INP,temp3
      close(11)
      if(INP.ne.ns+1) then
      open(11,file='sdimpar.inc',status='unknown')
      write(11,20)'      INTEGER sdim'
      write(11,10)'      parameter (sdim=',ns+1,')'
      close(11)
      endif

c-mb   gdim - grid sequencing something
c-mb   this is four for now
      open(11,file='gdimpar.inc',status='unknown')
      read(11,20)temp1
      read(11,10)temp2,INP,temp3
      close(11)
      if(INP.ne.4) then
      open(11,file='gdimpar.inc',status='unknown')
      write(11,20)'      INTEGER gdim'
      write(11,10)'      parameter (gdim=',4,')'
      close(11)
      endif
c-mb   jdim - something about the jacobian
c-mb   this was set to seven for some reason I think it needs to be
c-mb   well, in the ctrlsolve routine it needs +1
      open(11,file='jdimpar.inc',status='unknown')
      read(11,20)temp1
      read(11,10)temp2,INP,temp3
      close(11)
      if(INP.ne.(3*nr+1)) then
      open(11,file='jdimpar.inc',status='unknown')
      write(11,20)'      INTEGER jdim'
      write(11,10)'      parameter (jdim=',3*nr+1,')'
      close(11)
      endif

      write(*,*)'Set up the required dimensions for running FREEWAKE'

10    format(a22,I4,a1)
20    format(a18)
      stop
      end
