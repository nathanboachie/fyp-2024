c----------------------------------------------------------------------|
c    <Copyright (c) 1998 by Mahendra Bhagwat>                          |
c----------------------------------------------------------------------|
c
      subroutine farwake(rbc,zbc,avgzbc,rad,nr,nw,pw,np,nz,nzt,
     &       ppx,ppy,ppz,
     &       pox,poy,poz,rr_x,rr_y,rr_z,dp,dz,pif,zif,lami,Psi)
c----------------------------------------------------------------------|
c
c-mb  note that po[xyz] are useless as of now only extrapolate ppx
c-mb  for the interpolation to work correctly, we must interpolate
c-mb   all azimuths at each wake age i.e. z-loop outside p-loop
c-mb   i'm putting p=0 option so that I  will be able to use it with
c-mb   the transient code lateron
c
c----------------------------------------------------------------------|
c
      IMPLICIT none
c
      include 'adimpar.inc'
      include 'kdimpar.inc'
      include 'rdimpar.inc'
      include 'wdimpar.inc'
      include 'zdimpar.inc'
c
      INTEGER nr,np,nw,pw,nz,nzt,pif,zif,r,w,z,Psi,p,z1,nz1,nt,ntt,pm1
      INTEGER Psi1,Psi2
c
      DOUBLE PRECISION dp,dz,rr0,rr,rad,zr,sn,pi
      DOUBLE PRECISION a1,b1,a2,b2
c
      DOUBLE PRECISION rbc(rdim,wdim),zbc(rdim,wdim),avgzbc(rdim,wdim)
      DOUBLE PRECISION lami(rdim),rr_x(rdim),rr_y(rdim),rr_z(rdim)
c
      DOUBLE PRECISION ppx(rdim,adim,wdim,zdim),
     &       ppy(rdim,adim,wdim,zdim),
     &       ppz(rdim,adim,wdim,zdim),
     &       pox(rdim,adim,wdim,zdim),
     &       poy(rdim,adim,wdim,zdim),
     &       poz(rdim,adim,wdim,zdim)
c
      DOUBLE PRECISION eps
c----------------------------------------------------------------------|
c
      common/pidata/pi
      save/pidata/
c----------------------------------------------------------------------|
c
      eps=0.D+0
c
      Psi1=Psi
      Psi2=Psi
      if(Psi.eq.0) Psi1=1
      if(Psi.eq.0) Psi2=np
c
      nt=nint(2.D+0*pi/dz)*zif         !number of pts in one turn
      nz1=nz-nt                       !one turn before last free pt
c----------------------------------------------------------------------|
c
      do r=1,nr
       sn=(-1.D+0)**(r-1)
c-mb   sign of rotation is not used here
c----------------------------------------------------------------------|
c
       do w=1,nw-pw ! initial prescribed inboard sheet
c-mb-choice       do w=1,nw    ! solve for inboard sheet
c
        z1=0
c----------------------------------------------------------------------|
c
        do z=nz+1,nzt
        z1=z1+1
        ntt=((z1-1)/nt)*nt
c----------------------------------------------------------------------|
c
        do p=Psi1,Psi2
        pm1=p-1
        if(pm1.le.0) pm1=pm1+np
c
c-mb    x-interpolation
        a1=ppx(r,p,w,z-nt)  -ppx(r,pm1,w,z-nt-1)
        a2=ppx(r,p,w,z-nt-1)-ppx(r,pm1,w,z-nt-2)
        b2=ppx(r,p,w,z-1)-ppx(r,pm1,w,z-2)
c
        ppx(r,p,w,z)=ppx(r,pm1,w,z-1)
     &               +a1!+a2 - b2
c
c-mb    y-interpolation
        a1=ppy(r,p,w,z-nt)  -ppy(r,pm1,w,z-nt-1)
        a2=ppy(r,p,w,z-nt-1)-ppy(r,pm1,w,z-nt-2)
        b2=ppy(r,p,w,z-1)-ppy(r,pm1,w,z-2)
c
        ppy(r,p,w,z)=ppy(r,pm1,w,z-1)
     &               +a1!+a2 - b2
c
c-mb    z-interpolation
        a1=ppz(r,p,w,z-nt)  -ppz(r,pm1,w,z-nt-1)
        a2=ppz(r,p,w,z-nt-1)-ppz(r,pm1,w,z-nt-2)
        b2=ppz(r,p,w,z-1)-ppz(r,pm1,w,z-2)
c
        ppz(r,p,w,z)=ppz(r,pm1,w,z-1)
     &               +a1!+a2 - b2
c
        end do       !p-loop
c----------------------------------------------------------------------|
c
        end do       !z-loop
c----------------------------------------------------------------------|
c
       end do        !w-loop
c----------------------------------------------------------------------|
c
      end do         !r-loop
c----------------------------------------------------------------------|
c
      return
      end
