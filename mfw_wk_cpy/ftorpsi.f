c----------------------------------------------------------------------|
c    <Copyright (c) 1999 by MrB>
c----------------------------------------------------------------------|
c
      subroutine ftorpsi(r,p,s,nr,np,ns,dp,ip,asr,lsr,beta,
     &                vfx,vfy,vfz,vrx,vry,vrz,pif)
c     Fixed to rotating frame transformation for velocity
c----------------------------------------------------------------------|
c
      IMPLICIT none
c
        include 'adimpar.inc'
        include 'kdimpar.inc'
        include 'rdimpar.inc'
        include 'sdimpar.inc'
c
      INTEGER nr,np,ns,r,p,s
c
      DOUBLE PRECISION sn,aza,ip(rdim),dp,
     &       fr11,fr12,fr13,
     &       fr21,fr22,fr23,
     &       fr31,fr32,fr33,
     &       b0
c
      DOUBLE PRECISION asr(rdim),beta(rdim,adim),lsr(rdim)
c
      DOUBLE PRECISION vrx(rdim,adim,sdim),
     &       vry(rdim,adim,sdim),
     &       vrz(rdim,adim,sdim),
     &       vfx(rdim,adim,sdim),
     &       vfy(rdim,adim,sdim),
     &       vfz(rdim,adim,sdim)
c
      INTEGER pif
      INTEGER rotgeo
c
      common/rotdata/rotgeo
      save/rotdata/

c----------------------------------------------------------------------
c
      sn=1.D+0
      if(rotgeo.eq.1 .or. rotgeo.eq.2) sn=(-1.D+0)**(r-1)

c-mb VFF correctn
      aza=( ip(r)+DBLE(p-1)*dp/DBLE(pif) ) - ip(r)
c-
      b0=beta(r,p)
c-
      aza=(aza + ip(r))*sn
c-mb sign of rotation for coord transform
c
C/MR  Adding lateral shaft tilt angle      
      fr11= cos(-asr(r))*cos(aza)*cos(b0)
     $     -sn*sin(-asr(r))*sin(lsr(r))*sin(aza)*cos(b0) 
     $     +sin(-asr(r))*cos(lsr(r))*sin(b0)
      fr12= sn*cos(lsr(r))*sin(aza)*cos(b0)+sin(lsr(r))*sin(b0)
      fr13=-sin(-asr(r))*cos(aza)*cos(b0)
     $     -sn*cos(-asr(r))*sin(lsr(r))*sin(aza)*cos(b0)
     $     +cos(-asr(r))*cos(lsr(r))*sin(b0)

      fr21=-cos(-asr(r))*sin(aza)-sn*sin(asr(r))*sin(lsr(r))*cos(aza)
      fr22= cos(lsr(r))*cos(aza)*sn
      fr23= sin(-asr(r))*sin(aza)-sn*cos(asr(r))*sin(lsr(r))*cos(aza)

      fr31=-cos(-asr(r))*cos(aza)*sin(b0)
     $      +sn*sin(-asr(r))*sin(lsr(r))*sin(aza)*sin(b0) 
     $      +sin(-asr(r))*cos(lsr(r))*cos(b0)
      fr32=-sn*cos(lsr(r))*sin(aza)*sin(b0)+sin(lsr(r))*cos(b0)
      fr33= sin(-asr(r))*cos(aza)*sin(b0)
     $      +sn*cos(asr(r))*sin(lsr(r))*sin(aza)*sin(b0) 
     $      +cos(-asr(r))*cos(lsr(r))*cos(b0)

      vrx(r,p,s)=   fr11*vfx(r,p,s)+fr12*vfy(r,p,s)+fr13*vfz(r,p,s)
      vry(r,p,s)=+( fr21*vfx(r,p,s)+fr22*vfy(r,p,s)+fr23*vfz(r,p,s) )
      vrz(r,p,s)=   fr31*vfx(r,p,s)+fr32*vfy(r,p,s)+fr33*vfz(r,p,s)
c----------------------------------------------------------------------|
c
      return
      end
