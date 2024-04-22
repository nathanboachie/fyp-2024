C-$	$Id: ftor.f,v 1.1.1.1.6.1 2005/03/10 18:49:33 shreyas Exp $	

      subroutine ftor(r,p,s,dp,ip,b0r,asr,lsr,b1c,b1s,vfx,vfy,vfz 
     $     ,vrx,vry,vrz,pif)

c >>  Fixed to rotating frame transformation. <<


      IMPLICIT none

        include 'adimpar.inc'
        include 'kdimpar.inc'
        include 'rdimpar.inc'
        include 'sdimpar.inc'

      INTEGER r,p,s

      DOUBLE PRECISION sn,aza,ip(rdim),dp, fr11,fr12,fr13, fr21,fr22
     $     ,fr23, fr31,fr32,fr33, b0

      DOUBLE PRECISION asr(rdim),lsr(rdim),b0r(rdim),b1c(rdim),b1s(rdim)

      DOUBLE PRECISION vrx(rdim,adim,sdim), vry(rdim,adim,sdim),
     $     vrz(rdim,adim,sdim), vfx(rdim,adim,sdim), vfy(rdim,adim,sdim)
     $     , vfz(rdim,adim,sdim)

      INTEGER pif
      INTEGER rotgeo

      common/rotdata/rotgeo
      save/rotdata/

      sn=1.D+0
      if(rotgeo.eq.1 .or. rotgeo.eq.2) sn=(-1.D+0)**(r-1)

      aza=DBLE(p-1)*dp/DBLE(pif)

      b0=b0r(r)+b1c(r)*cos(aza)+b1s(r)*sin(aza)
      aza=(aza + ip(r))*sn

c      fr11= cos(-asr(r))*cos(aza)*cos(b0)+ sin(-asr(r))*sin(b0)
c      fr12= sin(aza)*cos(b0)
c      fr13=-sin(-asr(r))*cos(aza)*cos(b0)+ cos(-asr(r))*sin(b0)

c      fr21=-cos(-asr(r))*sin(aza)*sn
c      fr22= cos(aza)*sn
c      fr23= sin(-asr(r))*sin(aza)*sn

c      fr31=-cos(-asr(r))*cos(aza)*sin(b0)+ sin(-asr(r))*cos(b0)
c      fr32=-sin(aza)*sin(b0)
c      fr33= sin(-asr(r))*cos(aza)*sin(b0)+ cos(-asr(r))*cos(b0)


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


      return
      end
