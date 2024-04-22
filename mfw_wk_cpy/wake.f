C-$	$Id: wake.f,v 1.1.1.1.6.3 2005/06/11 15:57:37 shreyas Exp $	

      subroutine wake(ft,bct,dp,dz,mu,muc,nb,zif,pif, nt,ip,nr,np,nz,nzt
     $     ,anb)

      IMPLICIT none

      include 'adimpar.inc'
      include 'gdimpar.inc'
      include 'kdimpar.inc'
      include 'rdimpar.inc'
      include 'sdimpar.inc'
      include 'wdimpar.inc'
      include 'zdimpar.inc'

      INTEGER pif,zif,nb,nr,np,nz,nzt

      DOUBLE PRECISION dp,dz,pi,anb,ft,bct,nt,mu,muc,ip(rdim),vn


      common/pidata/pi
      save/pidata/

      common/vndata/vn
      save/vndata/


      dp=dp*pi/180.D+0
      dz=dz*pi/180.D+0
      zif=1
      pif=1
      if (dz .gt. dp) zif=nint(dz/dp)
      if (dp .gt. dz) pif=nint(dp/dz)

      anb=2.D+0*pi/DBLE(nb)
      np=nint(2.D+0*pi/dp)*pif

      nt=ft+bct
      ip(1)=nt*2.D+0*pi
      ip(1)=pi/2.D+0
      ip(1)=0.D+0
      if(nr.eq.2) then
         ip(nr)=0.0D0
      end if


      nz=nint(ft*2.D+0*pi/dz)*zif+1
      nzt=nint(nt*2.D+0*pi/dz)*zif+1


      return
      end
