C-$	$Id: rtof.f,v 1.1.1.1.6.1 2005/03/10 18:49:33 shreyas Exp $	

      subroutine rtof(nr,np,ns,pif,dp,ip,b0r,b1c,b1s,asr,lsr 
     $     ,rr_x,rr_y,rr_z,spx,spy,spz, bvx,bvy,bvz)

c     Rotating to fixed frame transformation.


      IMPLICIT none

      include 'adimpar.inc'
      include 'kdimpar.inc'
      include 'rdimpar.inc'
      include 'sdimpar.inc'

      INTEGER nr,np,ns,pif,r,p,s
      INTEGER rotgeo

      DOUBLE PRECISION sn,aza,ip(rdim),dp, tt11,tt12,tt13, tt21,tt22
     $     ,tt23, tt31,tt32,tt33, b0

      DOUBLE PRECISION spx(sdim),spy(sdim),spz(sdim), rr_x(rdim)
     $     ,rr_y(rdim),rr_z(rdim), asr(rdim),b0r(rdim),b1c(rdim)
     $     ,b1s(rdim), lsr(rdim)

      DOUBLE PRECISION bvx(rdim,adim,sdim), bvy(rdim,adim,sdim),
     $     bvz(rdim,adim,sdim)
     
      common/rotdata/rotgeo
      save/rotdata/


      sn=1.D+0
      do r=1,nr
         if(rotgeo.eq.1 .or. rotgeo.eq.2) sn=(-1.D+0)**(r-1)
         do p=1,np              !,pif
            aza=DBLE(p-1)*dp/DBLE(pif)
            b0=b0r(r)+b1c(r)*cos(aza)+b1s(r)*sin(aza)
            aza=(aza + ip(r))*sn

c            tt11= cos(-asr(r))*cos(aza)*cos(b0)+ sin(-asr(r))*sin(b0)
c            tt12=-sn*cos(-asr(r))*sin(aza)
c            tt13=-cos(-asr(r))*cos(aza)*sin(b0)+ sin(-asr(r))*cos(b0)

c            tt21= sin(aza)*cos(b0)
c            tt22= sn*cos(aza)
c            tt23= sin(-aza)*sin(b0)

c            tt31=-sin(-asr(r))*cos(aza)*cos(b0)+ cos(-asr(r))*sin(b0)
c            tt32= sn*sin(-asr(r))*sin(aza)
c            tt33= sin(-asr(r))*cos(aza)*sin(b0)+ cos(-asr(r))*cos(b0)

C/MR    Update with lateral shaft tilt
            tt11= cos(-asr(r))*cos(aza)*cos(b0)
     $           -sn*sin(-asr(r))*sin(lsr(r))*sin(aza)*cos(b0) 
     $           +sin(-asr(r))*cos(lsr(r))*sin(b0)
            tt12=-cos(-asr(r))*sin(aza)
     $       -sn*sin(asr(r))*sin(lsr(r))*cos(aza)
            tt13=-cos(-asr(r))*cos(aza)*sin(b0)
     $           +sn*sin(-asr(r))*sin(lsr(r))*sin(aza)*sin(b0) 
     $           +sin(-asr(r))*cos(lsr(r))*cos(b0)

            tt21= sn*cos(lsr(r))*sin(aza)*cos(b0)+sin(lsr(r))*sin(b0)
            tt22= cos(lsr(r))*cos(aza)*sn
            tt23= -sn*cos(lsr(r))*sin(aza)*sin(b0)+sin(lsr(r))*cos(b0)
            
            tt31=-sin(-asr(r))*cos(aza)*cos(b0)
     $           -sn*cos(-asr(r))*sin(lsr(r))*sin(aza)*cos(b0)
     $           +cos(-asr(r))*cos(lsr(r))*sin(b0)
            tt32= sin(-asr(r))*sin(aza)
     $           -sn*cos(asr(r))*sin(lsr(r))*cos(aza)
            tt33= sin(-asr(r))*cos(aza)*sin(b0)
     $           +sn*cos(asr(r))*sin(lsr(r))*sin(aza)*sin(b0) 
     $           +cos(-asr(r))*cos(lsr(r))*cos(b0)

            do s=1,ns+1
               bvx(r,p,s)=rr_x(r)+tt11*spx(s)+ tt12*spy(s)+ tt13
     $              *spz(s)
               bvy(r,p,s)=rr_y(r)+tt21*spx(s)+ tt22*spy(s)+ tt23
     $              *spz(s)
               bvz(r,p,s)=rr_z(r)+tt31*spx(s)+ tt32*spy(s)+ tt33
     $              *spz(s)
            end do

         end do
      end do

      return
      end
