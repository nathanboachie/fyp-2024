C-$	$Id: bcppsi.f,v 1.1.1.1.6.1 2005/03/10 18:49:33 shreyas Exp $	

      subroutine bcppsi(np,nr,ns,pif,ip,dp,asr,lsr,b0r,b1c,b1s,rad,om
     $     ,t0r,t1c,t1s,twt,tw,beta,theta, rr_x,rr_y,rr_z,bcx0,bcy0
     $     ,bcz0,cpx,cpy,cpz,p)

C-$   compute blade control points for each azimuth location


      IMPLICIT none

      include 'adimpar.inc'
      include 'rdimpar.inc'
      include 'sdimpar.inc'

      INTEGER r,nr,s,ns,p,np,pif

      DOUBLE PRECISION az,ip(rdim),dp,sn,pi,b0,rad

      DOUBLE PRECISION asr(rdim),b0r(rdim),b1c(rdim),b1s(rdim),
     $     rr_x(rdim),rr_y(rdim),rr_z(rdim),lsr(rdim)

      DOUBLE PRECISION bcx(rdim,sdim),bcy(rdim,sdim),bcz(rdim,sdim)

      DOUBLE PRECISION cpx(rdim,adim,sdim), cpy(rdim,adim,sdim),
     $     cpz(rdim,adim,sdim)

      DOUBLE PRECISION bcx0(rdim,sdim),bcy0(rdim,sdim),bcz0(rdim,sdim)

      DOUBLE PRECISION t0r(rdim),t1c(rdim),t1s(rdim)
      DOUBLE PRECISION tw(rdim,sdim)

      DOUBLE PRECISION om,twt,aza,aoac,pitch

      DOUBLE PRECISION beta(rdim,adim)
      DOUBLE PRECISION theta(rdim,adim)

      INTEGER rotgeo

      common/pidata/pi
      save/pidata/

      common/rotdata/rotgeo
      save/rotdata/

      do r=1,nr
         aza=DBLE(p-1)*dp/DBLE(pif)
         aoac=theta(r,p)
         
         do s=1,ns
            if(abs(twt).gt.100) aoac=t0r(r)/(bcx0(r,s)/rad)+t1c(r)
     $           *cos(aza)+t1s(r)*sin(aza)

            pitch=aoac+tw(r,s)

            bcx(r,s)=bcx0(r,s)
            bcy(r,s)= bcy0(r,s)*cos(pitch)+bcz0(r,s)*sin(pitch)
            bcz(r,s)=-bcy0(r,s)*sin(pitch)+bcz0(r,s)*cos(pitch)

         end do                 !s
      end do                    !r


      do r=1,nr
         sn=1.D+0
         if(rotgeo.eq.1 .or. rotgeo.eq.2) sn=(-1.D+0)**(r-1)
         az=(DBLE(p-1)*dp/DBLE(pif))

         b0=beta(r,p)

         az=(az + ip(r))*sn
         do s=1,ns
C           cpx(r,p,s)=rr_x(r)+ bcx(r,s)*(cos(b0)*cos(az)*cos(-asr(r))+
C    $           sin(b0)*sin(-asr(r)))+ (sn*bcy(r,s))*(-cos(-asr(r))
C    $           *sin(az))
C           cpy(r,p,s)=rr_y(r)+ bcx(r,s)*(sin(az)*cos(b0))+ (sn*bcy(r,s)
C    $           )*cos(az)
C           cpz(r,p,s)=rr_z(r)+ bcx(r,s)*(-sin(-asr(r))*cos(az)*cos(b0)+
C    $           cos(-asr(r))*sin(b0))+ (sn*bcy(r,s))*(sin(-asr(r))
C    $           *sin(az))

C/MR  Adding the lateral shaft tilt.
C/MR  Why is the cpz term not included?      
               cpx(r,p,s)=rr_x(r)+ bcx(r,s)*(cos(b0)*cos(az)*cos(-asr(r)
     $              )+sin(b0)*sin(-asr(r))*cos(lsr(r))-sn*sin(-asr(r))*
     $              sin(lsr(r))*sin(az)*cos(b0)) 
     $              +bcy(r,s)*(-cos(-asr(r))*sin(az)-sn*sin(-asr(r))*
     $              sin(lsr(r))*cos(az))
     $              +bcz(r,s)*(-cos(-asr(r))*cos(az)*sin(b0)
     $              +sn*sin(-asr(r))*sin(lsr(r))*sin(az)*sin(b0)
     $              +sin(-asr(r))*cos(lsr(r))*cos(b0))
               cpy(r,p,s)=rr_y(r)+ bcx(r,s)*(sn*cos(lsr(r))*sin(az)
     $                   *cos(b0)+sin(lsr(r))*sin(b0))
     $                   +bcy(r,s)*sn*cos(lsr(r))*cos(az)
     $                   +bcz(r,s)*(-sn*cos(lsr(r))*sin(az)*sin(b0)
     $                   +sin(lsr(r))*cos(b0))
               cpz(r,p,s)=rr_z(r)+ bcx(r,s)*(-sin(-asr(r))*cos(az)
     $              *cos(b0)+cos(-asr(r))*cos(lsr(r))*sin(b0)
     $              -sn*cos(-asr(r))*sin(lsr(r))*sin(az)*cos(b0))
     $              +bcy(r,s)*(sin(-asr(r))*sin(az)-sn*cos(-asr(r))
     $              *sin(lsr(r))*cos(az))
     $              +bcz(r,s)*(sin(-asr(r))*cos(az)*sin(b0)
     $              +sn*cos(-asr(r))*sin(lsr(r))*sin(az)*sin(b0)
     $              +cos(-asr(r))*cos(lsr(r))*cos(b0))

         end do                 !s
      end do                    !r


      return
      end
