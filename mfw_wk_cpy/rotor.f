C-$	$Id: rotor.f,v 1.3.2.3 2005/06/11 15:57:37 shreyas Exp $	

      subroutine rotor(np,nr,ns,pif,rad,flph,rcout,crd,twt,rada,ip,
     $     dp,den,ibeta,rr_x,rr_y,rr_z,asr0,asr,b0r,b1c,b1s,t0r,t1c,t1s
     $     ,spx ,spy,spz,seg,tw,taperst,taper,lock,bcx,bcy ,bcz,cpx,cpy
     $     ,cpz,lsr0,lsr)


      IMPLICIT none

      include 'adimpar.inc'
      include 'rdimpar.inc'
      include 'sdimpar.inc'

      INTEGER r,s,np,nr,ns,pif,rotgeo
      INTEGER nbtr

      DOUBLE PRECISION rad,flph0,flph,rcout,crd,twt,rada, ds,pi
     $     ,ip(rdim),dp,asr0,lsr0, lock,den,ibeta
     
      DOUBLE PRECISION trasr0,trlsr0,trrad,trflph,trrcout,trcrd,
     $     trtwt,trom,trtaperst,trtaper 

      DOUBLE PRECISION rr_x(rdim),rr_y(rdim),rr_z(rdim), asr(rdim)
     $     ,b0r(rdim), b1c(rdim),b1s(rdim),spx(sdim),spy(sdim),spz(sdim)
     $     , seg(sdim), lsr(rdim)

      DOUBLE PRECISION tw(rdim,sdim), bcx(rdim,sdim),bcy(rdim,sdim)
     $     ,bcz(rdim,sdim)

      DOUBLE PRECISION cpx(rdim,adim,sdim), cpy(rdim,adim,sdim),
     $     cpz(rdim,adim,sdim)

      DOUBLE PRECISION planf
      DOUBLE PRECISION taper,taperst, temp1,temp2,temp3

      external planf

      DOUBLE PRECISION t0r(rdim),t1c(rdim),t1s(rdim)
      
C/MR >>  Tail Rotor geometry: <<
      namelist/tailr/nbtr,trasr0,trlsr0,trrad,trflph, trrcout,trcrd,
     $     trtwt,trom,trtaperst,trtaper

      common/pidata/pi
      save/pidata/

      common/rotdata/rotgeo
      save/rotdata/

      flph0=flph
      if(flph.lt.0.D+0)flph=0.D+0
      flph=flph*rad
      rcout=rcout*rad
      lock=(crd*den*2.D+0*pi*rad**4)/ibeta !note Clalpha=2pi

      rr_x(1)= 0.0D+0*rad
      rr_y(1)= 0.0D+0*rad
      rr_z(1)= 0.1864D+0*rad

      if(nr.eq.2) then
         rr_x(nr)= 0.0D+0*rad
         rr_y(nr)= 0.0D+0*rad
         rr_z(nr)= 0.0D+0*rad
      endif

      if(rotgeo.eq.1) then     !Simulate ground effect
         rr_z(1)=1.5D+0*rad
         rr_z(nr)=-1.50D+0*rad
      endif
      
C/MR  For main rotor/tail rotor configuration      
      if(nr.eq.2 .and. rotgeo.eq.3) then
      
C/MR  Read the tail rotor geometry/configuration data
        open(54,file='tailrotor.input',status='old')
        read(54,tailr)
        close(54)  
    
     
C/MR  To start set reference in the main rotor hub      
         rr_x(1)= 0.0D+0*rad
         rr_y(1)= 0.0D+0*rad
         rr_z(1)= 0.0D+0*rad
C/MR  If the c.g. is used instead (for UH60):
c         rr_x(1)= -0.059D+0*rad
c         rr_y(1)= 0.0D+0*rad
c         rr_z(1)= 0.229D+0*rad

C/MR  Tail rotor coords relative to hub (hard-coded for UH60):
         rr_x(nr)= 1.241D+0*rad
         rr_y(nr)= -0.0445D+0*rad
         rr_z(nr)= 0.0308D+0*rad
C/MR  If the reference is the c.g. (for UH60)                  
c         rr_x(nr)= 1.182D+0*rad
c         rr_y(nr)= -0.0445D+0*rad
c         rr_z(nr)= 0.259D+0*rad      
      
      endif

      rada=rad-rcout
      ds=(rad-rcout)/DBLE(ns)
c     Blade segment end points
         spx(1)=rcout
         spy(1)=0.D+0
         spz(1)=0.D+0
      do s=2,ns+1
         spx(s)=spx(s-1)+ds
c         spx(s)=0.5D0*(rad+rcout)+0.5D0*(rad-rcout)*cos(pi/DBLE(ns
c     $        )*DBLE(ns+1-s))
         spy(s)=0.D+0
         spz(s)=0.D+0
      end do
      

c     Shaft tilt angles
      do r=1,nr
         asr(r)=asr0*pi/180.D+0
         lsr(r)=lsr0*pi/180.D+0
         if (r .eq. 2 .and. rotgeo .eq. 1) asr(r)=pi-asr0*pi/180.0D0
         if (r .eq. 2 .and. rotgeo .eq. 3) then
            asr(r)=trasr0*pi/180.D+0
            lsr(r)=trlsr0*pi/180.D+0
         endif
         tw(r,1)=0.D+0
         do s=1,ns+1
            if(abs(twt).gt.100) then
               tw(r,s)=0        !ideal
            else
               tw(r,s)=(twt*pi/180.D+0)*DBLE(s-1)/DBLE(ns) !linear
            endif
         end do
      end do
      
      if (rotgeo .eq. 2) then   !MTR coaxial rotor twist fix
         open(33,file='mtr_twist.data',status='old')
         do s=1,ns
            read(33,*) temp1, temp2,temp3
            tw(nr-1,s)=temp2
            tw(nr,s)=temp3
c            write(*,*) s, tw(nr-1,s), tw(nr,s)
         enddo
         close(33)
      endif

c     Blade segment lengths
      do s=1,ns
         seg(s)=spx(s+1)-spx(s)
      end do
      
c     Blade control points
      do r=1,nr
         do s=1,ns
            bcx(r,s)=(spx(s)+spx(s+1))*0.5D+0
            bcy(r,s)=-crd*planf(bcx(r,s)/rad,taperst,taper)*0.5D+0
            bcz(r,s)=0.D+0
         end do
      end do                    !r

c     Blade control points, fixed frame coordinates
      call bcp(np,nr,ns,pif,ip,dp,asr,lsr,b0r,b1c,b1s,rad,t0r,t1c,t1s
     $     ,twt,tw, rr_x,rr_y,rr_z,bcx,bcy,bcz, cpx,cpy,cpz)

      flph=flph0

      return
      end
