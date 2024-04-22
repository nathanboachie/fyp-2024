C-$	$Id: rsp.f,v 1.12.2.2 2005/03/23 21:24:38 shreyas Exp $	

      subroutine rsp(r,np,ns,nb,ip,dp,den,flph,om,rad,crd,clairfoil
     $     ,claoa, nbe,bmass,ibe,bep,kbe,seg,bcx,vr,vry,vrz,gb,Cnni,Cnqi
     $     ,sonic, cdairfoil,aoae,aoai,aoain,compress,BICM,NWICM,mux,
     $     Mach,Reynolds, taperst,taper,b0r,b1c,b1s,ct,cqi,cqp,pif,zif
     $     ,itern)

C-$   Compute blade flapping response


      IMPLICIT none


      include 'adimpar.inc'
      include 'gdimpar.inc'
      include 'kdimpar.inc'
      include 'rdimpar.inc'
      include 'sdimpar.inc'
      include 'wdimpar.inc'
      include 'zdimpar.inc'


      INTEGER r,p,s,np,ns,nb,dnp

      DOUBLE PRECISION thrust,amom,am0,am1c,am1s,aza,ip(rdim),dp, lift
     $     ,den,flph0,flph,pi,om,rad,nbe,nbe2,bmass,ibe, bep,kbe,crd
     $     ,pbar,qbar,sonic,mux

      DOUBLE PRECISION seg(sdim),b0r(rdim),b1c(rdim),b1s(rdim),ct(rdim)

      DOUBLE PRECISION lock,det
      DOUBLE PRECISION lock1,lock2,lock3,lock4

      DOUBLE PRECISION bcx(rdim,sdim)

      DOUBLE PRECISION vr(rdim,adim,sdim),gb(rdim,adim,sdim), Cnni(rdim
     $     ,adim,sdim),Cnqi(rdim,adim,sdim)

      DOUBLE PRECISION P_thrust(rdim,adim)
      DOUBLE PRECISION P_torque(rdim,adim)
      DOUBLE PRECISION P_moment(rdim,adim)
      DOUBLE PRECISION P_inflow(rdim,adim)

      INTEGER pif,zif

      DOUBLE PRECISION itorque,ptorque,cdprofile,drag
      DOUBLE PRECISION aoae(rdim,adim,sdim),aoai(rdim,adim,sdim)
      DOUBLE PRECISION aoain(rdim,adim,sdim)
      DOUBLE PRECISION cqi(rdim),cqp(rdim),vry(rdim,adim,sdim)
      DOUBLE PRECISION vrz(rdim,adim,sdim)

      DOUBLE PRECISION planf
      DOUBLE PRECISION taper,taperst

      DOUBLE PRECISION BICM(sdim,sdim),NWICM(sdim,sdim)
      DOUBLE PRECISION phi,Vinf

      DOUBLE PRECISION clairfoil(rdim,adim,sdim), cdairfoil(rdim,adim
     $     ,sdim), compress(rdim,adim,sdim),claoa(rdim,adim,sdim)
      DOUBLE PRECISION Mach(rdim,adim,sdim),Reynolds(rdim,adim,sdim)

      DOUBLE PRECISION eR

      INTEGER itern

      DOUBLE PRECISION hforce,yforce,rollmom,pitchmom

      INTEGER pnb,nn
      DOUBLE PRECISION psib
      DOUBLE PRECISION flmom(rdim,adim),refmom(rdim,adim)
      DOUBLE PRECISION Cl,autorot,autorot1

      external planf


      common/pidata/pi
      save/pidata/

      common/ssmdata/pbar,qbar
      save/ssmdata/


      flph0=flph
      if(flph.lt.0.D+0)flph=0.D+0

      dnp=np/nb

      thrust=0.D+0
      itorque=0.D+0
      ptorque=0.D+0
      hforce=0.D+0
      yforce=0.D+0
      rollmom=0.D+0
      pitchmom=0.D+0

      am0=0.D+0
      am1c=0.D+0
      am1s=0.D+0

      open(99,file='lift.dat',status='unknown')
      do p=1,np
         P_thrust(r,p)=0.d0
         if(p.gt.dnp) P_thrust(r,p)=P_thrust(r,p-dnp)
         P_torque(r,p)=0.d0
         if(p.gt.dnp) P_torque(r,p)=P_torque(r,p-dnp)
         P_moment(r,p)=0.d0
         if(p.gt.dnp) P_moment(r,p)=P_moment(r,p-dnp)
         P_inflow(r,p)=0.d0
         if(p.gt.dnp) P_inflow(r,p)=P_inflow(r,p-dnp)

         amom=0.D+0
         aza=DBLE(p-1)*dp/DBLE(pif)


         do s=1,ns
            Vinf=vr(r,p,s)
            phi=-aoain(r,p,s)*compress(r,p,s)
            Cl=clairfoil(r,p,s)+Cnni(r,p,s)+Cnqi(r,p,s)
            cdprofile=cdairfoil(r,p,s)

            lift=0.5D+0*den*(Vinf**2)*crd*planf(bcx(r,s)/rad,taperst
     $           ,taper)*Cl*seg(s)
            drag=0.5D+0*den*(Vinf**2)*crd*planf(bcx(r,s)/rad,taperst
     $           ,taper)*cdprofile*seg(s)

            autorot=Cl*sin(phi)+cdprofile*cos(phi)
            autorot1=Cl*phi+cdprofile

            thrust=thrust+lift*cos(phi)-drag*sin(phi)
            itorque=itorque+lift*sin(phi)*bcx(r,s)
            ptorque=ptorque+drag*cos(phi)*bcx(r,s)
            amom=amom+(lift*cos(phi)-drag*sin(phi))*(bcx(r,s)-flph)
            hforce=hforce+drag*cos(phi)*sin(aza)
            yforce=yforce-drag*cos(phi)*cos(aza)
            rollmom=rollmom+lift*(bcx(r,s)-flph)*sin(aza)
            pitchmom=pitchmom+lift*(bcx(r,s)-flph)*cos(aza)

            write(99,98)bcx(r,s),aza,bcx(r,s)*cos(aza),-bcx(r,s)*sin(aza
     $           ),aoae(r,p,s)*180.0D+0/pi,2.D+0*pi*aoae(r,p,s)*Vinf**2
     $           /sonic**2,clairfoil(r,p,s),clairfoil(r,p,s)*Vinf**2
     $           /sonic**2,Mach(r,p,s),compress(r,p,s)


            P_thrust(r,p)=P_thrust(r,p)+lift*cos(phi)-drag*sin(phi)
            
            P_torque(r,p)=P_torque(r,p)+lift*sin(phi)*bcx(r,s)+drag
     $           *cos(phi)*bcx(r,s)
            P_moment(r,p)=P_moment(r,p)+lift*(bcx(r,s)-flph)
            
            P_inflow(r,p)=P_inflow(r,p)-vrz(r,p,s)/om/rad
         end do                 !s
         
         flmom(r,p)=amom
         write(99,*)
         write(99,*)
         
         am0 =am0 +0.5D+0*amom*(dp/DBLE(pif))/pi
         am1c=am1c+amom*cos(aza)*(dp/DBLE(pif))/pi
         am1s=am1s+amom*sin(aza)*(dp/DBLE(pif))/pi
      end do                    !p
      flph=flph0

      close(99)
98    format(10(E13.7,2X))
122   format(5(E13.7,2X))


      if(flph.lt.0.D+0) then

         do p=1,np
            refmom(r,p)=0.D+0
            do nn=1,nb
               pnb=p+(nn-1)*dnp
               if(pnb.gt.np)pnb=pnb-np
               psib=DBLE(pnb-1)*dp/DBLE(pif) 
               psib=DBLE(nn-1)*2.D+0*pi/DBLE(nb)
               refmom(r,p)=refmom(r,p)+flmom(r,pnb)*cos(psib)
            end do              !nn
            refmom(r,p)=refmom(r,p)*2.D+0/DBLE(nb)
            if(nb.eq.2)refmom(r,p)=refmom(r,p)/2.D+0
         end do                 !p
         
         am0=0.D+0
         am1c=0.D+0
         am1s=0.D+0
         do p=1,np
            aza=DBLE(p-1)*dp/DBLE(pif)
            am0 =am0 +refmom(r,p)*(dp/DBLE(pif))/2.D+0/pi
            am1c=am1c+refmom(r,p)*cos(aza)*(dp/DBLE(pif))/pi
            am1s=am1s+refmom(r,p)*sin(aza)*(dp/DBLE(pif))/pi
         end do                 !p
      end if                    !flph<0 for teeter/gimball


      thrust=thrust*DBLE(nb)/DBLE(np)
      ct(r)=thrust/(den*pi*(rad**4)*om**2)

      itorque=itorque*DBLE(nb)/DBLE(np)
      cqi(r)=itorque/(den*(pi*rad**2)*om**2*rad**2*rad)

      ptorque=ptorque*DBLE(nb)/DBLE(np)
      cqp(r)=ptorque/(den*(pi*rad**2)*om**2*rad**2*rad)

      hforce=hforce*DBLE(nb)/DBLE(np)
      hforce=hforce/(den*pi*(rad**4)*om**2)

      yforce=yforce*DBLE(nb)/DBLE(np)
      yforce=yforce/(den*pi*(rad**4)*om**2)

      rollmom=rollmom*DBLE(nb)/DBLE(np)
      rollmom=rollmom/(den*(pi*rad**2)*om**2*rad**2*rad)

      pitchmom=pitchmom*DBLE(nb)/DBLE(np)
      pitchmom=pitchmom/(den*(pi*rad**2)*om**2*rad**2*rad)


      b0r(r)=(1.D+0/nbe**2)*(am0/(ibe*om**2)+kbe*bep/ibe/om**2)

      if(flph.lt.0.D+0) flph=0.D+0
      nbe2=nbe**2.D+0
      eR=flph/rad
      lock=(crd*den*2.D+0*pi*rad**4)/ibe !note Clalpha=2pi
      lock=lock/sqrt(1.D+0-(0.75D+0*om*rad/sonic)**2.D+0) !general
      lock1=lock*( 1 - 3.D+0/2.D+0*eR + 1.D+0/2.D+0*eR**3.D+0
     $     )
      lock2=lock*(1-8.D+0/3.D+0*eR+2.D+0*eR**2.D+0-1.D+0/3.D+0*eR**4.D
     $     +0)
      lock3=lock*( 1 - 3.D+0/1.D+0*eR + 3.D+0*eR**2.D+0 - 1.D+0/1.D+0*eR
     $     **4.D+0            )
      lock4=lock*( 1 - 2.D+0*eR + eR**2.D+0)

      if(flph0.lt.0.D+0) then
        lock1=0.D+0
        lock3=0.D+0
        lock4=0.D+0
      end if

      det=(nbe2-1.D+0)**2.D+0 + (lock2/8.D+0)**2.D+0 - (mux**2.D+0*lock4
     $     /16.D+0)**2.D+0 - mux**2.D+0/nbe2*(nbe2-1.D+0)*lock1/72.D+0
     $     *(lock3-lock1)

      b1c(r)=(nbe2-1.D+0)*(am1c - lock1/6.D+0*mux/nbe2*am0) - (lock2/8.D
     $     +0 + lock4/16.D+0*mux**2.D+0)*am1s

      b1s(r)=(lock2/8.D+0 - lock4/16.D+0*mux**2.D+0)* (am1c - lock1/6.D
     $     +0*mux/nbe2*am0) + ( (nbe2-1.D+0) + mux**2.D+0/nbe2*lock1/72
     $     .D+0*(lock3-lock1) )*am1s

      b1c(r)=b1c(r)/(ibe*om**2.D+0)/det
      b1s(r)=b1s(r)/(ibe*om**2.D+0)/det

      b0r(r)=( (am0+kbe*bep)/(ibe*om**2.D+0) - mux/12.D+0*(lock1-lock3)
     $     *b1c(r) )/nbe2

      flph=flph0
      if(flph.lt.0.D+0)b0r(r)=bep

      return
      end
