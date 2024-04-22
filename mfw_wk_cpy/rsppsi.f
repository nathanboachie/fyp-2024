C-$	$Id: rsppsi.f,v 1.11.2.2 2005/03/23 21:24:38 shreyas Exp $	

      subroutine rsppsi(r,np,ns,nb,ip,dp,den,flph,om,rad,crd,clairfoil
     $     ,claoa, nbe,bmass,ibe,bep,kbe,seg,bcx,vr,vry,vrz,gb,Cnni,Cnqi
     $     ,sonic, cdairfoil,aoae,aoai,aoain,compress,BICM,NWICM, Mach
     $     ,Reynolds, inflow,aoanc,beta,betastar, taperst,taper,ct,cqi
     $     ,cqp,pif,zif,amom,thrust,torque,p, Cldata, ClMsqr, liftx,
     $     lifty, liftz,cpitch)

C-$   Compute the blade flapping response for each azimuth position

      IMPLICIT none


      include 'adimpar.inc'
      include 'gdimpar.inc'
      include 'kdimpar.inc'
      include 'rdimpar.inc'
      include 'sdimpar.inc'
      include 'wdimpar.inc'
      include 'zdimpar.inc'


      INTEGER r,p,s,np,ns,nb,pm1

      DOUBLE PRECISION am0,am1c,am1s,aza,ip(rdim),dp, lift,den,flph,pi
     $     ,om,rad,nbe,bmass,ibe, bep,kbe,crd,pbar,qbar,sonic

      DOUBLE PRECISION seg(sdim),ct(rdim)

      DOUBLE PRECISION bcx(rdim,sdim)

      DOUBLE PRECISION vr(rdim,adim,sdim),gb(rdim,adim,sdim), Cnni(rdim
     $     ,adim,sdim),Cnqi(rdim,adim,sdim), clprm(rdim,adim,sdim)
      DOUBLE PRECISION Mach(rdim,adim,sdim),Reynolds(rdim,adim,sdim)

      INTEGER pif,zif

      DOUBLE PRECISION cdprofile,drag
      DOUBLE PRECISION aoae(rdim,adim,sdim),aoai(rdim,adim,sdim)
      DOUBLE PRECISION aoain(rdim,adim,sdim), aoanc(rdim,adim,sdim)
      DOUBLE PRECISION cqi(rdim),cqp(rdim),vry(rdim,adim,sdim)
      DOUBLE PRECISION vrz(rdim,adim,sdim)

      DOUBLE PRECISION planf
      DOUBLE PRECISION taper,taperst

      DOUBLE PRECISION BICM(sdim,sdim),NWICM(sdim,sdim)
      DOUBLE PRECISION phi,Vinf

      DOUBLE PRECISION clairfoil(rdim,adim,sdim), cdairfoil(rdim,adim
     $     ,sdim), claoa(rdim,adim,sdim),compress(rdim,adim,sdim),
     $     Cldata(rdim,adim,sdim), ClMsqr(rdim,adim,sdim), liftx(rdim
     $     ,adim,sdim), lifty(rdim,adim,sdim), liftz(rdim,adim,sdim)

      DOUBLE PRECISION amom(rdim,adim),thrust(rdim,adim) ,inflow(rdim
     $     ,adim) ,torque(rdim,adim,2)

      DOUBLE PRECISION beta(rdim,adim),betastar(rdim,adim)
      DOUBLE PRECISION Cl
      DOUBLE PRECISION normlt,tanglt
      double precision cpitch(rdim,adim), thmphi
      double precision cnplag(rdim,adim,sdim), deltalift(rdim,adim,sdim)
     $     , epressure, tpressure,deltat, cldash

      external planf


      common/pidata/pi
      save/pidata/

      common/ssmdata/pbar,qbar
      save/ssmdata/
      
      common/dynstall/cnplag,deltalift
      save/dynstall/

      Amom(r,p)=0.D+0
      deltat=dp/om
      tpressure=2.5D+0
      
      pm1=p-1
      if(pm1 .le. 0) pm1=pm1+np
      
      amom(r,p)=0.D+0
      thrust(r,p)=0.D+0
      inflow(r,p)=0.D+0
      torque(r,p,1)=0.D+0
      torque(r,p,2)=0.D+0

      do s=1,ns
         Vinf=vr(r,p,s)
         phi=-aoain(r,p,s)*compress(r,p,s)
         
         epressure=exp(-deltat*2*Vinf/(crd*tpressure))
         cnplag(r,p,s)=clairfoil(r,p,s)+Cnni(r,p,s)+Cnqi(r,p,s)
         deltalift(r,p,s)=deltalift(r,pm1,s)*epressure+
     $        (cnplag(r,p,s)-cnplag(r,pm1,s))*sqrt(epressure)
         Cl=cnplag(r,p,s)-deltalift(r,p,s)
         cdprofile=cdairfoil(r,p,s)
         
         lift=0.5D+0*den*(Vinf**2)*crd*planf(bcx(r,s)/rad,taperst,
     $        taper)*Cl*seg(s)
         drag=0.5D+0*den*(Vinf**2)*crd*planf(bcx(r,s)/rad,taperst,
     $        taper)*cdprofile*seg(s)
         
         thmphi=cpitch(r,p)-phi

         normlt=lift*cos(thmphi)+drag*sin(thmphi)
         tanglt=-lift*sin(thmphi)+drag*cos(thmphi)
         liftx(r,p,s)=tanglt
         lifty(r,p,s)=0.0D+0
         liftz(r,p,s)=normlt
         
         thrust(r,p)=thrust(r,p)+lift*cos(phi)-drag*sin(phi)
         torque(r,p,1)=torque(r,p,1)+lift*sin(phi)*bcx(r,s)
         torque(r,p,2)=torque(r,p,2)+drag*cos(phi)*bcx(r,s)
         inflow(r,p)=inflow(r,p)-vrz(r,p,s)/om/rad
         
         Cldata(r,p,s)=Cl
c         ClMsqr(r,p,s)=Cl*Vinf**2/sonic**2
         ClMsqr(r,p,s)=cdairfoil(r,p,s)
         amom(r,p)=amom(r,p)+lift*(bcx(r,s)-flph)
      end do                    !s

      inflow(r,p)=inflow(r,p)/DBLE(ns)
      amom(r,p)=amom(r,p)/(ibe*om**2)

      return
      end
