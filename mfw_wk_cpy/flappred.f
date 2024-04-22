C-$	$Id: flappred.f,v 1.9.2.3 2005/06/11 15:57:37 shreyas Exp $	

      subroutine flappred(nr,np,ns,nw,pw,nz,nk,nb,ip,dp,dz,den,flph,om
     $     ,rad, compress,nubeta,bmass,ibeta,betap,kbeta,crd,clairfoil
     $     ,claoa, cttol,fltol,seg,tw,bicm,nwicm,bcx,bcy,twt,cdairfoil
     $     ,kinv, sonic,aoai,aoain,Cnni,Cnqi,vr,vry,gb,anb,rv,dcy,gv,asr
     $     ,lsr,rc0, sco,sce,gsc,csc,pcx,pcy,pcz,ppx,ppy,ppz,amom,thrust
     $     ,torque, inflow, cpx,cpy,cpz,taperst,vinfx,vinfy,vinfz,rr_x
     $     ,rr_y,rr_z, taper,t0r,t1c,t1s,b0r,b1c,b1s,ct,cqi,cqp,pif,zif
     $     ,theta, bvx,bvy,bvz,rcb,lengthnw,psi,beta,betastar, eqtim,
     $     Cldata, ClMsqr,liftx,lifty,liftz,alpha0)

C-$   Solve blade flapping response using time-accurate predictor
C-$   corrector integration. This subroutine performs the predictor step


      IMPLICIT none


      include 'adimpar.inc'
      include 'gdimpar.inc'
      include 'jdimpar.inc'
      include 'kdimpar.inc'
      include 'rdimpar.inc'
      include 'sdimpar.inc'
      include 'wdimpar.inc'
      include 'zdimpar.inc'


      INTEGER r,nr,p,np,s,ns,nb,j,c,nw,pw,nz,nk,w,z

      DOUBLE PRECISION aza,ip(rdim),dp,dz,aoac,den,flph,om,rad,nubeta,
     $     sonic,kinv,lin,gav,dcy, bmass,ibeta,betap,kbeta,crd,cttol
     $     ,fltol,pi,rtd,anb

      DOUBLE PRECISION t0r(rdim),t1c(rdim),t1s(rdim), b0r(rdim),b1c(rdim
     $     ),b1s(rdim), b0r_0(rdim),b1c_0(rdim),b1s_0(rdim), ct(rdim)
     $     ,ct_0(rdim),ect(rdim),efl(rdim), err(jdim),crr(jdim),seg(sdim
     $     )

      DOUBLE PRECISION er0(rdim,jdim),cr0(rdim,jdim), tw(rdim,sdim)
     $     ,nwicm(sdim,sdim),bcx(rdim,sdim), bcy(rdim,sdim)

      DOUBLE PRECISION bicm(sdim,sdim)

      DOUBLE PRECISION aoae(rdim,adim,sdim),aoai(rdim,adim,sdim),
     $     wlfv(rdim,adim,sdim),vr(rdim,adim,sdim), Cnni(rdim,adim,sdim)
     $     ,Cnqi(rdim,adim,sdim), gb(rdim,adim,sdim),cpx(rdim,adim,sdim)
     $     , cpy(rdim,adim,sdim),cpz(rdim,adim,sdim), aoanc(rdim,adim
     $     ,sdim)

      CHARACTER*1 trimmed,WantTrim

c-mb to include ideal twist
      DOUBLE PRECISION twt,alpha0

      INTEGER pif,zif

      DOUBLE PRECISION vry(rdim,adim,sdim),vrz(rdim,adim,sdim)
      DOUBLE PRECISION cqi(rdim)
      DOUBLE PRECISION cqp(rdim)
      DOUBLE PRECISION aflp(rdim,adim,sdim),apit(rdim,adim,sdim)
      DOUBLE PRECISION aoain(rdim,adim,sdim)

      INTEGER s75

      DOUBLE PRECISION clairfoil(rdim,adim,sdim), cdairfoil(rdim,adim
     $     ,sdim), claoa(rdim,adim,sdim),compress(rdim,adim,sdim)
      DOUBLE PRECISION Mach(rdim,adim,sdim),Reynolds(rdim,adim,sdim)

      DOUBLE PRECISION planf
      DOUBLE PRECISION taper,taperst

      DOUBLE PRECISION Timex,maxChange

      DOUBLE PRECISION gbold(rdim,adim,sdim)
      DOUBLE PRECISION t0rold(rdim),t1cold(rdim),t1sold(rdim)

      INTEGER n,itmp,dnp,psi
      
      DOUBLE PRECISION bvx(rdim,adim,sdim), bvy(rdim,adim,sdim),
     $     bvz(rdim,adim,sdim)
      DOUBLE PRECISION rcb

      INTEGER lengthnw(sdim)

      DOUBLE PRECISION gv(rdim,adim,wdim,zdim,kdim), rc0(rdim,adim,wdim
     $     ,kdim), pcx(rdim,adim,wdim,zdim), pcy(rdim,adim,wdim,zdim),
     $     pcz(rdim,adim,wdim,zdim), ppx(rdim,adim,wdim,zdim), ppy(rdim
     $     ,adim,wdim,zdim), ppz(rdim,adim,wdim,zdim), eqtim(rdim,adim
     $     ,wdim,zdim,kdim)

      DOUBLE PRECISION asr(rdim),rr_x(rdim),rr_y(rdim),rr_z(rdim),
     $                 lsr(rdim) 

      DOUBLE PRECISION vibx(rdim,adim,sdim), viby(rdim,adim,sdim),
     $     vibz(rdim,adim,sdim), vibxi(rdim,adim,sdim), vibyi(rdim,adim
     $     ,sdim), vibzi(rdim,adim,sdim), Cldata(rdim,adim,sdim),
     $     ClMsqr(rdim,adim,sdim), liftx(rdim,adim,sdim),lifty(rdim,adim
     $     ,sdim), liftz(rdim,adim,sdim),cpitch(rdim,adim)

      CHARACTER*1 sco

      DOUBLE PRECISION gsc(rdim,adim,wdim,zdim,kdim),
     &       csc(rdim,adim,wdim,zdim,kdim),sce

      DOUBLE PRECISION vinfx,vinfy,vinfz

      DOUBLE PRECISION amom(rdim,adim),thrust(rdim,adim)
     &                                ,inflow(rdim,adim)
     &                                ,torque(rdim,adim,2)

      INTEGER usaero
      DOUBLE PRECISION ub1,ub2,ua1,ua2,uf1,uf2,ug1,ug2

      INTEGER pm1
      DOUBLE PRECISION beta(rdim,adim),betastar(rdim,adim)
      DOUBLE PRECISION theta(rdim,adim)

      DOUBLE PRECISION sn,az,ze,wage,tt11,tt21,tt31,b0
      DOUBLE PRECISION rv(rdim,adim,wdim)

      DOUBLE PRECISION aflptmp,apittmp,det

      INTEGER pm2,pm3
      DOUBLE PRECISION betadot,thetadot

      DOUBLE PRECISION b1psi,b1stpsi,b1pm1,b1stpm1
      DOUBLE PRECISION am1pm1
      INTEGER nn,pm1dnp
      DOUBLE PRECISION b1cpsi,b1spsi,b1cstpsi,b1sstpsi
      DOUBLE PRECISION b1cpm1,b1spm1,b1cstpm1,b1sstpm1
      DOUBLE PRECISION am1cpm1,am1spm1
      INTEGER pdnp
      INTEGER rotgeo

      double precision ahdot(rdim,adim,sdim),alpdot(rdim,adim,sdim)
      double precision ualpha(rdim,adim,sdim),delalp(rdim,adim,sdim)
      double precision xcoeff(rdim,adim,sdim),ycoeff(rdim,adim,sdim)
      double precision uacirc(rdim,adim,sdim),uancir(rdim,adim,sdim)
      double precision secaoa(rdim,adim,sdim)

      external planf


      common/usadata/ub1,ub2,ua1,ua2,uf1,uf2,ug1,ug2,usaero
      save/usadata/


      common/pidata/pi
      save/pidata/

      common/unsdata/ualpha,delalp,ahdot,alpdot,xcoeff,ycoeff,uacirc,
     $     secaoa
      save/unsdata/

      common/rotdata/rotgeo
      save/rotdata/



      dnp=np/nb
      p=psi
      pm1=psi-1
      if(pm1.le.0)pm1=pm1+np

      call rindvtpsi(nr,np,nw,pw,nz,ns,nk,dp,dz,ip,anb,crd,rad, dcy
     $     ,kinv,om,lin,gav,gv,rc0,asr,beta,rr_x,rr_y,rr_z,pcx,pcy,pcz
     $     ,cpx,cpy,cpz, sco,sce,gsc,csc,vibx,viby,vibz, vibxi,vibyi
     $     ,vibzi,compress,lengthnw, pif,zif,bvx,bvy,bvz,rcb,gb,BICM
     $     ,NWICM,pm1, eqtim,lsr)

      do r=1,nr
         aza=DBLE(pm1-1)*dp/DBLE(pif)
         aoac=t0r(r)+t1c(r)*cos(aza)+t1s(r)*sin(aza)
         do s=1,ns
            if(abs(twt).gt.100) aoac=t0r(r)/(bcx(r,s)/rad)+t1c(r)
     $           *cos(aza)+t1s(r)*sin(aza)
            call aoapsi(r,pm1,s,cpx,cpy,cpz,vinfx,vinfy,vinfz, vibx,viby
     $           ,vibz,vibxi,vibyi,vibzi, nr,np,ns,dp,ip,b0r,b1c,b1s,t1c
     $           ,t1s,asr,om,bcx,bcy, aoac,tw,aza,pif,beta,theta,crd,
     $           aoai,aoain,aflp,apit,vr,vry,vrz,aoae,aoanc,lsr)
            aoae(r,pm1,s)=aoae(r,pm1,s)-aoai(r,pm1,s)+aoain(r,pm1,s)
            aoae(r,pm1,s)=aoae(r,pm1,s)-xcoeff(r,pm1,s)-ycoeff(r,pm1,s)
            
            Mach(r,pm1,s)=abs(vry(r,pm1,s))/sonic
            Reynolds(r,pm1,s)=vr(r,pm1,s) *crd*planf(bcx(r,s)/rad
     $           ,taperst,taper)/kinv

         end do                 !s
      end do                    !r

      do r=1,nr
         call airfoilpsi(bcx,aoae,Mach,Reynolds,nr,np,ns, clairfoil
     $        ,cdairfoil,claoa,compress,pm1,0,alpha0)

         call  rsppsi(r,np,ns,nb,ip,dp,den,flph,om,rad,crd,clairfoil,
     $        claoa, nubeta,bmass,ibeta,betap,kbeta,seg,bcx,vr,vry,vrz,
     $        gb,Cnni,Cnqi, sonic,cdairfoil,aoae,aoai,aoain,compress,
     $        BICM,NWICM, Mach,Reynolds, inflow,aoanc,beta,betastar,
     $        taperst,taper,ct,cqi,cqp,pif,zif,amom,thrust,torque,pm1,
     $        Cldata, ClMsqr, liftx, lifty, liftz,cpitch) 


c-mb explicit (Fwd Euler) predictor
         if(flph.ge.0.0D+0*rad) then
            betastar(r,p)=betastar(r,pm1) + dp*( amom(r,pm1)-nubeta**2.D
     $           +0*beta(r,pm1) )
            beta(r,p)=beta(r,pm1) + dp*betastar(r,pm1)
         end if
         
         
c-    mb  this part handles the teetering/gimballed rotor
         if(flph.lt.0.D+0) then
            if(nb.eq.2) then
               if(p.gt.np-dnp) then
                  
                  pm1dnp=pm1-dnp
                  if(pm1dnp.le.0) pm1dnp=pm1dnp+np
                  b1stpm1=(betastar(r,pm1)-betastar(r,pm1dnp))/DBLE(nb)
                  b1pm1=  (beta(r,pm1)-beta(r,pm1dnp))/DBLE(nb)
                  am1pm1 =(amom(r,pm1)-amom(r,pm1dnp))/DBLE(nb)
                  
                  b1stpsi=b1stpm1 + dp*( am1pm1 - nubeta**2.D+0*b1pm1)
                  b1psi  =b1pm1   + dp*b1stpm1
                  
                  beta(r,p)    = b1psi+betap
                  beta(r,p-dnp)=-b1psi+betap
                  betastar(r,p)    = b1stpsi
                  betastar(r,p-dnp)=-b1stpsi
                  
               end if           ! p > Np- dNp (second blade)
            end if              ! nb=2
            
            if(nb.ge.3) then
               if(p.gt.np-dnp) then
                  
                  b1cpm1  =0.D+0
                  b1spm1  =0.D+0
                  b1cstpm1=0.D+0
                  b1sstpm1=0.D+0
                  am1cpm1 =0.D+0
                  am1spm1 =0.D+0
                  do nn=1,nb
                     pm1=p-(nb-nn)*dnp-1
                     if(pm1.le.0.D+0)pm1=pm1+np
                     aza=ip(r)+DBLE(pm1-1)*dp/DBLE(pif) - ip(r)
                     b1cpm1  =b1cpm1  +cos(aza)*beta(r,pm1)
                     b1cstpm1=b1cstpm1+cos(aza)*betastar(r,pm1)
                     am1cpm1 =am1cpm1 +cos(aza)*amom(r,pm1)
                     b1spm1  =b1spm1  +sin(aza)*beta(r,pm1)
                     b1sstpm1=b1sstpm1+sin(aza)*betastar(r,pm1)
                     am1spm1 =am1spm1 +sin(aza)*amom(r,pm1)
                  end do
                  b1cpm1  =b1cpm1  *2.D+0/DBLE(nb)
                  b1cstpm1=b1cstpm1*2.D+0/DBLE(nb)
                  am1cpm1 =am1cpm1 *2.D+0/DBLE(nb)
                  b1spm1  =b1spm1  *2.D+0/DBLE(nb)
                  b1sstpm1=b1sstpm1*2.D+0/DBLE(nb)
                  am1spm1 =am1spm1 *2.D+0/DBLE(nb)
                  
                  b1cstpm1=b1cstpm1 - b1spm1
                  b1sstpm1=b1sstpm1 + b1cpm1
                  
                  b1cstpsi=b1cstpm1+dp*(am1cpm1 - 2.D+0*b1sstpm1
     &                 -(nubeta**2.D+0 -1.D+0)*b1cpm1)
                  b1sstpsi=b1sstpm1+dp*(am1spm1 + 2.D+0*b1cstpm1
     &                 -(nubeta**2.D+0 -1.D+0)*b1spm1)
                  b1cpsi  =b1cpm1  +dp*b1cstpm1
                  b1spsi  =b1spm1  +dp*b1sstpm1
                  
                  
                  do nn=1,nb
                     pdnp=p-(nb-nn)*dnp
                     aza=ip(r)+DBLE(pdnp-1)*dp/DBLE(pif) - ip(r)
                     betastar(r,pdnp)=cos(aza)*(b1cstpsi+b1spsi)
     &                    +sin(aza)*(b1sstpsi-b1cpsi)
                     beta(r,pdnp)=cos(aza)*b1cpsi+sin(aza)*b1spsi+betap
                  end do 
                  
               end if           ! p > Np- dNp (last blade)
            end if              ! nb>=3
         end if                 ! flph < 0.
      end do                    ! r
      
      do r=1,nr
         sn=1.0D0
         if(rotgeo.eq.1 .or. rotgeo.eq.2)  sn=(-1.D+0)**(r-1)
         az=ip(r)+DBLE(p-1)*dp/DBLE(pif) - ip(r)
         b0=beta(r,p)
         do w=1,nw
            wage=sn*(az + ip(r))
            tt11=cos(-asr(r))*cos(wage)*cos(b0)
     $           -sn*sin(-asr(r))*sin(lsr(r))*sin(wage)*cos(b0) 
     $           +sin(-asr(r))*cos(lsr(r))*sin(b0)
            tt21= sn*cos(lsr(r))*sin(wage)*cos(b0)+sin(lsr(r))*sin(b0)
            tt31=-sin(-asr(r))*cos(wage)*cos(b0)
     $           -sn*cos(-asr(r))*sin(lsr(r))*sin(wage)*cos(b0)
     $           +cos(-asr(r))*cos(lsr(r))*sin(b0)
     
            ppx(r,p,w,1)=rr_x(r)+rv(r,p,w)*tt11
            ppy(r,p,w,1)=rr_y(r)+rv(r,p,w)*tt21
            ppz(r,p,w,1)=rr_z(r)+rv(r,p,w)*tt31
         end do                 !w
      end do                    !r

      return
      end
