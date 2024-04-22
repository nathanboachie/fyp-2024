C-$	$Id: flapcorr.f,v 1.15.2.3 2005/06/11 15:57:37 shreyas Exp $	

      subroutine flapcorr(nr,np,ns,nw,pw,nz,nk,nb,ip,dp,dz,den,flph,om
     $     ,rad, compress,nubeta,bmass,ibeta,betap,kbeta,crd,clairfoil
     $     ,claoa, cttol,fltol,seg,tw,bicm,nwicm,bcx,bcy,twt,cdairfoil
     $     ,kinv, sonic,aoai,aoain,Cnni,Cnqi,vr,vry,gb,anb,rv,dcy,gv,asr
     $     ,lsr,rc0, sco,sce,gsc,csc,pcx,pcy,pcz,ppx,ppy,ppz,amom,thrust
     $     ,torque, inflow, cpx,cpy,cpz,taperst,vinfx,vinfy,vinfz,rr_x
     $     ,rr_y,rr_z, taper,t0r,t1c,t1s,b0r,b1c,b1s,ct,cqi,cqp,pif,zif
     $     ,theta, bvx,bvy,bvz,rcb,lengthnw,psi,beta,betastar,Iteratn,
     $     eqtim, vibxi, vibyi, vibzi, Cldata, ClMsqr, liftx, lifty,
     $     liftz,alpha0)

C-$   Corrector step for the time integration of flap equations PC2B.


      IMPLICIT none


      include 'adimpar.inc'
      include 'gdimpar.inc'     
      include 'jdimpar.inc'
      include 'kdimpar.inc'
      include 'rdimpar.inc'
      include 'sdimpar.inc'
      include 'wdimpar.inc'     
      include 'zdimpar.inc'


      INTEGER r,nr,p,np,s,ns,nb,j,c,nw,pw,nz,nk,w,z,i

      DOUBLE PRECISION aza,ip(rdim),dp,dz,aoac,den,flph,om,rad,nubeta,
     &       sonic,kinv,lin,gav,dcy,
     &       bmass,ibeta,betap,kbeta,crd,cttol,fltol,pi,rtd,anb

      DOUBLE PRECISION t0r(rdim),t1c(rdim),t1s(rdim),
     &       b0r(rdim),b1c(rdim),b1s(rdim),
     &       b0r_0(rdim),b1c_0(rdim),b1s_0(rdim),
     &       ct(rdim),ct_0(rdim),ect(rdim),efl(rdim),
     &       err(jdim),crr(jdim),seg(sdim)

      DOUBLE PRECISION er0(rdim,jdim),cr0(rdim,jdim),
     &       tw(rdim,sdim),nwicm(sdim,sdim),bcx(rdim,sdim),
     &       bcy(rdim,sdim)

      DOUBLE PRECISION bicm(sdim,sdim)

      DOUBLE PRECISION aoae(rdim,adim,sdim),aoai(rdim,adim,sdim),
     &       wlfv(rdim,adim,sdim),vr(rdim,adim,sdim),
     &       Cnni(rdim,adim,sdim),Cnqi(rdim,adim,sdim),
     &       gb(rdim,adim,sdim),cpx(rdim,adim,sdim),
     &       cpy(rdim,adim,sdim),cpz(rdim,adim,sdim),
     &       aoanc(rdim,adim,sdim),alpha0

      CHARACTER*1 trimmed,WantTrim

c-mb to include ideal twist
      DOUBLE PRECISION twt

      INTEGER pif,zif

      DOUBLE PRECISION vry(rdim,adim,sdim),vrz(rdim,adim,sdim)
      DOUBLE PRECISION cqi(rdim)
      DOUBLE PRECISION cqp(rdim)
      DOUBLE PRECISION aflp(rdim,adim,sdim),apit(rdim,adim,sdim)
      DOUBLE PRECISION aoain(rdim,adim,sdim)

      INTEGER s75

      DOUBLE PRECISION clairfoil(rdim,adim,sdim),
     &       cdairfoil(rdim,adim,sdim),
     &       claoa(rdim,adim,sdim),compress(rdim,adim,sdim)
      DOUBLE PRECISION Mach(rdim,adim,sdim),Reynolds(rdim,adim,sdim)

      DOUBLE PRECISION planf
      DOUBLE PRECISION taper,taperst

      DOUBLE PRECISION Timex,maxChange

      DOUBLE PRECISION gbold(rdim,adim,sdim)
      DOUBLE PRECISION t0rold(rdim),t1cold(rdim),t1sold(rdim)

      INTEGER n,itmp,dnp,psi
      
      DOUBLE PRECISION bvx(rdim,adim,sdim),
     &       bvy(rdim,adim,sdim),
     &       bvz(rdim,adim,sdim)
      DOUBLE PRECISION rcb

      INTEGER lengthnw(sdim)

      DOUBLE PRECISION gv(rdim,adim,wdim,zdim,kdim),
     &       rc0(rdim,adim,wdim,kdim),
     &       pcx(rdim,adim,wdim,zdim),
     &       pcy(rdim,adim,wdim,zdim),
     &       pcz(rdim,adim,wdim,zdim),
     &       ppx(rdim,adim,wdim,zdim),
     &       ppy(rdim,adim,wdim,zdim),
     &       ppz(rdim,adim,wdim,zdim),
     &       eqtim(rdim,adim,wdim,zdim,kdim)

      DOUBLE PRECISION asr(rdim),rr_x(rdim),rr_y(rdim),rr_z(rdim),
     &                 lsr(rdim) 

      DOUBLE PRECISION vibx(rdim,adim,sdim),
     &       viby(rdim,adim,sdim),
     &       vibz(rdim,adim,sdim),
     &       vibxi(rdim,adim,sdim),
     &       vibyi(rdim,adim,sdim),
     &       vibzi(rdim,adim,sdim),
     $     Cldata(rdim,adim,sdim), ClMsqr(rdim,adim,sdim),
     $     liftx(rdim,adim,sdim),lifty(rdim,adim,sdim),
     $     liftz(rdim,adim,sdim)

      CHARACTER*1 sco

      DOUBLE PRECISION gsc(rdim,adim,wdim,zdim,kdim),
     &       csc(rdim,adim,wdim,zdim,kdim),sce

      DOUBLE PRECISION vinfx,vinfy,vinfz

      DOUBLE PRECISION amom(rdim,adim),thrust(rdim,adim)
     &                                ,inflow(rdim,adim)
     &                                ,torque(rdim,adim,2)

      INTEGER usaero
      DOUBLE PRECISION ub1,ub2,ua1,ua2,uf1,uf2,ug1,ug2

      INTEGER pm1,pm2,pm3
      DOUBLE PRECISION beta(rdim,adim),betastar(rdim,adim)
      DOUBLE PRECISION theta(rdim,adim)

      DOUBLE PRECISION sn,az,ze,wage,tt11,tt21,tt31,b0
      DOUBLE PRECISION rv(rdim,adim,wdim)

      DOUBLE PRECISION rhs1,rhs2

      DOUBLE PRECISION aflptmp,apittmp

      DOUBLE PRECISION CONSTANT,det

      DOUBLE PRECISION betadot,thetadot

      INTEGER Iteratn,p1,nblade
      DOUBLE PRECISION lam0,lam1c,lam1s,lam1cc,lam1ss

      INTEGER stmp

      DOUBLE PRECISION b1psi,b1stpsi,b1pm1,b1stpm1
      DOUBLE PRECISION am1pm1
      DOUBLE PRECISION am1psi
      INTEGER nn,pm1dnp
      DOUBLE PRECISION b1cpsi,b1spsi,b1cstpsi,b1sstpsi
      DOUBLE PRECISION b1cpm1,b1spm1,b1cstpm1,b1sstpm1
      DOUBLE PRECISION am1cpm1,am1spm1
      DOUBLE PRECISION am1cpsi,am1spsi
      INTEGER pdnp
      DOUBLE PRECISION rhs3,rhs4,pbar,qbar,pdot,qdot
      INTEGER rotgeo

      double precision ahdot(rdim,adim,sdim),alpdot(rdim,adim,sdim)
      double precision ualpha(rdim,adim,sdim),delalp(rdim,adim,sdim)
      double precision xcoeff(rdim,adim,sdim),ycoeff(rdim,adim,sdim)
      double precision uacirc(rdim,adim,sdim),uancir(rdim,adim,sdim)
      double precision secaoa(rdim,adim,sdim)

      double precision cpitch(rdim,adim)

      external planf


      common/usadata/ub1,ub2,ua1,ua2,uf1,uf2,ug1,ug2,usaero
      save/usadata/


      common/pidata/pi
      save/pidata/

      common/ssmdata/pbar,qbar
      save/ssmdata/

      common/unsdata/ualpha,delalp,ahdot,alpdot,xcoeff,ycoeff,uacirc,
     $     secaoa
      save/unsdata/

      common/mandata/pdot,qdot
      save/mandata/

      common/rotdata/rotgeo
      save/rotdata/



      dnp=np/nb
      p=psi
      pm1=psi-1
      if(pm1.le.0)pm1=pm1+np

      call rindvtpsi(nr,np,nw,pw,nz,ns,nk,dp,dz,ip,anb,crd,rad, dcy,kinv
     $     ,om,lin,gav,gv,rc0,asr,beta, rr_x,rr_y,rr_z,ppx,ppy,ppz,cpx
     $     ,cpy,cpz, sco,sce,gsc,csc,vibx,viby,vibz, vibxi,vibyi,vibzi
     $     ,compress,lengthnw, pif,zif,bvx,bvy,bvz,rcb,gb,BICM,NWICM,p,
     $     eqtim,lsr)
      

      do r=1,nr
         aza=DBLE(p-1)*dp/DBLE(pif)
         aoac=t0r(r)+t1c(r)*cos(aza)+t1s(r)*sin(aza)
         
         do s=1,ns
            if(abs(twt).gt.100) aoac=t0r(r)/(bcx(r,s)/rad)+t1c(r)
     $           *cos(aza)+t1s(r)*sin(aza)
            call aoapsi(r,p,s,cpx,cpy,cpz,vinfx,vinfy,vinfz, vibx,viby
     $           ,vibz,vibxi,vibyi,vibzi, nr,np,ns,dp,ip,b0r,b1c,b1s,t1c
     $           ,t1s,asr,om,bcx,bcy, aoac,tw,aza,pif,beta,theta,crd,
     $           aoai,aoain,aflp,apit,vr,vry,vrz,aoae,aoanc,lsr)
            aoae(r,p,s)=aoae(r,p,s)-aoai(r,p,s)
            aoae(r,p,s)=aoae(r,p,s)+aoain(r,p,s)
            aoae(r,p,s)=aoae(r,p,s)-xcoeff(r,p,s)-ycoeff(r,p,s)
            Mach(r,p,s)=abs(vry(r,p,s))/sonic
            Reynolds(r,p,s)=vr(r,p,s) *crd*planf(bcx(r,s)/rad,taperst
     $           ,taper)/kinv
         end do                 !s 
      end do                    !r

      call airfoilpsi(bcx,aoae,Mach,Reynolds,nr,np,ns, clairfoil
     $     ,cdairfoil,claoa,compress,p,0,alpha0)

      do r=1,nr
         call rsppsi(r,np,ns,nb,ip,dp,den,flph,om,rad,crd,clairfoil,
     $        claoa,nubeta,bmass,ibeta,betap,kbeta,seg,bcx,vr,vry,vrz,
     $        gb,Cnni,Cnqi, sonic,cdairfoil,aoae,aoai,aoain,compress,
     $        BICM,NWICM, Mach,Reynolds, inflow,aoanc,beta,betastar,
     $        taperst,taper,ct,cqi,cqp,pif,zif,amom,thrust,torque,p,
     $        Cldata, ClMsqr, liftx, lifty, liftz,cpitch) 

         if(flph.ge.0.0D+0) then
            CONSTANT=0.D+0
            pm2=psi-2
            if(pm2.le.0)pm2=pm2+np
            pm3=psi-3
            if(pm3.le.0)pm3=pm3+np
            
            rhs1=dp/2.D+0*( amom(r,pm1)+amom(r,p) - nubeta**2.D+0
     $           *(beta(r,pm1)+beta(r,p))-2*qbar*sin(aza)/om-2*pbar
     $           *cos(aza)/om+(qdot/om)*cos(aza)+(pdot/om)*cos(aza))
            rhs2=dp/2.D+0*( betastar(r,pm1)+betastar(r,p) )
            
            betastar(r,p)=betastar(r,pm1)+rhs1+CONSTANT/16.D+0
     $           *(1.D+0*betastar(r,p)-3.D+0*betastar(r,pm1)
     $           +3.D+0*betastar(r,pm2)-1.D+0*betastar(r,pm3) )
            
            beta(r,p)=beta(r,pm1) + rhs2
     &           +CONSTANT/16.D+0*(1.D+0*beta(r,p)-3.D+0*beta(r,pm1)
     &           +3.D+0*beta(r,pm2)-1.D+0*beta(r,pm3) )
            
         end if

c-mb  this part handles the teetering/gimballed rotor
         if(flph.lt.0.D+0) then
            if(nb.eq.2) then
               if(p.gt.np-dnp) then
                  
                  pm1dnp=pm1-dnp
                  if(pm1dnp.le.0) pm1dnp=pm1dnp+np
                  pdnp=p-dnp
                  b1stpm1=(betastar(r,pm1)-betastar(r,pm1dnp))/DBLE(nb)
                  b1pm1  =(beta(r,pm1)-beta(r,pm1dnp))/DBLE(nb)
                  am1pm1 =(amom(r,pm1)-amom(r,pm1dnp))/DBLE(nb)
                  b1stpsi=(betastar(r,p)-betastar(r,pdnp))/DBLE(nb)
                  b1psi  =(beta(r,p)-beta(r,pdnp))/DBLE(nb)
                  am1psi =(amom(r,p)-amom(r,pdnp))/DBLE(nb)
                  
                  rhs1=dp/2.D+0*( am1pm1 - nubeta**2.D+0*b1pm1
     &                 +am1psi - nubeta**2.D+0*b1psi )
                  rhs2=dp/2.D+0*( b1stpm1 + b1stpsi )

                  b1stpsi=b1stpm1 + rhs1
                  b1psi  =b1pm1   + rhs2
                  
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
                  b1cpsi  =0.D+0
                  b1spsi  =0.D+0
                  b1cstpsi=0.D+0
                  b1sstpsi=0.D+0
                  am1cpsi =0.D+0
                  am1spsi =0.D+0
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
                     
                     pdnp=p-(nb-nn)*dnp
                     aza=ip(r)+DBLE(pdnp-1)*dp/DBLE(pif) - ip(r)
                     b1cpsi  =b1cpsi  +cos(aza)*beta(r,pdnp)
                     b1cstpsi=b1cstpsi+cos(aza)*betastar(r,pdnp)
                     am1cpsi =am1cpsi +cos(aza)*amom(r,pdnp)
                     b1spsi  =b1spsi  +sin(aza)*beta(r,pdnp)
                     b1sstpsi=b1sstpsi+sin(aza)*betastar(r,pdnp)
                     am1spsi =am1spsi +sin(aza)*amom(r,pdnp)
                  end do
                  b1cpm1  =b1cpm1  *2.D+0/DBLE(nb)
                  b1cstpm1=b1cstpm1*2.D+0/DBLE(nb)
                  am1cpm1 =am1cpm1 *2.D+0/DBLE(nb)
                  b1spm1  =b1spm1  *2.D+0/DBLE(nb)
                  b1sstpm1=b1sstpm1*2.D+0/DBLE(nb)
                  am1spm1 =am1spm1 *2.D+0/DBLE(nb)
                  b1cpsi  =b1cpsi  *2.D+0/DBLE(nb)
                  b1cstpsi=b1cstpsi*2.D+0/DBLE(nb)
                  am1cpsi =am1cpsi *2.D+0/DBLE(nb)
                  b1spsi  =b1spsi  *2.D+0/DBLE(nb)
                  b1sstpsi=b1sstpsi*2.D+0/DBLE(nb)
                  am1spsi =am1spsi *2.D+0/DBLE(nb)
                  
                  b1cstpm1=b1cstpm1 - b1spm1
                  b1sstpm1=b1sstpm1 + b1cpm1
                  b1cstpsi=b1cstpsi - b1spsi
                  b1sstpsi=b1sstpsi + b1cpsi
                  
                  rhs1=dp/2.D+0*( am1cpm1 - 2.D+0*b1sstpm1
     &                 -(nubeta**2.D+0 -1.D+0)*b1cpm1
     &                 +am1cpsi - 2.D+0*b1sstpsi
     &                 -(nubeta**2.D+0 -1.D+0)*b1cpsi )
                  rhs2=dp/2.D+0*( am1spm1 + 2.D+0*b1cstpm1
     &                 -(nubeta**2.D+0 -1.D+0)*b1spm1
     &                 +am1spsi + 2.D+0*b1cstpsi
     &                 -(nubeta**2.D+0 -1.D+0)*b1spsi )
                  rhs3=dp/2.D+0*( b1cstpm1 + b1cstpsi )
                  rhs4=dp/2.D+0*( b1sstpm1 + b1sstpsi )
                  
                  b1cstpsi=b1cstpm1+rhs1
                  b1sstpsi=b1sstpm1+rhs2
                  b1cpsi  =b1cpm1  +rhs3
                  b1spsi  =b1spm1  +rhs4
                  
                  
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
         sn=1.D+0
         if(rotgeo.eq.1 .or. rotgeo.eq.2) sn=(-1.D+0)**(r-1)
         az=DBLE(p-1)*dp/DBLE(pif)
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

            pcx(r,p,w,1)=rr_x(r)+rv(r,p,w)*tt11
            pcy(r,p,w,1)=rr_y(r)+rv(r,p,w)*tt21
            pcz(r,p,w,1)=rr_z(r)+rv(r,p,w)*tt31
         end do                 !w
      end do                    !r

c$$$      call rindvt(nr,np,nw,pw,nz,ns,nk,dp,dz,ip,anb,crd,rad, dcy,kinv,om
c$$$     $     ,lin,gav,gv,rc0,asr,b0r,b1c,b1s, rr_x,rr_y,rr_z,pcx,pcy,pcz
c$$$     $     ,cpx,cpy,cpz, sco,sce,gsc,csc,vibx,viby,vibz, vibxi,vibyi
c$$$     $     ,vibzi,compress,lengthnw, pif,zif,bvx,bvy,bvz,rcb,gb,bicm
c$$$     $     ,nwicm,.false., eqtim)
c$$$      
c$$$      if(p.eq.np) then
c$$$         open(71,file='VBX.dat',status='unknown')
c$$$         open(72,file='VBY.dat',status='unknown')
c$$$         open(73,file='VBZ.dat',status='unknown')
c$$$         
c$$$         do r=1,nr
c$$$            
c$$$            do p=1,1
c$$$               do s=ns,1,-1
c$$$                  write(71,*)(cpx(r,p,s)-rr_x(r))/rad*cos(-asr(r))
c$$$     $                 -(cpz(r,p,s)-rr_z(r))/rad*sin(-asr(r)), vibx(r,p
c$$$     $                 ,s),vibxi(r,p,s)
c$$$                  write(72,*)(cpx(r,p,s)-rr_x(r))/rad*cos(-asr(r))
c$$$     $                 -(cpz(r,p,s)-rr_z(r))/rad*sin(-asr(r)), viby(r,p
c$$$     $                 ,s),vibyi(r,p,s)
c$$$                  write(73,314)(cpx(r,p,s)-rr_x(r))/rad*cos(-asr(r))
c$$$     $                 -(cpz(r,p,s)-rr_z(r))/rad*sin(-asr(r)), vibz(r,p
c$$$     $                 ,s),vibzi(r,p,s)
c$$$               end do
c$$$            end do
c$$$            
c$$$            write(71,*)
c$$$            write(72,*)
c$$$            write(73,*)
c$$$            
c$$$            do p=1+int(pi/(dp/DBLE(pif))),1+int(pi/(dp/DBLE(pif)))
c$$$               do s=1,ns
c$$$                  write(71,*)(cpx(r,p,s)-rr_x(r))/rad*cos(-asr(r))
c$$$     $                 -(cpz(r,p,s)-rr_z(r))/rad*sin(-asr(r)), vibx(r,p
c$$$     $                 ,s),vibxi(r,p,s)
c$$$                  write(72,*)(cpx(r,p,s)-rr_x(r))/rad*cos(-asr(r))
c$$$     $                 -(cpz(r,p,s)-rr_z(r))/rad*sin(-asr(r)), viby(r,p
c$$$     $                 ,s),vibyi(r,p,s)
c$$$                  write(73,314)(cpx(r,p,s)-rr_x(r))/rad*cos(-asr(r))
c$$$     $                 -(cpz(r,p,s)-rr_z(r))/rad*sin(-asr(r)), vibz(r,p
c$$$     $                 ,s),vibzi(r,p,s)
c$$$               end do
c$$$            end do
c$$$
c$$$            write(71,*)
c$$$            write(72,*)
c$$$            write(73,*)
c$$$            
c$$$
c$$$            do p=int(pi/(2.D+0*(dp/DBLE(pif))))+1,
c$$$     &           int(pi/(2.D+0*(dp/DBLE(pif))))+1
c$$$
c$$$               do s=ns,1,-1
c$$$                  write(71,*)(cpy(r,p,s)-rr_y(r))/rad, vibx(r,p,s)
c$$$     $                 ,vibxi(r,p,s)
c$$$                  write(72,*)(cpy(r,p,s)-rr_y(r))/rad,
c$$$     &                 viby(r,p,s),vibyi(r,p,s)
c$$$                  write(73,314)(cpy(r,p,s)-rr_y(r))/rad, vibz(r,p,s)
c$$$     $                 ,vibzi(r,p,s)
c$$$               end do
c$$$            end do
c$$$            write(71,*)
c$$$            write(72,*)
c$$$            write(73,*)
c$$$
c$$$            do p=int(pi/(2.D+0*(dp/DBLE(pif))))+1 +int(pi/(dp/DBLE(pif))
c$$$     $           ),int(pi/(2.D+0*(dp/DBLE(pif)))) +1+int(pi/(dp/DBLE(pif
c$$$     $           )))
c$$$
c$$$               do s=1,ns
c$$$                  write(71,*)(cpy(r,p,s)-rr_y(r))/rad, vibx(r,p,s)
c$$$     $                 ,vibxi(r,p,s)
c$$$                  write(72,*)(cpy(r,p,s)-rr_y(r))/rad, viby(r,p,s)
c$$$     $                 ,vibyi(r,p,s)
c$$$                  write(73,314)(cpy(r,p,s)-rr_y(r))/rad, vibz(r,p,s)
c$$$     $                 ,vibzi(r,p,s)
c$$$               end do
c$$$            end do
c$$$            
c$$$            write(71,*)
c$$$            write(72,*)
c$$$            write(73,*)
c$$$            write(71,*)
c$$$            write(72,*)
c$$$            write(73,*)
c$$$         end do
c$$$         
c$$$         close(71)
c$$$         close(72)
c$$$         close(73)
c$$$         
c$$$ 314     format(3(E12.6,1X))
c$$$         open(74,file='VBZ_PSI.dat',status='unknown')
c$$$         do r=1,nr
c$$$            do p=1,np
c$$$               do s=1,ns
c$$$                  write(74,315)cpx(r,p,s),cpy(r,p,s),cpz(r,p,s),
c$$$     $                 vibxi(r,p,s),vibyi(r,p,s),vibzi(r,p,s),
c$$$     $                 vibz(r,p,s)
c$$$                  
c$$$               end do
c$$$               write(74,*)
c$$$               write(74,*)
c$$$            end do
c$$$            write(74,*)
c$$$         end do
c$$$
c$$$         close(74)
c$$$ 315     format(7(E12.6,1X))
c$$$
c$$$         r=1                    !for now ...
c$$$         do p=1,dnp,1
c$$$            lam0 =0.d0
c$$$            lam1c=0.d0
c$$$            lam1s=0.d0
c$$$            lam1cc=0.d0
c$$$            lam1ss=0.d0
c$$$            do nblade=1,nb
c$$$               p1=p+(nblade-1)*dnp
c$$$               do s=1,ns
c$$$                  lam0 =lam0 +vibzi(r,p1,s)/om/rad
c$$$                  lam1c=lam1c+vibzi(r,p1,s)/om/rad*cos(dp*DBLE(p1-1))
c$$$                  lam1s=lam1s+vibzi(r,p1,s)/om/rad*sin(dp*DBLE(p1-1))
c$$$                  lam1cc=lam1cc+vibzi(r,p1,s)/om/rad*cos(dp*DBLE(p1-1))
c$$$     &                 *bcx(r,s)/rad
c$$$                  lam1ss=lam1ss+vibzi(r,p1,s)/om/rad*sin(dp*DBLE(p1-1))
c$$$     &                 *bcx(r,s)/rad
c$$$               end do           !s
c$$$            end do              !nblade
c$$$            lam0 =lam0       /DBLE(nb)/DBLE(ns)
c$$$            lam1c=lam1c*2.D+0/DBLE(nb)/DBLE(ns)
c$$$            lam1s=lam1s*2.D+0/DBLE(nb)/DBLE(ns)
c$$$            lam1cc=lam1cc*2.D+0/DBLE(nb)/DBLE(ns)
c$$$            lam1ss=lam1ss*2.D+0/DBLE(nb)/DBLE(ns)
c$$$            write(77,316)r,dp*DBLE(p-1)+2.D+0*pi*DBLE(Iteratn-1)
c$$$     $           /DBLE(nb),lam0,lam1c,lam1s,lam1cc,lam1ss
c$$$         end do                 !p1
c$$$ 316     format(I2,1X,6(E13.6,1X))
c$$$         p=np                   !to restore p to np
c$$$      endif
c$$$
      return
      end
