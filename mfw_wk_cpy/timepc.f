C-$	$Id: timepc.f,v 1.14.2.1 2005/03/10 18:49:33 shreyas Exp $	

      subroutine timepc(nr,np,nw,pw,nz,nzt,nk,ns,nb,pif,zif,ip,dp,dz,
     $     sonic,vinfx,vinfy,vinfz,om,asr,b0r,anb,rad,rcout,rcb,gb, t0r
     $     ,t1c,t1s,flph,nubeta,kbeta,bmass,ibeta,betap, lin,gv,gvold
     $     ,rc0,rv,kinv,dcy,gav,sco,sce,gsc,csc, rpar,rrn,zrn,rbc,zbc
     $     ,den,compress,lengthnw, clairfoil,cdairfoil,claoa,Mach
     $     ,Reynolds,crd, avgzbc,lami,rr_x,rr_y,rr_z,spx,spy,spz,seg, tw
     $     ,twt,aoai,aoain,bcx,bcy,bcz,bicm,nwicm, vr,vry,taper,taperst
     $     ,cnni,cnqi,ct,b1s,b1c, cpx,cpy,cpz,pcx,pcy,pcz,beta,betastar
     $     ,theta, cqi,cqp,gammaRMS,betaRMS, pox,poy,poz,ppx,ppy,ppz,bvx
     $     ,bvy,bvz,Iteratn,outplot, eqtim, Cldata, ClMsqr, liftx, lifty
     $     , liftz, trm,alpha0,lsr)

C-$   MrB's PC2B scheme

      IMPLICIT none

      include 'adimpar.inc'
      include 'gdimpar.inc'
      include 'jdimpar.inc'
      include 'kdimpar.inc'
      include 'rdimpar.inc'
      include 'sdimpar.inc'
      include 'wdimpar.inc'
      include 'zdimpar.inc'

      INTEGER nr,np,nw,pw,nz,nzt,nk,ns,nb,pif,zif, r,p,w,z,n,s,dnp
     $     ,pblade,nblade,Psi,p0,k,Iteratn,i

      DOUBLE PRECISION ip(rdim), dp,dz,vinfx,vinfy,vinfz,om,anb,rad
     $     ,rcout, rcb,lin,kinv,dcy,gav,rpar,rrn,zrn,sce, ze,wage,tt11
     $     ,tt21,tt31,b0,sonic,pbar,qbar,pi

      DOUBLE PRECISION rbc(rdim,wdim),zbc(rdim,wdim), avgzbc(rdim,wdim)

      DOUBLE PRECISION asr(rdim),b0r(rdim),lami(rdim), rr_x(rdim)
     $     ,rr_y(rdim),rr_z(rdim), spx(sdim),spy(sdim),spz(sdim),
     $     lsr(rdim)

      DOUBLE PRECISION rv(rdim,adim,wdim), gb(rdim,adim,sdim), cpx(rdim
     $     ,adim,sdim), cpy(rdim,adim,sdim), cpz(rdim,adim,sdim),
     $     vibx(rdim,adim,sdim), viby(rdim,adim,sdim), vibz(rdim,adim
     $     ,sdim), vibxi(rdim,adim,sdim), vibyi(rdim,adim,sdim),
     $     vibzi(rdim,adim,sdim), bvx(rdim,adim,sdim), bvy(rdim,adim
     $     ,sdim), bvz(rdim,adim,sdim), Cldata(rdim,adim,sdim)
     $     ,ClMsqr(rdim,adim,sdim), liftx(rdim,adim,sdim),lifty(rdim
     $     ,adim,sdim), liftz(rdim,adim,sdim)

      DOUBLE PRECISION gv(rdim,adim,wdim,zdim,kdim), gvold(rdim,adim
     $     ,wdim,kdim), gvnew(rdim,adim,wdim,kdim), rc0(rdim,adim,wdim
     $     ,kdim), ppx(rdim,adim,wdim,zdim), ppy(rdim,adim,wdim,zdim),
     $     ppz(rdim,adim,wdim,zdim), pcx(rdim,adim,wdim,zdim), pcy(rdim
     $     ,adim,wdim,zdim), pcz(rdim,adim,wdim,zdim), pox(rdim,adim
     $     ,wdim,zdim), poy(rdim,adim,wdim,zdim), poz(rdim,adim,wdim
     $     ,zdim), vipx(rdim,adim,wdim,zdim), vipy(rdim,adim,wdim,zdim),
     $     vipz(rdim,adim,wdim,zdim), vicx(rdim,adim,wdim,zdim),
     $     vicy(rdim,adim,wdim,zdim), vicz(rdim,adim,wdim,zdim),
     $     eqtim(rdim,adim,wdim,zdim,kdim)

      CHARACTER*1 sco,trm

      DOUBLE PRECISION gsc(rdim,adim,wdim,zdim,kdim), csc(rdim,adim,wdim
     $     ,zdim,kdim)

      DOUBLE PRECISION t0r(rdim),t1c(rdim),t1s(rdim)
      DOUBLE PRECISION t0r_0(rdim),t1c_0(rdim),t1s_0(rdim)

      INTEGER nn,zz
      LOGICAL outplot
      INTEGER pminus

      INTEGER zen

      DOUBLE PRECISION flph,nubeta,kbeta,bmass,ibeta,betap

      DOUBLE PRECISION den,crd,cttol,fltol
      DOUBLE PRECISION cqi(rdim)
      DOUBLE PRECISION cqp(rdim)


      INTEGER lengthnw(sdim)

      DOUBLE PRECISION clairfoil(rdim,adim,sdim),
     &       cdairfoil(rdim,adim,sdim),
     &       claoa(rdim,adim,sdim),compress(rdim,adim,sdim)
      DOUBLE PRECISION Mach(rdim,adim,sdim),Reynolds(rdim,adim,sdim)
      DOUBLE PRECISION ct(rdim),b1c(rdim),b1s(rdim),seg(sdim),
     &       bcx(rdim,sdim),bcy(rdim,sdim),bcz(rdim,sdim)
      DOUBLE PRECISION bicm(sdim,sdim),nwicm(sdim,sdim)
      DOUBLE PRECISION aoai(rdim,adim,sdim),aoain(rdim,adim,sdim),
     &       vr(rdim,adim,sdim),vry(rdim,adim,sdim),
     &       Cnni(rdim,adim,sdim),Cnqi(rdim,adim,sdim)

      DOUBLE PRECISION taper,taperst
      DOUBLE PRECISION tw(rdim,sdim),twt

      DOUBLE PRECISION beta(rdim,adim),betastar(rdim,adim)
      DOUBLE PRECISION oldbeta(rdim,adim)
      DOUBLE PRECISION betaRMS(rdim)
      DOUBLE PRECISION theta(rdim,adim)
      DOUBLE PRECISION aza

      INTEGER pout1,pout2,pout3,pout4

      INTEGER pnp
      DOUBLE PRECISION tmp(rdim,adim),tmpstar(rdim,adim)

      DOUBLE PRECISION sumgamma,gammaRMS,OLDgammaRMS

      DOUBLE PRECISION amom(rdim,adim),thrust(rdim,adim) ,inflow(rdim
     $     ,adim) ,torque(rdim,adim,2)
      DOUBLE PRECISION nthrust(rdim,adim),ntorque(rdim,adim)
     $     ,ninflow(rdim,adim) ,inertia(rdim,adim)

      INTEGER pm1,count,t0steps
      DOUBLE PRECISION dthetac,thetac

      DOUBLE PRECISION hf(rdim,adim),yf(rdim,adim),thincr
      DOUBLE PRECISION timet,blarea

      double precision ahdot(rdim,adim,sdim),alpdot(rdim,adim,sdim)
      double precision ualpha(rdim,adim,sdim),delalp(rdim,adim,sdim)
      double precision xcoeff(rdim,adim,sdim),ycoeff(rdim,adim,sdim)
      double precision uacirc(rdim,adim,sdim),uancir(rdim,adim,sdim)
      double precision secaoa(rdim,adim,sdim)

      double precision vrz(rdim,adim,sdim), aflp(rdim,adim,sdim),
     $     apit(rdim,adim,sdim), aoae(rdim,adim,sdim), aoanc(rdim,adim
     $     ,sdim),aoac,alpha0

      double precision px_ml(rdim*adim*zdim*wdim)
      double precision py_ml(rdim*adim*zdim*wdim) 
      double precision pz_ml(rdim*adim*zdim*wdim)
      double precision vcx_ml(rdim*adim*zdim*wdim)
      double precision vcy_ml(rdim*adim*zdim*wdim)
      double precision vcz_ml(rdim*adim*zdim*wdim)



      common/pidata/pi
      save/pidata/

      common/ssmdata/pbar,qbar
      save/ssmdata/

      common/changec/count
      save/changec/

      common/origthet/t0r_0,t1c_0,t1s_0
      save/origthet/

      common/unsdata/ualpha,delalp,ahdot,alpdot,xcoeff,ycoeff,uacirc,
     $     secaoa
      save/unsdata/

      common/indv_ml/vcx_ml,vcy_ml,vcz_ml
      save/indv_ml/


      do r=1,nr
         do p=1,np
            oldbeta(r,p)=beta(r,p)
         end do
      end do

      call ml_coord(pcx,pcy,pcz,adim*rdim*zdim*wdim,px_ml,py_ml,pz_ml)
      call in2out(px_ml,py_ml,pz_ml,adim*rdim*zdim*wdim)
      call wake_pertub(px_ml, py_ml, pz_ml, rdim*adim*zdim*wdim,
     $  vcx_ml, vcy_ml, vcz_ml)
      call ml_inv(vcx_ml,vcy_ml,vcz_ml,rdim*adim*zdim*wdim)
 
      
      
      zen=1                     !zen=1 is PC2B 
      dnp=np/nb

      do pblade=1,dnp

C-$   Read in changes in control inputs or flight conditions... 
         call controll(Iteratn,nb,pblade,dp,pif,om,dnp,np,nr,theta,t0r
     $        ,t1c,t1s,ip,asr,pbar,qbar,tw,vinfx,vinfy,vinfz,rad)
C-$   write out flight conditions for each timestep...
         r=1
         timet=2.D+0*pi*DBLE(Iteratn-1)/DBLE(nb)+dp*DBLE(pblade-1)
         write(69,10) timet, t0r(r), t1c(r), t1s(r), vinfx, vinfy,
     $        vinfz, pbar, qbar
 10      format(9(E15.7,2x))


         do nblade=1,nb           
            Psi=pblade+(nblade-1)*dnp
            do r=1,nr
               do w=1,nw
                  do z=1,nz
                     vipx(r,Psi,w,z)=0.D+0
                     vipy(r,Psi,w,z)=0.D+0
                     vipz(r,Psi,w,z)=0.D+0
                     vicx(r,Psi,w,z)=0.D+0
                     vicy(r,Psi,w,z)=0.D+0
                     vicz(r,Psi,w,z)=0.D+0
                  end do
               end do
            end do

            p0=Psi-1
            if(p0 .le. 0) p0=p0+np !will occur only once
            call indv(nr,np,nw,pw,nz,nzt,nk,ns,pif,zif,dp,dz,rad, om,anb
     $           ,lin,dcy,gav,kinv,gv,rc0,sco,sce,gsc,csc, rcb,gb,ip,b0r
     $           ,asr,lsr,rr_x,rr_y,rr_z, t0r,t1c,t1s,beta, spx,spy,spz
     $           ,pcx,pcy,pcz,bvx,bvy,bvz, vipx,vipy,vipz,p0,cpx,cpy,cpz
     $           ,eqtim)


            p=Psi
            if(zif.gt.1) call vinterpe(nr,np,nw,pw,nz,pif,zif,vipx,vipy
     $           ,vipz,p)

            do r=1,nr
               aza=DBLE(p-1)*dp/DBLE(pif)
               aoac=t0r(r)+t1c(r)*cos(aza)+t1s(r)*sin(aza)
               do s=1,ns
                  call aoapsi(r,p,s,cpx,cpy,cpz,vinfx,vinfy,vinfz, vibx
     $                 ,viby,vibz,vibxi,vibyi,vibzi, nr,np,ns,dp,ip,b0r
     $                 ,b1c,b1s,t1c,t1s,asr,om,bcx,bcy, aoac,tw,aza,pif
     $                 ,beta,theta,crd, aoai,aoain,aflp,apit,vr,vry,vrz
     $                 ,aoae,aoanc,lsr)
                  secaoa(r,p,s)=aoac
               end do           !s
            end do              !r
         end do                 !nblade
         
         if(trm .ne. 'y') then
            call unsteady_aero(nr,np,ns,nb,dp,ip,tw,crd,cpx,cpy,cpz,
     $           vinfx, vinfy,vinfz,vibx,viby,vibz, t0r,t1c,t1s,theta
     $           ,b0r,b1c,b1s,beta,betastar,asr,pif, om,bcx,bcy,bcz
     $           ,sonic,Mach,Reynolds,den,taper,taperst, kinv,rad,vr,vry
     $           ,rcout,Cnni,Cnqi,pblade,aoain)
         endif

         do nblade=1,nb           
            Psi=pblade+(nblade-1)*dnp
            p=Psi
            p0=Psi-1
            if(p0 .le. 0) p0=p0+np

            call flappred(nr,np,ns,nw,pw,nz,nk,nb,ip,dp,dz,den,flph,om
     $           ,rad,compress,nubeta,bmass,ibeta,betap,kbeta,crd
     $           ,clairfoil,claoa,cttol,fltol,seg,tw,bicm,nwicm,bcx,bcy
     $           ,twt,cdairfoil,kinv,sonic,aoai,aoain,Cnni,Cnqi,vr,vry
     $           ,gb,anb,rv,dcy,gv,asr,lsr,rc0,sco,sce,gsc,csc,pcx,pcy
     $           ,pcz,ppx,ppy,ppz,amom,thrust,torque,inflow,cpx,cpy,cpz
     $           ,taperst,vinfx,vinfy,vinfz,rr_x,rr_y,rr_z,taper,t0r,t1c
     $           ,t1s,b0r,b1c,b1s,ct,cqi,cqp,pif,zif,theta,bvx,bvy,bvz
     $           ,rcb,lengthnw,p,beta,betastar,eqtim, Cldata, ClMsqr,
     $           liftx, lifty, liftz,alpha0)


         end do                 !nblade

         do nblade=1,nb           
            Psi=pblade+(nblade-1)*dnp
            p=Psi
            p0=Psi-1
            if(p0 .le. 0) p0=p0+np

c-mb  update blade control points with the predicted flap!
            call bcppsi(np,nr,ns,pif,ip,dp,asr,lsr,b0r,b1c,b1s,rad,om
     $           ,t0r,t1c,t1s,twt,tw,beta,theta,rr_x,rr_y,rr_z,bcx,bcy
     $           ,bcz,cpx,cpy,cpz,p)

c-mb  update vortex strength

            call vortexpsi(nr,np,nw,ns,nzt,nk,dp,anb,rad,cpx,spx,gb,rv,
     $           gv,p)


c-mb needs vipx,pcx etc at p-1 uses pcx @ p as its = pox

            call pred(nr,np,nw,pw,nz,nzt,nk,dp,dz,pif,zif,rad,om,rbc,
     $           zbc,avgzbc,lami,gv,rr_x,rr_y,rr_z,vinfx,vinfy,vinfz,
     $           pcx,pcy,pcz,vipx,vipy,vipz,ppx,ppy,ppz, p)

C-$   Apply strain model for the vortex filament
            call strainmb(nr,np,nw,pw,nz,nzt,nk,ns,nb,pif,zif,ip,dp,dz,
     $           gv,rc0,kinv,dcy,om,gav, ppx,ppy,ppz, pcx,pcy,pcz,p,'p',
     $           eqtim)


c-mb     calculating velocities using ppx@p so update gv NOW

            p=Psi
            p0=Psi-1
            if(p0 .le. 0) p0=p0+np !will occur only once
            do r=1,nr
               do w=1,nw
                  do z=2,nzt
                     do k=1,nk
                        gv(r,p,w,z,k)=gv(r,p0,w,z-1,k)
                     end do
                  end do
               end do
            end do

         end do                 !nblade

         do nblade=1,nb           
            Psi=pblade+(nblade-1)*dnp

            p=Psi
            call indv(nr,np,nw,pw,nz,nzt,nk,ns,pif,zif,dp,dz,rad, om,anb
     $           ,lin,dcy,gav,kinv,gv,rc0,sco,sce,gsc,csc, rcb,gb,ip,b0r
     $           ,asr,lsr,rr_x,rr_y,rr_z, t0r,t1c,t1s,beta, spx,spy,spz
     $           ,ppx,ppy,ppz,bvx,bvy,bvz, vicx,vicy,vicz,p,cpx,cpy,cpz,
     $           eqtim)


            if(zif.gt.1) call vinterpe(nr,np,nw,pw,nz,pif,zif,vipx,vipy
     $           ,vipz,p)

         end do                 !nblade

         do nblade=1,nb           
            Psi=pblade+(nblade-1)*dnp
            p=Psi
            p0=Psi-1
            if(p0 .le. 0) p0=p0+np

            call flapcorr(nr,np,ns,nw,pw,nz,nk,nb,ip,dp,dz,den,flph,om
     $           ,rad,compress,nubeta,bmass,ibeta,betap,kbeta,crd
     $           ,clairfoil,claoa,cttol,fltol,seg,tw,bicm,nwicm,bcx,bcy
     $           ,twt,cdairfoil,kinv,sonic,aoai,aoain,Cnni,Cnqi,vr,vry
     $           ,gb,anb,rv,dcy,gv,asr,lsr,rc0,sco,sce,gsc,csc,pcx,pcy
     $           ,pcz,ppx,ppy,ppz,amom,thrust,torque,inflow,cpx,cpy,cpz
     $           ,taperst,vinfx,vinfy,vinfz,rr_x,rr_y,rr_z,taper,t0r,t1c
     $           ,t1s,b0r,b1c,b1s,ct,cqi,cqp,pif,zif,theta,bvx,bvy,bvz
     $           ,rcb,lengthnw,p,beta,betastar,Iteratn,eqtim, vibxi,
     $           vibyi, vibzi, Cldata, ClMsqr, liftx, lifty,liftz,alpha0
     $           )

         end do                 !nblade

C-$   Write out inflow history. Feb. 19, 2003
C-$   Earlier in flapcorr.f
C-$   This is a better location imo. 
         Psi=pblade+mod(Iteratn-1,nb)*dnp
         do r=1,nr
            write(78,317)r,(dp*DBLE(pblade-1)*180.0D+0/pi+360.0D+0
     $           *DBLE(Iteratn-1)/DBLE(nb)),(vibzi(r,Psi,i),i=1,ns)
            do nblade=1,nb
               p=Psi+(nblade-1)*dnp
               if(p.gt.np) p=p-np
               write(56,318)r,nblade, (dp*DBLE(pblade-1)*180.0D+0/pi+360
     $              .0D+0 *DBLE(Iteratn-1)/DBLE(nb)),(Cldata(r,p,i),i=1
     $              ,ns)
               write(57,318)r,nblade, (dp*DBLE(pblade-1)*180.0D+0/pi+360
     $              .0D+0 *DBLE(Iteratn-1)/DBLE(nb)),(ClMsqr(r,p,i),i=1
     $              ,ns)
               write(60,318)r,nblade, (dp*DBLE(pblade-1)*180.0D+0/pi+360
     $              .0D+0 *DBLE(Iteratn-1)/DBLE(nb)),(liftx(r,p,i),i=1
     $              ,ns)
               write(64,318)r,nblade, (dp*DBLE(pblade-1)*180.0D+0/pi+360
     $              .0D+0 *DBLE(Iteratn-1)/DBLE(nb)),(lifty(r,p,i),i=1
     $              ,ns)
               write(67,318)r,nblade, (dp*DBLE(pblade-1)*180.0D+0/pi+360
     $              .0D+0 *DBLE(Iteratn-1)/DBLE(nb)),(liftz(r,p,i),i=1
     $              ,ns)
            end do
         end do
 317     format(i4,E13.6,1x,100(E13.6,1x))
 318     format(i4,i4,E13.6,1x,100(E13.6,1x))

         do nblade=1,nb           
            Psi=pblade+(nblade-1)*dnp
            p=Psi
            p0=Psi-1
            if(p0 .le. 0) p0=p0+np

            call bcppsi(np,nr,ns,pif,ip,dp,asr,lsr,b0r,b1c,b1s,rad,om
     $           ,t0r,t1c,t1s,twt,tw,beta,theta, rr_x,rr_y,rr_z,bcx,bcy 
     $           ,bcz,cpx,cpy,cpz,p)
            call vortexpsi(nr,np,nw,ns,nzt,nk,dp,anb,rad,cpx,spx,gb,rv,
     $           gv,p)

            p=Psi
            call corr(nr,np,nw,pw,nz,nzt,nk,dp,dz,pif,zif,rad,om,rbc,
     $           zbc,avgzbc,lami,gv,rr_x,rr_y,rr_z,vinfx,vinfy,vinfz,
     $           rpar,vipx,vipy,vipz,ppx,ppy,ppz,vicx,vicy, vicz,pcx,pcy
     $           ,pcz,p,zen)

            call strainmb(nr,np,nw,pw,nz,nzt,nk,ns,nb,pif,zif,ip,dp,dz,
     $           gv,rc0,kinv,dcy,om,gav, ppx,ppy,ppz, pcx,pcy,pcz,p,'c',
     $           eqtim)

            p=Psi
            p0=Psi-1
            if(p0 .le. 0) p0=p0+np !will occur only once
            do r=1,nr
               do w=1,nw
                  do z=2,nzt
                     do k=1,nk
                        gv(r,p,w,z,k)=gv(r,p0,w,z-1,k)
                     end do
                  end do
               end do
            end do


         end do                 !nblade 

         do r=1,nr
            pout1=pblade+dnp*mod(Iteratn-1,nb)
            if(pout1.gt.np)pout1=pout1-np
            pout2=pout1+dnp
            if(pout2.gt.np)pout2=pout2-np
            pout3=pout2+dnp
            if(pout3.gt.np)pout3=pout3-np
            pout4=pout3+dnp
            if(pout4.gt.np)pout4=pout4-np
            write(44,'(I2,1x,9(E13.6,1X))')r, 360.D+0*DBLE(Iteratn-1)
     $           /DBLE(nb)+dp*180.D+0/pi*DBLE(pblade-1), beta(r,pout1)
     $           ,beta(r,pout2),beta(r,pout3),beta(r,pout4), amom(r
     $           ,pout1),amom(r,pout2),amom(r,pout3),amom(r,pout4)
         end do                 !r

      end do                    !pblade

c     calculate the rotor thrust !

      do r=1,nr
         ct(r)=0.D+0
         cqi(r)=0.D+0
         cqp(r)=0.D+0
         do p=1,np
            aza=ip(r)+DBLE(p-1)*dp/DBLE(pif) - ip(r)
            ct(r)=ct(r)+thrust(r,p)

            inertia(r,p)=(amom(r,p)-nubeta**2.D+0*beta(r,p))*om**2.D+0
     &           *bmass*(rad**2.D+0-rcout**2.D+0)/2.D+0

            if(flph.lt.0.D+0)inertia(r,p)=0.D+0

            ct(r)=ct(r)-inertia(r,p)
            if(p.gt.dnp)inertia(r,p)=inertia(r,p)+inertia(r,p-dnp)

            cqi(r)=cqi(r)+torque(r,p,1)
            cqp(r)=cqp(r)+torque(r,p,2)
            nthrust(r,p)=thrust(r,p)
            ntorque(r,p)=torque(r,p,1)+torque(r,p,2)
            ninflow(r,p)=inflow(r,p)

c-mb note that I am using beta to calculate drag/side force but for
c-mb thrust I assume beta=0, i.e., cos(beta)=1 ... not that it would ...

            hf(r,p)=-sin(beta(r,p))*thrust(r,p)*cos(aza)
            yf(r,p)=-sin(beta(r,p))*thrust(r,p)*sin(aza)
            if(p.gt.dnp)nthrust(r,p)=nthrust(r,p)+nthrust(r,p-dnp)
            if(p.gt.dnp)ninflow(r,p)=ninflow(r,p)+ninflow(r,p-dnp)
            if(p.gt.dnp)ntorque(r,p)=ntorque(r,p)+ntorque(r,p-dnp)
            if(p.gt.dnp)hf(r,p)=hf(r,p)+hf(r,p-dnp)
            if(p.gt.dnp)yf(r,p)=yf(r,p)+yf(r,p-dnp)
            if(p.gt.(np-dnp)) write(55,'(I2,1X,7(E13.6,1X))')r, 360.D+0
     $           *DBLE(Iteratn-1)/DBLE(nb) +180.D+0*dp/pi*DBLE((p-(np
     $           -dnp))-1), (nthrust(r,p)-inertia(r,p))/(den*pi*(rad**4)
     $           *om**2), ntorque(r,p)/(den*pi*(rad**4)*om**2*rad),
     $           ninflow(r,p)/DBLE(nb), inertia(r,p)/(den*pi*(rad**4)*om
     $           **2) ,hf(r,p)/(den*pi*(rad**4)*om**2) ,yf(r,p)/(den*pi
     $           *(rad**4)*om**2)
         end do                 !p
         ct(r)=ct(r)*DBLE(nb)/DBLE(np)
         ct(r)=ct(r)/(den*pi*(rad**4)*om**2)
         cqi(r)=cqi(r)*DBLE(nb)/DBLE(np)
         cqi(r)=cqi(r)/(den*pi*(rad**4)*om**2*rad)
         cqp(r)=cqp(r)*DBLE(nb)/DBLE(np)
         cqp(r)=cqp(r)/(den*pi*(rad**4)*om**2*rad)
      end do                    !r

c$$$      do r=1,nr
c$$$       b0r(r)=0.D+0
c$$$       b1c(r)=0.D+0
c$$$       b1s(r)=0.D+0
c$$$       do p=1,np
c$$$        aza=ip(r)+DBLE(p-1)*dp/DBLE(pif) - ip(r)
c$$$        b0r(r)=b0r(r)+beta(r,p)
c$$$        b1c(r)=b1c(r)+beta(r,p)*cos(aza)
c$$$        b1s(r)=b1s(r)+beta(r,p)*sin(aza)
c$$$       end do !p
c$$$       b0r(r)=b0r(r)*dp/2.D+0/pi
c$$$       b1c(r)=b1c(r)*dp/pi
c$$$       b1s(r)=b1s(r)*dp/pi
c$$$       end do !r

c     need to cycle beta to maintain continuity of blade(s)
c     also calculate RMS change

      do r=1,nr
         betaRMS(r)=0.D+0
         do p=1,np
            pnp=p+dnp
            if(pnp.le.0)pnp=pnp+np
            if(pnp.gt.np)pnp=pnp-np
            tmp(r,p)=beta(r,pnp)
            tmpstar(r,p)=betastar(r,pnp)
            betaRMS(r)=betaRMS(r)+(beta(r,p)-oldbeta(r,p))**2.D+0
         end do                 ! p
         betaRMS(r)=sqrt(betaRMS(r)/DBLE(NP))
      end do                    ! r

      if(Iteratn.gt.1) OLDgammaRMS=gammaRMS
      
      call wrtgbmb(anb,dp,nr,np,ns,rad,cpx,cpz,rr_x,rr_z,asr,lsr,
     &             gb,pif,zif)

      sumgamma=0.D+0
      do r=1,nr
         do p=1,np
            do w=1,nw
               do k=1,nk
                  sumgamma=sumgamma+(gv(r,p,w,1,k)-gvold(r,p,w,k))**2.D
     $                 +0
               end do
            end do
         end do
      end do

      gammaRMS=sqrt(sumgamma/DBLE(nr*np*nw*nk))
      if(gammaRMS.eq.0) gammaRMS=OLDgammaRMS
      write(*,*)'Change in bound circn: ',gammaRMS
      write(*,*)'Change in blade flap angles: '
      do r=1,nr
         write(*,'(A14,I2,A2,E13.6)')'BetaRMS (Rotor',r,') ',betaRMS(r)
      end do

      return
      end

