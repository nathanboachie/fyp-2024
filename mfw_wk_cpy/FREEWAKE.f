C-$	$Id: FREEWAKE.f,v 1.10.2.4 2005/06/11 15:57:36 shreyas Exp $	

C-$   FREE VORTEX WAKE METHOD
C-$   DEPARTMENT OF AEROSPACE ENGINEERING 
C-$   UNIVERSITY OF MARYLAND COLLEGE PARK

      IMPLICIT none


      include 'adimpar.inc'
      include 'gdimpar.inc'
      include 'jdimpar.inc'
      include 'kdimpar.inc'
      include 'rdimpar.inc'
      include 'sdimpar.inc'
      include 'wdimpar.inc'
      include 'zdimpar.inc'


      INTEGER nb,zif,pif,z,nz,nzt,p,np,r,nr,s,ns,w,nw,pw,nk,li,gs,k,
     $     nctrl,usaero

      DOUBLE PRECISION pi,ft,dp,dz,mu,muc,nt,ip(rdim),rad,flph0,flph,
     $     crd,twt,om,rada,anb, rcout,sce,vn,rcb,dcy,rpar,rrn,zrn,cnv
     $     ,cvt,lin,ct0,asr0, b1c0,b1s0,t0_0,t1c0,t1s0,b0_0, cttol,fltol
     $     ,lock,den,bmass,ibeta,sbeta,nubeta, kbeta,betap ,vinfx ,vinfy
     $     ,vinfz,sonic,sn,kinv, ub1,ub2,ua1,ua2,uf1,uf2,ug1 ,ug2 ,pbar
     $     ,qbar,lsr0
      DOUBLE PRECISION bct

      DOUBLE PRECISION clairfoil(rdim,adim,sdim), cdairfoil(rdim,adim
     $     ,sdim), claoa(rdim,adim,sdim),compress(rdim,adim,sdim)
      DOUBLE PRECISION Mach(rdim,adim,sdim),Reynolds(rdim,adim,sdim)

      DOUBLE PRECISION rr_x(rdim),rr_y(rdim),rr_z(rdim), asr(rdim)
     $     ,lami(rdim), b0r(rdim),b1c(rdim),b1s(rdim), t0r(rdim)
     $     ,t1c(rdim),t1s(rdim), ct(rdim),spx(sdim),spy(sdim),spz(sdim)
     $     ,seg(sdim),lsr(rdim)

      DOUBLE PRECISION tw(rdim,sdim), bcx(rdim,sdim),bcy(rdim,sdim)
     $     ,bcz(rdim,sdim), nwicm(sdim,sdim),jac(jdim,jdim)

      DOUBLE PRECISION bicm(sdim,sdim)

      DOUBLE PRECISION rv(rdim,adim,wdim), gb(rdim,adim,sdim), cpx(rdim
     $     ,adim,sdim), cpy(rdim,adim,sdim), cpz(rdim,adim,sdim),
     $     vr(rdim,adim,sdim), aoai(rdim,adim,sdim), Cnni(rdim,adim,sdim
     $     ), Cnqi(rdim,adim,sdim)

      DOUBLE PRECISION aoain(rdim,adim,sdim),alpha0

      DOUBLE PRECISION gv(rdim,adim,wdim,zdim,kdim), rc0(rdim,adim,wdim
     $     ,kdim), pcx(rdim,adim,wdim,zdim), pcy(rdim,adim,wdim,zdim),
     $     pcz(rdim,adim,wdim,zdim), eqtim(rdim,adim,wdim,zdim,kdim)

      CHARACTER*1 sco,trm,cvr,method,initial

      DOUBLE PRECISION gsc(rdim,adim,wdim,zdim,kdim),
     &       csc(rdim,adim,wdim,zdim,kdim)

c-mb these are for pretrim.f which was included here
      DOUBLE PRECISION aza,aoac, vfx(rdim,adim,sdim), vfy(rdim,adim,sdim
     $     ), vfz(rdim,adim,sdim), vrx(rdim,adim,sdim), vry(rdim,adim
     $     ,sdim), vrz(rdim,adim,sdim), aoae(rdim,adim,sdim), aflp(rdim
     $     ,adim,sdim), apit(rdim,adim,sdim), wlfv(rdim,adim,sdim),
     $     b0r_0(rdim),b1c_0(rdim),b1s_0(rdim), ct_0(rdim)
  

c-mb these are for the mfw.f never knew why it was a routine
      DOUBLE PRECISION gav,gsv(gdim,gdim)
      INTEGER j1,j2,gsn

      INTEGER n1,dn1

      INTEGER nfw

      DOUBLE PRECISION cqi(rdim)
      DOUBLE PRECISION cqp(rdim)

      INTEGER lengthnw(sdim)

      DOUBLE PRECISION planf
      DOUBLE PRECISION taper,taperst,altitude,isaplus
      INTEGER rotgeo
      LOGICAL outplot

      DOUBLE PRECISION rcore(rdim,adim,wdim,zdim) ,squire(rdim,adim,wdim
     $     ,zdim) ,osquire(rdim,adim,wdim,zdim)
      DOUBLE PRECISION trb,zet0,wage, factor


      external planf


c >>  DATA INPUT FILES <<
c >>  User defined inputs: <<
      namelist/user/nw,pw,nfw,nk,ft,bct,dp,dz,sco,sce,vn,rcb,dcy,rpar,
     $     rrn,zrn,cnv,trm,cvr,cvt,lin,li,gs,method,initial, outplot

c >>  Rotor flight and operating conditions: <<
      namelist/flight/mu,muc,pbar,qbar,ct0,t0_0,t1c0,t1s0, b0_0,b1c0
     $     ,b1s0,cttol,fltol,altitude,isaplus

c >>  Rotor geometry: <<
      namelist/geometry/nr,nb,ns,asr0,lsr0,rad,flph, rcout,crd,twt,om
     $     ,taperst,taper,rotgeo

c >>  Rotor mass and inertial properties: <<
      namelist/rotprop/bmass,ibeta,sbeta,nubeta,kbeta,betap

c >>  Unsteady aerodynamics coefficients: <<
      namelist/usa/usaero,ub1,ub2,ua1,ua2,uf1,uf2,ug1,ug2,alpha0


      common/pidata/pi
      save/pidata/

      common/vndata/vn
      save/vndata/

      common/usadata/ub1,ub2,ua1,ua2,uf1,uf2,ug1,ug2,usaero
      save/usadata/

      common/ssmdata/pbar,qbar
      save/ssmdata/

      common/coredata/rcore,squire,osquire
      save/coredata/

      common/rotdata/rotgeo
      save/rotdata/

      common/ctdata/ct0
      save/ctdata/


      open(33,file='user.input',status='old')
      read(33,user)
      close(33)
 
      open(44,file='flight.input',status='old')
      read(44,flight)
      close(44)

C/MR  Define lsr0 for older geometry input files
      lsr0=0.D+0
      open(55,file='geometry.input',status='old')
      read(55,geometry)
      close(55)

      open(66,file='rotprop.input',status='old')
      read(66,rotprop)
      close(66)

      open(77,file='usa.input',status='old')
      read(77,usa)
      close(77)

      pi=4.D+0*atan(1.D+0)
      betap=betap*pi/180.D+0 !read in beta_precone in degrees


      flph0=flph

      if(flph.lt.0.D+0) then
         write(*,*)'Negative flap hinge : setting to zero'
         write(*,*)'                    : teetering/gimballed hub'
         flph=0.D+0
      end if
      sbeta=(bmass/2.D+0)*(rad**2-2.D+0*flph*rad+flph**2)
      ibeta=(bmass/1.D+0)*(rad**3-3.D+0*flph*rad**2+ 3.D+0*rad*flph**2
     $     -flph**3)
      ibeta=ibeta/3.D+0
      nubeta=sqrt(1.D+0+flph*(sbeta/ibeta)+(kbeta/ibeta)/(om**2))
      flph=flph0

c >>  Some Necessary Initializations: <<

      call isa(den,sonic,altitude,isaplus, kinv)
      nctrl=3

      pbar=pbar*om
      qbar=qbar*om
      alpha0=alpha0*pi/180.0D+0

      if (pw .ge. nw) then
          write(*,*)' nw =',nw,',    pw =',pw
          write(*,*)' pw must be less than nw [user.input]'
          write(*,*)' Enter new value for nw             :'
          read(*,*)nw
          write(*,*)' Enter new value for pw             :'
          read(*,*)pw
      endif

c >>  Wake Discretized Model <<
      call wake(ft,bct,dp,dz,mu,muc,nb,zif,pif,nt,ip,nr,np,nz,nzt,anb)
c      write(*,*) 'nzt',nzt
      do r=1,nr

         ct(r)=ct0

         t0r(r)=t0_0*pi/180.D+0
         t1c(r)=t1c0*pi/180.D+0
         t1s(r)=t1s0*pi/180.D+0

         b0r(r)=b0_0*pi/180.D+0
         b1c(r)=b1c0*pi/180.D+0
         b1s(r)=b1s0*pi/180.D+0

         b0r_0(r)=b0r(r)
         b1c_0(r)=b1c(r)
         b1s_0(r)=b1s(r)

      end do

      vinfx=mu*om*rad
      vinfy=0.D+0
      vinfz=-muc*om*rad

c >>  Rotor/Blade Discretized Model <<
      call rotor(np,nr,ns,pif,rad,flph,rcout,crd,twt,rada,ip, dp,den
     $     ,ibeta,rr_x,rr_y,rr_z,asr0,asr,b0r,b1c,b1s,t0r,t1c,t1s,spx
     $     ,spy,spz,seg,tw,taperst,taper,lock,bcx,bcy ,bcz,cpx,cpy,cpz
     $     ,lsr0,lsr)


c >>  Blade Aerodynamic Model/Near-wake Influence Coefficients <<
      call blade(nr,ns,crd,rcb,bcx,bcy,spx,spy,spz,rad,taperst,taper
     $     ,dp,lengthnw,bicm,nwicm,tw,twt,t0r)

c >>  Screen output: <<
      call wrtscr(nr,nb,nw,pw,ns,nz,np,pif,zif, mu,muc,ft,bct,dp,dz,rad
     $     ,crd,om,ct,asr,lsr,sco, rr_x,rr_y,rr_z)

      rcb=rcb*crd
      do r=1,nr
         do p=1,np
            do w=1,nw
               do k=1,nk
                  rc0(r,p,w,k)=rcb*planf(spx(ns+2-w)/rad,taperst,taper)
               end do
            end do
         end do
      end do



      do r=1,nr
         sn=(-1.D+0)**(r-1)
         do p=1,np
            do w=1,nw-1
               do z=1,nz
                  do k=1,nk
                     gsc(r,p,w,z,k)=-1.D+0*sn*1.5D+0
                     csc(r,p,w,z,k)=0.1D+0*crd*planf(spx(ns+2-w)/rad
     $                    ,taperst,taper)
                  end do
               end do
            end do
         end do
      end do


      do r=1,nr
         do p=1,np
            do w=1,nw
               do k=1,nk
                  gsc(r,p,w,1,k)=0.0D+0
                  gsc(r,p,w,nz,k)=0.0D+0
               end do
            end do
         end do
      end do


c     >>  Prescribed Wake Geometry <<
      do r=1,nr
         do p=1,np
            do s=1,ns
               Cnni(r,p,s)=0.D+0
               Cnqi(r,p,s)=0.D+0
            end do
            do w=1,nw
               do z=1,nzt
                  pcx(r,p,w,z)=0.D+0
                  pcy(r,p,w,z)=0.D+0
                  pcz(r,p,w,z)=0.D+0
               end do
c$$$               dn1=ns
c$$$               if(nw.gt.1) dn1=ns/(nw-1)
c$$$               n1=ns+1 - (w-1)*dn1
               n1=ns+1-4*(w-1)
               rv(r,p,w)=spx(n1)

            end do
         end do
      end do


      call pwake(nr,np,nw,nzt,pif,zif,ip,dp,dz,mu,muc,om, vinfx,vinfy
     $     ,vinfz,asr,lsr,ct,b0r,b1c,b1s,anb,rad,crd,rr_x,rr_y,rr_z,rv
     $     ,twt,taperst,taper, lami,pcx,pcy,pcz)

      call pwake(nr,np,nw,nzt,pif,zif,ip,dp,dz,mu,muc,om, vinfx,vinfy
     $     ,vinfz,asr,lsr,ct,b0r,b1c,b1s,anb,rad,crd,rr_x,rr_y,rr_z,rv
     $     ,twt,taperst,taper, lami,pcx,pcy,pcz)


      do r=1,nr
         do p=1,np
            aza=DBLE(p-1)*dp/DBLE(pif)
            
            aoac=t0r(r)+t1c(r)*cos(aza)+t1s(r)*sin(aza)

            do s=1,ns

               vfx(r,p,s)=vinfx
               vfy(r,p,s)=vinfy
               vfz(r,p,s)=vinfz

               call ftor(r,p,s,dp,ip,b0r,asr,lsr,b1c,b1s,vfx,vfy,vfz
     $              ,vrx,vry,vrz,pif)

               vry(r,p,s)=vry(r,p,s)-om*bcx(r,s)
               vr(r,p,s)=sqrt(vrx(r,p,s)**2+vry(r,p,s)**2+vrz(r,p,s)**2)
               aoai(r,p,s)=atan((vrz(r,p,s)-lami(r)*om*rad)/(-vry(r,p,s)
     $              ))

               aflp(r,p,s)=-bcx(r,s)*om*( -b1c(r)*sin(aza) +
     &              b1s(r)*cos(aza) )/(-vry(r,p,s))
               aflp(r,p,s)=atan(aflp(r,p,s))

               apit(r,p,s)=-om*(-t1c(r)*sin(aza)+t1s(r)*cos(aza))
               apit(r,p,s)=atan(apit(r,p,s)*bcy(r,s)/(-vry(r,p,s)))


               if(abs(twt).gt.100) aoac=t0r(r)/(bcx(r,s)/rad)+t1c(r)
     $              *cos(aza)+t1s(r)*sin(aza)

               aoae(r,p,s)=aoac+aoai(r,p,s)+tw(r,s)+ aflp(r,p,s)+apit(r
     $              ,p,s)

               
               Mach(r,p,s)=abs(vry(r,p,s))/sonic
               Reynolds(r,p,s)=den*vr(r,p,s)* crd*planf(bcx(r,s)/rad
     $              ,taperst,taper)/kinv
               compress(r,p,s)=(1.0D0-Mach(r,p,s)**2.0D+0)
               wlfv(r,p,s)=vr(r,p,s)*sin(aoae(r,p,s)-alpha0)
            end do
            call wlsolve(r,p,ns,compress,bicm,nwicm,wlfv, gb)
         end do                 !p
      end do                    !r


      if(initial.eq.'b'.OR.trm.eq.'r') then
         open(45,file="controls.dat",status='unknown')
         do r=1,nr
            write(45,*)ct(r)    !ct0
            write(45,*)t0r(r)
            write(45,*)t1c(r)
            write(45,*)t1s(r)
            write(45,*)b0r(r)
            write(45,*)b1c(r)
            write(45,*)b1s(r)
         end do
         close(45)
      endif

      call airfoil(bcx,aoae,Mach,Reynolds,nr,np,ns,clairfoil, cdairfoil
     $     ,claoa,compress,0,alpha0)

      do r=1,nr
         call rsp(r,np,ns,nb,ip,dp,den,flph,om,rad,crd,clairfoil,claoa,
     $        nubeta,bmass,ibeta,betap,kbeta,seg,bcx,vr,vry,vrz,gb,Cnni
     $        ,Cnqi, sonic,cdairfoil,aoae,aoai,aoai,compress,bicm,nwicm
     $        ,mu, Mach,Reynolds, taperst,taper,b0r_0,b1c_0,b1s_0,ct_0
     $        ,cqi,cqp,pif,zif,0)
      end do

c-mb prej also uses AOAi based on uniform unflow

      if (nr .eq. 1) then 
         call prej(t0r,t1c,t1s,b0r_0,b1c_0,b1s_0,ct_0,clairfoil,claoa
     $        ,crd,sonic,nr,np,ns,nb,ip,dp,den,flph,om,rad,bicm,nwicm
     $        ,nubeta,bmass,ibeta,betap,kbeta,seg,bcx,tw,vr,aoai,nctrl
     $        ,anb,cpx,twt,taperst,taper,mu,cpz,asr,rr_x,rr_z, jac,pif
     $        ,zif,cdairfoil,vry,vrz,compress,kinv,alpha0)
         call trimtest(nr,np,ns,nb,nctrl,ip,dp,den,flph,om,rad,compress,
     $        nubeta,bmass,ibeta,betap,kbeta,crd,clairfoil,claoa,cttol
     $        ,fltol, seg,tw,bicm,nwicm,bcx,bcy,twt,cdairfoil,kinv,
     $        sonic ,aoai,aoai,jac,Cnni,Cnqi,vr,vry,vrz,gb,anb,cpx,
     $        taperst,taper ,mu,cpz,asr,rr_x,rr_z, t0r,t1c,t1s,b0r,b1c
     $        ,b1s,ct,cqi,cqp,trm ,pif,zif,0,alpha0)
      endif

c >>  Vortex Strengths/Release Points: <<
      call vortex(nr,np,nw,ns,nzt,nk,dp,anb,rad,cpx,spx,gb,rv, gv,1)

      gav=0.D+0
      do r=1,nr
         do p=1,np
            do k=1,nk
               gav=gav+abs(gv(r,p,1,1,k))
            end do
         end do
      end do
      gav=gav/DBLE(nr*np)
      factor=1.0D+0
      do r=1,nr
         do p=1,np
            do w=1,nw
               do z=1,nzt
                  do k=1,nk
                     wage=DBLE(z-1)*dz/DBLE(zif)
                     trb=1.D+0+dcy*factor*abs(gv(r,p,w,z,1)/kinv)
                     trb=1.D+0+dcy*factor*gav/kinv
                     zet0=(rc0(r,p,w,k)/.00855D+0)**2*om/trb/Nb
                     eqtim(r,p,w,z,k)=zet0+wage
                  end do
               end do
            end do
         end do
      end do

      do r=1,nr
         do p=1,np
            do w=1,nw
               do z=1,nzt
                  squire(r,p,w,z)=1.D+0
                  osquire(r,p,w,z)=1.D+0
                  wage=DBLE(z-1)*dz/DBLE(zif)
                  trb=1.D+0+dcy*gav/kinv
                  zet0=(rc0(r,p,w,1)/.00855D+0)**2*om/trb/nb
                  rcore(r,p,w,z)=0.00855D+0*sqrt(trb*(wage+zet0)*nb/om)
               end do
            end do
         end do
      end do


      if (trm .eq. 'y') then
         call jacobian(nr,np,nw,pw,nz,ns,nk,nb,nctrl,dp,dz,ip,anb, rad
     $        ,crd,clairfoil,claoa,dcy,kinv,om,lin,gv,rc0,asr,b0r, b1c
     $        ,b1s,rr_x,rr_y,rr_z,pcx,pcy,pcz,cpx,cpy,cpz, t0r,t1c,t1s
     $        ,vinfx,vinfy,vinfz,bcx,bcy,lengthnw, tw,bicm,nwicm,den
     $        ,flph,nubeta,bmass,ibeta,betap,sonic, kbeta,seg,spx,spy
     $        ,spz,sco,sce,gsc,csc,vr,aoai,aoain,twt,nzt, taperst,taper,
     $        jac,pif,zif,lsr,cdairfoil,rcb,gb,compress,eqtim ,alpha0)

c     >>  Rotor Retrim: <<

         if(initial.eq.'n') then
            call trimtest(nr,np,ns,nb,nctrl,ip,dp,den,flph,om,rad
     $           ,compress,nubeta,bmass,ibeta,betap,kbeta,crd,clairfoil
     $           ,claoa,cttol,fltol,seg,tw,bicm,nwicm,bcx,bcy,twt
     $           ,cdairfoil,kinv,sonic,aoai,aoain,jac,Cnni,Cnqi,vr,vry
     $           ,vrz,gb,anb,cpx,taperst,taper,mu,cpz,asr,rr_x,rr_z,t0r
     $           ,t1c,t1s,b0r,b1c,b1s,ct,cqi,cqp,trm,pif,zif,0,alpha0)


            do r=1,nr
               do p=1,np
                  aza=DBLE(p-1)*dp/DBLE(pif)
                  aoac=t0r(r)+t1c(r)*cos(aza)+t1s(r)*sin(aza)
                  do s=1,ns
                     
                     if(abs(twt).gt.100) aoac=t0r(r)/(bcx(r,s)/rad)
     $                    +t1c(r)*cos(aza)+t1s(r)*sin(aza)
                     
                     aoae(r,p,s)=aoac+tw(r,s)+aoai(r,p,s)
                     wlfv(r,p,s)=vr(r,p,s)*sin(aoae(r,p,s)-alpha0)
                  enddo
                  call wlsolve(r,p,ns,compress,bicm,nwicm,wlfv,gb)  
               enddo
            enddo
            
            call bcp(np,nr,ns,pif,ip,dp,asr,lsr,b0r,b1c,b1s,rad,t0r,t1c
     $           ,t1s,twt,tw, rr_x,rr_y,rr_z,bcx,bcy,bcz, cpx,cpy,cpz)

            call vortex(nr,np,nw,ns,nzt,nk,dp,anb,rad,cpx,spx,gb,rv, gv
     $           ,1)

            call pwake(nr,np,nw,nzt,pif,zif,ip,dp,dz,mu,muc,om, vinfx
     $           ,vinfy,vinfz,asr,lsr,ct,b0r,b1c,b1s,anb,rad,crd,rr_x
     $           ,rr_y,rr_z,rv,twt,taperst,taper, lami,pcx,pcy,pcz)

         endif                  !initial

      endif                     !trm


      if(initial.eq.'b'.OR.trm.eq.'r') then
         open(45,file="controls.dat",status='old')
         do r=1,nr
            read(45,*)ct(r)     !ct0
            read(45,*)t0r(r)
            read(45,*)t1c(r)
            read(45,*)t1s(r)
            read(45,*)b0r(r)
            read(45,*)b1c(r)
            read(45,*)b1s(r)
            write(*,*)'Re initialized trim state: rotor',r
            write(*,*)' C_T',ct(r)
            write(*,*)' t0 ',t0r(r)*180/pi
            write(*,*)' t1c',t1c(r)*180/pi
            write(*,*)' t1s',t1s(r)*180/pi
            write(*,*)' b0 ',b0r(r)*180/pi
            write(*,*)' b1c',b1c(r)*180/pi
            write(*,*)' b1s',b1s(r)*180/pi
         end do
         close(45)
      endif

      if(initial.ne.'n') then
         call readwake(anb,sco,rad,dp,dz,nr,np,nw,pw,nz,nzt,pif,zif, pcx
     $        ,pcy,pcz, eqtim, om, kinv, gav, dcy, Nb)

         call readgb(anb,dp,nr,np,ns,rad,gb,pif,zif)
         call bcp(np,nr,ns,pif,ip,dp,asr,lsr,b0r,b1c,b1s,rad,t0r,t1c,t1s
     $        ,twt,tw, rr_x,rr_y,rr_z,bcx,bcy,bcz, cpx,cpy,cpz)

         call vortex(nr,np,nw,ns,nzt,nk,dp,anb,rad,cpx,spx,gb,rv, gv,1)
      endif


      gav=0.D+0
      do r=1,nr
         do p=1,np
            do k=1,nk
               gav=gav+abs(gv(r,p,1,1,k))
            end do
         end do
      end do
      gav=gav/DBLE(nr*np)


      do j1=1,gdim,1
         do j2=1,gdim,1
            gsv(j1,j2)=0.D+0
         end do
      end do
      gsn=0

      call gridseq(gs,li,dp,gsv,gsn)

      call march(nr,np,nw,pw,nz,nzt,nk,ns,nb,pif,zif,ip,dp,dz, sonic
     $     ,vinfx,vinfy,vinfz,om,asr,b0r,anb,rad,rcout,rcb,gb, cvr,
     $     cvt,trm,cnv,lin,gv,rc0,rv,gs,li,gsv, gsn,kinv,dcy,gav,sco,sce
     $     ,gsc,csc,lengthnw, flph,crd,rpar,rrn,zrn,lami,rr_x,rr_y,rr_z,
     $     spx,spy,spz,cpx,cpy,cpz,pcx,pcy,pcz,sbeta, nctrl,den,nubeta
     $     ,bmass,ibeta,betap,kbeta,clairfoil,claoa, cdairfoil,cttol
     $     ,fltol,t0r,t1c,t1s,ct,cqi,cqp,b1c,b1s,seg, tw,bicm,nwicm,bcx
     $     ,bcy,bcz,jac,method,nfw,twt, taperst,taper,initial, aoai
     $     ,aoain,vr,vry,outplot, eqtim,alpha0,lsr)


      write(*,*)
      write(*,*) "Success: simulation complete!"
      write(*,*)

      stop
      end
