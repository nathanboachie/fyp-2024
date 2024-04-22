C-$	$Id: march.f,v 1.18.2.4 2005/06/11 15:57:37 shreyas Exp $	

      subroutine march(nr,np,nw,pw,nz,nzt,nk,ns,nb,pif,zif,ip,dp,dz,
     $     sonic,vinfx,vinfy,vinfz,om,asr,b0r,anb,rad,rcout,rcb,gb, cvr
     $     ,cvt,trm,cnv,lin,gv,rc0,rv,gs,li,gsv, gsn,kinv,dcy,gav,sco
     $     ,sce,gsc,csc,lengthnw, flph,crd,rpar,rrn,zrn,lami,rr_x,rr_y
     $     ,rr_z, spx,spy,spz,cpx,cpy,cpz,pcx,pcy,pcz,sbeta, nctrl,den
     $     ,nubeta,bmass,ibeta,betap,kbeta,clairfoil,claoa, cdairfoil
     $     ,cttol,fltol,t0r,t1c,t1s,ct,cqi,cqp,b1c,b1s,seg, tw,bicm
     $     ,nwicm,bcx,bcy,bcz,jac,method,nfw,twt, taperst,taper,initial,
     $     aoai,aoain,vr,vry,outplot, eqtim,alpha0,lsr)

C-$   The engine. Calls the numerical schemes based on 'method' variable
C-$   in user.input. Runs the code for a specified number of iterations
C-$   or until convergence is achieved (in relaxation methods). 



      IMPLICIT none


      include 'adimpar.inc'
      include 'gdimpar.inc'
      include 'jdimpar.inc'
      include 'kdimpar.inc'
      include 'rdimpar.inc'
      include 'sdimpar.inc'
      include 'wdimpar.inc'
      include 'zdimpar.inc'


      INTEGER nr,np,nw,pw,nz,nzt,nk,ns,nb,pif,zif,gs,li, r,p,w,z,n0,nfw
     $     ,gsno,n1,n,ii,gsn,nctrl,s,k, periodic

      DOUBLE PRECISION ip(rdim), dp,dz,vinfx,vinfy,vinfz,om,anb,rad
     $     ,rcout, rcb,lin,kinv,dcy,gav,flph0,flph,crd,rpar,rrn, zrn,cnv
     $     ,cvt,sce,den,nubeta,bmass,ibeta,betap,sbeta, kbeta,cttol
     $     ,fltol,aza,aoac,sn,az, ze,wage,tt11,tt21,tt31,b0,sonic,pbar
     $     ,qbar,pi,mu

      DOUBLE PRECISION rbc(rdim,wdim),zbc(rdim,wdim), avgzbc(rdim,wdim)

      DOUBLE PRECISION asr(rdim),ct(rdim),b0r(rdim),lami(rdim),
     $     rr_x(rdim),rr_y(rdim),rr_z(rdim), spx(sdim),spy(sdim)
     $     ,spz(sdim), t0r(rdim),t1c(rdim),t1s(rdim), b1c(rdim),b1s(rdim
     $     ),seg(sdim),lsr(rdim)

      DOUBLE PRECISION cqi(rdim)
      DOUBLE PRECISION cqp(rdim)

      DOUBLE PRECISION gsv(gdim,gdim),rms(rdim,wdim),rms0(rdim,wdim),
     $     tw(rdim,sdim),nwicm(sdim,sdim), bcx(rdim,sdim),bcy(rdim,sdim)
     $     ,bcz(rdim,sdim), jac(jdim,jdim)

      DOUBLE PRECISION bicm(sdim,sdim)

      DOUBLE PRECISION rv(rdim,adim,wdim), gb(rdim,adim,sdim), cpx(rdim
     $     ,adim,sdim), cpy(rdim,adim,sdim), cpz(rdim,adim,sdim),
     $     vibx(rdim,adim,sdim), viby(rdim,adim,sdim), vibz(rdim,adim
     $     ,sdim), vibxin(rdim,adim,sdim), vibyin(rdim,adim,sdim),
     $     vibzin(rdim,adim,sdim), bvx(rdim,adim,sdim), bvy(rdim,adim
     $     ,sdim), bvz(rdim,adim,sdim), aoai(rdim,adim,sdim), aoain(rdim
     $     ,adim,sdim), aflp(rdim,adim,sdim), apit(rdim,adim,sdim),
     $     aoae(rdim,adim,sdim), aoanc(rdim,adim,sdim), Cnni(rdim,adim
     $     ,sdim), Cnqi(rdim,adim,sdim), vry(rdim,adim,sdim), vrz(rdim
     $     ,adim,sdim), vr(rdim,adim,sdim)

      DOUBLE PRECISION gv(rdim,adim,wdim,zdim,kdim), gvold(rdim,adim
     $     ,wdim,kdim), rc0(rdim,adim,wdim,kdim), ppx(rdim,adim,wdim
     $     ,zdim), ppy(rdim,adim,wdim,zdim), ppz(rdim,adim,wdim,zdim),
     $     pcx(rdim,adim,wdim,zdim), pcy(rdim,adim,wdim,zdim), pcz(rdim
     $     ,adim,wdim,zdim), pox(rdim,adim,wdim,zdim), poy(rdim,adim
     $     ,wdim,zdim), poz(rdim,adim,wdim,zdim), eqtim(rdim,adim,wdim
     $     ,zdim,kdim), Cldata(rdim,adim,sdim), ClMsqr(rdim,adim,sdim),
     $     liftx(rdim,adim,sdim), lifty(rdim,adim,sdim), liftz(rdim,adim
     $     ,sdim),wlfv(rdim,adim,sdim)

      CHARACTER*1 sco,trm,cvr,method

      DOUBLE PRECISION gsc(rdim,adim,wdim,zdim,kdim), csc(rdim,adim,wdim
     $     ,zdim,kdim)

      DOUBLE PRECISION twt,alpha0

      DOUBLE PRECISION gammaRMS


      DOUBLE PRECISION clairfoil(rdim,adim,sdim), cdairfoil(rdim,adim
     $     ,sdim), claoa(rdim,adim,sdim),compress(rdim,adim,sdim)
      DOUBLE PRECISION Mach(rdim,adim,sdim),Reynolds(rdim,adim,sdim)

      INTEGER lengthnw(sdim)
      DOUBLE PRECISION taper,taperst
      LOGICAL outplot

      DOUBLE PRECISION beta(rdim,adim),betastar(rdim,adim)
      DOUBLE PRECISION betaRMS(rdim)
      DOUBLE PRECISION theta(rdim,adim)

      DOUBLE PRECISION amom(rdim,adim),thrust(rdim,adim) ,inflow(rdim
     $     ,adim) ,torque(rdim,adim,2)

      INTEGER itmp
      DOUBLE PRECISION  eR ,lock,tmp
      INTEGER rotgeo


      double precision ahdot(rdim,adim,sdim),alpdot(rdim,adim,sdim)
      double precision ualpha(rdim,adim,sdim),delalp(rdim,adim,sdim)
      double precision xcoeff(rdim,adim,sdim),ycoeff(rdim,adim,sdim)
      double precision uacirc(rdim,adim,sdim),uancir(rdim,adim,sdim)
      double precision secaoa(rdim,adim,sdim)
      double precision ncdalp(rdim,adim,sdim)
      double precision ncdq(rdim,adim,sdim), delqdt(rdim,adim,sdim)
      double precision qval(rdim,adim,sdim)
      double precision cpitch(rdim,adim)

      double precision cnplag(rdim,adim,sdim),deltalift(rdim,adim,sdim)

      CHARACTER*1 initial

      common/pidata/pi
      save/pidata/

      common/ssmdata/pbar,qbar
      save/ssmdata/

      common/unsdata/ualpha,delalp,ahdot,alpdot,xcoeff,ycoeff,uacirc,
     $     secaoa
      save/unsdata/

      common/unncdata/qval,delqdt,ncdalp,ncdq
      save/unncdata/

      common/dynstall/cnplag,deltalift
      save/dynstall/

      common/rotdata/rotgeo
      save/rotdata/


      mu=vinfx/om/rad

      gammaRMS=0.D+0
        do r=1,nr
        do p=1,np
        aza=DBLE(p-1)*dp/DBLE(pif)
        beta(r,p)    =b0r(r)+b1c(r)*cos(aza)+b1s(r)*sin(aza)
        betastar(r,p)=      -b1c(r)*sin(aza)+b1s(r)*cos(aza)

        do s=1,ns
           compress(r,p,s)=1.D+0
           ualpha(r,p,s)=t0r(r)+t1c(r)*cos(aza)+t1s(r)*sin(aza)+tw(r,s)
           delalp(r,p,s)=0.0D+0
           xcoeff(r,p,s)=0.0D+0
           ycoeff(r,p,s)=0.0D+0
           uacirc(r,p,s)=0.0D+0
           secaoa(r,p,s)=0.0D+0
           ahdot(r,p,s)=0.0D+0
           alpdot(r,p,s)=0.0D+0
           qval(r,p,s)=0.0D+0
           delqdt(r,p,s)=0.0D+0
           ncdq(r,p,s)=0.0D+0
           ncdalp(r,p,s)=0.0D+0
           cnplag(r,p,s)=0.0D+0
           deltalift(r,p,s)=0.0D+0
        end do
      end do
      end do

      if(initial.ne.'n') then
         open(66,file='flap.data',status='old')
         do r=1,nr
            do p=1,np
               read(66,'(I3,1X,F7.3,2(1X,E13.6))')itmp,tmp,
     &              beta(r,p),betastar(r,p)
            end do
         end do
         close(66)
      end if

c-mb initializing flap angle and velocity

      call rindvt(nr,np,nw,pw,nz,ns,nk,dp,dz,ip,anb,crd,rad, dcy,kinv,om
     $     ,lin,gav,gv,rc0,asr,b0r,b1c,b1s, rr_x,rr_y,rr_z,pcx,pcy,pcz
     $     ,cpx,cpy,cpz, sco,sce,gsc,csc,vibx,viby,vibz, vibxin,vibyin
     $     ,vibzin,compress,lengthnw, pif,zif,bvx,bvy,bvz,rcb,gb,bicm
     $     ,nwicm,.false.,eqtim,lsr)

      do r=1,nr
         do p=1,np
            aza=DBLE(p-1)*dp/DBLE(pif)
            aoac=t0r(r)+t1c(r)*cos(aza)+t1s(r)*sin(aza)
            do s=1,ns
               if(abs(twt).gt.100) aoac=t0r(r)/(bcx(r,s)/rad)+t1c(r)
     $              *cos(aza)+t1s(r)*sin(aza)
               call aoa(r,p,s,cpx,cpy,cpz,vinfx,vinfy,vinfz,vibx,viby
     $              ,vibz,vibxin,vibyin,vibzin,nr,np,ns,dp,ip,b0r,b1c
     $              ,b1s,t1c,t1s,asr,om,bcx,bcy,aoac,tw,aza,pif,aoai
     $              ,aoain,aflp,apit,vr,vry,vrz,aoae,aoanc,lsr)
               aoae(r,p,s)=aoae(r,p,s)-aoai(r,p,s)+aoain(r,p,s)
            end do
         end do
      end do
      
      call airfoil(bcx,aoae,Mach,Reynolds,nr,np,ns,clairfoil, cdairfoil
     $     ,claoa,compress,0,alpha0)

      do r=1,nr
         do p=1,np
            call  rsppsi(r,np,ns,nb,ip,dp,den,flph,om,rad,crd,clairfoil
     $           ,claoa,nubeta,bmass,ibeta,betap,kbeta,seg,bcx,vr,vry
     $           ,vrz,gb,Cnni,Cnqi,sonic,cdairfoil,aoae,aoai,aoain
     $           ,compress,BICM,NWICM,Mach,Reynolds,inflow,aoanc,beta
     $           ,betastar,taperst,taper,ct,cqi,cqp,pif,zif,amom,thrust
     $           ,torque,p,Cldata, ClMsqr, liftx, lifty, liftz,cpitch) 
         end do
      end do

      eR=flph/rad
      if(eR.lt.0.D+0) eR=0.D+0
      lock=(crd*den*2.D+0*pi*rad**4)/ibeta !note Clalpha=2pi
      lock=lock/sqrt(1.D+0-(0.75D+0*om*rad/sonic)**2.D+0) !general

      do r=1,nr
         do p=1,np
            aza=ip(r)+DBLE(p-1)*dp/DBLE(pif) - ip(r)
            theta(r,p)=t0r(r)+t1c(r)*cos(aza)+t1s(r)*sin(aza)
         end do
      end do


      call rtof(nr,np,ns,pif,dp,ip,b0r,b1c,b1s,asr,lsr, rr_x,rr_y,rr_z
     $     ,spx,spy,spz, bvx,bvy,bvz)



      do r=1,nr
         do w=1,nw
            rms(r,w)=0.D+0
            rms0(r,w)=0.D+0
         end do
      end do

      n0=1

      gsno=0
      if(gs .eq. 1) then
         gsno=1

         write(*,*)
         write(*,22)' Grid Sequencing Block      .......... : ',gsno
         write(*,*)

         call gsassign(gsn,gsno,gsv,gs,li,dp,dz,pif,zif,n1)
         nfw=n1
      endif

      periodic=1
      if(method.eq.'t') periodic=0

 11   do n=n0,nfw

         write(*,*)
         if(method.eq.'r') write(*,22)
     &        ' Free-Wake Relaxation Iteration ------- : ',n
         if(method.eq.'t') write(*,22)
     &        ' Free-Wake Transient PC Iteration ----- : ',n
         if(method.eq.'T') write(*,22)
     &        ' Free-Wake Transient PIPC Iteration --- : ',n

         call udgeom(nr,np,nw,nz,nzt,pox,poy,poz,ppx,ppy,ppz, pcx,pcy
     $        ,pcz,n,rbc,zbc,avgzbc,dz,zif)      


         if(method.eq.'r') then 
            call pipcab(nr,np,nw,pw,nz,nzt,nk,ns,nb,pif,zif,ip,dp,dz,
     $           sonic,vinfx,vinfy,vinfz,om,asr,b0r,anb,rad,rcb,gb, t0r
     $           ,t1c,t1s,b1c,b1s, lin,gv,rc0,rv,kinv,dcy,gav,sco,sce
     $           ,gsc,csc, rpar,rrn,zrn,rbc,zbc, avgzbc,lami,rr_x,rr_y
     $           ,rr_z,spx,spy,spz, cpx,cpy,cpz,pcx,pcy,pcz, pox,poy,poz
     $           ,ppx,ppy,ppz,bvx,bvy,bvz, eqtim)

         elseif(method.eq.'t') then 
            if(n.eq.1) then
               open(44,file='timeflap.dat',status='new')
               open(55,file='timectcq.dat',status='new')
               open(77,file='fcinflow.dat',status='new')
               open(78,file='inflowhist.dat',status='new')
               open(79,file='ctrlhist.dat',status='new')
               open(56,file='clhist.dat', status='new')
               open(57,file='lifthist.dat', status='new')
               open(58,file='ctrlinp.data', status='old')
               open(59,file='velocity.data', status='old')
               open(68,file='orient.data', status='old')
               open(60,file='liftxcomp.dat', status='new')
               open(64,file='liftycomp.dat',status='new')
               open(67,file='liftzcomp.dat',status='new')
               open(69,file='timeflight.dat',status='new')
               open(70,file='rates.data',status='old')
            end if
            call timepc(nr,np,nw,pw,nz,nzt,nk,ns,nb,pif,zif,ip,dp,dz,
     $           sonic,vinfx,vinfy,vinfz,om,asr,b0r,anb,rad,rcout,rcb,gb
     $           , t0r,t1c,t1s,flph,nubeta,kbeta,bmass,ibeta,betap, lin
     $           ,gv,gvold,rc0,rv,kinv,dcy,gav,sco,sce,gsc,csc, rpar,rrn
     $           ,zrn,rbc,zbc,den,compress,lengthnw, clairfoil,cdairfoil
     $           ,claoa,Mach,Reynolds,crd, avgzbc,lami,rr_x,rr_y,rr_z
     $           ,spx,spy,spz,seg, tw,twt,aoai,aoain,bcx,bcy,bcz,bicm
     $           ,nwicm, vr,vry,taper,taperst,cnni,cnqi,ct,b1s,b1c, cpx
     $           ,cpy,cpz,pcx,pcy,pcz,beta,betastar,theta, cqi,cqp
     $           ,gammaRMS,betaRMS, pox,poy,poz,ppx,ppy,ppz,bvx,bvy,bvz
     $           ,n,outplot, eqtim, Cldata, ClMsqr, liftx, lifty, liftz
     $           ,trm,alpha0,lsr)

         else
            call prescr(nr,np,nw,pw,nz,nzt,nk,ns,nb,pif,zif,ip,dp,dz,
     $           sonic,vinfx,vinfy,vinfz,om,asr,b0r,b1c,b1s,ct, anb,rad
     $           ,rcb,gb,crd,twt, lin,gv,rc0,rv,kinv,dcy,gav,sco,sce,gsc
     $           ,csc, rpar,rrn,zrn,rbc,zbc, avgzbc,lami,rr_x,rr_y,rr_z
     $           ,spx,spy,spz, cpx,cpy,cpz,pcx,pcy,pcz, taperst,taper,
     $           pox,poy,poz,ppx,ppy,ppz,bvx,bvy,bvz,n)

         endif
         
         do r=1,nr
            do p=1,np
               do w=1,nw
                  do k=1,nk
                     gvold(r,p,w,k)=gv(r,p,w,1,k)
                  end do
               end do
            end do
         end do

         call rindvt(nr,np,nw,pw,nz,ns,nk,dp,dz,ip,anb,crd,rad, dcy
     $        ,kinv,om,lin,gav,gv,rc0,asr,b0r,b1c,b1s, rr_x,rr_y,rr_z
     $        ,pcx,pcy,pcz,cpx,cpy,cpz, sco,sce,gsc,csc,vibx,viby ,vibz,
     $        vibxin,vibyin,vibzin,compress,lengthnw, pif,zif ,bvx,bvy
     $        ,bvz,rcb,gb,bicm,nwicm,.false.,eqtim,lsr)

         do r=1,nr
            do p=1,np
               aza=DBLE(p-1)*dp/DBLE(pif)
               aoac=t0r(r)+t1c(r)*cos(aza)+t1s(r)*sin(aza)
               do s=1,ns
                  if(abs(twt).gt.100) aoac=t0r(r)/(bcx(r,s)/rad)
     $                 +t1c(r)*cos(aza)+t1s(r)*sin(aza)
                  call aoa(r,p,s,cpx,cpy,cpz,vinfx,vinfy,vinfz,vibx
     $                 ,viby,vibz,vibxin,vibyin,vibzin,nr,np,ns,dp,ip
     $                 ,b0r,b1c,b1s,t1c,t1s,asr,om,bcx,bcy,aoac,tw ,aza
     $                 ,pif,aoai,aoain,aflp,apit,vr,vry,vrz,aoae ,aoanc
     $                 ,lsr)
                  wlfv(r,p,s)=vr(r,p,s)*sin((aoae(r,p,s)-alpha0))
               end do
               call wlsolve(r,p,ns,compress,bicm,nwicm,wlfv,gb)  
            end do
         end do
         call vortex(nr,np,nw,ns,nzt,nk,dp,anb,rad,cpx,spx,gb,rv, gv
     $        ,periodic)

         call wrtwke(anb,sco,rad,dp,dz,nr,np,nw,pw,nz,nzt,pif,zif, pcx
     $        ,pcy,pcz)

         call wrtwkmb(anb,sco,rad,dp,dz,nr,np,nw,nz,nzt,pif,zif, pcx,pcy
     $        ,pcz,gv, eqtim, nb, om, gav, kinv, dcy)

         call wrtgbmb(anb,dp,nr,np,ns,rad,cpx,cpz,rr_x,rr_z, asr,lsr 
     $        ,gb,pif,zif)

         call l2norm(nr,np,nw,nz,pox,poy,poz,pcx,pcy,pcz, gs,gsv,gsno
     $        ,rms)

         if(n .eq. 1) then
            do r=1,nr
               do w=1,nw
                  rms0(r,w)=rms(r,w)
                  if(rms0(r,w) .eq. 0.D+0) rms0(r,w)=1.D+0
               end do
            end do
         endif


         ii=0
         do r=1,nr
            do w=1,nw
               rms(r,w)=rms(r,w)/rms0(r,w)
               write(*,'(a15,a26,e11.6)') ' L_2 Norm [RMS]'
     $              ,' ---------------------- : ',rms(r,w)
               if(rms(r,w)*rms0(r,w) .le. cnv) then
                  if(n.gt.1) ii=ii+1
               endif
            end do
         end do
         write(*,*)

         if((mod(n,nb) .eq. 0) .and. trm .ne. 'n') then
            call rindvt(nr,np,nw,pw,nz,ns,nk,dp,dz,ip,anb,crd,rad,
     $           dcy,kinv,om,lin,gav,gv,rc0,asr,b0r,b1c,b1s, rr_x
     $           ,rr_y,rr_z,pcx,pcy,pcz,cpx,cpy,cpz, sco,sce,gsc,csc
     $           ,vibx,viby,vibz, vibxin,vibyin,vibzin,compress
     $           ,lengthnw, pif,zif,bvx,bvy,bvz,rcb,gb,bicm,nwicm,
     $           .false., eqtim,lsr)

            do r=1,nr
               do p=1,np
                  aza=DBLE(p-1)*dp/DBLE(pif)
                  aoac=t0r(r)+t1c(r)*cos(aza)+t1s(r)*sin(aza)
                  do s=1,ns
                     if(abs(twt).gt.100) aoac=t0r(r)/(bcx(r,s)/rad)
     $                    +t1c(r)*cos(aza)+t1s(r)*sin(aza)
                     call aoa(r,p,s,cpx,cpy,cpz,vinfx,vinfy,vinfz
     $                    ,vibx,viby,vibz,vibxin,vibyin,vibzin,nr,np
     $                    ,ns,dp,ip,b0r,b1c,b1s,t1c,t1s,asr,om,bcx
     $                    ,bcy,aoac,tw,aza,pif,aoai,aoain,aflp,apit
     $                    ,vr,vry,vrz,aoae,aoanc,lsr)
                  end do
               end do
            end do

c-    mb here the Gamma_bound will change
c-    mb      so store the old vortex strengths
            do r=1,nr
               do p=1,np
                  do w=1,nw
                     do k=1,nk
                        gvold(r,p,w,k)=gv(r,p,w,1,k)
                     end do
                  end do
               end do
            end do

            call jacobian(nr,np,nw,pw,nz,ns,nk,nb,nctrl,dp,dz,ip,anb,
     $           rad,crd,clairfoil,claoa,dcy,kinv,om,lin,gv,rc0,asr,b0r,
     $           b1c,b1s,rr_x,rr_y,rr_z,pcx,pcy,pcz,cpx,cpy,cpz, t0r,t1c
     $           ,t1s,vinfx,vinfy,vinfz,bcx,bcy,lengthnw, tw,bicm,nwicm
     $           ,den,flph,nubeta,bmass,ibeta,betap,sonic, kbeta,seg,spx
     $           ,spy,spz,sco,sce,gsc,csc,vr,aoai,aoain,twt,nzt, taperst
     $           ,taper,jac,pif,zif,lsr,cdairfoil,rcb,gb,compress,eqtim
     $           ,alpha0)
            call trimtest(nr,np,ns,nb,nctrl,ip,dp,den,flph,om,rad
     $           ,compress,nubeta,bmass,ibeta,betap,kbeta,crd
     $           ,clairfoil,claoa,cttol,fltol,seg,tw,bicm,nwicm,bcx
     $           ,bcy,twt,cdairfoil,kinv,sonic,aoai,aoain,jac,Cnni
     $           ,Cnqi,vr,vry,vrz,gb,anb,cpx,taperst,taper,mu,cpz,asr
     $           ,rr_x,rr_z,t0r,t1c,t1s,b0r,b1c,b1s,ct,cqi,cqp,trm
     $           ,pif,zif,n,alpha0)


            call bcp(np,nr,ns,pif,ip,dp,asr,lsr,b0r,b1c,b1s,rad,t0r,t1c
     $           ,t1s,twt,tw,rr_x,rr_y,rr_z,bcx,bcy,bcz,cpx,cpy,cpz)


            do r=1,nr
               sn=1.D+0
               if(rotgeo.eq.1 .or. rotgeo.eq.2) sn=(-1.D+0)**(r-1)
               do p=1,np
                  az=DBLE(p-1)*dp/DBLE(pif)
                  b0=b0r(r)+b1c(r)*cos(az)+b1s(r)*sin(az)
                  do w=1,nw
                     do z=1,1   !nzt
                        ze=DBLE(z-1)*dz/DBLE(zif)
                        wage=sn*(az + ip(r) -ze)

                        tt11=cos(-asr(r))*cos(wage)*cos(b0)
     $                   -sn*sin(-asr(r))*sin(lsr(r))*sin(wage)*cos(b0) 
     $                   +sin(-asr(r))*cos(lsr(r))*sin(b0)
                        tt21=sn*cos(lsr(r))*sin(wage)*cos(b0)
     $                     +sin(lsr(r))*sin(b0)
                        tt31=-sin(-asr(r))*cos(wage)*cos(b0)
     $                   -sn*cos(-asr(r))*sin(lsr(r))*sin(wage)*cos(b0)
     $                   +cos(-asr(r))*cos(lsr(r))*sin(b0)
     


                        pox(r,p,w,z)=pcx(r,p,w,z)
                        poy(r,p,w,z)=pcy(r,p,w,z)
                        poz(r,p,w,z)=pcz(r,p,w,z)


                        pcx(r,p,w,z)=rr_x(r)+rv(r,p,w)*tt11
                        pcy(r,p,w,z)=rr_y(r)+rv(r,p,w)*tt21
                        pcz(r,p,w,z)=rr_z(r)+rv(r,p,w)*tt31

                     end do
                  end do
               end do
            end do
         endif

         open(66,file='flap.dat',status='unknown')
         do r=1,nr
            do p=1,np
               aza=DBLE(p-1)*dp/DBLE(pif)
               if(method.eq.'r') then
                  beta(r,p)    =b0r(r)+b1c(r)*cos(aza)+b1s(r)*sin(aza)
                  betastar(r,p)=      -b1c(r)*sin(aza)+b1s(r)*cos(aza)
               end if
               write(66,'(I3,1X,F7.3,4(1X,E13.6))')p,dp*180.D+0/pi
     $              *DBLE(p-1),beta(r,p),betastar(r,p),b0r(r)+b1c(r)
     $              *cos(aza)+b1s(r)*sin(aza),-b1c(r)*sin(aza)+b1s(r)
     $              *cos(aza)
            end do
         end do
         close(66)

         if(ii .eq. nr*nw .AND. gammaRMS.le.cnv*10.D+0) then
            n0=n+1
            go to 17
         endif
         if(ii .eq. nr*nw) then
            n0=n+1
            go to 17
         endif
      enddo


      n0=n

 17   if(gs .eq. 1) then
         if(gsno .le. gsn) then

            write(*,*)
            write(*,22)' Grid Sequencing Block      .......... : ',gsno
            write(*,*)

            call gsassign(gsn,gsno,gsv,gs,li,dp,dz,pif,zif,n1)
            nfw=n+n1-1
            go to 11
         endif
      endif


      call wrtwkmb(anb,sco,rad,dp,dz,nr,np,nw,nz,nzt,pif,zif, pcx,pcy
     $     ,pcz,gv, eqtim, nb, om, gav, kinv, dcy)

      call wrtgbmb(anb,dp,nr,np,ns,rad,cpx,cpz,rr_x,rr_z,asr,lsr
     $     ,gb,pif,zif)

      open(66,file='flap.dat',status='unknown')
      do r=1,nr
         do p=1,np
            aza=ip(r)+DBLE(p-1)*dp/DBLE(pif) - ip(r)
            if(method.eq.'r') then
               beta(r,p)    =b0r(r)+b1c(r)*cos(aza)+b1s(r)*sin(aza)
               betastar(r,p)=      -b1c(r)*sin(aza)+b1s(r)*cos(aza)
            end if
            write(66,'(I3,1X,F7.3,4(1X,E13.6))')p,dp*180.D+0/pi*DBLE(p-1
     $           ),beta(r,p),betastar(r,p),b0r(r)+b1c(r)*cos(aza)+b1s(r)
     $           *sin(aza),-b1c(r)*sin(aza)+b1s(r)*cos(aza)
         end do
      end do
      close(66)


22    format(a41,i10)

      if(method.eq.'t') then
         close(44)              !timeflap.dat
         close(55)              !timectcq.dat
         close(77)              !fcinflow.dat
         close(78)              !inflowhist.dat
         close(79)              !ctrlhist.dat
         close(56)              !clhist.dat
         close(57)              !lifthist.dat
         close(58)              !ctrlinp.data
         close(59)              !velocity.data
         close(68)              !orient.data
         close(69)              !timeflight.dat
         close(70)              !rates.data
      end if

      return
      end
