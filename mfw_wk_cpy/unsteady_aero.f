C-$	$Id: unsteady_aero.f,v 1.9.4.2 2005/03/23 21:24:38 shreyas Exp $	

      subroutine unsteady_aero(nr,np,ns,nb,dp,ip,tw,crd,cpx,cpy,cpz,
     $     vinfx, vinfy,vinfz,vibx,viby,vibz,t0r,t1c,t1s,theta, b0r,b1c
     $     ,b1s,beta,betastar,asr,pif, om,bcx,bcy,bcz,sonic,Mach
     $     ,Reynolds, den,taper,taperst,kinv,rad,vr,vry,rcout,Cnni,Cnqi
     $     ,pblade, aoain)

C-$   Subroutine for computing unsteady aerodynamics

      implicit none

      include 'adimpar.inc'
      include 'rdimpar.inc'
      include 'sdimpar.inc'

      integer nr,np,ns,nb,pif,usaero

      integer r,p,s,b,pm1,pm2,dnp,pblade,nblade,pm3

      double precision ip(rdim),dp,tw(rdim,sdim), crd, vinfx,
     $     vinfy, vinfz,asr(rdim),om,pi,sonic,rcout

      double precision t0r(rdim),t1c(rdim), t1s(rdim),
     $     b0r(rdim), b1c(rdim), b1s(rdim),
     $     theta(rdim,adim), beta(rdim,adim), betastar(rdim,adim)

      double precision bcx(rdim,sdim),bcy(rdim,sdim),bcz(rdim,sdim)
      double precision taper,taperst,den,rad,kinv

      double precision cpx(rdim,adim,sdim), cpy(rdim,adim,sdim),
     $     cpz(rdim,adim,sdim), vibx(rdim,adim,sdim),
     $     viby(rdim,adim,sdim), vibz(rdim,adim,sdim),
     $     vr(rdim,adim,sdim),vrx(rdim,adim,sdim),vry(rdim,adim,sdim),
     $     vrz(rdim,adim,sdim),vfx(rdim,adim,sdim),vfy(rdim,adim,sdim),
     $     vfz(rdim,adim,sdim),Mach(rdim,adim,sdim),
     $     Reynolds(rdim,adim,sdim), aoain(rdim,adim,sdim)

      double precision dbetadt,dalpdt,psi,machno,pgcf,stime
      double precision ahdot(rdim,adim,sdim),alpdot(rdim,adim,sdim)
      double precision ualpha(rdim,adim,sdim),delalp(rdim,adim,sdim)
      double precision xcoeff(rdim,adim,sdim),ycoeff(rdim,adim,sdim)
      double precision uacirc(rdim,adim,sdim),uancir(rdim,adim,sdim)
      double precision secaoa(rdim,adim,sdim)

      double precision tinc, temp1, temp2,deltim, ednca, edncq,
     $     kapalp, kappaq
      double precision Cnni(rdim,adim,sdim),Cnqi(rdim,adim,sdim)
      double precision ncdalp(rdim,adim,sdim)
      double precision ncdq(rdim,adim,sdim), delqdt(rdim,adim,sdim)
      double precision qval(rdim,adim,sdim)

      double precision ub1,ub2,ua1,ua2,uf1,uf2,ug1,ug2
      double precision velox,veloy,veloz,planf
      external velox,veloy,veloz,planf

      common/pidata/pi
      save/pidata/

C-$   from usa.input coefficients for Duhamel integration
      common/usadata/ub1,ub2,ua1,ua2,uf1,uf2,ug1,ug2,usaero
      save/usadata/

C-$   unsteady aerodynamics variables - needed to modify aoae and Cn
      common/unsdata/ualpha,delalp,ahdot,alpdot,xcoeff,ycoeff,uacirc,
     $     secaoa
      save/unsdata/

      common/unncdata/qval,delqdt,ncdalp,ncdq
      save/unncdata/
      
      dnp=np/nb
      tinc=crd/sonic
      deltim=dp/om
      open(111,file='noncirc.dat',status='unknown')
      do r=1,nr
         do nblade=1,nb
            p=pblade+(nblade-1)*dnp
            pm1=p-1
            pm2=p-2
            pm3=p-3
            if(pm1 .le. 0) pm1=pm1+np
            if(pm2 .le. 0) pm2=pm2+np
            if(pm3 .le. 0) pm3=pm3+np
            dbetadt=(3.0D+0*beta(r,pm1)-4.0D+0*beta(r,pm2)+beta(r,pm3))
     $           *0.5D+0*om/dp
            dalpdt=(3.0D+0*theta(r,p)-4.0D+0*theta(r,pm1)+
     $           theta(r,pm2))*0.5D+0*om/dp
            
            do s=1,ns
               dalpdt=(3.0D+0*secaoa(r,pm1,s)-4.0D+0*secaoa(r,pm2,s)+
     $              secaoa(r,pm3,s))*0.5D+0*om/dp
               
               Mach(r,p,s)=abs(vry(r,p,s))/sonic
               Reynolds(r,p,s)=vr(r,p,s)*crd*planf(bcx(r,s)/rad, taperst
     $              ,taper)/kinv
               
C-$   AoA contribution from blade flapping (hdot effect)
               ahdot(r,p,s)=atan(-dbetadt*(bcx(r,s))
     $              /(-vry(r,p,s)))
               
C-$   AoA contribution from pitch (alpha and alpha dot effect)
               alpdot(r,p,s)=atan(dalpdt*bcy(r,s)/(-vry(r,p,s)))
               
C-$   unsteady instantaneous angle of attack 
               ualpha(r,p,s)=secaoa(r,p,s)+tw(r,s)+alpdot(r,p,s)
     $              +ahdot(r,p,s)
               delalp(r,p,s)=ualpha(r,p,s)-ualpha(r,pm1,s)
     $              +aoain(r,pm1,s)-aoain(r,pm2,s)
               
               machno=vr(r,p,s)/sonic
               pgcf=(1.0D+0-machno**2)
               stime=2.0D+0*vr(r,p,s)*dp/(om*crd)
               
C-$   Circulatory effect corrections
               xcoeff(r,p,s)=xcoeff(r,pm1,s)*exp(-ub1*pgcf*stime)+
     $              delalp(r,p,s)*ua1*exp(-0.5D+0*ub1*pgcf*stime)
               ycoeff(r,p,s)=ycoeff(r,pm1,s)*exp(-ub2*pgcf*stime)+
     $              delalp(r,p,s)*ua2*exp(-0.5D+0*ub2*pgcf*stime)
               
               uacirc(r,p,s)=ualpha(r,p,s)-xcoeff(r,p,s)
     $              -ycoeff(r,p,s)
               
C-$   Non-circulatory effect contributions
               qval(r,p,s)=dalpdt
               temp1=(1.0D+0-machno)
               temp2=sqrt(pgcf)*machno**2.0D+0*(ua1*ub1+ua2*ub2)*pi
               kapalp=1.0D+0/(temp1+temp2)
               kappaq=1.0D+0/(temp1+2.0D+0*temp2)
               
               ednca=exp(-deltim/(kapalp*tinc))
               edncq=exp(-deltim/(kappaq*tinc))
               
               ncdalp(r,p,s)=ncdalp(r,pm1,s)*ednca+(delalp(r,p,s)
     $              -delalp(r,pm1,s))*sqrt(ednca)/deltim
               Cnni(r,p,s)=(4.0D+0*kapalp*tinc/machno)*
     $              (delalp(r,p,s)/deltim -ncdalp(r,p,s))
               
               delqdt(r,p,s)=qval(r,p,s)-qval(r,pm1,s)
               ncdq(r,p,s)=ncdq(r,pm1,s)*edncq+(delqdt(r,p,s)
     $              -delqdt(r,pm1,s))*sqrt(edncq)/deltim
               Cnqi(r,p,s)=(kappaq*tinc/machno)*(delqdt(r,p,s)
     $              /deltim-ncdq(r,p,s))
               write(111,10) p,s,xcoeff(r,p,s)*180.0D+0/pi,
     $              ycoeff(r,p,s)*180.0D+0/pi,
     $              Cnni(r,p,s),Cnqi(r,p,s)
            end do
         end do
      end do
      close(111)

 10   format(2(i3,1x),4(E15.8,1x))
      return
      end
