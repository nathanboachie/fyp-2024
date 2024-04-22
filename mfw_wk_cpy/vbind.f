C-$	$Id: vbind.f,v 1.1.1.1.6.2 2005/06/11 15:57:37 shreyas Exp $	

      subroutine vbind(r0,p0,s0,r1,p1,w1,z1,np,nz,nk,dz, gv,rc,lin,trb
     $     ,om,rc0,cpx,cpy,cpz, pcx,pcy,pcz,vbx,vby,vbz,zif,Nb)

c     Vortex/blade induced velocity calculations.

      IMPLICIT none

      include 'adimpar.inc'
      include 'kdimpar.inc'
      include 'rdimpar.inc'
      include 'sdimpar.inc'
      include 'wdimpar.inc'
      include 'zdimpar.inc'

      INTEGER r0,p0,s0,r1,p1,w1,z1,np,nk,k,p3,p4,chi,chi2

      DOUBLE PRECISION dz,lin,trb,om,vbx,vby,vbz,wage,
     &     r1x,r1y,r1z,mr1,r2x,r2y,r2z,mr2,
     &     lx,ly,lz,ml,theta1,theta2,h,ux,uy,uz,
     &     umag,ca,cb,cot1,cot2,gam1,gam2,rc,vvp,
     &     icv,icv2,icvx,icvy,icvz,icvx2,icvy2,icvz2,
     &     vorp,sn

      DOUBLE PRECISION cpx(rdim,adim,sdim),
     &     cpy(rdim,adim,sdim),
     &     cpz(rdim,adim,sdim)

      DOUBLE PRECISION gv(rdim,adim,wdim,zdim,kdim),
     &     rc0(rdim,adim,wdim,kdim),
     &     pcx(rdim,adim,wdim,zdim),
     &     pcy(rdim,adim,wdim,zdim),
     &     pcz(rdim,adim,wdim,zdim)

      DOUBLE PRECISION zet0,deltacg,umaginv,ang1,ang2,pi
      INTEGER nz,Nb

      INTEGER zif

      external vorp

      pi=4.0D0*atan(1.0D0)
      sn=1.D+0
      wage=DBLE(z1-1)*dz/DBLE(zif)
      if(z1.gt.nz) wage=DBLE(nz-1)*dz/DBLE(zif)

      r1x=cpx(r0,p0,s0)-pcx(r1,p1,w1,z1-1)
      r1y=cpy(r0,p0,s0)-pcy(r1,p1,w1,z1-1)
      r1z=cpz(r0,p0,s0)-pcz(r1,p1,w1,z1-1)
      mr1=sqrt(r1x**2+r1y**2+r1z**2)

      r2x=cpx(r0,p0,s0)-pcx(r1,p1,w1,z1)
      r2y=cpy(r0,p0,s0)-pcy(r1,p1,w1,z1)
      r2z=cpz(r0,p0,s0)-pcz(r1,p1,w1,z1)
      mr2=sqrt(r2x**2+r2y**2+r2z**2)

      lx=pcx(r1,p1,w1,z1)-pcx(r1,p1,w1,z1-1)
      ly=pcy(r1,p1,w1,z1)-pcy(r1,p1,w1,z1-1)
      lz=pcz(r1,p1,w1,z1)-pcz(r1,p1,w1,z1-1)
      ml=sqrt(lx**2+ly**2+lz**2)

      if ((mr1 .le. 1.0D-4) .or. (mr2 .le. 1.0D-4)) then
         ang1=0.0D+0
         ang2=0.0D+0
         h=0.0D+0
      else
         ang1=acos((lx*r1x+ly*r1y+lz*r1z)/(ml*mr1))
         ang2=acos((lx*r2x+ly*r2y+lz*r2z)/(ml*mr2))
         h=mr1*sin(ang1)
      endif
     
      if((mr1-3.0D+0*rc).lt.2.0D-15) h=h*00.D+0
      if((mr2-3.0D+0*rc).lt.2.0D-15) h=h*00.D+0

      ux=ly*r1z-lz*r1y
      uy=lz*r1x-lx*r1z
      uz=lx*r1y-ly*r1x
      umag=sqrt(ux**2+uy**2+uz**2)
      
      if ( (umag .eq. 0.D+0) .or. (h .eq. 0.D+0) ) then
         ux=0.D+0
         uy=0.D+0
         uz=0.D+0
         vbx=0.D+0
         vby=0.D+0
         vbz=0.D+0
      else
         ux=ux/umag
         uy=uy/umag
         uz=uz/umag
         vvp=gv(r1,p1,w1,(z1-1),1)*h/(4.0D+0*pi*((h**4+rc**4)**(0.5D+0))
     $        )
         icv=vvp*(cos(ang1)-cos(ang2))
         vbx=icv*ux
         vby=icv*uy
         vbz=icv*uz
      endif

c$$$      ux=ly*r1z-lz*r1y
c$$$      uy=lz*r1x-lx*r1z
c$$$      uz=lx*r1y-ly*r1x
c$$$      umag=sqrt(ux**2+uy**2+uz**2)
c$$$
c$$$      if(umag .gt. 0.0D0) then
c$$$         umaginv=1.0D0/umag
c$$$         h=umag/ml
c$$$         theta1=(lx*r1x+ly*r1y+lz*r1z)/(ml*mr1)
c$$$         theta2=(lx*r2x+ly*r2y+lz*r2z)/(ml*mr2)
c$$$
c$$$         vvp=vorp(h,rc)*gv(r1,p1,w1,(z1-1),1)*(theta1-theta2)*umaginv
c$$$         vbx=ux*vvp
c$$$         vby=uy*vvp
c$$$         vbz=uz*vvp
c$$$      else
c$$$         vbx=0.0D+0
c$$$         vby=0.0D+0
c$$$         vbz=0.0D+0
c$$$      endif
c$$$      if((s0 .eq. 15) .and. (p1 .eq. 1) .and. (p0.eq. 1)) then
c$$$         write(333,*) z1, vbz
c$$$      endif
c$$$                           
c$$$
c$$$      if(.false.)then
c$$$      if (umag .eq. 0.0) then
c$$$         ux=0.D+0
c$$$         uy=0.D+0
c$$$         uz=0.D+0
c$$$         ca=0.D+0
c$$$         cb=0.D+0
c$$$         vbx=0.D+0
c$$$         vby=0.D+0
c$$$         vbz=0.D+0
c$$$         h=0.D+0
c$$$      else
c$$$         if ( (mr1 .le. 0.1D-6) .or. (mr2 .le. 0.1D-6) ) then
c$$$            theta1=0.D+0
c$$$            theta2=0.D+0
c$$$            h=0.D+0
c$$$         else
c$$$            theta1=acos((lx*r1x+ly*r1y+lz*r1z)/(ml*mr1))
c$$$            theta2=acos((lx*r2x+ly*r2y+lz*r2z)/(ml*mr2))
c$$$            h=mr1*sin(theta1)
c$$$         endif
c$$$         ux=ux/umag
c$$$         uy=uy/umag
c$$$         uz=uz/umag
c$$$
c$$$         chi=int(abs(DBLE(np-(p1-z1+2))/DBLE(np)))
c$$$         p3=(1+chi)*np+p1-z1+2-np
c$$$
c$$$         chi2=int(abs(DBLE(np-(p1-z1+1))/DBLE(np)))
c$$$         p4=(1+chi2)*np+p1-z1+1-np
c$$$         cot1=1.D+0/tan(theta1)
c$$$         cot2=1.D+0/tan(theta2)
c$$$         do k=1,nk
c$$$            gam1=gv(r1,p1,w1,(z1-1),k)
c$$$            gam2=lin*gv(r1,p1,w1,(z1-0),k)+ (1.D+0-lin)*gv(r1,p1,w1,(z1
c$$$     $           -1),k)
c$$$
c$$$            gam1=sn*gam1
c$$$            gam2=sn*gam2
c$$$            ca=(gam1*cot2-gam2*cot1)/(cot2-cot1)
c$$$            cb=(gam2-gam1)/(cot2-cot1)
c$$$            vvp=vorp(h,rc)
c$$$            icv=vvp*ca*(cos(theta1)-cos(theta2))
c$$$            icv2=vvp*cb*(sin(theta2)-sin(theta1))
c$$$            icvx=icv*ux
c$$$            icvy=icv*uy
c$$$            icvz=icv*uz
c$$$            icvx2=icv2*ux
c$$$            icvy2=icv2*uy
c$$$            icvz2=icv2*uz
c$$$            vbx=icvx+icvx2
c$$$            vby=icvy+icvy2
c$$$            vbz=icvz+icvz2
c$$$         end do
c$$$      endif
c$$$      endif
c$$$
      return
      end
