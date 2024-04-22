C-$	$Id: sbind.f,v 1.1.1.1.6.1 2005/03/10 18:49:33 shreyas Exp $	

      subroutine sbind(r0,p0,s0,r1,p1,w1,z1,np,nk,dz,om,rc0, gv,rc,lin
     $     ,trb,sce,gsc,csc,cpx,cpy,cpz, pcx,pcy,pcz,vbx,vby,vbz,zif)

c     Shed-cirulation/blade induced velocity calculations.


      IMPLICIT none

      include 'adimpar.inc'
      include 'kdimpar.inc'
      include 'rdimpar.inc'
      include 'sdimpar.inc'
      include 'wdimpar.inc'
      include 'zdimpar.inc'

      INTEGER r0,p0,s0,r1,p1,w1,z1,np,nk,k,p3,p4,chi,chi2

      DOUBLE PRECISION dz,lin,trb,om,vbx,vby,vbz,wage, r1x,r1y,r1z,mr1
     $     ,r2x,r2y,r2z,mr2, lx,ly,lz,ml,theta1,theta2,h,ux,uy,uz, umag
     $     ,ca,cb,cot1,cot2,g1,g2,rc,vvp, icv,icv2,icvx,icvy,icvz,icvx2
     $     ,icvy2,icvz2, vorp,sce,sn

      DOUBLE PRECISION cpx(rdim,adim,sdim), cpy(rdim,adim,sdim),
     $     cpz(rdim,adim,sdim)

      DOUBLE PRECISION gv(rdim,adim,wdim,zdim,kdim), rc0(rdim,adim,wdim
     $     ,kdim), pcx(rdim,adim,wdim,zdim), pcy(rdim,adim,wdim,zdim),
     $     pcz(rdim,adim,wdim,zdim)

      DOUBLE PRECISION gsc(rdim,adim,wdim,zdim,kdim), csc(rdim,adim,wdim
     $     ,zdim,kdim)

      INTEGER zif

      external vorp


      sn=(-1.D+0)**(r1-1)
      wage=DBLE(z1-1)*dz/DBLE(zif)

      r1x=cpx(r0,p0,s0)-pcx(r1,p1,w1,z1)
      r1y=cpy(r0,p0,s0)-pcy(r1,p1,w1,z1)
      r1z=cpz(r0,p0,s0)-pcz(r1,p1,w1,z1)
      mr1=sqrt(r1x**2+r1y**2+r1z**2)

      r2x=cpx(r0,p0,s0)-pcx(r1,p1,w1-1,z1)
      r2y=cpy(r0,p0,s0)-pcy(r1,p1,w1-1,z1)
      r2z=cpz(r0,p0,s0)-pcz(r1,p1,w1-1,z1)
      mr2=sqrt(r2x**2+r2y**2+r2z**2)

      lx=pcx(r1,p1,w1-1,z1)-pcx(r1,p1,w1,z1)
      ly=pcy(r1,p1,w1-1,z1)-pcy(r1,p1,w1,z1)
      lz=pcz(r1,p1,w1-1,z1)-pcz(r1,p1,w1,z1)
      ml=sqrt(lx**2+ly**2+lz**2)

      ux=ly*r1z-lz*r1y
      uy=lz*r1x-lx*r1z
      uz=lx*r1y-ly*r1x
      umag=sqrt(ux**2+uy**2+uz**2)

      if (umag .eq. 0.0) then
         ux=0.D+0
         uy=0.D+0
         uz=0.D+0

         ca=0.D+0
         cb=0.D+0

         vbx=0.D+0
         vby=0.D+0
         vbz=0.D+0

         h=0.D+0

      else
         if ( (mr1 .le. 0.1D-6) .or. (mr2 .le. 0.1D-6) ) then
            theta1=0.D+0
            theta2=0.D+0
            h=0.D+0
         else
            theta1=acos((lx*r1x+ly*r1y+lz*r1z)/(ml*mr1))
            theta2=acos((lx*r2x+ly*r2y+lz*r2z)/(ml*mr2))
            h=mr1*sin(theta1)
         endif

         ux=ux/umag
         uy=uy/umag
         uz=uz/umag
         chi=int(abs(DBLE(np-(p1-z1+2))/DBLE(np)))
         p3=(1+chi)*np+p1-z1+2-np
         chi2=int(abs( DBLE(np-(p1-z1+1))/DBLE(np) ))
         p4=(1+chi2)*np+p1-z1+1-np

         cot1=1.D+0/tan(theta1)
         cot2=1.D+0/tan(theta2)

         do k=1,nk

            g1=gv(r1,p1,w1-1,z1-0,k)-gv(r1,p1,w1-1,z1-1,k)
            g2=gv(r1,p1,w1,z1-0,k)-gv(r1,p1,w1,z1-1,k)

            ca=lin*(g1*cot2-g2*cot1)/(cot2-cot1)+(1.D+0-lin)*(g1+g2)
            cb=lin*(g2-g1)/(cot2-cot1)

            ca=sce*gsc(r1,p1,w1-1,z1,k)+(1.D+0-sce)*ca
            cb=(1.D+0-sce)*cb

            rc=sce*csc(r1,p1,w1-1,z1,k)+(1.D+0-sce)*rc

            vvp=vorp(h,rc)

            icv=vvp*ca*(cos(theta1)-cos(theta2))
            icv2=vvp*cb*(sin(theta2)-sin(theta1))

            icvx=icv*ux
            icvy=icv*uy
            icvz=icv*uz

            icvx2=icv2*ux
            icvy2=icv2*uy
            icvz2=icv2*uz

            vbx=icvx+icvx2
            vby=icvy+icvy2
            vbz=icvz+icvz2

         end do

      endif

      return
      end
