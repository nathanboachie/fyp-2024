C-$	$Id: svind.f,v 1.1.1.1.6.1 2005/03/10 18:49:33 shreyas Exp $	

      subroutine svind(r0,p0,w0,z0,r1,p1,w1,z1,np,nk, pif,zif,dz,om,anb
     $     ,gv,rc,rc0,lin, trb,sce,gsc,csc,px,py,pz,vix,viy,viz)

c     Shed-circulation/free-vortex induced velocity calculations.


      IMPLICIT none

      include 'adimpar.inc'
      include 'kdimpar.inc'
      include 'rdimpar.inc'
      include 'wdimpar.inc'
      include 'zdimpar.inc'

      INTEGER r0,p0,w0,z0,r1,p1,w1,z1,np, nk,pif,zif,chi1,chi2,p3,p4,k

      DOUBLE PRECISION dz,om,anb,lin,trb,vix,viy,viz,wage, r1x,r1y,r1z
     $     ,mr1,r2x,r2y,r2z,mr2, lx,ly,lz,ml,ang1,ang2,h,ux,uy,uz, umag
     $     ,ca,cb,cot1,cot2,g1,g2,rc,vvp, ic,ic2,icx,icy,icz,icx2,icy2
     $     ,icz2, vorp,ag1,ag2,sce

      DOUBLE PRECISION gv(rdim,adim,wdim,zdim,kdim), rc0(rdim,adim,wdim
     $     ,kdim), px(rdim,adim,wdim,zdim), py(rdim,adim,wdim,zdim),
     $     pz(rdim,adim,wdim,zdim)

      DOUBLE PRECISION gsc(rdim,adim,wdim,zdim,kdim), csc(rdim,adim,wdim
     $     ,zdim,kdim)

      external vorp


      wage=DBLE(z1-1)*dz/DBLE(zif)

      r1x=px(r0,p0,w0,z0)-px(r1,p1,w1,z1)
      r1y=py(r0,p0,w0,z0)-py(r1,p1,w1,z1)
      r1z=pz(r0,p0,w0,z0)-pz(r1,p1,w1,z1)
      mr1=sqrt(r1x**2+r1y**2+r1z**2)

      r2x=px(r0,p0,w0,z0)-px(r1,p1,w1-1,z1)
      r2y=py(r0,p0,w0,z0)-py(r1,p1,w1-1,z1)
      r2z=pz(r0,p0,w0,z0)-pz(r1,p1,w1-1,z1)
      mr2=sqrt(r2x**2+r2y**2+r2z**2)

      lx=px(r1,p1,w1-1,z1)-px(r1,p1,w1,z1)
      ly=py(r1,p1,w1-1,z1)-py(r1,p1,w1,z1)
      lz=pz(r1,p1,w1-1,z1)-pz(r1,p1,w1,z1)
      ml=sqrt(lx**2+ly**2+lz**2)

      if ( (mr1 .le. 0.1D-6)  .or.
     &     (mr2 .le. 0.1D-6)  .or.
     &     (ml  .le. 0.1D-6) ) then
         ag1=0.D+0
         ag2=0.D+0
         h=0.D+0
      else
         ag1=(lx*r1x+ly*r1y+lz*r1z)/(ml*mr1)
         ag2=(lx*r2x+ly*r2y+lz*r2z)/(ml*mr2)
      endif

      if ( (abs(ag1) .ge. 1.D+0)   .or. 
     &     (abs(ag2) .ge. 1.D+0) )  then
         ag1=ag1-1.D+0*(ag1/abs(ag1))*1.0D-7
         ag2=ag2-1.D+0*(ag2/abs(ag2))*1.0D-8
      endif
      
      ang1=acos(ag1)
      ang2=acos(ag2)
      h=mr1*sin(ang1)

      if (ang1 .eq. ang2) h=0.D+0

      ux=ly*r1z-lz*r1y
      uy=lz*r1x-lx*r1z
      uz=lx*r1y-ly*r1x
      umag=sqrt(ux**2+uy**2+uz**2)

      if ( (umag .eq. 0.D+0) .or. (h .eq. 0.D+0) ) then

         ux=0.D+0
         uy=0.D+0
         uz=0.D+0
         ca=0.D+0
         cb=0.D+0
         vix=0.D+0
         viy=0.D+0
         viz=0.D+0

      else

         ux=ux/umag
         uy=uy/umag
         uz=uz/umag

         chi1=int(abs(DBLE(np-(p1-z1+2))/DBLE(np)))
         p3=(1+chi1)*np+p1-z1+2-np

         chi2=int(abs(DBLE(np-(p1-z1+1))/DBLE(np)))
         p4=(1+chi2)*np+p1-z1+1-np

         cot1=1.D+0/tan(ang1)
         cot2=1.D+0/tan(ang2)

         do k=1,nk

            g1=gv(r1,p1,w1-1,z1-0,k)-gv(r1,p1,w1-1,z1-1,k)
            g2=gv(r1,p1,w1,z1-0,k)-gv(r1,p1,w1,z1-1,k)

            ca=lin*(g1*cot2-g2*cot1)/(cot2-cot1) +(1.D+0-lin)*(g1-g2)/2
     $           .D+0
            cb=lin*(g2-g1)/(cot2-cot1)

            ca=sce*gsc(r1,p1,w1-1,z1,k)+(1.D+0-sce)*ca
            cb=(1.D+0-sce)*cb

            rc=sce*csc(r1,p1,w1-1,z1,k)+(1.D+0-sce)*rc

            vvp=vorp(h,rc)

            ic=vvp*ca*(cos(ang1)-cos(ang2))
            ic2=vvp*cb*(sin(ang2)-sin(ang1))

            icx=ic*ux
            icy=ic*uy
            icz=ic*uz

            icx2=ic2*ux
            icy2=ic2*uy
            icz2=ic2*uz

            vix=icx+icx2
            viy=icy+icy2
            viz=icz+icz2

         end do

      endif


      return
      end
