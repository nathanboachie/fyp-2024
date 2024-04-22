C-$	$Id: bvind.f,v 1.1.1.1.6.2 2005/06/11 15:57:36 shreyas Exp $	

      subroutine bvind(r0,p0,w0,z0,r1,p1,s1,rcb,gb, px,py,pz,bvx,bvy,bvz
     $     ,vix,viy,viz)

c     Bound-vortex/free-vortex induced velocity calculations.


      IMPLICIT none

      include 'adimpar.inc'
      include 'kdimpar.inc'
      include 'rdimpar.inc'
      include 'sdimpar.inc'
      include 'wdimpar.inc'
      include 'zdimpar.inc'

      INTEGER r0,p0,w0,z0,r1,p1,s1


      DOUBLE PRECISION rcb,vix,viy,viz,vorp, r1x,r1y,r1z,mr1,r2x,r2y,r2z
     $     ,mr2, lx,ly,lz,ml,ang1,ang2,h,ux,uy,uz, umag,vvp,ic,icx,icy
     $     ,icz,umaginv

      DOUBLE PRECISION gb(rdim,adim,sdim), bvx(rdim,adim,sdim), bvy(rdim
     $     ,adim,sdim), bvz(rdim,adim,sdim)

      DOUBLE PRECISION px(rdim,adim,wdim,zdim), py(rdim,adim,wdim,zdim),
     $     pz(rdim,adim,wdim,zdim)

      external vorp

      common/velcomp/r1x,r1y,r1z,mr1,r2x,r2y,r2z,mr2
      save/velcomp/

      if (s1 .gt. 1) then
         r1x=r2x
         r1y=r2y
         r1z=r2z
         mr1=mr2
      else
         r1x=px(r0,p0,w0,z0)-bvx(r1,p1,s1)
         r1y=py(r0,p0,w0,z0)-bvy(r1,p1,s1)
         r1z=pz(r0,p0,w0,z0)-bvz(r1,p1,s1)
         mr1=sqrt(r1x*r1x+r1y*r1y+r1z*r1z)
      endif

      r2x=px(r0,p0,w0,z0)-bvx(r1,p1,s1+1)
      r2y=py(r0,p0,w0,z0)-bvy(r1,p1,s1+1)
      r2z=pz(r0,p0,w0,z0)-bvz(r1,p1,s1+1)
      mr2=sqrt(r2x*r2x+r2y*r2y+r2z*r2z)

      lx=bvx(r1,p1,s1+1)-bvx(r1,p1,s1)
      ly=bvy(r1,p1,s1+1)-bvy(r1,p1,s1)
      lz=bvz(r1,p1,s1+1)-bvz(r1,p1,s1)
      ml=sqrt(lx*lx+ly*ly+lz*lz)

      ux=ly*r1z-lz*r1y
      uy=lz*r1x-lx*r1z
      uz=lx*r1y-ly*r1x
      umag=sqrt(ux*ux+uy*uy+uz*uz)
      h=umag/ml


      if (umag .gt. 1.0D-8) then 
         umaginv=1.0D0/umag
         ang1=(lx*r1x+ly*r1y+lz*r1z)/(ml*mr1)
         ang2=(lx*r2x+ly*r2y+lz*r2z)/(ml*mr2)

         vvp=vorp(h,rcb)*gb(r1,p1,s1)*(ang1-ang2)*umaginv
         vix=ux*vvp
         viy=uy*vvp
         viz=uz*vvp
      else
         vix=0.0D0
         viy=0.0D0
         viz=0.0D0
      endif

      if(.false.)then
      if ( (mr1 .le. 0.1D-4) .or. (mr2 .le. 0.1D-4) ) then
         ang1=0.D+0
         ang2=0.D+0
         h=0.D+0
      else
         ang1=acos((lx*r1x+ly*r1y+lz*r1z)/(ml*mr1))
         ang2=acos((lx*r2x+ly*r2y+lz*r2z)/(ml*mr2))
         h=mr1*sin(ang1)
      endif

      ux=ly*r1z-lz*r1y
      uy=lz*r1x-lx*r1z
      uz=lx*r1y-ly*r1x
      umag=sqrt(ux**2+uy**2+uz**2)

      if ( (umag .le. 0.1D-6) .or. (h .le. 0.1D-6) ) then
         ux=0.D+0
         uy=0.D+0
         uz=0.D+0
      else
         ux=ux/umag
         uy=uy/umag
         uz=uz/umag

         vvp=vorp(h,rcb)
      endif

      ic=vvp*(cos(ang1)-cos(ang2))
      icx=ic*ux
      icy=ic*uy
      icz=ic*uz

      vix=icx*gb(r1,p1,s1)
      viy=icy*gb(r1,p1,s1)
      viz=icz*gb(r1,p1,s1)

      endif

      return
      end
