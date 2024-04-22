
      subroutine nwind(r0,p0,w0,z0,r1,p1,s1,rcb,gb,px,py,pz,bvx,bvy,bvz
     $     ,nwx,nwy,nwz,vix,viy,viz)

      IMPLICIT none

      include 'adimpar.inc'
      include 'kdimpar.inc'
      include 'rdimpar.inc'
      include 'sdimpar.inc'
      include 'wdimpar.inc'
      include 'zdimpar.inc'

      integer r0,p0,w0,z0,r1,p1,s1

      double precision rcb, pi, vorp

      double precision gb(rdim,adim,sdim), bvx(rdim,adim,sdim), bvy(rdim
     $     ,adim,sdim), bvz(rdim,adim,sdim)

      double precision px(rdim,adim,wdim,zdim), py(rdim,adim,wdim,zdim),
     $     pz(rdim,adim,wdim,zdim)
      
      double precision nwx(rdim,adim,sdim), nwy(rdim,adim,sdim),
     $     nwz(rdim,adim,sdim)


      double precision r1x, r1y, r1z, mr1, r2x, r2y, r2z, mr2, r3x,r3y
     $     ,r3z, mr3, r4x,r4y,r4z,mr4, l1x,l1y,l1z, m1l, l2x,l2y,l2z,m2l

      double precision cosa1, cosa2, h, ux,uy,uz, umag,umaginv, vvp

      double precision i1x, i1y, i1z, i2x, i2y, i2z,vix,viy,viz

      external vorp

      common/nwvelo/r1x,r1y,r1z,mr1,r2x,r2y,r2z,mr2,r3x,r3y ,r3z, mr3,
     $     r4x,r4y,r4z,mr4
      save/nwvelo/

      r1x=px(r0,p0,w0,z0)-nwx(r1,p1,s1)
      r1y=py(r0,p0,w0,z0)-nwy(r1,p1,s1)
      r1z=pz(r0,p0,w0,z0)-nwz(r1,p1,s1)
      mr1=sqrt(r1x*r1x+r1y*r1y+r1z*r1z)
         
      r2x=px(r0,p0,w0,z0)-bvx(r1,p1,s1)
      r2y=py(r0,p0,w0,z0)-bvy(r1,p1,s1)
      r2z=pz(r0,p0,w0,z0)-bvz(r1,p1,s1)
      mr2=sqrt(r2x*r2x+r2y*r2y+r2z*r2z)

      l1x=-nwx(r1,p1,s1)+bvx(r1,p1,s1)
      l1y=-nwy(r1,p1,s1)+bvy(r1,p1,s1)
      l1z=-nwz(r1,p1,s1)+bvz(r1,p1,s1)
      m1l=sqrt(l1x*l1x+l1y*l1y+l1z*l1z)

      if ( (mr1 .le. 1.0D-5) .or. (mr2 .le. 1.0D-5)) then
         i1x=0.0D0
         i1y=0.0D0
         i1z=0.0D0
      else
         cosa1=(l1x*r1x+l1y*r1y+l1z*r1z)/(m1l*mr1)
         cosa2=(l1x*r2x+l1y*r2y+l1z*r2z)/(m1l*mr2)

         ux=l1y*r1z-l1z*r1y
         uy=l1z*r1x-l1x*r1z
         uz=l1x*r1y-l1y*r1x
         umag=sqrt(ux*ux+uy*uy+uz*uz)
         umaginv=1.0D0/umag
         h=umag/m1l

         ux=ux*umaginv
         uy=uy*umaginv
         uz=uz*umaginv

         vvp=vorp(h,rcb)*(cosa1-cosa2)
         i1x=vvp*ux
         i1y=vvp*uy
         i1z=vvp*uz
      endif

      r3x=px(r0,p0,w0,z0)-bvx(r1,p1,s1+1)
      r3y=py(r0,p0,w0,z0)-bvy(r1,p1,s1+1)
      r3z=pz(r0,p0,w0,z0)-bvz(r1,p1,s1+1)
      mr3=sqrt(r3x*r3x+r3y*r3y+r3z*r3z)
      
      r4x=px(r0,p0,w0,z0)-nwx(r1,p1,s1+1)
      r4y=py(r0,p0,w0,z0)-nwy(r1,p1,s1+1)
      r4z=pz(r0,p0,w0,z0)-nwz(r1,p1,s1+1)
      mr4=sqrt(r4x*r4x+r4y*r4y+r4z*r4z)

      l2x=nwx(r1,p1,s1+1)-bvx(r1,p1,s1+1)
      l2y=nwy(r1,p1,s1+1)-bvy(r1,p1,s1+1)
      l2z=nwz(r1,p1,s1+1)-bvz(r1,p1,s1+1)
      m2l=sqrt(l2x*l2x+l2y*l2y+l2z*l2z)

      if ( (mr1 .le. 1.0D-5) .or. (mr2 .le. 1.0D-5)) then
         i2x=0.0D0
         i2y=0.0D0
         i2z=0.0D0
      else
         cosa1=(l2x*r3x+l2y*r3y+l2z*r3z)/(m2l*mr3)
         cosa2=(l2x*r4x+l2y*r4y+l2z*r4z)/(m2l*mr4)

         ux=l2y*r3z-l2z*r3y
         uy=l2z*r3x-l2x*r3z
         uz=l2x*r3y-l2y*r3x
         umag=sqrt(ux*ux+uy*uy+uz*uz)
         umaginv=1.0D0/umag
         h=umag/m2l

         ux=ux*umaginv
         uy=uy*umaginv
         uz=uz*umaginv

         vvp=vorp(h,rcb)*(cosa1-cosa2)
         i2x=vvp*ux
         i2y=vvp*uy
         i2z=vvp*uz
      endif
      
      vix=(i1x+i2x)*gb(r1,p1,s1)
      viy=(i1y+i2y)*gb(r1,p1,s1)
      viz=(i1z+i2z)*gb(r1,p1,s1)

c$$$      if (z0 .lt. 3) then 
c$$$         write(*,*) p1, s1, i1z,i2z
c$$$      endif

      return
      end
