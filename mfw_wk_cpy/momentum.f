c----------------------------------------------------------------------|
c     Copyright violation -MrB
c----------------------------------------------------------------------|
c-mb inflow velo induced by bound vortex alone
c-mb this is assumed to be the momentum added by the blade
c----------------------------------------------------------------------|
c
      subroutine momentum(r0,p0,s0,ns,r1,p1,s1,rcb,gb,rad,aza,
     &     t0r,t1c,t1s,cpx,cpy,cpz,bvx,bvy,bvz,vix,viy,viz)
c----------------------------------------------------------------------|
c
      IMPLICIT none
c
        include 'adimpar.inc'
        include 'kdimpar.inc'
        include 'rdimpar.inc'
        include 'sdimpar.inc'
        include 'wdimpar.inc'
        include 'zdimpar.inc'
c
      INTEGER r0,p0,s0,z0,r1,p1,s1,ns

c
      DOUBLE PRECISION rcb,vix,viy,viz,vorp,
     &       r1x,r1y,r1z,mr1,r2x,r2y,r2z,mr2,
     &       lx,ly,lz,ml,ang1,ang2,h,ux,uy,uz,
     &       umag,vvp,ic,icx,icy,icz
c
      DOUBLE PRECISION gb(rdim,adim,sdim),
     &       bvx(rdim,adim,sdim),
     &       bvy(rdim,adim,sdim),
     &       bvz(rdim,adim,sdim)
c
      DOUBLE PRECISION cpx(rdim,adim,sdim),
     &       cpy(rdim,adim,sdim),
     &       cpz(rdim,adim,sdim)
c
      DOUBLE PRECISION cx,cy,cz,rad,ratio,qx,qy,qz
c
      DOUBLE PRECISION vix1,viy1,viz1,pitch,aza
      DOUBLE PRECISION t0r(rdim),t1c(rdim),t1s(rdim)
c
      external vorp
c----------------------------------------------------------------------|
c-mb the blade control points are at segment  midpoints
c-mb I want to calculate inflow at the tip so extrapolate
c
c-mb   segment end-points corresponding to vortex trailer(s)
c-mb    at 3/4 chords though (just like the CP)
      if(s0.eq.ns) then
       cx=cpx(r0,p0,ns-1)+1.5D+0*(cpx(r0,p0,ns)-cpx(r0,p0,ns-1))
       cy=cpy(r0,p0,ns-1)+1.5D+0*(cpy(r0,p0,ns)-cpy(r0,p0,ns-1))
       cz=cpz(r0,p0,ns-1)+1.5D+0*(cpz(r0,p0,ns)-cpz(r0,p0,ns-1))
      elseif(s0.eq.0) then
       cx=cpx(r0,p0,2)+1.5D+0*(cpx(r0,p0,1)-cpx(r0,p0,2))
       cy=cpy(r0,p0,2)+1.5D+0*(cpy(r0,p0,1)-cpy(r0,p0,2))
       cz=cpz(r0,p0,2)+1.5D+0*(cpz(r0,p0,1)-cpz(r0,p0,2))
      else
       cx=(cpx(r0,p0,s0)+cpx(r0,p0,s0+1))/2.D+0
       cy=(cpy(r0,p0,s0)+cpy(r0,p0,s0+1))/2.D+0
       cz=(cpz(r0,p0,s0)+cpz(r0,p0,s0+1))/2.D+0
      endif

c-mb first get the quarter chord point using the bvx
c-mb note the +1 so we can have s0 becoming 0
c-mb    qx=(bvx(r0,p0,s0+1)+bvx(r0,p0,s0+1))/2.D+0
c-mb    qy=(bvy(r0,p0,s0+1)+bvy(r0,p0,s0+1))/2.D+0
c-mb    qz=(bvz(r0,p0,s0+1)+bvz(r0,p0,s0+1))/2.D+0
c
c-mb now interpolate towards TE
c-mb       cx=qx+1.0D+0*(cx-qx)
c-mb       cy=qy+1.0D+0*(cy-qy)
c-mb       cz=qz+1.0D+0*(cz-qz)
c
c-mb now we are calculating the induced velos at the TE
          

          r1x=cx-bvx(r1,p1,s1)
          r1y=cy-bvy(r1,p1,s1)
          r1z=cz-bvz(r1,p1,s1)
          mr1=sqrt(r1x**2+r1y**2+r1z**2)
c
          r2x=cx-bvx(r1,p1,s1+1)
          r2y=cy-bvy(r1,p1,s1+1)
          r2z=cz-bvz(r1,p1,s1+1)
          mr2=sqrt(r2x**2+r2y**2+r2z**2)
c
          lx=bvx(r1,p1,s1+1)-bvx(r1,p1,s1)
          ly=bvy(r1,p1,s1+1)-bvy(r1,p1,s1)
          lz=bvz(r1,p1,s1+1)-bvz(r1,p1,s1)
          ml=sqrt(lx**2+ly**2+lz**2)
c
          if ( (mr1 .le. 0.1D-4) .or. (mr2 .le. 0.1D-4) ) then
           ang1=0.D+0
           ang2=0.D+0
           h=0.D+0
          else
           ang1=acos((lx*r1x+ly*r1y+lz*r1z)/(ml*mr1))
           ang2=acos((lx*r2x+ly*r2y+lz*r2z)/(ml*mr2))
           h=mr1*sin(ang1)
          endif
c
          ux=ly*r1z-lz*r1y
          uy=lz*r1x-lx*r1z
          uz=lx*r1y-ly*r1x
          umag=sqrt(ux**2+uy**2+uz**2)
c
          if ( (umag .le. 0.1D-4) .or. (h .le. 0.1D-4) ) then
           ux=0.D+0
           uy=0.D+0
           uz=0.D+0
          else
           ux=ux/umag
           uy=uy/umag
           uz=uz/umag
c
           vvp=vorp(h,rcb*0.D+0)
c-mb                      Note core size set to zero for bound
          endif
c
          ic=vvp*(cos(ang1)-cos(ang2))
          icx=ic*ux
          icy=ic*uy
          icz=ic*uz
c
          vix=icx*gb(r1,p1,s1)
          viy=icy*gb(r1,p1,s1)
          viz=icz*gb(r1,p1,s1)
c----------------------------------------------------------------------|
c-mb  note that this velocity is nomal to blade surface!
c-mb  will need to include effects of both flap & pitch
c-mb  the effect of flap is already included because using bvx/y/z
c-mb  however, in calculating bcp, we dont include pitch!
      pitch=t0r(r0)+t1c(r0)*cos(aza)+t1s(r0)*sin(aza)
      vix1=vix
      viy1=viy
      viz1=viz
c
c-mb  this will slightly decrease the increased momentum (:-)
      viz=viz1*cos(pitch)
c
c>mb  just chek if I've already included pitch effects in bcp.f
c----------------------------------------------------------------------|
c
      return
      end
