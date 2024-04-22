C-$	$Id: corrab.f,v 1.1.1.1.6.2 2005/03/23 21:24:38 shreyas Exp $	
c----------------------------------------------------------------------|
c    <Copyright (c) 1997 by Ashish Bagai>
c    <Copyright (c) 2000 by Mahendra Bhagwat> 
c----------------------------------------------------------------------|

      subroutine corrab(nr,np,nw,pw,nz,nzt,nk,dp,dz,pif,zif,rad,om,rbc,
     &             zbc,avgzbc,lami,gv,rr_x,rr_y,rr_z,vinfx,vinfy,vinfz,
     &             rpar,vipx,vipy,vipz,ppx,ppy,ppz,vicx,vicy,
     &             vicz,pcx,pcy,pcz)

c     Pseudo-IMPLICIT corrector -- free-vortex position vectors.
c----------------------------------------------------------------------|

      IMPLICIT none
c----------------------------------------------------------------------|
c
        include 'adimpar.inc'
        include 'kdimpar.inc'
        include 'rdimpar.inc'
        include 'wdimpar.inc'
        include 'zdimpar.inc'
c----------------------------------------------------------------------|
c
      INTEGER nr,np,nw,pw,nz,nzt,pif,zif,r,w,z,p,p1,p2,z1,k,nk,p3,chi
c
      DOUBLE PRECISION dp,dz,rad,om,rbc(rdim,wdim),zbc(rdim,wdim),
     &       avgzbc(rdim,wdim),
     &       vinfx,vinfy,vinfz,dp1,dz1,
     &       vx,vy,vz,velox,veloy,veloz,rr,zr,rpar,gt
c
      DOUBLE PRECISION lami(rdim),rr_x(rdim),rr_y(rdim),rr_z(rdim)
c
      DOUBLE PRECISION ppx(rdim,adim,wdim,zdim),
     &       ppy(rdim,adim,wdim,zdim),
     &       ppz(rdim,adim,wdim,zdim),
     &       pcx(rdim,adim,wdim,zdim),
     &       pcy(rdim,adim,wdim,zdim),
     &       pcz(rdim,adim,wdim,zdim),
     &       vipx(rdim,adim,wdim,zdim),
     &       vipy(rdim,adim,wdim,zdim),
     &       vipz(rdim,adim,wdim,zdim),
     &       vicx(rdim,adim,wdim,zdim),
     &       vicy(rdim,adim,wdim,zdim),
     &       vicz(rdim,adim,wdim,zdim),
     &       gv(rdim,adim,wdim,zdim,kdim)
c
      external velox
      external veloy
      external veloz
c----------------------------------------------------------------------|
c     Pseudo-Implicit Corrector using Implicit Boundary Conditions:
c     Velocity relaxation -- predictor/corrector:
c
      do r=1,nr
      do p=1,np
      do w=1,nw    ! solve for inboard sheet
c-mb-choice      do w=1,nw-pw ! initial prescribed inboard sheet
      do z=1,nz

         vicx(r,p,w,z)=rpar*vicx(r,p,w,z)+(1.D+0-rpar)*vipx(r,p,w,z)
         vicy(r,p,w,z)=rpar*vicy(r,p,w,z)+(1.D+0-rpar)*vipy(r,p,w,z)
         vicz(r,p,w,z)=rpar*vicz(r,p,w,z)+(1.D+0-rpar)*vipz(r,p,w,z)

      end do
      end do
      end do
      end do
c----------------------------------------------------------------------|
c
      dp1=dp/DBLE(pif)
      dz1=dz/DBLE(zif)
c----------------------------------------------------------------------|
c
      do r=1,nr
      do w=1,nw    ! solve for inboard sheet
c-mb-choice      do w=1,nw-pw ! initial prescribed inboard sheet
      do z=2,nz
      do p=1,np

         p1=p-1
         p2=p-2
         p3=p-3
         if (p1 .le. 0) p1=p1+np
         if (p2 .le. 0) p2=p2+np
         if (p3 .le. 0) p3=p3+np

c
c        Velocity averaging and concatonate free-stream velocities 
c        to induced velocities at all free and ersatz-free points:
c
         vx=0.25D+0*(vicx(r,p1,w,z-1)+vicx(r,p1,w,z)+
     &              vicx(r,p,w,z-1) +vicx(r,p,w,z) )
         vx=vx+vinfx+velox(ppx(r,p,w,z),ppy(r,p,w,z),ppz(r,p,w,z),z)

         vy=0.25D+0*(vicy(r,p1,w,z-1)+vicy(r,p1,w,z)+
     &              vicy(r,p,w,z-1) +vicy(r,p,w,z) )
         vy=vy+vinfy+veloy(ppx(r,p,w,z),ppy(r,p,w,z),ppz(r,p,w,z),z)

         vz=0.25D+0*(vicz(r,p1,w,z-1)+vicz(r,p1,w,z)+
     &              vicz(r,p,w,z-1) +vicz(r,p,w,z) )
         vz=vz+vinfz+veloz(ppx(r,p,w,z),ppy(r,p,w,z),ppz(r,p,w,z),z)

c        Corrector:

         pcx(r,p,w,z)=pcx(r,p1,w,z-1)+
     &                (2.D+0/om)*(dp1*dz1/(dp1+dz1))*vx
     &   +0.d0/16.d0*( pcx(r,p,w,z)-3.d0*pcx(r,p1,w,z)
     &               +3.d0*pcx(r,p2,w,z)-pcx(r,p3,w,z)
     &               +pcx(r,p,w,z-1)-3.d0*pcx(r,p1,w,z-1)
     &               +3.d0*pcx(r,p2,w,z-1)-pcx(r,p3,w,z-1) )

         pcy(r,p,w,z)=pcy(r,p1,w,z-1)+
     &                (2.D+0/om)*(dp1*dz1/(dp1+dz1))*vy
     &   +0.d0/16.d0*( pcy(r,p,w,z)-3.d0*pcy(r,p1,w,z)
     &               +3.d0*pcy(r,p2,w,z)-pcy(r,p3,w,z)
     &               +pcy(r,p,w,z-1)-3.d0*pcy(r,p1,w,z-1)
     &               +3.d0*pcy(r,p2,w,z-1)-pcy(r,p3,w,z-1) )

         pcz(r,p,w,z)=pcz(r,p1,w,z-1)+
     &                (2.D+0/om)*(dp1*dz1/(dp1+dz1))*vz
     &   +0.d0/16.d0*( pcz(r,p,w,z)-3.d0*pcz(r,p1,w,z)
     &               +3.d0*pcz(r,p2,w,z)-pcz(r,p3,w,z)
     &               +pcz(r,p,w,z-1)-3.d0*pcz(r,p1,w,z-1)
     &               +3.d0*pcz(r,p2,w,z-1)-pcz(r,p3,w,z-1) )

      end do
      end do
      end do
      end do
c----------------------------------------------------------------------|
c        Hover boundary condition wake geometry:
c
c-mb    this is for the NEW farwake
         p=0

         call farwake(rbc,zbc,avgzbc,rad,nr,nw,pw,np,nz,nzt,
     &    pcx,pcy,pcz,
     &    ppx,ppy,ppz,rr_x,rr_y,rr_z,dp,dz,pif,zif,lami,p)
c-mb
c----------------------------------------------------------------------|
c
      return
      end
