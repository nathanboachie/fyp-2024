
      subroutine corr(nr,np,nw,pw,nz,nzt,nk,dp,dz,pif,zif,rad,om,rbc,
     &             zbc,avgzbc,lami,gv,rr_x,rr_y,rr_z,vinfx,vinfy,vinfz,
     &                rpar,vipx,vipy,vipz,ppx,ppy,ppz,vicx,vicy,
     &                vicz,pcx,pcy,pcz,p,zen)
c-mb     explicit corrector we know the pisitions at previous Psi
c-mb     and predicted at the next axim Psi+Delta Psi i.e
c-mb     time march by Delta t :)
c----------------------------------------------------------------------
c
      IMPLICIT none
c
        include 'adimpar.inc'
        include 'kdimpar.inc'
        include 'rdimpar.inc'
        include 'wdimpar.inc'
        include 'zdimpar.inc'
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
      INTEGER zen
c
      DOUBLE PRECISION pcxx,pcyy,pczz
c
      external velox
      external veloy
      external veloz


c----------------------------------------------------------------------
c
c     Explicit Corrector
c
      dp1=dp/DBLE(pif)
      dz1=dz/DBLE(zif)

      do r=1,nr
      do w=1,nw-pw
      do z=2,nz

         p1=p-1
         p2=p-1
         if (p1 .le. 0) p1=p1+np
         if (p2 .le. 0) p2=p2+np

c        Velocity averaging and concatonate free-stream velocities 

c-mb note velocities at p-1 changed to vipx
         vx=0.25D+0*(vipx(r,p2,w,z-1)+vipx(r,p2,w,z)+
     &              vicx(r,p,w,z-1) +vicx(r,p,w,z) )
         vx=vx+vinfx+velox(ppx(r,p,w,z),ppy(r,p,w,z),ppz(r,p,w,z),z)

         vy=0.25D+0*(vipy(r,p2,w,z-1)+vipy(r,p2,w,z)+
     &              vicy(r,p,w,z-1) +vicy(r,p,w,z) )
         vy=vy+vinfy+veloy(ppx(r,p,w,z),ppy(r,p,w,z),ppz(r,p,w,z),z)

         vz=0.25D+0*(vipz(r,p2,w,z-1)+vipz(r,p2,w,z)+
     &              vicz(r,p,w,z-1) +vicz(r,p,w,z) )
         vz=vz+vinfz+veloz(ppx(r,p,w,z),ppy(r,p,w,z),ppz(r,p,w,z),z)

c        Corrector:

         pcxx=pcx(r,p,w,z)
         pcyy=pcy(r,p,w,z)
         pczz=pcz(r,p,w,z)

         pcx(r,p,w,z)=pcx(r,p1,w,z-1)+
     &                (2.D+0/om)*(dp1*dz1/(dp1+dz1))*vx

         pcy(r,p,w,z)=pcy(r,p1,w,z-1)+
     &                (2.D+0/om)*(dp1*dz1/(dp1+dz1))*vy

         pcz(r,p,w,z)=pcz(r,p1,w,z-1)+
     &                (2.D+0/om)*(dp1*dz1/(dp1+dz1))*vz

         if(zen.eq.0) then
         pcx(r,p,w,z)=0.5D+0*(pcx(r,p,w,z)+pcxx)
         pcy(r,p,w,z)=0.5D+0*(pcy(r,p,w,z)+pcyy)
         pcz(r,p,w,z)=0.5D+0*(pcz(r,p,w,z)+pczz)
         endif

      end do
      end do
      end do
c
c        Hover boundary condition wake geometry:

         call farwake(rbc,zbc,avgzbc,rad,nr,nw,pw,np,nz,nzt,
     &      pcx,pcy,pcz,
     &      ppx,ppy,ppz,rr_x,rr_y,rr_z,dp,dz,pif,zif,lami,p)

c----------------------------------------------------------------------
c

      return
      end
