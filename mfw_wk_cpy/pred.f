C-$	$Id: pred.f,v 1.1.1.1.6.1 2005/03/10 18:49:33 shreyas Exp $	

      subroutine pred(nr,np,nw,pw,nz,nzt,nk,dp,dz,pif,zif,rad,om,rbc,
     $     zbc,avgzbc,lami,gv,rr_x,rr_y,rr_z,vinfx,vinfy,vinfz, pcx,pcy
     $     ,pcz,vipx,vipy,vipz,ppx,ppy,ppz, p)

C-$   Predictor step for wake control points (PC2B scheme)


      IMPLICIT none

      include 'adimpar.inc'
      include 'kdimpar.inc'
      include 'rdimpar.inc'
      include 'wdimpar.inc'
      include 'zdimpar.inc'

      INTEGER nr,np,nw,pw,nz,nzt,pif,zif,r,w,z,p,p1,p2,z1,k,nk,p3,chi

      DOUBLE PRECISION dp,dz,rad,om,rbc(rdim,wdim),zbc(rdim,wdim),
     $     avgzbc(rdim,wdim), vinfx,vinfy,vinfz,dp1,dz1, vx,vy,vz,velox
     $     ,veloy,veloz,rr,zr,gt

      DOUBLE PRECISION lami(rdim),rr_x(rdim),rr_y(rdim),rr_z(rdim)
      
      DOUBLE PRECISION pcx(rdim,adim,wdim,zdim), pcy(rdim,adim,wdim,zdim
     $     ), pcz(rdim,adim,wdim,zdim), ppx(rdim,adim,wdim,zdim),
     $     ppy(rdim,adim,wdim,zdim), ppz(rdim,adim,wdim,zdim), vipx(rdim
     $     ,adim,wdim,zdim), vipy(rdim,adim,wdim,zdim), vipz(rdim,adim
     $     ,wdim,zdim), gv(rdim,adim,wdim,zdim,kdim)

      external velox
      external veloy
      external veloz


      dp1=dp/DBLE(pif)
      dz1=dz/DBLE(zif)

      do r=1,nr
         do w=1,nw              ! solve for inboard sheet
            do z=2,nz

               p1=p-1
               p2=p-1
               if (p1 .le. 0) p1=p1+np
               if (p2 .le. 0) p2=p2+np
               
               vx=vipx(r,p2,w,z-1)
               vx=vx+vinfx+velox(pcx(r,p,w,z),pcy(r,p,w,z),pcz(r,p,w,z)
     $              ,z)

               vy=vipy(r,p2,w,z-1)
               vy=vy+vinfy+veloy(pcx(r,p,w,z),pcy(r,p,w,z),pcz(r,p,w,z)
     $              ,z)
               
               vz=vipz(r,p2,w,z-1)
               vz=vz+vinfz+veloz(pcx(r,p,w,z),pcy(r,p,w,z),pcz(r,p,w,z)
     $              ,z)

c     Predictor:
               ppx(r,p,w,z)=pcx(r,p1,w,z-1)+ (2.D+0/om)*(dp1*dz1/(dp1
     $              +dz1))*vx

               ppy(r,p,w,z)=pcy(r,p1,w,z-1)+ (2.D+0/om)*(dp1*dz1/(dp1
     $              +dz1))*vy

               ppz(r,p,w,z)=pcz(r,p1,w,z-1)+ (2.D+0/om)*(dp1*dz1/(dp1
     $              +dz1))*vz
            end do
         end do
      end do

c        Hover boundary condition wake geometry:

      call farwake(rbc,zbc,avgzbc,rad,nr,nw,pw,np,nz,nzt, ppx,ppy,ppz,
     $     pcx,pcy,pcz,rr_x,rr_y,rr_z,dp,dz,pif,zif,lami,p)


      return
      end
