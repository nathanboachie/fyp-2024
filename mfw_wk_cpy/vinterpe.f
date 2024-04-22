c----------------------------------------------------------------------|
c     Copyright violation -MrB small changes
c-mb  modification of vinterp but only do one Psi at a time
c----------------------------------------------------------------------|
c
      subroutine vinterpe(nr,np,nw,pw,nz,pif,zif,vx,vy,vz,Pin)
c     Velocity interpolation for ersatz-free
c     collocation point convection velocities.
c----------------------------------------------------------------------|
c
      IMPLICIT none
c
        include 'adimpar.inc'
        include 'rdimpar.inc'
        include 'wdimpar.inc'
        include 'zdimpar.inc'
c
      INTEGER nr,np,nw,pw,nz,pif,zif,pf,zf,pf1,
     &        zin,r,p,w,z,p1,z1,p2,z2,p3
c
      DOUBLE PRECISION vx(rdim,adim,wdim,zdim),
     &       vy(rdim,adim,wdim,zdim),
     &       vz(rdim,adim,wdim,zdim)
c
      INTEGER Pin
c----------------------------------------------------------------------|
c
      if(pif.gt.1) then
c     For boundary release points -- z=1:
      z1=1
      z2=1
      z=1
         do r=1,nr
         do w=1,nw-pw
c-mb this interpolates the points between p1 and p1+pif
c-mb  i.e. the points p1+1,...,p1+pif-1 insteps of 1
         p1=Pin-pif
         if(p1.le.0) p1=p1+np
c-mb
         p2=p1+pif
         if (p2 .gt. np) p2=p2-np
         pf1=0
         do p=p1+1,p1+pif-1,1
         pf1=pf1+1

         vx(r,p,w,z)=vx(r,p1,w,z1)+
     &   DBLE(pf1)*(vx(r,p2,w,z2)-vx(r,p1,w,z1))/DBLE(pif)
         vy(r,p,w,z)=vy(r,p1,w,z1)+
     &   DBLE(pf1)*(vy(r,p2,w,z2)-vy(r,p1,w,z1))/DBLE(pif)
         vz(r,p,w,z)=vz(r,p1,w,z1)+
     &   DBLE(pf1)*(vz(r,p2,w,z2)-vz(r,p1,w,z1))/DBLE(pif)

         end do
         end do
         end do
      endif
c
c --- Interpolation on zeta -- Dz .gt. Dp:
c-mb at least for the explicit time marching, we interpolate
c-mb  along zeta at all points at once so no Zin!

      if (zif .gt. 1) then

      do r=1,nr
      do w=1,nw-pw
      do p=1,np,pif
      do z=1,nz-zif,zif
      zf=1
      do zin=z+1,z+zif-1,1

         vx(r,p,w,zin)=vx(r,p,w,z)-
     &   DBLE(zf)*(vx(r,p,w,z)-vx(r,p,w,z+zif))/DBLE(zif)
         vy(r,p,w,zin)=vy(r,p,w,z)-
     &   DBLE(zf)*(vy(r,p,w,z)-vy(r,p,w,z+zif))/DBLE(zif)
         vz(r,p,w,zin)=vz(r,p,w,z)-
     &   DBLE(zf)*(vz(r,p,w,z)-vz(r,p,w,z+zif))/DBLE(zif)
         zf=zf+1

      end do
      end do
      end do
      end do
      end do

      endif
c
c     Interpolation on psi -- Dp .gt. Dz:
      if (pif .gt. 1) then

c     Standard points -
      do r=1,nr
      do w=1,nw-pw
c-mb  this interpolates the points between (p3-1) and (p3+pif-1)
c-mb  remember p3 started with 2!!! that's why the -1
c-mb  otherwise essentially the same as above
c-mb  i.e. Pin=p3+pif-1
c-mb   -->
      p3=Pin-pif+1
      if(p3.le.0) p3=p3+np

      pf=0
      zf=0
      do p=p3,(p3-1)+(pif-1),1
      pf=pf+1
      zf=zf+1
      do z=1+pf,nz-(pif-pf)   

         p1=p3-1
         if(p1.le.0) p1=p1+np
c-mb
         z1=z-zf
         p2=(p3-1)+pif
         if (p2 .gt. np) p2=p2-np
         z2=z+(pif-zf)

         vx(r,p,w,z)=vx(r,p1,w,z1)+
     &   DBLE(pf)*(vx(r,p2,w,z2)-vx(r,p1,w,z1))/DBLE(pif)
         vy(r,p,w,z)=vy(r,p1,w,z1)+
     &   DBLE(pf)*(vy(r,p2,w,z2)-vy(r,p1,w,z1))/DBLE(pif)
         vz(r,p,w,z)=vz(r,p1,w,z1)+
     &   DBLE(pf)*(vz(r,p2,w,z2)-vz(r,p1,w,z1))/DBLE(pif)

      end do
      end do
      end do
      end do
 
c     Vortex filament tails (Approximate interpolation)
      do r=1,nr
      do w=1,nw-pw
c-mb see above
      p3=Pin-pif+1
      if(p3.le.0) p3=p3+np
      pf=0
      zf=0
      do p=p3,(p3-1)+(pif-1),1
      pf=pf+1
      zf=zf+1
      do z=nz-(pif-pf)+1,nz

         p1=p3-1
         if(p1.le.0) p1=p1+np
c
         z1=z-zf
         p2=(p3-1)+pif
         if (p2 .gt. np) p2=p2-np
         z2=nz

         vx(r,p,w,z)=vx(r,p1,w,z1)+
     &   DBLE(pf)*(vx(r,p2,w,z2)-vx(r,p1,w,z1))/DBLE(pif)
         vy(r,p,w,z)=vy(r,p1,w,z1)+
     &   DBLE(pf)*(vy(r,p2,w,z2)-vy(r,p1,w,z1))/DBLE(pif)
         vz(r,p,w,z)=vz(r,p1,w,z1)+
     &   DBLE(pf)*(vz(r,p2,w,z2)-vz(r,p1,w,z1))/DBLE(pif)

      end do
      end do
      end do
      end do
 
c     Special points -
      do r=1,nr
      do w=1,nw-pw
c-mb see above
      p3=Pin-pif+1
      if(p3.le.0) p3=p3+np
      z1=1
      p2=p3-1+pif
      if (p2 .gt. np) p2=p2-np
      zf=0
      do p=p3+1,(p3-1)+(pif-1),1
      pf=1
      do z=2,2+zf,1

         p1=p-pf
         if(p1.le.0) p1=p1+np
c-mb
         z2=(pif-zf+pf-1)

         vx(r,p,w,z)=vx(r,p1,w,z1)+
     &   DBLE(pf)*(vx(r,p2,w,z2)-vx(r,p1,w,z1))/DBLE(pif)
         vy(r,p,w,z)=vy(r,p1,w,z1)+
     &   DBLE(pf)*(vy(r,p2,w,z2)-vy(r,p1,w,z1))/DBLE(pif)
         vz(r,p,w,z)=vz(r,p1,w,z1)+
     &   DBLE(pf)*(vz(r,p2,w,z2)-vz(r,p1,w,z1))/DBLE(pif)

      pf=pf+1

      end do
      zf=zf+1
      end do
      end do
      end do

      endif
c ----------------------------------------------------------------------
c
      return
      end
