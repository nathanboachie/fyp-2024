c----------------------------------------------------------------------|
c    <Copyright (c) 1997 by Ashish Bagai>
c     Copyright violation -MrB small corr
c----------------------------------------------------------------------|
c-mb I personally do not believe in this myth of centroid of vorticity
c-mb but well I will need to find something cleverer or else just ignore
c----------------------------------------------------------------------|
c
      subroutine centroid(nr,np,ns,dz,cvt,bvx,bvy,bvz,pcx,pcy,pcz,zif)
c     Centroid of vorticity shift -- applied only to tip vortex.
c----------------------------------------------------------------------|
c
      IMPLICIT none
c
        include 'adimpar.inc'
        include 'rdimpar.inc'
        include 'sdimpar.inc'
        include 'wdimpar.inc'
        include 'zdimpar.inc'
c
      INTEGER ncz,w,z,r,p,nr,np,ns,s 
c
      DOUBLE PRECISION pi,cvt,dz,Dx0,Dy0,Dz0,mx,my,mz,snx,sny,snz
c
      DOUBLE PRECISION bvx(rdim,adim,sdim),
     &       bvy(rdim,adim,sdim),
     &       bvz(rdim,adim,sdim)
c
      DOUBLE PRECISION Rx(rdim,adim),
     &       Ry(rdim,adim),
     &       Rz(rdim,adim)
c
      DOUBLE PRECISION pcx(rdim,adim,wdim,zdim),
     &       pcy(rdim,adim,wdim,zdim),
     &       pcz(rdim,adim,wdim,zdim)
c
      INTEGER zif
c
c----------------------------------------------------------------------
c
      common/pidata/pi
      save/pidata/
c----------------------------------------------------------------------
c
      s=ns+1

      do r=1,nr
      do p=1,np
         Rx(r,p)=bvx(r,p,s)
         Ry(r,p)=bvy(r,p,s)
         Rz(r,p)=bvz(r,p,s)
      end do
      end do

c-mb VFF correctn            __________
      ncz=int(cvt*2.D+0*pi/(dz/DBLE(zif)))

      do r=1,nr
      do p=1,np

         w=1
         z=1

         Dx0=abs(Rx(r,p))-abs(pcx(r,p,w,z))
         Dy0=abs(Ry(r,p))-abs(pcy(r,p,w,z))
         Dz0=abs(Rz(r,p))-abs(pcz(r,p,w,z))
 
      do z=1,ncz-1

         mx=Dx0/DBLE(ncz-1)
         my=Dy0/DBLE(ncz-1)
         mz=Dz0/DBLE(ncz-1)

         snx=pcx(r,p,w,z)/abs(pcx(r,p,w,z)+1.D-8)
         sny=pcy(r,p,w,z)/abs(pcy(r,p,w,z)+1.D-8)
         snz=pcz(r,p,w,z)/abs(pcz(r,p,w,z)+1.D-8)

         pcx(r,p,w,z)=snx*DBLE(ncz-z)*mx+pcx(r,p,w,z)
         pcy(r,p,w,z)=sny*DBLE(ncz-z)*my+pcy(r,p,w,z)
         pcz(r,p,w,z)=snz*DBLE(ncz-z)*mz+pcz(r,p,w,z)

       end do

       end do
       end do
c----------------------------------------------------------------------|
c
      return
      end
