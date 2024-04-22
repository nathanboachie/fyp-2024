C-$	$Id: wrtwke.f,v 1.1.1.1.6.1 2005/03/10 18:49:33 shreyas Exp $	

      subroutine wrtwke(anb,sco,rad,dp,dz,nr,np,nw,pw,nz,nzt,pif,zif,
     $     pcx,pcy,pcz)

C-$   Write out the top,side and rear views of wake geometry

      IMPLICIT none

      include 'adimpar.inc'
      include 'rdimpar.inc'
      include 'wdimpar.inc'
      include 'zdimpar.inc'

      INTEGER nr,np,nw,pw,nz,nzt,pif,zif,r,p,w,z,dnp

      DOUBLE PRECISION anb,rad,dp,dz,pi,xp,yp

      DOUBLE PRECISION pcx(rdim,adim,wdim,zdim), pcy(rdim,adim,wdim,zdim
     $     ), pcz(rdim,adim,wdim,zdim)

      DOUBLE PRECISION gt

      CHARACTER*1 sco 


      common/pidata/pi
      save/pidata/


      dnp=nint(anb/(dp/DBLE(pif)))

      open(81,file='FWT.dat',status='unknown')
      open(82,file='FWS.dat',status='unknown')
      open(83,file='FWR.dat',status='unknown')
      do r=1,nr
         do p=1,np,dnp
            do w=1,nw-pw
               do z=1,nzt
                  write(81,*)pcx(r,p,w,z)/rad,pcy(r,p,w,z)/rad
                  write(82,*)pcx(r,p,w,z)/rad,pcz(r,p,w,z)/rad
                  write(83,*)pcy(r,p,w,z)/rad,pcz(r,p,w,z)/rad
               end do
               write(81,*)
               write(82,*)
               write(83,*)
            end do
            write(81,*)
            write(82,*)
            write(83,*)
         end do
         write(81,*)
         write(82,*)
         write(83,*)
      end do

      if (sco .eq. 'y') then
         do r=1,nr
            do p=1,np,dnp
               do z=1,nz-1,2
                  do w=1,nw-pw
                     write(81,*)pcx(r,p,w,z)/rad,pcy(r,p,w,z)/rad
                     write(82,*)pcx(r,p,w,z)/rad,pcz(r,p,w,z)/rad
                     write(83,*)pcy(r,p,w,z)/rad,pcz(r,p,w,z)/rad
                  end do
                  do w=nw-pw,1,-1
                     write(81,*)pcx(r,p,w,z+1)/rad,pcy(r,p,w,z+1)/rad
                     write(82,*)pcx(r,p,w,z+1)/rad,pcz(r,p,w,z+1)/rad
                     write(83,*)pcy(r,p,w,z+1)/rad,pcz(r,p,w,z+1)/rad
                  end do
               end do
               write(81,*)
               write(82,*)
               write(83,*)
            end do
            write(81,*)
            write(82,*)
            write(83,*)
         end do
      endif

      close(81)
      close(82)
      close(83)

      open(84,file='FWISO.dat',status='unknown')
      do r=1,nr
         do z=1,nz
            do w=1,nw-pw
               write(84,'(8(f10.5,1x))')(
     &              (pcx(r,p,w,z)/rad)*cos(30.D+0*pi/180.D+0)+
     &              (pcy(r,p,w,z)/rad)*sin(30.D+0*pi/180.D+0),
     &              (pcz(r,p,w,z)/rad)*cos(30.D+0*pi/180.D+0)+
     &              (pcy(r,p,w,z)/rad)*sin(30.D+0*pi/180.D+0)
     &              ,p=1,np,dnp)
            end do
         end do
      end do
      close(84)

      return
      end 
