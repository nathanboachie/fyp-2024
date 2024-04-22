C-$	$Id: wrtwkmb.f,v 1.1.1.1.6.2 2005/05/10 16:14:30 shreyas Exp $	

      subroutine wrtwkmb(anb,sco,rad,dp,dz,nr,np,nw,nz,nzt,pif,zif, pcx
     $     ,pcy,pcz,gv, eqtim, nb, om, gav, kinv, dcy)

C-$   Complete wake geometry information output


      IMPLICIT none

      include 'adimpar.inc'
      include 'rdimpar.inc'
      include 'wdimpar.inc'
      include 'zdimpar.inc'
      include 'kdimpar.inc'

      INTEGER nr,np,nw,nz,nzt,pif,zif,r,p,w,z,dnp,nb

      DOUBLE PRECISION anb,rad,dp,dz,pi, om,gav, kinv, dcy, factor,crad,
     $     crad1

      DOUBLE PRECISION pcx(rdim,adim,wdim,zdim), pcy(rdim,adim,wdim,zdim
     $     ), pcz(rdim,adim,wdim,zdim), eqtim(rdim,adim,wdim,zdim,kdim)

      DOUBLE PRECISION gv(rdim,adim,wdim,zdim,kdim)

      CHARACTER*1 sco 

      INTEGER now,itmp




      common/pidata/pi
      save/pidata/

      dnp=nint(anb/(dp/DBLE(pif)))
      open(86,file='FWGEOM.dat',status='unknown')
      do r=1,nr
         do p=1,np
            do w=1,nw
               do z=1,nzt

                  now=w
                  factor=1.0D+0
                  crad1=1.D+0+1.0D+0*dcy*factor*gav/kinv
                  crad=0.00855D+0*sqrt(crad1*eqtim(r,p,w,z,1)*nb/om)
                  write(86,10)((dp/DBLE(pif))*(p-1)*180.D+0/pi), ((dz
     $                 /DBLE(zif))*(z-1)*180.D+0/pi), (pcx(r,p,now,z)
     $                 /rad), (pcy(r,p,now,z)/rad), (pcz(r,p,now,z)/rad)
     $                 , (gv(r,p,now,z,1)*1.D+0), (crad)
               end do
               write(86,*)
            end do
            write(86,*)
         end do
         write(86,*)
      end do

      close(86)

10    format(6(f15.8,1x),E15.6)


      return
      end 
