C-$   $Id: readwake.f,v 1.2.6.2 2005/05/10 16:14:30 shreyas Exp $	

      subroutine readwake(anb,sco,rad,dp,dz,nr,np,nw,pw,nz,nzt,pif,zif,
     $     pcx,pcy,pcz, eqtim, om, kinv, gav, dcy, Nb)

C-$   Read initial wake geometry


      IMPLICIT none

      include 'adimpar.inc'
      include 'rdimpar.inc'
      include 'wdimpar.inc'
      include 'zdimpar.inc'
      include 'kdimpar.inc'

      INTEGER nr,np,nw,pw,nz,nzt,pif,zif,r,p,w,z,dnp, Nb

      DOUBLE PRECISION anb,rad,dp,dz

      DOUBLE PRECISION pcx(rdim,adim,wdim,zdim),
     &     pcy(rdim,adim,wdim,zdim),
     &     pcz(rdim,adim,wdim,zdim),
     &     eqtim(rdim, adim, wdim, zdim, kdim)

      DOUBLE PRECISION temp1,temp2,temp3,temp4,temp5,temp6, temp7

      DOUBLE PRECISION kinv, gav, om, trb, factor, dcy, zet0

      CHARACTER*1 sco 


      dnp=nint(anb/(dp/DBLE(pif)))

      factor=1.0D+0
      open(86,file='IWGEOM.data',status='old')
      do r=1,nr
         do p=1,np
            do w=1,nw-pw
               do z=1,nzt

                  read(86,10)temp1,temp2,
     &                 temp3,temp4,temp5,temp6,temp7

                  pcx(r,p,w,z)=temp3*rad
                  pcy(r,p,w,z)=temp4*rad
                  pcz(r,p,w,z)=temp5*rad
                  trb=1.D+0+dcy*factor*gav/kinv
                  zet0=(temp7/.00855D+0)**2*om/trb/Nb
                  eqtim(r,p,w,z,1)=zet0
               end do
               read(86,*)
            end do
            read(86,*)
         end do
         read(86,*)
      end do

      close(86)

10    format(6(f15.8,1x),E15.6)


      return
      end 
