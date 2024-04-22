C-$   $Id: readgb.f,v 1.2.6.1 2005/03/10 18:49:33 shreyas Exp $	

      subroutine readgb(anb,dp,nr,np,ns,rad,gb,pif,zif)


      IMPLICIT none

      include 'adimpar.inc'
      include 'rdimpar.inc'
      include 'sdimpar.inc'

      INTEGER nr,np,ns,r,p,s,dnp

      DOUBLE PRECISION anb,rad,dp

      DOUBLE PRECISION gb(rdim,adim,sdim)

      DOUBLE PRECISION tmp1,tmp2
      INTEGER pif,zif


      dnp=nint(anb*100000.D+0/(dp/DBLE(pif))/100000.D+0)
      open(44,file='IWG_b.data',status='old')
      do r=1,nr
         do p=1,np              !,dnp  need gb at all azimuths
            do s=1,ns
               read(44,'(f15.8,x,f15.8)')tmp1,tmp2
               gb(r,p,s)=tmp2
            end do
            read(44,*)
         end do
         read(44,*)
         read(44,*)
      end do
      close(44)


      return
      end 
