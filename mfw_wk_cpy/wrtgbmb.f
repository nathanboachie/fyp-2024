C-$	$Id: wrtgbmb.f,v 1.1.1.1.6.1 2005/03/10 18:49:33 shreyas Exp $	

      subroutine wrtgbmb(anb,dp,nr,np,ns,rad,cpx,cpz,rr_x,rr_z, asr,lsr
     $     ,gb,pif,zif)

C-$   Output bound circulation information


      IMPLICIT none

      include 'adimpar.inc'
      include 'rdimpar.inc'
      include 'sdimpar.inc'

      INTEGER nr,np,ns,r,p,s,dnp

      DOUBLE PRECISION anb,rad,dp

      DOUBLE PRECISION cpx(rdim,adim,sdim),cpz(rdim,adim,sdim),
     &     gb(rdim,adim,sdim)

      DOUBLE PRECISION rr_x(rdim),rr_z(rdim),asr(rdim),lsr(rdim)

      INTEGER pif,zif


      dnp=nint(anb/(dp/DBLE(pif)))
      open(43,file='FWG_b.dat',status='unknown')
      do r=1,nr
         do p=1,np              !,dnp need to write for all azimuths
            do s=1,ns

C/MR  This transformation is still valid without lsr for main rotor
C/MR  but need to look closer into other orientations. 
               write(43,'(f15.8,x,f15.8)') ( (cpx(r,1,s)-rr_x(r))*cos(
     $              -asr(r)) -(cpz(r,1,s)-rr_z(r))*sin(-asr(r)) )/rad
     $              ,gb(r,p,s)

            end do
            write(43,*)
         end do
         write(43,*)
         write(43,*)
      end do
      close(43)


      return
      end 
