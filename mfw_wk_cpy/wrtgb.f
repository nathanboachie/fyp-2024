c----------------------------------------------------------------------|
c    <Copyright (c) 1997 by Ashish Bagai>
c     Copyright violation -MrB small changes
c----------------------------------------------------------------------|
c
      subroutine wrtgb(anb,dp,nr,np,ns,rad,cpx,cpz,rr_x,rr_z,
     &           asr,gb,pif,zif)
c----------------------------------------------------------------------|
c
      IMPLICIT none
c
         include 'adimpar.inc'
         include 'rdimpar.inc'
         include 'sdimpar.inc'
c
      INTEGER nr,np,ns,r,p,s,dnp
c
      DOUBLE PRECISION anb,rad,dp
c
      DOUBLE PRECISION cpx(rdim,adim,sdim),cpz(rdim,adim,sdim),
     &       gb(rdim,adim,sdim)
c
      DOUBLE PRECISION rr_x(rdim),rr_z(rdim),asr(rdim)
c
      INTEGER pif,zif
c----------------------------------------------------------------------|
c
c-mb VFF correctn
      dnp=nint(anb*100000.D+0/(dp/DBLE(pif))/100000.D+0)
      open(43,file='G_b.dat',status='unknown')
      do r=1,nr
      do p=1,np,dnp
      do s=1,ns

         write(43,'(31f10.4)')
     & ( (cpx(r,1,s)-rr_x(r))*cos(-asr(r))
     &  -(cpz(r,1,s)-rr_z(r))*sin(-asr(r)) )/rad,gb(r,p,s)

      end do
         write(43,*)
      end do
         write(43,*)
         write(43,*)
      end do
      close(43)
c----------------------------------------------------------------------|
c
      return
      end 
