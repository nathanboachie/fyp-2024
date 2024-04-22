c----------------------------------------------------------------------|
c    <Copyright (c) 1998 by Mahendra Bhagwat>
c----------------------------------------------------------------------|
c
      subroutine
     &     wrtrms(nr,nw,pw,n,rms,rms0,gammaRMS,betaRMS,time1,tarray)
c----------------------------------------------------------------------|
c>mb note that this will give error if nw-pw > 5 and nr > 2
c     Wake convergence history.
c----------------------------------------------------------------------|
c
      IMPLICIT none
c
         include 'rdimpar.inc'
         include 'wdimpar.inc'
c
      INTEGER nr,nw,pw,n,r,w
c
      DOUBLE PRECISION rms(rdim,wdim),rms0(rdim,wdim)
      DOUBLE PRECISION gammaRMS
      DOUBLE PRECISION betaRMS(rdim)
      DOUBLE PRECISION time1,dtime,tarray(2)
      external dtime
c----------------------------------------------------------------------|
c
         time1=time1+dtime(tarray)
         write(21,10)n,time1,(rms(r,1),rms0(r,1)*rms(r,1),r=1,nr)
     &               ,gammaRMS,(betaRMS(r),r=1,nr)
10     format(i3,1x,e11.6,8(1x,e14.8))
c-mb                      ^ this should be 3*nr+2
c----------------------------------------------------------------------|
c
      return
      end 

