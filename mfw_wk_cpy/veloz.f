c----------------------------------------------------------------------|
      DOUBLE PRECISION function veloz(x,y,z,wa)
c----------------------------------------------------------------------|
c Non-uniform freestream velocity in the z direction as a function 
c of x,y,z; i.e. veloz=f(x,y,z) -- fixed frame (global) coordinates.
c
c-mb  here the signs are right!
c
c----------------------------------------------------------------------|
c
      IMPLICIT none
c
      DOUBLE PRECISION x,y,z,pbar,qbar
      INTEGER wa
c----------------------------------------------------------------------|
      common/ssmdata/pbar,qbar
      save/ssmdata/
c----------------------------------------------------------------------|
      veloz=1.0D+0*x*cos(1.0D+0*DBLE(wa-1)*15.D+0/(57.3D+0*220.D+0))
      veloz=0.0D+0+qbar*x-pbar*y
c----------------------------------------------------------------------|
c
      return
      end
