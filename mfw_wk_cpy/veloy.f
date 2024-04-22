c----------------------------------------------------------------------|
      DOUBLE PRECISION function veloy(x,y,z,wa)
c----------------------------------------------------------------------|
c Non-uniform freestream velocity in the y direction as a function 
c of x,y,z; i.e. veloy=f(x,y,z) -- fixed frame (global) coordinates.
c
c-mb for whatever reason, the code I was given had the wrong sign
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
      veloy=0.D+0+pbar*z
c----------------------------------------------------------------------|
c
      return
      end
