c----------------------------------------------------------------------|
c    <Copyright (c) 1997 by Ashish Bagai>
c     Copyright violation -MrB never used never read no clue!
c----------------------------------------------------------------------|

      subroutine gridseq(gs,li,dp,gsv,gsn)
c     Grid-sequencing discretization blocks
c----------------------------------------------------------------------
c
      IMPLICIT none
c
      include 'gdimpar.inc'
c
      INTEGER gs,li,gsn,gs1,gs2,gs3,gs0
c
      DOUBLE PRECISION pi,dp,dp1,dp2,dp3,dp0,dz1,dz2,dz3,dz0,
     &       r00,r10,r20,r30
c
      DOUBLE PRECISION gsv(gdim,gdim)
c----------------------------------------------------------------------
c
      common/pidata/pi
      save/pidata/
c----------------------------------------------------------------------
c
      namelist/gseq/dp1,dp2,dp3,gs1,gs2,gs3,gs0
      namelist/gsue/dp1,dp2,dp3,dz1,dz2,dz3,gs1,gs2,gs3,gs0
c----------------------------------------------------------------------
c
c     Grid-sequencing = true
      if (gs .eq. 1) then
c
c      Using equal step size discretization resolutions for
c      adaptive grid-sequencing iteration blocks:
       if (li .ne. 1) then
c       Grid sequencing resolution blocks: 4 iteration blocks at
c       four equal resolutions. Read step sizes for 3 GS blocks
c       and No. of iterations for each of the 4 GS blocks -- 4th.
c       is for nominal resolution:
        open(13,file='gseq.input',status='old')
        read(13,gseq)
        close(13)
         dp1=dp1*pi/180.D+0
         dp2=dp2*pi/180.D+0
         dp3=dp3*pi/180.D+0
c        Nominal resolution:
         dp0=dp
c
c        Equal step sizes => dz=dp
         dz1=dp1
         dz2=dp2
         dz3=dp3
         dz0=dp0
c
c       Successive discretization ratios
        r10=dp1/dp0
        r20=dp2/dp0
        r30=dp3/dp0
        r00=dp0/dp0
c
c       Number of grid-sequenced iteration blocks
        gsn=4
c
c       Store all 16 variables for all 4 iteration blocks in a single
c       2-D array:
        gsv(1,1)=dp1
        gsv(1,2)=dz1
        gsv(1,3)=r10
        gsv(1,4)=DBLE(gs1)
c
        gsv(2,1)=dp2
        gsv(2,2)=dz2
        gsv(2,3)=r20
        gsv(2,4)=DBLE(gs2)
c
        gsv(3,1)=dp3
        gsv(3,2)=dz3
        gsv(3,3)=r30
        gsv(3,4)=DBLE(gs3)
c
        gsv(4,1)=dp0
        gsv(4,2)=dz0
        gsv(4,3)=r00
        gsv(4,4)=DBLE(gs0)
c
       endif
c
c      Using unequal step size discretization resolutions for
c      adaptive grid-sequencing iteration blocks:
       if (li .eq. 1) then
c
c       Grid-sequencing resolution blocks: 4 iteration blocks at
c       3 unequal resolutions. 
c       Unequal step sizes => dz > dp (max. advantage)
c       Unequal step sizes => dz < dp (max. accuracy)
        open(14,file='gsue.input',status='old')
        read(14,gsue)
        close(14)
         dp1=dp1*pi/180.D+0
         dp2=dp2*pi/180.D+0
         dp3=dp3*pi/180.D+0
         dz1=dz1*pi/180.D+0
         dz2=dz2*pi/180.D+0
         dz3=dz3*pi/180.D+0
c
        dp0=dp
        dz0=dp0
c
        r10=dp1/dp0
        r20=dp2/dp0
        r30=dp3/dp0
        r00=dp0/dp0
c
        gsn=4
c
        gsv(1,1)=dp1
        gsv(1,2)=dz1
        gsv(1,3)=r10
        gsv(1,4)=DBLE(gs1)
c
        gsv(2,1)=dp2
        gsv(2,2)=dz2
        gsv(2,3)=r20
        gsv(2,4)=DBLE(gs2)
c
        gsv(3,1)=dp3
        gsv(3,2)=dz3
        gsv(3,3)=r30
        gsv(3,4)=DBLE(gs3)
c
        gsv(4,1)=dp0
        gsv(4,2)=dz0
        gsv(4,3)=r00
        gsv(4,4)=DBLE(gs0)
c
       endif
c
      endif
c----------------------------------------------------------------------
c
      return
      end
