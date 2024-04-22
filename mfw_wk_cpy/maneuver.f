c----------------------------------------------------------------------|
c    <Copyright (c) 1997 by Ashish Bagai>
c     Copyright violation -MrB no changes
c----------------------------------------------------------------------|
c
      subroutine maneuver()
c----------------------------------------------------------------------|
c
      IMPLICIT none
c
      DOUBLE PRECISION pbar,qbar
c----------------------------------------------------------------------|
c
      common/ssmdata/pbar,qbar
      save/ssmdata/
c----------------------------------------------------------------------|
c
c     Set pbar/qbar in <ssm.input> to first trim-converge to 
c     steady state non-maneuver and then enter maneuver. Set pbar/qbar
c     in <flight.input> to always be in maneuver.
c
      namelist/ssm/pbar,qbar

      open(29,file='ssm.input',status='old')
      read(29,ssm)
      close(29)
c----------------------------------------------------------------------|
c

      if (pbar .ne. 0.D+0 .or. qbar .ne. 0.D+0) then
      write(*,*)
      write(*,*)'   *  **   ***    ****    *****    ****    ***   **  * 
     &  '
      write(*,*)'*                      MANEUVERING                     
     & *'
      write(*,*)'   *  **   ***    ****    *****    ****    ***   **  * 
     &  '
      write(*,'(a45,f10.5)')'    Steady Roll  Rate: P_bar =           '
     &,pbar
      write(*,'(a45,f10.5)')'    Steady Pitch Rate: Q_bar =           '
     &,qbar
      endif
c----------------------------------------------------------------------|
c
      return
      end
