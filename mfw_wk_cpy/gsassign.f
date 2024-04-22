C-$	$Id: gsassign.f,v 1.1.1.1.6.1 2005/03/10 18:49:33 shreyas Exp $	

      subroutine gsassign(gsn,gsno,gsv,gs,li,dp,dz,pif,zif,n1)

c     Grid-sequencing iteration resolution and index
c     increment assignment between iteration blocks. 


      IMPLICIT none

      include 'gdimpar.inc'

      INTEGER gsn,gsno,n1,gs,li,pif,zif

      DOUBLE PRECISION dp,dz

      DOUBLE PRECISION gsv(gdim,gdim)


c     Number of variables = 4
 
c     Equal step sizes
      if (li .ne. 1) then
         dp=gsv(gsno,1) 
         dz=gsv(gsno,2) 
         pif=nint( gsv(gsno,3) )
         zif=pif
         n1=nint( gsv(gsno,4) )
      endif

c     Unequal step sizes
      if (li .eq. 1) then
         dp=gsv(gsno,1) 
         dz=gsv(gsno,2) 
         pif=nint( gsv(gsno,3) )
         zif=nint(gsv(gsno,2)/gsv(gsn,2))
         n1=nint(gsv(gsno,4))
      endif

      gsno=gsno+1


      return
      end
