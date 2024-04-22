C-$	$Id: velox.f,v 1.1.1.1.6.1 2005/03/10 18:49:33 shreyas Exp $	


      DOUBLE PRECISION function velox(x,y,z,wa)

      IMPLICIT none

      DOUBLE PRECISION x,y,z,pbar,qbar
      INTEGER wa

      common/ssmdata/pbar,qbar
      save/ssmdata/

      velox=-0.0D+0-qbar*z


      return
      end
