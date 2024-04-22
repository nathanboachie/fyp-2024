C-$	$Id: myerf.f,v 1.1.1.1.6.1 2005/03/23 21:24:38 shreyas Exp $	

      DOUBLE PRECISION function myerf(x)

c-mb  (approximate error funciton with error lt 1.2E-7 everywhere :)
c-mb   stolen from numerical recipes

      IMPLICIT none

      DOUBLE PRECISION x,z,t,erfcc

      myerf=1.D+0
      z=abs(x)
      t=1.D+0/(1.D+0+0.5D+0*z)
      erfcc=t*exp( -z*z -1.26551223D+0 + t*( 1.00002368D+0 + t*(
     $     .37409196D+0 + t*( .09678418D+0 + t*( -.18628806D+0 + t*(
     $     .27886807D+0 + t*( -1.13520398D+0 + t*( 1.48851587D+0+t*(-
     $     .82215223D+0+t*.17087277D+0)))))))))
      if (x.lt.0.D+0) erfcc=2.D+0-erfcc
      myerf=1.D+0-erfcc

      return
      end
