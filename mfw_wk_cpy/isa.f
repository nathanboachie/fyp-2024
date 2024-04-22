C-$	$Id: isa.f,v 1.1.4.2 2005/03/23 21:24:38 shreyas Exp $	

      subroutine isa(den,sonic,altitude,isaplus, kinv)

      implicit none

      double precision den, sonic, altitude, isaplus, kinv

      double precision tsl, talt, rhosl, rhoalt, prsl, pralt

      double precision gravity, lapse, gasconst, gbyar, visc


      gravity=9.81D0              !acceleration due to gravity
      lapse=-0.0065D0           !lapse rate
      gasconst=287.05D0           !Gas constant for air
      gbyar=gravity/(lapse*gasconst)
      

      tsl=288.15D0                !ISA SL temperature in Kelvin
      prsl=1.013D+5             !Pressure at sea level
      rhosl=1.2256D+0           !Density at sea level
      
      talt=tsl+isaplus+lapse*altitude !ISA+ temp at req. altitude
      pralt=prsl*(tsl/talt)**gbyar !pressure 
      rhoalt=pralt/(gasconst*talt) !density

      !Sutherland's law
      visc= 1.458D-6*sqrt(talt)/(1.0D+0+110.0D0/talt)
      kinv=visc/rhoalt
      sonic=sqrt(1.4D0*pralt/rhoalt)
      den=rhoalt
      
      return
      end
      
      
