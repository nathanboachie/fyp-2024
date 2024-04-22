C-$	$Id: beddoesa.f,v 1.4.2.2 2005/03/23 21:24:38 shreyas Exp $	

         subroutine beddoesa(alpha,Mno,Cl,Cd,clsl,alpha0)

C-$   Beddoes non-linear airfoil model


      DOUBLE PRECISION alpha,Mno,Cl,Cd
      DOUBLE PRECISION pi,effalpha,clsl,effalp1,absaoa
      DOUBLE PRECISION alpha0,alpha1,alphadd,S1,S2,ff,df,cd0,KD,Cn
      DOUBLE PRECISION myerf
      DOUBLE PRECISION Mno0
      external myerf


      pi=4.D+0*atan(1.D+0)

      Mno0=Mno

      if(Mno.lt.0.3D+0)Mno=0.3D+0
      if(Mno.gt.0.85D+0)Mno=0.85D+0
      alpha1=21.5D+0-25.D+0*Mno +2.D+0*exp(-((Mno-0.65D+0)/0.125D+0)**2
     $     .D+0)
      alpha1=alpha1*pi/180.D+0
      S1=1.8D+0*exp(-((Mno-0.45D+0)/0.3D+0)**2.D+0)
      S2=3.6D+0*exp(-((Mno-0.525D+0)/0.25D+0)**2.D+0)

      effalpha=alpha-alpha0
      absaoa=abs(alpha)
      if(absaoa.le.alpha1) then
         ff=1.D+0-0.3D+0*exp((absaoa-alpha1)*180.D+0/pi/S1)
      else
         ff=.04D+0+0.66D+0*exp((alpha1-absaoa)*180.D+0/pi/S2)
      end if
      
      Cn=2.D+0*pi/sqrt(1.D+0-Mno**2.D+0)* ( (1.D+0+sqrt(ff))/2.D+0 )**2
     $     .D+0*effalpha

c      cd0=0.01D+0 + 0.002D+0*myerf( 50.D+0*(Mno-0.8D+0) )
      cd0=0.011D+0! + 0.002D+0*myerf( 50.D+0*(Mno-0.8D+0) )
      df=6.1D+0 - 7.D+0*Mno + .5D+0*exp(-((Mno-.6D+0)/.125D+0)**2.D+0)
      alphadd=16.0D+0 - 20.D+0*Mno + .5D+0*exp(-((Mno-.6D+0)/.125D+0)**2
     $     .D+0)
      alphadd=alphadd*pi/180.D+0
      KD=0.D+0
      if(abs(alpha).ge.alphadd) KD=2.7D+0*exp(-df*ff)

      if(alpha .ge. 0.0D0) then
         effalp1=effalpha-alphadd
      else
         effalp1=alphadd-effalpha
      endif

      Cl=Cn*cos(effalpha)
      Cd=cd0 + 0.035D+0*Cn*sin(effalpha) + KD*Cn*sin(effalp1)
      clsl=2.D+0*pi/sqrt(1.D+0-Mno**2.D+0)* ( (1.D+0 +sqrt(ff))/2.D+0 )
     $     **2.D+0
      Mno=Mno0

      return
      end
