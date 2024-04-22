C-$	$Id: airfoilpsi.f,v 1.6.2.2 2005/03/23 21:24:38 shreyas Exp $	

      subroutine airfoilpsi(bcx,aoae,Mach,Reynolds,nr,np,ns, clairfoil
     $     ,cdairfoil,claoa,compress,p,lin,alpha0)

C-$   Compute Cl and Cd for a given airfoil section 


      include 'rdimpar.inc'
      include 'adimpar.inc'
      include 'sdimpar.inc'

      INTEGER r,p,s,nr,np,ns
      DOUBLE PRECISION bcx(rdim,sdim),aoae(rdim,adim,sdim), Mach(rdim
     $     ,adim,sdim),Reynolds(rdim,adim,sdim)

      DOUBLE PRECISION clairfoil(rdim,adim,sdim), cdairfoil(rdim,adim
     $     ,sdim), claoa(rdim,adim,sdim),compress(rdim,adim,sdim)

      DOUBLE PRECISION pi,effalpha,lmt12
      DOUBLE PRECISION alpha0,alpha1,alphadd,S1,S2,ff,df,cd0,KD,Cn
      DOUBLE PRECISION myerf
      DOUBLE PRECISION Mno
      INTEGER lin
      external myerf


      pi=4.D+0*atan(1.D+0)

      do r=1,nr
         do s=1,ns
            effalpha=aoae(r,p,s)-alpha0
            Mno=Mach(r,p,s)
            if(lin.eq.1) then
               if(Mno.lt.0.3D+0)Mno=0.3D+0
               if(Mno.gt.0.85D+0)Mno=0.85D+0
               claoa(r,p,s)=2.D+0*pi/sqrt(1.D+0-Mno**2.D+0)
               clairfoil(r,p,s)=claoa(r,p,s)*(effalpha)
               cdairfoil(r,p,s)=0.008D+0 +0.400D+0*effalpha**2/sqrt(1.D
     $              +0-Mno**2.D+0)
               compress(r,p,s)=claoa(r,p,s)/2.D+0/pi
               cdairfoil(r,p,s)=0.008D+0
            else
               call beddoesa(aoae(r,p,s),Mno,clairfoil(r,p,s)
     $              ,cdairfoil(r,p,s),claoa(r,p,s),alpha0)
               compress(r,p,s)=1.0D0/sqrt(1.0D0-Mach(r,p,s)**2.0D0)
            end if              !lin.eq.1
         end do
      end do


      return
      end
