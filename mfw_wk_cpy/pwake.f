C-$	$Id: pwake.f,v 1.1.1.1.6.2 2005/06/11 15:57:37 shreyas Exp $	

      subroutine pwake(nr,np,nw,nzt,pif,zif,ip,dp,dz,mu,muc,om, vinfx
     $     ,vinfy,vinfz,asr,lsr,ct,b0r,b1c,b1s,anb,rad,crd,rr_x,rr_y 
     $     ,rr_z,rv,twt,taperst,taper,lami,pcx,pcy,pcz)

c >>  Prescribed, initial, epicycloidal, undistorted wake geometry. <<


      IMPLICIT none

      include 'adimpar.inc'
      include 'rdimpar.inc'
      include 'sdimpar.inc'
      include 'wdimpar.inc'
      include 'zdimpar.inc'

      INTEGER nr,np,nzt,pif,zif,r,p,z,nw,w

      DOUBLE PRECISION ip(rdim), dp,dz,mu,muc,mux,muz,lam1,lam0,rerr,om
     $     ,rad, mui,vinfx,muj,vinfy,muk,vinfz,az,tt11,tt21,tt31, wage
     $     ,sn,ze,anb,velox,veloy,veloz,b0

      DOUBLE PRECISION asr(rdim),ct(rdim), b0r(rdim),b1c(rdim),b1s(rdim)
     $     , uiv(rdim),lami(rdim),rr_x(rdim),rr_y(rdim),rr_z(rdim)
     $     ,lsr(rdim)

      DOUBLE PRECISION rv(rdim,adim,wdim)

      DOUBLE PRECISION pcx(rdim,adim,wdim,zdim), pcy(rdim,adim,wdim,zdim
     $     ), pcz(rdim,adim,wdim,zdim)

      DOUBLE PRECISION Aland,Lland,k1land,k2land,crd,pi,twt,
     $     K1r1land,K2r1land,K2r0land,zs0,zs1

      DOUBLE PRECISION planf
      DOUBLE PRECISION taper,taperst

      INTEGER lambdaitn
      INTEGER rotgeo

      external planf

      external velox
      external veloy
      external veloz

      common/rotdata/rotgeo
      save/rotdata/


      pi=4.D+0*atan(1.D+0)

      do r=1,nr

         Aland=1.D+0/sqrt(2.D+0)
         Aland=.78D+0
         if(mu.ne.0.D+0)Aland=1.0D+0
         if(muc.ne.0.D+0)Aland=1.0D+0
         Lland=0.145D+0+27.D+0*ct(r)
         k1land=0.25D+0*ct(r)/ (2.D+0/anb*crd*planf(-1.D+0,taperst,taper
     $        )/rad)
         k2land=sqrt(ct(r))
         k1land=0.25D+0*( ct(r)/ (2.D+0/anb*crd*planf(-1.D+0,taperst
     $        ,taper)/rad) + 0.001D+0*twt )
         k2land=sqrt(ct(r))*(1.D+0+0.01D+0*twt)
         K1r1land=-2.2D+0*sqrt(ct(r)/2.D+0) 
         K2r1land=-2.7D+0*sqrt(ct(r)/2.D+0) 
         K2r0land=(twt/128.D+0*(0.45D+0*twt+18.D+0)) *sqrt(ct(r)/2.D+0) 

c >>     Uniform inflow distribution: <<
c         mux=mu*cos(-asr(r))+muc*sin(-asr(r))
c         muz=mu*sin(-asr(r))-muc*cos(-asr(r))
C/MR     Correcting signs and adding lateral angle
C/MR     asr is positive back, so this neg sign must not be canceled
C/MR     in input file.
         mux=mu*cos(-asr(r))-muc*sin(-asr(r))
         muz=mu*sin(-asr(r))*cos(lsr(r))+muc*cos(-asr(r))*cos(lsr(r))
         
         lam1=sqrt(ct(r)*0.5D+0)
         if(mux.eq.0.) then
            lam0=-muz+lam1
         else
            lambdaitn=0
            rerr=1.0D+0
            do while (lambdaitn.le.100 .OR. rerr .ge. 0.01)
               lam0=ct(r)/(2.D+0*sqrt(mux**2+(muz-lam1)**2))
               rerr=abs(lam0-lam1)/lam0
               lam1=lam0
               lambdaitn=lambdaitn+1
            enddo
         endif
         uiv(r)=lam0*om*rad
         lami(r)=uiv(r)/(om*rad)
         uiv(r)=-uiv(r)

c >>     Freestream velocity components: <<
         mui=vinfx/(om*rad)
         muj=vinfy/(om*rad)
         muk=vinfz/(om*rad)
         
         sn=1.D+0
         if(rotgeo.eq.1 .or. rotgeo.eq.2) sn=(-1.D+0)**(r-1)
         do p=1,np
            az=DBLE(p-1)*dp/DBLE(pif)
            b0=b0r(r)+b1c(r)*cos(az)+b1s(r)*sin(az)
            do w=1,nw
               do z=1,nzt
                  ze=DBLE(z-1)*dz/DBLE(zif)
                  wage=sn*(az + ip(r) -ze)
                  
c                  tt11=cos(-asr(r))*cos(wage)*cos(b0)+ sin(-asr(r))
c     $                 *sin(b0)
c                  tt21=sin(wage)*cos(b0)
c                  tt31=-sin(-asr(r))*cos(wage)*cos(b0)+ cos(-asr(r))
c     $                 *sin(b0)
     
C/MR    Adding lateral tilt 
                  tt11=cos(-asr(r))*cos(wage)*cos(b0) - 
     $             sin(-asr(r))*sin(lsr(r))*sin(wage)*cos(b0) +
     $             + sin(-asr(r))*cos(lsr(r))*sin(b0)
                  tt21=sin(wage)*cos(b0)
                  tt31=-sin(-asr(r))*cos(wage)*cos(b0)+ cos(-asr(r))
     $                 *sin(b0)


                  pcx(r,p,w,z)=rr_x(r)+rad*(mui+ velox(pcx(r,p,w,z)
     $                 ,pcy(r,p,w,z),pcz(r,p,w,z),z)/ (om*rad)-lami(r)
     $                 *sin(-asr(r)))*ze +rv(r,p,w)*tt11 *(Aland+(1.D+0
     $                 -Aland)*exp(-Lland*ze))

                  pcy(r,p,w,z)=rr_y(r)+rad*(muj+ veloy(pcx(r,p,w,z)
     $                 ,pcy(r,p,w,z),pcz(r,p,w,z),z)/ (om*rad))*ze+rv(r
     $                 ,p,w)*tt21 *(Aland+(1.D+0-Aland)*exp(-Lland*ze))

                  pcz(r,p,w,z)=rr_z(r)+rad*(muk+ veloz(pcx(r,p,w,z)
     $                 ,pcy(r,p,w,z),pcz(r,p,w,z),z)/ (om*rad)-lami(r)
     $                 *cos(-asr(r)))*ze +rv(r,p,w)*tt31 *(Aland+(1.D+0
     $                 -Aland)*exp(-Lland*ze))

                  if(mu.eq.0.D+0 .AND. muc.eq.0.D+0) then
                     if(w.eq.1) then

                        pcz(r,p,w,z)=rr_z(r)-rad*k1land*ze
                        pcz(r,p,w,z)=rr_z(r)+rad*(muk+ veloz(pcx(r,p,w,z
     $                       ),pcy(r,p,w,z),pcz(r,p,w,z),z)/ (om*rad)
     $                       -k1land*cos(-asr(r)))*ze +rv(r,p,w)*tt31
     $                       *(Aland+(1.D+0-Aland)*exp(-Lland*ze))

                        if(ze.gt.anb) then
                           pcz(r,p,w,z)=rr_z(r)-rad*(k1land*anb+k2land
     $                          *(ze-anb))
                           pcz(r,p,w,z)=rr_z(r)+rad*(muk+ veloz(pcx(r,p
     $                          ,w,z),pcy(r,p,w,z),pcz(r,p,w,z),z)/ (om
     $                          *rad) -(k1land*anb+k2land*(ze-anb))/ze
     $                          *cos(-asr(r)))*ze +rv(r,p,w)*tt31
     $                          *(Aland+(1.D+0-Aland)*exp(-Lland*ze))
                        endif
                     else
                        zs0=0.D+0
                        if(ze.gt.(pi/2.D+0)) zs0=K2r0land*(ze-pi/2.D+0)

                        zs1=K1r1land*ze
                        if(ze.gt.anb) zs1=K1r1land*anb + K2r1land*(ze
     $                       -anb)

                        pcz(r,p,w,z)=rr_z(r)+rad*(zs0+rv(r,p,w)/rad*(zs1
     $                       -zs0))

c$$$                        pcz(r,p,w,z)=rr_z(r)+rad*(muk+ veloz(pcx(r,p,w,z
c$$$     $                       ),pcy(r,p,w,z),pcz(r,p,w,z),z)/ (om*rad)
c$$$     $                       -(zs0+rv(r,p,w)/rad*(zs1-zs0))/ze*cos(
c$$$     $                       -asr(r)))*ze +rv(r,p,w)*tt31 *(Aland+(1.D+0
c$$$     $                       -Aland)*exp(-Lland*ze))

                     endif
                  endif         !if(mu.eq.0)...
               end do           !z
            end do              !w
         end do                 !p
      end do                    !r

      return
      end
