c----------------------------------------------------------------------|
c    <Copyright (c) 1998 by Mahendra Bhagwat>
c----------------------------------------------------------------------|
c     prescribed wake based on Landgrebe's model
c
      subroutine prescr(nr,np,nw,pw,nz,nzt,nk,ns,nb,pif,zif,ip,dp,dz,
     &                sonic,vinfx,vinfy,vinfz,om,asr,b0r,b1c,b1s,ct,
     &                anb,rad,rcb,gb,crd,twt,
     &                lin,gv,rc0,rv,kinv,dcy,gav,sco,sce,gsc,csc,
     &                rpar,rrn,zrn,rbc,zbc,
     &                avgzbc,lami,rr_x,rr_y,rr_z,spx,spy,spz,
     &                cpx,cpy,cpz,pcx,pcy,pcz,
     &                taperst,taper,
     &                pox,poy,poz,ppx,ppy,ppz,bvx,bvy,bvz,itn)
c----------------------------------------------------------------------|
c
      IMPLICIT none
c
        include 'adimpar.inc'
        include 'gdimpar.inc'
        include 'jdimpar.inc'
        include 'kdimpar.inc'
        include 'rdimpar.inc'
        include 'sdimpar.inc'
        include 'wdimpar.inc'
        include 'zdimpar.inc'
c
      INTEGER nr,np,nw,pw,nz,nzt,nk,ns,nb,pif,zif,
     &        r,p,w,z,n,s
c
      DOUBLE PRECISION ip(rdim),dp,dz,vinfx,vinfy,vinfz,om,anb,rad,
     &       rcb,lin,kinv,dcy,gav,rpar,rrn,zrn,sce,
     &       ze,wage,tt11,tt21,tt31,b0,sonic,pbar,qbar,pi
c
      DOUBLE PRECISION rbc(rdim,wdim),zbc(rdim,wdim),
     &       avgzbc(rdim,wdim)
c
      DOUBLE PRECISION asr(rdim),b0r(rdim),lami(rdim),
     &       rr_x(rdim),rr_y(rdim),rr_z(rdim),
     &       spx(sdim),spy(sdim),spz(sdim),
     &       b1c(rdim),b1s(rdim),ct(rdim)
c
      DOUBLE PRECISION rv(rdim,adim,wdim),
     &       gb(rdim,adim,sdim),
     &       cpx(rdim,adim,sdim),
     &       cpy(rdim,adim,sdim),
     &       cpz(rdim,adim,sdim),
     &       vibx(rdim,adim,sdim),
     &       viby(rdim,adim,sdim),
     &       vibz(rdim,adim,sdim),
     &       bvx(rdim,adim,sdim),
     &       bvy(rdim,adim,sdim),
     &       bvz(rdim,adim,sdim)
c
      DOUBLE PRECISION gv(rdim,adim,wdim,zdim,kdim),
     &       rc0(rdim,adim,wdim,kdim),
     &       ppx(rdim,adim,wdim,zdim),
     &       ppy(rdim,adim,wdim,zdim),
     &       ppz(rdim,adim,wdim,zdim),
     &       pcx(rdim,adim,wdim,zdim),
     &       pcy(rdim,adim,wdim,zdim),
     &       pcz(rdim,adim,wdim,zdim),
     &       pox(rdim,adim,wdim,zdim),
     &       poy(rdim,adim,wdim,zdim),
     &       poz(rdim,adim,wdim,zdim),
     &       vipx(rdim,adim,wdim,zdim),
     &       vipy(rdim,adim,wdim,zdim),
     &       vipz(rdim,adim,wdim,zdim),
     &       vicx(rdim,adim,wdim,zdim),
     &       vicy(rdim,adim,wdim,zdim),
     &       vicz(rdim,adim,wdim,zdim)
c
      CHARACTER*1 sco
c
      DOUBLE PRECISION gsc(rdim,adim,wdim,zdim,kdim),
     &       csc(rdim,adim,wdim,zdim,kdim)
c
      DOUBLE PRECISION zec,Aland,Lland,k1land,k2land,crd,twt,
     &       K1r1land,K2r1land,K2r0land,zs0,zs1,sn,az
c
      INTEGER itn
      DOUBLE PRECISION BKT,CKT,mKT,nKT,CT0KT
c
c
      DOUBLE PRECISION planf
      DOUBLE PRECISION taperst,taper
c
      external planf
c
      external velox
      external veloy
      external veloz
c----------------------------------------------------------------------|
c
      common/pidata/pi
      save/pidata/
c
      common/ssmdata/pbar,qbar
      save/ssmdata/
c----------------------------------------------------------------------|
c
      do r=1,nr

         pi=4.D+0*atan(1.D+0)
         Aland=0.78D+0    !1.D+0/sqrt(2.D+0)
         Lland=0.145D+0+27.D+0*ct(r)
         k1land=0.25D+0*( ct(r)/
     &                  (2.D+0/anb*crd*planf(-1.D+0,taperst,taper)/rad) 
     &                   + 0.001D+0*twt )
         k2land=sqrt(ct(r))*(1.D+0+0.01D+0*twt)
         K1r1land=-2.2D+0*sqrt(ct(r)/2.D+0)
         K2r1land=-2.7D+0*sqrt(ct(r)/2.D+0)
         K2r0land=(twt/128.D+0*(0.45D+0*twt+18.D+0))*sqrt(ct(r)/2.D+0)

c-mb Kocurek-Tangler

      Aland=0.78D+0
      Lland=4.D+0*sqrt(ct(r))
      BKT=-0.000729D+0*twt
      CKT=-2.3D+0+0.206D+0*twt
      mKT=1.D+0-0.25D+0*exp(0.04D+0*twt)
      nKT=0.5D+0-0.0172D+0*twt
      CT0KT=Nb**nKT*(-BKT/CKT)**(1/mKT)

      k1land=BKT + CKT*(ct(r)/Nb**nKT)**mKT
      k2land=-sqrt(ct(r)-CT0KT)

c-mb note that the sign convention is opposite to that of the paper
c-mb for the slopes because I was following Landgrebe earlier
      k1land=-k1land
      k2land=-k2land

      sn=(-1.D+0)**(r-1)
      do p=1,np
      az=ip(r)+DBLE(p-1)*dp/DBLE(pif) - ip(r)
      b0=b0r(r)+b1c(r)*cos(az)+b1s(r)*sin(az)
      do w=1,nw
      do z=1,nzt
      ze=DBLE(z-1)*dz/DBLE(zif)
      zec=DBLE(z-1)*dz/DBLE(zif)
      if(z.gt.nz) zec=DBLE(nz-1)*dz/DBLE(zif)
      wage=sn*(az + ip(r) -ze)

      tt11=cos(-asr(r))*cos(wage)*cos(b0)+
     &     sin(-asr(r))*sin(b0)
      tt21=sin(wage)*cos(b0)
      tt31=-sin(-asr(r))*cos(wage)*cos(b0)+
     &      cos(-asr(r))*sin(b0)

      pcx(r,p,w,z)=rr_x(r)+rv(r,p,w)*tt11
     &             *(Aland+(1.D+0-Aland)*exp(-Lland*ze))

      pcy(r,p,w,z)=rr_y(r)+rv(r,p,w)*tt21
     &             *(Aland+(1.D+0-Aland)*exp(-Lland*ze))

c-mb  Tip vortex:
      if(w.eq.1) then

      pcz(r,p,w,z)=rr_z(r)-rad*k1land*ze
      if(ze.gt.anb) 
     & pcz(r,p,w,z)=rr_z(r)-rad*(k1land*anb+k2land*(ze-anb))

      else
c-mb inboard sheet

      zs0=0.D+0
      if(ze.gt.(pi/2.D+0)) zs0=K2r0land*(ze-pi/2.D+0)

      zs1=K1r1land*ze
      if(ze.gt.anb) zs1=K1r1land*anb + K2r1land*(ze-anb)

      pcz(r,p,w,z)=rr_z(r)+rad*(zs0+rv(r,p,w)/rad*(zs1-zs0))
c-mb   note        -->  _^_

      endif

      end do
      end do
      end do

      end do
c
      write(*,*)'Prescribed wake iteration: ',itn
c----------------------------------------------------------------------|
c
      return
      end
